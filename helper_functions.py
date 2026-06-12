# helper_functions.py
import os
import re
import numpy as np
import pandas as pd
import scipy.ndimage
import pickle
import gzip
import lzma
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pathlib import Path
from scipy import fftpack
from scipy.stats import kurtosis, skew, entropy
from scipy.ndimage import label as label_region, gaussian_filter
from scipy.ndimage import convolve1d, map_coordinates
from astropy.io import fits
from astropy.time import Time
from sunpy.net import Fido
from sunpy.net import attrs as a
from typing import Optional, Union, List, Tuple
from datetime import datetime, timedelta, timezone
from tqdm import tqdm
from collections import defaultdict
    
    


def pad_3d(a, n):
    """
    Pad edges of a 3D (or N-D) NumPy array `a` with `n` zeros.
    
    Example:
        If a.shape = (30, 30, 30) and n = 5,
        returns an array of shape (40, 40, 40)
        with `a` placed at [5:35, 5:35, 5:35].
    
    Parameters
    ----------
    a : np.ndarray
        Input N-dimensional array.
    n : int
        Number of zeros to pad on each edge.

    Returns
    -------
    np.ndarray
        Zero-padded array.
    """
    a = np.asarray(a)
    ndim = a.ndim
    pad_width = [(n, n)] * ndim
    a2 = np.pad(a, pad_width, mode='constant', constant_values=0)
    return a2




def unpad_3d(a, n):
    """
    Remove a margin of n pixels from each edge of a 3D (or N-D) NumPy array.
    
    Example:
        If a.shape = (40, 40, 40) and n = 5,
        returns an array of shape (30, 30, 30)
        corresponding to a[5:35, 5:35, 5:35].
    
    Parameters
    ----------
    a : np.ndarray
        Input N-dimensional array (e.g., 3D image cube).
    n : int
        Number of pixels to remove from each edge.

    Returns
    -------
    np.ndarray
        Cropped array with margins removed.
    """
    a = np.asarray(a)
    ndim = a.ndim

    # Build slicing tuple dynamically for arbitrary dimensions
    slices = tuple(slice(n, -n if n != 0 else None) for _ in range(ndim))

    return a[slices]




def unpad_2d(a, n):
    """
    Remove a margin of n pixels from each edge of a 2D NumPy array.
    
    Example:
        If a.shape = (40, 40) and n = 5,
        returns an array of shape (30, 30)
        corresponding to a[5:35, 5:35].
    
    Parameters
    ----------
    a : np.ndarray
        Input 2D array.
    n : int
        Number of pixels to remove from each edge.
    
    Returns
    -------
    np.ndarray
        Cropped array with margins removed.
    """
    a = np.asarray(a)
    if a.ndim != 2:
        raise ValueError("unpad_2d expects a 2D array")

    return a[n:-n, n:-n]




def pad_2d(a, n):
    """
    Pad edges of a 2D NumPy array `a` with `n` zeros.
    
    Example:
        If a.shape = (30, 30) and n = 5,
        returns an array of shape (40, 40)
        with `a` placed at [5:35, 5:35].
    
    Parameters
    ----------
    a : np.ndarray
        Input 2D array.
    n : int
        Number of zeros to pad on each edge.
    
    Returns
    -------
    np.ndarray
        Zero-padded array.
    """
    a = np.asarray(a)
    if a.ndim != 2:
        raise ValueError("pad_2d expects a 2D array")

    pad_width = ((n, n), (n, n))
    return np.pad(a, pad_width, mode='constant', constant_values=0)




def one2n(index: Union[int, List[int], np.ndarray],
          array: Optional[np.ndarray] = None,
          minmax: bool = False,
          dimensions: Optional[List[int]] = None
         ) -> Tuple[np.ndarray, ...]:
    """
    Convert 1D flat indices to multi-dimensional indices, mimicking IDL's one2n.
    
    Returns exactly as many arrays as the number of dimensions of `array` or `dimensions`.
    
    Parameters
    ----------
    index : int or array_like
        Flattened indices.
    array : np.ndarray, optional
        The array from which indices are derived.
    minmax : bool, optional
        If True, returns [min, max] for each dimension.
    dimensions : list[int], optional
        Shape of the array, if `array` is not provided.
    
    Returns
    -------
    tuple of np.ndarray
        One array per dimension, e.g., (i0, i1, i2) for a 3D array.
    """
    index = np.atleast_1d(index)

    # Determine array shape
    if dimensions is not None:
        shape = tuple(dimensions)
    elif array is not None:
        shape = array.shape
    else:
        raise ValueError("Must provide either `array` or `dimensions`.")

    ndim = len(shape)
    multi_idx = np.unravel_index(index, shape)  # tuple of arrays per dimension

    # Apply minmax if requested
    if minmax:
        multi_idx = tuple(np.array([np.min(arr), np.max(arr)]) for arr in multi_idx)

    # Return exactly ndim arrays
    return multi_idx



def region_size_filter(image, threshold, npix,
                       label=None, mask=None,
                       max_region_size=None, all_neighbors=1,
                       return_mask=False, smooth=None, npixinitial=None):
    """
    Filter image regions by size with optional multi-dimensional Gaussian smoothing.

    Parameters
    ----------
    image : ndarray
        Input 2D or 3D image.
    threshold : float
        Minimum value for a pixel/voxel to be considered part of a region.
    npix : int
        Minimum number of pixels/voxels in a region to keep.
    label : ndarray, optional
        Label array for connected components (will be computed if not provided).
    mask : ndarray, optional
        Boolean mask of pixels/voxels above threshold (will be computed if not provided).
    max_region_size : int, optional
        Maximum number of pixels/voxels in a region to keep (if provided).
    all_neighbors : int, default 1
        Connectivity for labeling (1=faces, 2=faces+edges, 3=faces+edges+corners in 3D).
    return_mask : bool, default False
        If True, return the boolean mask instead of masked image.
    smooth : float or list of floats, optional
        Standard deviation(s) for Gaussian smoothing along each axis.
    npixinitial : int, optional
        Remove regions smaller than this before smoothing.

    Returns
    -------
    ndarray
        Masked image (same shape as input) or boolean mask if return_mask=True.

    Examples
    --------
    2D Example:

    >>> import numpy as np
    >>> image = np.zeros((6,6))
    >>> image[0:2, 0:2] = 5     # small region
    >>> image[3:6, 3:6] = 10    # large region
    >>> filtered = region_size_filter(image, threshold=1, npix=5)
    >>> print(filtered)
    [[0. 0. 0. 0. 0. 0.]
     [0. 0. 0. 0. 0. 0.]
     [0. 0. 0. 0. 0. 0.]
     [0. 0. 0. 10. 10. 10.]
     [0. 0. 0. 10. 10. 10.]
     [0. 0. 0. 10. 10. 10.]]

    Explanation:
    - The small 2x2 region (value 5) is removed because it has fewer than 5 pixels.
    - The large 3x3 region (value 10) is kept because it has more than 5 pixels.
    - Padding/unpadding ensures edge effects do not influence the labeling.
    """

    ndim = image.ndim
    if ndim > 3:
        print("Region size filter only implemented for 2 or 3 dimensions, returning -1")
        return -1

    # Initial threshold mask with padding
    if ndim == 2:
        mask = pad_2d(image, 1) > threshold
    elif ndim == 3:
        mask = pad_3d(image, 1) > threshold

    # Remove small regions initially if requested
    if npixinitial is not None:
        lbl, n_features = label_region(mask, structure=np.ones((3,)*ndim))
        counts = np.bincount(lbl.ravel())
        too_small = np.where(counts < npixinitial)[0]
        if too_small.size > 0:
            mask[np.isin(lbl, too_small)] = 0

    # Optional Gaussian smoothing
    if smooth is not None:
        sigma = np.full(ndim, smooth) if np.ndim(smooth) == 0 else np.array(smooth)
        mask_float = mask.astype(float)
        masksmo = gaussian_filter(mask_float, sigma=sigma, mode='nearest')
        masksmo = masksmo > 1e-2
    else:
        masksmo = mask

    if not np.any(masksmo):
        return np.zeros_like(image)

    # Label connected regions
    lbl, n_features = label_region(masksmo, structure=np.ones((3,)*ndim))
    lbl *= mask  # keep original mask

    counts = np.bincount(lbl.ravel())
    # Keep only regions within size limits
    if max_region_size is not None:
        reject = np.where((counts < npix) | (counts > max_region_size))[0]
    else:
        reject = np.where(counts < npix)[0]

    if reject.size > 0:
        mask[np.isin(lbl, reject)] = 0

    # Recompute labels for final mask
    lbl, _ = label_region(mask, structure=np.ones((3,)*ndim))

    # Unpad
    if ndim == 2:
        mask = unpad_2d(mask, 1)
        lbl = unpad_2d(lbl, 1)
    elif ndim == 3:
        mask = unpad_3d(mask, 1)
        lbl = unpad_3d(lbl, 1)

    return mask if return_mask else image * mask




def gamma_transform(imin, gamma=3.5, a0=None, a1=None):
    """Perform gamma transformation on input array."""
    imin = imin.astype(float)
    if a0 is None:
        a0 = np.nanmin(imin)
    if a1 is None:
        a1 = np.nanmax(imin)
    return ((imin - a0) / (a1 - a0)) ** (1.0 / gamma)




def gaussian_function(sigma, normalize=True):
    """
    Return a 1D Gaussian kernel with standard deviation sigma.
    Normalized to sum=1 by default.
    """
    radius = int(np.ceil(sigma * 3))
    if radius % 2 == 0:
        radius += 1
    x = np.arange(-radius//2 + 1, radius//2 + 1)
    g = np.exp(-(x**2)/(2*sigma**2))
    if normalize:
        g /= g.sum()
    return g




def mgn(imin, a0=None, a1=None, gamma=3.5, h=0.9, k=1.0):
    """
    Multi-scale Gaussian Normalization (MGN) in Python.

    Parameters
    ----------
    imin : ndarray
        2D input image array.
    a0, a1 : float
        Optional min/max clipping for gamma transform.
    gamma : float
        Gamma for global normalization (default 3.5).
    h : float
        Weighting for gamma-transformed image vs local MGN (default 0.9).
    k : float
        Contrast stretching factor for final arctangent transform (default 1.0).

    Returns
    -------
    imout : ndarray
        2D processed image array.
    """
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    
    imin2 = imin.astype(float)
    
    # Global gamma transform
    im = gamma_transform(imin2, gamma, a0, a1)

    # Gaussian kernel widths (different spatial scales)
    w = [2.5, 5, 10, 20, 40, 80]
    nw = len(w)

    s = imin2.shape
    imp = np.zeros(s, dtype=float)

    for iw in range(nw):
        ker = gaussian_function(w[iw] / 4.0, normalize=True)
        # Local mean via separable convolution
        m = scipy.ndimage.convolve1d(
            scipy.ndimage.convolve1d(imin2, ker, axis=0, mode='mirror'), 
            ker, axis=1, mode='mirror'
        )
        # Local standard deviation
        md = np.sqrt(
            scipy.ndimage.convolve1d(
                scipy.ndimage.convolve1d((imin2 - m)**2, ker, axis=0, mode='mirror'), 
                ker, axis=1, mode='mirror'
            )
        )
        const = np.nanmean(md) / 10.0
        imp += (imin2 - m) / (md + const)

    # Average over scales
    imp /= nw

    # Arctangent contrast stretch
    imp = np.arctan(k * imp)

    # Combine global gamma and MGN
    imout = h * im + (1 - h) * imp

    return imout
    
    
    

def gaussian_function(sigma, normalize=True):
    """Return a 1D gaussian kernel (same behaviour as earlier helper)."""
    radius = int(np.ceil(sigma * 3))
    if radius % 2 == 0:
        radius += 1
    x = np.arange(-radius//2 + 1, radius//2 + 1)
    g = np.exp(-(x**2) / (2.0 * sigma**2))
    if normalize:
        s = g.sum()
        if s != 0:
            g = g / s
    return g




def normalize_to_uint8(arr):
    """
    Normalize a numeric NumPy array to the 8-bit range [0, 255].

    This function rescales the values of a 2D or 3D array to the range [0, 1],
    handles NaNs and infinities safely, and then converts the result to unsigned
    8-bit integers (uint8). If the array has no variation (i.e. min == max),
    the output will be a zero array.

    Parameters
    ----------
    arr : np.ndarray
        Input numeric array (e.g. float or int). Can be 2D or 3D.

    Returns
    -------
    np.ndarray
        The normalized array, converted to dtype uint8 and scaled to [0, 255].

    Examples
    --------
    >>> import numpy as np
    >>> a = np.array([[0.0, 0.5], [1.0, 2.0]])
    >>> normalize_to_uint8(a)
    array([[  0,  63],
           [127, 255]], dtype=uint8)
    """
    arr = np.nan_to_num(arr)  # handle NaNs or infs
    amin, amax = arr.min(), arr.max()
    if amax > amin:  # avoid divide-by-zero
        arr = (arr - amin) / (amax - amin)
    else:
        arr = np.zeros_like(arr)
    return (arr * 255).astype(np.uint8)




def split_by_wavelength(files: list[str], status: list[bool]) -> dict[int, list[str]]:
    """
    Split a list of AIA FITS filenames into wavelength-based groups,
    sorted by timestamp, and filter using status flags.

    Parameters
    ----------
    files : list of str
        Full paths to FITS files.
    status : list of bool
        Whether each file is valid (True) or should be excluded (False).

    Returns
    -------
    grouped_sorted : dict[int, list[str]]
        Keys: integer wavelengths (94, 131, 171, ...)
        Values: sorted list of *filenames only*, not full paths.
    """

    pattern = re.compile(r"AIA(\d{8}_\d{4})_(\d{4})\.fits$")
    grouped = defaultdict(list)

    for f, ok in zip(files, status):
        if not ok:
            continue  # skip missing/failed downloads

        base = os.path.basename(f)
        m = pattern.search(base)
        if not m:
            continue

        timestamp_str, wavelength_str = m.groups()
        t = datetime.strptime(timestamp_str, "%Y%m%d_%H%M")
        wavelength = int(wavelength_str)

        # only store the filename; run_wave builds the full path later
        grouped[wavelength].append((t, base))

    # sort each wavelength list by timestamp
    grouped_sorted = {}
    for wl in sorted(grouped.keys()):
        grouped_sorted[wl] = [
            fname for _, fname in sorted(grouped[wl], key=lambda x: x[0])
        ]

    return grouped_sorted




def compute_time_strings(end_time: str, offset=8):
    """
    Given an end_time string 'YYYY-MM-DD HH:MM', compute:
      - start_time (8 hours before) as 'YYYY-MM-DD HH:MM'
      - out string for output directories as 'YYYY_MM_DD'
    Returns:
      start_time_str, out_str
    """
    dt = datetime.strptime(end_time, "%Y-%m-%d %H:%M")
    start_time = dt - timedelta(hours=offset)
    start_time_str = start_time.strftime("%Y-%m-%d %H:%M")
    out_str = dt.strftime("%Y_%m_%d")
    return start_time_str, out_str




def snap_to_cadence(dt: datetime, cadence_minutes: int = 6) -> datetime:
    total_minutes = dt.hour*60 + dt.minute
    snapped_minutes = (total_minutes // cadence_minutes) * cadence_minutes
    new_hour = snapped_minutes // 60
    new_minute = snapped_minutes % 60
    return dt.replace(hour=new_hour, minute=new_minute, second=0, microsecond=0)




def kurtosis_timeseries_cube(datacube, fisher=True, bias=False, nan_policy="omit"):
    nt = datacube.shape[2]
    kappa = np.empty(nt, dtype=np.float64)

    for t in range(nt):
        frame = datacube[:, :, t].ravel()
        kappa[t] = kurtosis(frame, fisher=fisher, bias=bias, nan_policy=nan_policy)

    return kappa




def skew_timeseries_cube(datacube, bias=False, nan_policy="omit"):
    nt = datacube.shape[2]
    s = np.empty(nt, dtype=np.float64)

    for t in range(nt):
        frame = datacube[:, :, t].ravel()
        s[t] = skew(frame, bias=bias, nan_policy=nan_policy)

    return s




def shannon_entropy_timeseries_cube(datacube, bins=256):
    nt = datacube.shape[2]
    H = np.empty(nt, dtype=np.float64)

    for t in range(nt):
        frame = datacube[:, :, t].ravel()
        frame = frame[np.isfinite(frame)]

        if frame.size == 0:
            H[t] = np.nan
            continue

        hist, _ = np.histogram(frame, bins=bins, density=True)
        hist = hist[hist > 0]

        if hist.size == 0:
            H[t] = np.nan
        else:
            H[t] = entropy(hist)

    return H




def save_daily_kurtosis(
    df_new,
    wave,
    outdir,
    running_window_points,
    min_periods
):

    stat_cols = [
        "kurtosis_intensity",
        "kurtosis_ratio",
        "skew_intensity",
        "skew_ratio",
        "entropy_intensity",
        "entropy_ratio",
    ]

    df_new["date"] = df_new["time"].dt.date

    for date, g in df_new.groupby("date"):

        outfile = os.path.join(
            outdir,
            f"kurtosis_{date}_wvlnth_{wave}.pkl.gz"
        )

        # ----------------------------
        # Load existing file if present
        # ----------------------------

        if os.path.exists(outfile):

            with gzip.open(outfile, "rb") as fp:
                df_old = pickle.load(fp)

            df_all = pd.concat([df_old, g], ignore_index=True)

        else:
            df_all = g.copy()

        # ----------------------------
        # Enforce schema
        # ----------------------------

        df_all = df_all[df_new.columns]

        # ----------------------------
        # Remove duplicate timestamps
        # ----------------------------

        df_all = (
            df_all
            .set_index("time")
            .sort_index()
            .loc[lambda x: ~x.index.duplicated(keep="last")]
            .reset_index()
        )

        # ----------------------------
        # Recompute rolling statistics
        # ----------------------------

        for col in stat_cols:

            mu = df_all[col].rolling(
                running_window_points,
                min_periods=min_periods
            ).mean()

            sig = df_all[col].rolling(
                running_window_points,
                min_periods=min_periods
            ).std()

            sig = sig.replace(0, np.nan)

            df_all[f"{col}_mu"] = mu
            df_all[f"{col}_sig"] = sig
            df_all[f"{col}_z"] = (df_all[col] - mu) / sig

        # ----------------------------
        # Save file
        # ----------------------------

        with gzip.open(outfile, "wb") as fp:
            pickle.dump(df_all.drop(columns="date"), fp)




def update_kurtosis_file(
    wave,
    im_cube,
    ratio_cube,
    hdrs,
    outdir="./kurtosis",
    running_window_points=10,
    min_periods=3
):

    os.makedirs(outdir, exist_ok=True)

    # ----------------------------
    # Observation times
    # ----------------------------

    times = pd.to_datetime(
        [h["T_OBS"] for h in hdrs],
        errors="coerce",
        utc=True
    )

    # ----------------------------
    # Raw statistics
    # ----------------------------

    k_im = kurtosis_timeseries_cube(im_cube)
    k_ratio = kurtosis_timeseries_cube(ratio_cube)

    s_im = skew_timeseries_cube(im_cube)
    s_ratio = skew_timeseries_cube(ratio_cube)

    H_im = shannon_entropy_timeseries_cube(im_cube)
    H_ratio = shannon_entropy_timeseries_cube(ratio_cube)

    # ----------------------------
    # Build dataframe
    # ----------------------------

    df_new = pd.DataFrame({

        "time": times,

        "kurtosis_intensity": k_im,
        "kurtosis_ratio": k_ratio,

        "skew_intensity": s_im,
        "skew_ratio": s_ratio,

        "entropy_intensity": H_im,
        "entropy_ratio": H_ratio

    })

    df_new = (
        df_new
        .dropna(subset=["time"])
        .sort_values("time")
        .reset_index(drop=True)
    )

    # optional safety
    df_new = df_new.dropna()

    # ----------------------------
    # Save daily files
    # ----------------------------

    save_daily_kurtosis(
        df_new,
        wave,
        outdir,
        running_window_points,
        min_periods
    )

    return df_new




def kurt_skew_entr(end_time=None,
                   start_time=8,
                   wavelengths=[94, 131, 171, 193, 211, 304, 335],
                   sharp_region=None,
                   plot_correlations=True,
                   kurt_range=None,
                   skew_range=None,
                   entr_range=None,
                   wind_range=None,
                   helc_range=None
                  ):
    """
    Generate multi-panel time series plots of kurtosis, skewness, and entropy
    (raw and smoothed) for multiple EUV wavelengths within a specified time window,
    together with ARTop-derived deltaL and deltaH signals and event markers from
    ALMANAC CME detections and HEK flare events.

    The function reads precomputed statistical outputs for each wavelength,
    aligns them within a common analysis window, applies optional smoothing,
    and produces publication-style PDF figures containing:

        • Raw intensity statistics
        • Smoothed intensity statistics
        • Raw ratio statistics
        • Smoothed ratio statistics
        • deltaL signal with running mean ±1σ/±2σ
        • deltaH signal with running mean ±1σ/±2σ

    If `plot_correlations=True`, the function additionally computes and plots
    rolling inter-wavelength correlations for each statistic (intensity and ratio),
    including both raw and Z-score variants.

    Parameters
    ----------
    end_time : str or datetime, optional
        End time of the analysis window. Can be any valid datetime string
        (e.g., '2011-02-15 23:59') or a datetime object.
        The analysis window extends backward by `start_time` hours.

    start_time : float, optional
        Duration in hours prior to `end_time` defining the start of the
        analysis window. Default is 8 hours.

    wavelengths : list of int, optional
        EUV wavelengths (Å) to process. Default is:
        [94, 131, 171, 193, 211, 304, 335].

    sharp_region : int, optional
        HMI SHARP active region number used to locate ARTop outputs
        (deltaL/deltaH maps).

    plot_correlations : bool, optional
        If True (default), compute and generate additional figures showing
        rolling inter-wavelength correlations for each statistic.

    Returns
    -------
    list
        A two-element list:
            [dfs, correlations]

        dfs : dict
            Dictionary of pandas DataFrames keyed by wavelength.
            Each DataFrame contains time-filtered statistics for
            kurtosis, skewness, and entropy (raw and Z-score).

        correlations : dict
            Dictionary of rolling inter-wavelength correlation time series
            (empty if `plot_correlations=False`).

    Expected Directory Structure
    ----------------------------
    ./kurtosis/kurtosis_<YYYY-MM-DD>_wvlnth_<wv>.pkl.gz
    ./ARTop/AR_<sharp_region>_Output/*metadata*
    ./ARTop/AR_<sharp_region>_Output/*deltaL*
    ./ARTop/AR_<sharp_region>_Output/*deltaH*
    ./detections/detection_<YYYY-MM-DD_HHMM>.pkl

    Notes
    -----
    • Time alignment across wavelengths is handled internally.
    • Z-score panels include optional ±1σ and ±2σ shading.
    • deltaL and deltaH signals are clipped to non-negative values.
    • Correlation plots use rolling window statistics with pairwise
      complete observations.
    • All figures are formatted with synchronized time axes and
      annotated event markers for CME and flare classifications (C, M, X).
    """
    
    '''
    READ IN ALL KURTOSIS FILES AND MAKE DATAFRAMES
    '''
    # Convert to datetime (this is a DatetimeIndex)
    cme_times = pd.to_datetime(end_time)
    
    # Make them Series so concat works
    times = pd.concat([
        pd.Series(cme_times),
        pd.Series(cme_times - pd.Timedelta(hours=start_time))
    ])
    
    # Extract unique dates
    dates = sorted(times.dt.strftime('%Y-%m-%d').unique().tolist())
    ftimes = sorted(times.dt.strftime('%Y-%m-%d %H:%M').unique().tolist())
    
    dfs = {}  # dictionary to hold combined dataframes per wavelength
    
    for wv in wavelengths:
        combined = []
        for d in dates:
            path = Path(f"./kurtosis/kurtosis_{d}_wvlnth_{wv}.pkl.gz")
            if path.exists():
                with gzip.open(path, "rb") as f:
                    df = pickle.load(f)
                    combined.append(df)
            else:
                print(f"WARNING: Missing file {path}")
        if combined:
            # concat along rows and sort by time
            dfs[wv] = pd.concat(combined,
                ignore_index=True
                ).sort_values("time").reset_index(drop=True)
        else:
            print(f"No files found for wavelength {wv}")
    
    
    '''
    FIND ALL CME AND FLARE TIMES FROM ALMANAC AND HEK
    '''
    # Define time window
    t_start = Time(f"{ftimes[0]}")
    t_end   = Time(f"{ftimes[-1]}")
    
    det_dir = Path(os.path.join(os.getcwd(),"detections"))
    t_dt,t_almanac,almanac_files = [],[],[]
    
    for file in det_dir.glob("detection_*.pkl"):
        # Remove prefix and suffix
        time_str = file.stem.replace("detection_", "")  # '2011-02-15_0520'
        
        # Parse into datetime
        dt = datetime.strptime(time_str, "%Y-%m-%d_%H%M")
        
        # Convert to Astropy Time
        t_det = Time(dt)
    
        # Keep only if within window
        if t_start <= t_det <= t_end:
            t_dt.append(t_det.to_datetime())
            t_almanac.append(time_str)
            almanac_files.append(file)
    
    # HEK flare search
    result = Fido.search(
        a.Time(t_start, t_end),
        a.hek.EventType("FL")  # FL = flare
    )
    
    hek_table = result['hek']
    X_times, M_times, C_times = [], [], []
    
    if len(hek_table) != 0:
        for i, time in enumerate(hek_table['event_starttime']):
            fclass = hek_table['fl_goescls'][i]
            try:
                if fclass[0] == 'C':
                    C_times.append(time.to_datetime())
                elif fclass[0] == 'M':
                    M_times.append(time.to_datetime())
                elif fclass[0] == 'X':
                    X_times.append(time.to_datetime())
            except Exception:
                pass
    
    # Determine date/time range across all wavelengths
    all_times = pd.concat([df['time'] for df in dfs.values()])
    start_str = all_times.min().strftime("%Y%m%d_%H%M")
    end_str   = all_times.max().strftime("%Y%m%d_%H%M")
    
    # Ensure output directory exists
    outdir = Path("./kurtosis")
    outdir.mkdir(exist_ok=True, parents=True)

   
    '''
    CREATE PLOTS
    '''
    # Colors for intensity and ratio (7 wavelengths)
    colors_int = ['navy', 'purple', 'blue', 'green', 'orange', 'brown', 'red']
    colors_ratio = ['cyan', 'magenta', 'lightblue', 'lime', 'gold', 'sandybrown', 'pink']
    
    cwd = os.getcwd()
    if sharp_region != None:
        sharp  = os.path.join(cwd,'ARTop',f'AR_{sharp_region}_Output')
    
        metadata_files = list(Path(sharp).glob("*metadata*"))
        
        if metadata_files:
            with lzma.open(str(list(Path(sharp).rglob("*metadata*"))[0]), 'rb') as f:
                meta = pickle.load(f)
            
            with lzma.open(str(list(Path(sharp).rglob("*deltaL*"))[0]), 'rb') as f:
                deltaL = pickle.load(f)
            
            with lzma.open(str(list(Path(sharp).rglob("*deltaH*"))[0]), 'rb') as f:
                deltaH = pickle.load(f)
    
            def maps_to_timeseries(maps, reducer=np.nanmean):
                """
                Convert list of SunPy maps to 1D time series + times.
                Skips missing / None maps safely.
                """
            
                times = []
                values = []
            
                for m in maps:
                    if m is None:
                        continue
            
                    if not hasattr(m, "data"):
                        continue
            
                    values.append(reducer(m.data))
                    times.append(pd.to_datetime(m.date.datetime, utc=True))
            
                return pd.DatetimeIndex(times), np.asarray(values)
            
            
            time, dL = maps_to_timeseries(deltaL)
            time, dH = maps_to_timeseries(deltaH)
        
        else:
            header_path = Path(os.path.join(sharp,"header.txt"))
            with header_path.open() as f:
                lines = [line.strip() for line in f.readlines()]
            
            header = {
                'rpx': lines[0],
                'rpy': lines[1],
                'observatory': lines[2],
                'instrument': lines[3],
                'detector': lines[4],
                'obsTime': lines[5],
                'rcx': lines[6],
                'rcy': lines[7],
                'scaleX': lines[8],
                'scaleY': lines[9],
            }
            
            spec_path = Path(os.path.join(sharp,"specifications.txt"))
            with spec_path.open() as f:
                lines = [line.strip() for line in f.readlines()]
            
            with spec_path.open() as f:
                lines = [line.strip() for line in f.readlines()]
            
            specification = {
                'sharp_number': lines[0],
                'nx': lines[1],
                'ny': lines[2],
                'nt': lines[3],
                'vs': lines[4],
                'co': lines[5],
                'd': lines[6],
            }
            
            
            import zarr
            root = zarr.open(os.path.join(sharp,f"topology_{sharp_region}.zarr"), mode="r")
            
            def cube_to_timeseries(cube, header, specification, reducer=np.nanmean):
                """
                Convert a 3D Zarr array (nx, ny, nt) into a 1D time series
                with automatically generated 720 s cadence timestamps.
            
                Parameters
                ----------
                cube : zarr.core.array.Array or ndarray
                    3D array with shape (nx, ny, nt)
                header : dict
                    Header dictionary containing 'obsTime'
                specification : dict
                    Specification dictionary containing 'nt'
                reducer : function
                    Reduction function applied over spatial axes (default: nanmean)
            
                Returns
                -------
                pd.DatetimeIndex
                    Time index of length nt
                np.ndarray
                    1D time series of length nt
                """
            
                # Reduce spatial dimensions (0,1), keep time (2)
                values = reducer(cube, axis=(0, 1))
            
                # Parse start time
                start_time = pd.to_datetime(header['obsTime'], utc=True)
            
                # Number of timesteps
                nt = int(specification['nt'])
            
                # Create datetime index with 720-second cadence
                times = pd.date_range(
                    start=start_time,
                    periods=nt,
                    freq="720s",
                    tz="UTC"
                )
            
                return times, np.asarray(values)
            
            time,dL = cube_to_timeseries(root['dLFlux'],header,specification)
            time,dH = cube_to_timeseries(root['dHFlux'],header,specification)
    
        dL = np.clip(dL, 0, None)
        dH = np.clip(dH, 0, None)

    # -----------------------------
    # PARAMETERS
    # -----------------------------
    STAT_CONFIGS = {
        "Kurtosis": {
            "int_z": "kurtosis_intensity",
            "rat_z": "kurtosis_ratio",
            "label1": "Kurtosis Intensity",
            "label2": "Kurtosis Ratio",
            "fname": "kurtosis",
            "shade": False
        },
        "Kurtosis Z-score": {
            "int_z": "kurtosis_intensity_z",
            "rat_z": "kurtosis_ratio_z",
            "label1": "Kurtosis Intensity Z-score",
            "label2": "Kurtosis Ratio Z-score",
            "fname": "kurtosis_z",
            "shade": True
        },
        "Skew": {
            "int_z": "skew_intensity",
            "rat_z": "skew_ratio",
            "label1": "Skewness Intensity",
            "label2": "Skewness Ratio",
            "fname": "skew",
            "shade": False
        },
        "Skew Z-score": {
            "int_z": "skew_intensity_z",
            "rat_z": "skew_ratio_z",
            "label1": "Skewness Intensity Z-score",
            "label2": "Skewness Ratio Z-score",
            "fname": "skew_z",
            "shade": True
        },
        "Entropy": {
            "int_z": "entropy_intensity",
            "rat_z": "entropy_ratio",
            "label1": "Entropy Intensity",
            "label2": "Entropy Ratio",
            "fname": "entropy",
            "shade": False
        },
        "Entropy Z-score": {
            "int_z": "entropy_intensity_z",
            "rat_z": "entropy_ratio_z",
            "label1": "Entropy Intensity Z-score",
            "label2": "Entropy Ratio Z-score",
            "fname": "entropy_z",
            "shade": True
        }
    }

    smooth_z = 5          # smoothing for z-score
    dl_dh_window = 15     # running mean window for dL / dH
    
    outdir = Path("figures")
    outdir.mkdir(exist_ok=True, parents=True)
    
    # -----------------------------
    # TIME RANGE
    # -----------------------------
    # Define window
    end_dt = pd.to_datetime(end_time, utc=True)
    start_dt = end_dt - pd.Timedelta(hours=start_time)
    
    # Filter each dataframe
    for wv in wavelengths:
        df = dfs[wv]
        mask = (df['time'] >= start_dt) & (df['time'] <= end_dt)
        dfs[wv] = df.loc[mask].copy()
    
    all_times = pd.concat([dfs[w]["time"] for w in wavelengths]).sort_values()
    xmin = all_times.min() - pd.Timedelta(minutes=20)
    xmax = all_times.max() + pd.Timedelta(minutes=20)
    
    start_str = all_times.min().strftime("%Y%m%d_%H%M")
    end_str   = all_times.max().strftime("%Y%m%d_%H%M")


    for stat_name, cfg in STAT_CONFIGS.items():
        int_z = cfg["int_z"]
        rat_z = cfg["rat_z"]
        slabel1 = cfg["label1"]
        slabel2 = cfg["label2"]
        fname = cfg["fname"]
        shade = cfg["shade"]
        
        # -----------------------------
        # FIGURE LAYOUT
        # -----------------------------
        if sharp_region != None:
            plot_num = 6
        else:
            plot_num = 4
            
        fig, axes = plt.subplots(
            plot_num, 1,
            figsize=(10, 12),
            sharex=True
        )

        if sharp_region != None:
            ax_int_raw, ax_int_smt, ax_rat_raw, ax_rat_smt, ax_dL, ax_dH = axes
        else:
            ax_int_raw, ax_int_smt, ax_rat_raw, ax_rat_smt = axes
        
        # -----------------------------
        # Z-SCORE PANELS
        # -----------------------------
        for wave in wavelengths:
            df = dfs[wave].sort_values("time")
        
            # Raw
            z_int_raw = df[int_z]
            z_rat_raw = df[rat_z]
            
            ax_int_raw.plot(df["time"], z_int_raw, lw=1.8, label=f"{wave} Å")
            ax_rat_raw.plot(df["time"], z_rat_raw, lw=1.8)

            if shade:
                # 2 sigma
                ax_int_raw.axhspan(-2, 2, color="grey", alpha=0.03)
                ax_rat_raw.axhspan(-2, 2, color="grey", alpha=0.03)
                # 1 sigma
                ax_int_raw.axhspan(-1, 1, color="grey", alpha=0.03)
                ax_rat_raw.axhspan(-1, 1, color="grey", alpha=0.03)
                # limit axes
                ax_int_raw.set_ylim(-3,3)
                ax_rat_raw.set_ylim(-3,3)
            
            # Smoothed
            z_int_smt = z_int_raw.rolling(smooth_z, center=False, min_periods=1).mean()
            z_rat_smt = z_rat_raw.rolling(smooth_z, center=False, min_periods=1).mean()
        
            ax_int_smt.plot(df["time"], z_int_smt, lw=1.8, label=f"{wave} Å")
            ax_rat_smt.plot(df["time"], z_rat_smt, lw=1.8)
            
            if stat_name == 'Kurtosis' and kurt_range:
                ax_int_raw.set_ylim(kurt_range[0],kurt_range[1])
                ax_int_smt.set_ylim(kurt_range[0],kurt_range[1]/3.)
            
            if stat_name == 'Skew' and skew_range:
                ax_int_raw.set_ylim(skew_range[0],skew_range[1])
                ax_int_smt.set_ylim(skew_range[0],skew_range[1])
            
            if stat_name == 'Entropy' and entr_range:
                ax_int_raw.set_ylim(entr_range[0],entr_range[1])
                ax_int_smt.set_ylim(entr_range[0],entr_range[1])
            
            if shade:
                # 2 sigma
                ax_int_smt.axhspan(-2, 2, color="grey", alpha=0.03)
                ax_rat_smt.axhspan(-2, 2, color="grey", alpha=0.03)
                #1 sigma
                ax_int_smt.axhspan(-1, 1, color="grey", alpha=0.03)
                ax_rat_smt.axhspan(-1, 1, color="grey", alpha=0.03)
                # limit axes
                ax_int_smt.set_ylim(-3,3)
                ax_rat_smt.set_ylim(-3,3)
                
        # Zero lines
        for ax in [ax_int_raw, ax_int_smt, ax_rat_raw, ax_rat_smt]:
            ax.axhline(0, color="k", lw=0.8, alpha=0.3)
        
        # -----------------------------
        # dL / dH PANELS
        # -----------------------------
        if sharp_region != None:
            for ax, y, label, color in [
                (ax_dL, dL, "dL", "tab:blue"),
                (ax_dH, dH, "dH", "tab:purple"),
            ]:
                y = pd.Series(y, index=time)
                mu  = y.rolling(dl_dh_window, center=False, min_periods=3).mean()
                sig = y.rolling(dl_dh_window, center=False, min_periods=3).std()
            
                # Raw
                ax.plot(y.index, y.values, lw=1.0, color=color, label=f"{label} (raw)")
                # Mean
                ax.plot(mu.index, mu.values, lw=2.0, color=color, label=f"{label} ⟨{dl_dh_window}⟩")
                # ±2σ envelope
                ax.fill_between(
                    mu.index, 
                    mu - 2*sig, 
                    mu + 2*sig, 
                    color=color, 
                    alpha=0.15, 
                    label=f"{label} ±2σ"
                )
                ax.fill_between(
                    mu.index, 
                    mu - 1*sig, 
                    mu + 1*sig, 
                    color=color, 
                    alpha=0.15, 
                    label=f"{label} ±1σ"
                )
        
                # Y-limits from visible x-range
                xmask = (y.index >= xmin) & (y.index <= xmax)
                visible = np.concatenate([
                y.loc[xmask].to_numpy(),
                (mu - 2 * sig).loc[xmask].to_numpy(),
                (mu + 2 * sig).loc[xmask].to_numpy()
                ])
            
                visible = visible[np.isfinite(visible)]
                
                if visible.size > 0:
                    ymin = visible.min()
                    ymax = visible.max()
                    ax.set_ylim(-10, ymax * 1.05)
        
        # -----------------------------
        # Event markers (all panels)
        # -----------------------------
        # Major ticks: every hour (labeled)
        major_locator = mdates.HourLocator(interval=2)
        major_formatter = mdates.DateFormatter('%m-%d %H')
        
        # Minor ticks: every 30 minutes (no labels)
        minor_locator = mdates.HourLocator(byhour=[1])
        for ax in axes:
            for t in t_dt:
                ax.axvline(t, color="black", lw=2, alpha=0.625, zorder=0, ls='--')
            for t in X_times:
                ax.axvline(t, color="red", lw=2, alpha=0.35, zorder=0, ls='--')
            for t in M_times:
                ax.axvline(t, color="green", lw=2, alpha=0.35, zorder=0, ls='--')
            for t in C_times:
                ax.axvline(t, color="blue", lw=2, alpha=0.25, zorder=0, ls='--')
            
            ax.set_xlim(xmin, xmax)
            ax.xaxis.set_major_locator(major_locator)
            ax.xaxis.set_major_formatter(major_formatter)
            ax.xaxis.set_minor_locator(minor_locator)
            ax.tick_params(axis='x', which='minor', length=4)
            #ax.title.set_fontsize(16)
            #ax.xaxis.label.set_size(15)
            #ax.yaxis.label.set_size(15)
            
        # -----------------------------
        # Titles & labels
        # -----------------------------
        ax_int_raw.set_title(f"{slabel1} (raw)")
        ax_int_smt.set_title(f"{slabel1} (smoothed, {smooth_z} pts)")
        ax_rat_raw.set_title(f"{slabel2} (raw)")
        ax_rat_smt.set_title(f"{slabel2} (smoothed, {smooth_z} pts)")
        if sharp_region != None:
            ax_dL.set_title("Magnetic Winding (running mean ±2σ)")
            ax_dH.set_title("Magnetic Helicity (running mean ±2σ)")
        
        ax_int_raw.set_ylabel(f"{stat_name}")
        ax_int_smt.set_ylabel(f"{stat_name}")
        ax_rat_raw.set_ylabel(f"{stat_name}")
        ax_rat_smt.set_ylabel(f"{stat_name}")
        if sharp_region != None:
            ax_dL.set_ylabel(r'$\delta \mathcal{L}_c^{\prime}\quad\left[\text{km}^2\,\text{s}^{-1} \right]$')
            ax_dH.set_ylabel(r'$\delta \mathcal{H}_c^{\prime}\quad\left[\text{G}^2\text{km}^2\,\text{s}^{-1} \right]$')
            ax_dH.set_xlabel("UTC Time [hours]")
        else:
            ax_rat_smt.set_xlabel("UTC Time [hours]")
        
        # -----------------------------
        # Legend (top, horizontal)
        # -----------------------------
        handles, labels = ax_int_raw.get_legend_handles_labels()
        fig.legend(
            handles, labels,
            loc="upper center",
            ncol=len(wavelengths),
            frameon=False,
            bbox_to_anchor=(0.5, 0.995)
        )
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        # -----------------------------
        # Save
        # -----------------------------
        if sharp_region != None:
            fname = outdir / f"{fname}_dLdH_{start_str}_to_{end_str}.pdf"
        else:
            fname = outdir / f"{fname}_{start_str}_to_{end_str}.pdf"
        fig.savefig(fname, format="pdf")
        print(f"Saved: {fname}")
        plt.show()


    
    # -----------------------------
    # Plot correlations for all stats
    # -----------------------------
    if plot_correlations:

        def generate_full_index(dfs, cadence=6):
            """Create a global timestamp index covering all wavelengths at fixed cadence."""
            all_times = pd.concat([df.index.to_series() for df in dfs.values()])
            start = all_times.min()
            end = all_times.max()
            return pd.date_range(start=start, end=end, freq=f'{cadence}min')
        
        
        
        def stack_quantity(dfs, wavelengths, quantity, cadence=6):
            """
            Stack a quantity across wavelengths, aligning timestamps.
            Missing timestamps are added as NaNs.
            """
            # Generate full index
            full_index = generate_full_index(dfs, cadence)
        
            # Reindex each wavelength
            aligned = {wl: dfs[wl].reindex(full_index) for wl in wavelengths}
        
            # Concatenate
            stacked = pd.concat({wl: aligned[wl][quantity] for wl in wavelengths},
                                axis=1).sort_index()
            return stacked
        
        
        
        def rolling_inter_wavelength_correlation(df, window, min_waves=2):
            out = np.full(len(df), np.nan)
        
            for i in range(window - 1, len(df)):
                block = df.iloc[i - window + 1 : i + 1]
        
                # Drop columns that are completely NaN in window
                block = block.dropna(axis=1, how='all')
        
                if block.shape[1] < min_waves:
                    continue
        
                # Compute correlation using pairwise complete observations
                C = block.corr().values
        
                iu = np.triu_indices_from(C, k=1)
                out[i] = np.nanmean(C[iu])
        
            return pd.Series(out, index=df.index, name="inter_wavelength_corr")
        
        
        
        def prepare_df(df, cadence=6):
        
            # Case 1: already datetime index
            if isinstance(df.index, pd.DatetimeIndex):
                pass
        
            # Case 2: has time column
            elif 'time' in df.columns:
                df['time'] = pd.to_datetime(df['time'], utc=True)
                df = df.set_index('time')
        
            else:
                raise ValueError("No datetime index or 'time' column found.")
        
            df = df.sort_index()
        
            # Snap to cadence grid
            df.index = df.index.floor(f'{cadence}min')
        
            # Remove duplicates after flooring
            df = df[~df.index.duplicated(keep='first')]
        
            return df


        for wl in wavelengths:
            dfs[wl] = prepare_df(dfs[wl])
        
        correlations = {}
        for stat_name, cfg in STAT_CONFIGS.items():
            int_z = cfg["int_z"]
            rat_z = cfg["rat_z"]
        
            # intensity quantity correlations
            df = stack_quantity(dfs,wavelengths,int_z)
            correlations[int_z] = rolling_inter_wavelength_correlation(
                df,
                window=10
            )
        
            # ratio quantity correlations
            df = stack_quantity(dfs,wavelengths,rat_z)
            correlations[rat_z] = rolling_inter_wavelength_correlation(
                df,
                window=10
            )
           
        # separate raw and z-score stats
        keys_z,keys_nz = [],[]
        for k in correlations.keys():
            if k.endswith('z'):
                keys_z.append(k)
            else:
                keys_nz.append(k)
        zcorr = {k: correlations[k] for k in keys_z}
        nzcorr = {k: correlations[k] for k in keys_nz}
        ctime = df.index
        
        
        # separate out the variables and do the plotting
        for it in range(2):
            if it == 0: 
                corr = zcorr
                klab = 'Kurtosis Z-Score'
                slab = 'Skew Z-Score'
                elab = 'Entropy Z-Score'
                fname = 'correlations_zscore'
            else:
                corr = nzcorr
                klab = 'Kurtosis'
                slab = 'Skew'
                elab = 'Entropy'
                fname = 'correlations'
        
            # intensity correlations
            kurt = pd.DataFrame({k: v for k, v in corr.items() if "kurtosis_int" in k})
            skew = pd.DataFrame({k: v for k, v in corr.items() if "skew_int" in k})
            entr = pd.DataFrame({k: v for k, v in corr.items() if "entropy_int" in k})
            # ratio correlations
            kurtR = pd.DataFrame({k: v for k, v in corr.items() if "kurtosis_rat" in k})
            skewR = pd.DataFrame({k: v for k, v in corr.items() if "skew_rat" in k})
            entrR = pd.DataFrame({k: v for k, v in corr.items() if "entropy_rat" in k})
            # smoothed correlations
            kurt_smth = pd.DataFrame({
                k: v.rolling(smooth_z, center=True, min_periods=1).mean()
                for k, v in kurt.items()
            })
            skew_smth = pd.DataFrame({
                k: v.rolling(smooth_z, center=True, min_periods=1).mean()
                for k, v in skew.items()
            })
            entr_smth = pd.DataFrame({
                k: v.rolling(smooth_z, center=True, min_periods=1).mean()
                for k, v in entr.items()
            })
            kurtR_smth = pd.DataFrame({
                k: v.rolling(smooth_z, center=True, min_periods=1).mean()
                for k, v in kurtR.items()
            })
            skewR_smth = pd.DataFrame({
                k: v.rolling(smooth_z, center=True, min_periods=1).mean()
                for k, v in skewR.items()
            })
            entrR_smth = pd.DataFrame({
                k: v.rolling(smooth_z, center=True, min_periods=1).mean()
                for k, v in entrR.items()
            })
            
            # -----------------------------
            # FIGURE LAYOUT
            # -----------------------------
            if sharp_region != None:
                plot_num = 5
            else:
                plot_num = 3
            fig, axes = plt.subplots(
                plot_num, 1,
                figsize=(10, 12),
                sharex=True
            )
            if sharp_region != None:
                ax_kur, ax_skw, ax_ent, ax_dL, ax_dH = axes
            else:
                ax_kur, ax_skw, ax_ent = axes
            
            # -----------------------------
            # CORRELATION PANELS
            # -----------------------------
            # kurtosis correlations
            ax_kur.plot(ctime, kurt, lw=1.8, c='tab:blue', label=f"Intensity Correlation")
            ax_kur.plot(ctime, kurt_smth, lw=1.8, c='tab:blue', linestyle='-.', 
                        alpha=0.8, label=f"Smoothed Intensity Correlation")
            ax_kur.plot(ctime, kurtR, lw=1.8, c='tab:orange', label=f"Ratio Correlation")
            ax_kur.plot(ctime, kurtR_smth, lw=1.8, c='tab:orange', linestyle='-.', 
                        alpha=0.8, label=f"Smoothed Ratio Correlation")
            ax_kur.legend()
            # skew correlations
            ax_skw.plot(ctime, skew, lw=1.8, c='tab:blue', label=f"Intensity Correlation")
            ax_skw.plot(ctime, skew_smth, lw=1.8, c='tab:blue', linestyle='-.', 
                        alpha=0.8, label=f"Smoothed Intensity Correlation")
            ax_skw.plot(ctime, skewR, lw=1.8, c='tab:orange', label=f"Ratio Correlation")
            ax_skw.plot(ctime, skewR_smth, lw=1.8, c='tab:orange', linestyle='-.', 
                        alpha=0.8, label=f"Smoothed Ratio Correlation")
            # entropy correlations
            ax_ent.plot(ctime, entr, lw=1.8, c='tab:blue', label=f"Intensity Correlation")
            ax_ent.plot(ctime, entr_smth, lw=1.8, c='tab:blue', linestyle='-.', 
                        alpha=0.8, label=f"Smoothed Intensity Correlation")
            ax_ent.plot(ctime, entrR, lw=1.8, c='tab:orange', label=f"Ratio Correlation")
            ax_ent.plot(ctime, entrR_smth, lw=1.8, c='tab:orange', linestyle='-.', 
                        alpha=0.8, label=f"Smoothed Ratio Correlation")
        
            '''
            if shade:
                # 2 sigma
                ax_int_raw.axhspan(-2, 2, color="grey", alpha=0.03)
                ax_rat_raw.axhspan(-2, 2, color="grey", alpha=0.03)
                # 1 sigma
                ax_int_raw.axhspan(-1, 1, color="grey", alpha=0.03)
                ax_rat_raw.axhspan(-1, 1, color="grey", alpha=0.03)
                # limit axes
                ax_int_raw.set_ylim(-3,3)
                ax_rat_raw.set_ylim(-3,3)
            '''
            
            # Zero lines
            for ax in [ax_kur, ax_skw, ax_ent]:
                ax.axhline(0, color="k", lw=0.8, alpha=0.3)
            
            # -----------------------------
            # dL / dH PANELS
            # -----------------------------
            if sharp_region != None:
                for ax, y, label, color in [
                    (ax_dL, dL, "dL", "tab:blue"),
                    (ax_dH, dH, "dH", "tab:purple"),
                ]:
                    y = pd.Series(y, index=time)
                    mu  = y.rolling(dl_dh_window, center=True, min_periods=3).mean()
                    sig = y.rolling(dl_dh_window, center=True, min_periods=3).std()
                
                    # Raw
                    ax.plot(y.index, y.values, lw=1.0, color=color, label=f"{label} (raw)")
                    # Mean
                    ax.plot(mu.index, mu.values, lw=2.0, color=color, label=f"{label} ⟨{dl_dh_window}⟩")
                    # ±2σ envelope
                    ax.fill_between(
                        mu.index, 
                        mu - 2*sig, 
                        mu + 2*sig, 
                        color=color, 
                        alpha=0.15, 
                        label=f"{label} ±2σ"
                    )
                    ax.fill_between(
                        mu.index, 
                        mu - 1*sig, 
                        mu + 1*sig, 
                        color=color, 
                        alpha=0.15, 
                        label=f"{label} ±1σ"
                    )
                
                    # Y-limits from visible x-range
                    xmask = (y.index >= xmin) & (y.index <= xmax)
                    visible = np.concatenate([
                    y.loc[xmask].to_numpy(),
                    (mu - 2 * sig).loc[xmask].to_numpy(),
                    (mu + 2 * sig).loc[xmask].to_numpy()
                    ])
                
                    visible = visible[np.isfinite(visible)]
                    
                    if visible.size > 0:
                        ymin = visible.min()
                        ymax = visible.max()
                        ax.set_ylim(ymin * 0.95, ymax * 1.05)
            
            # -----------------------------
            # Event markers (all panels)
            # -----------------------------
            # Major ticks: every hour (labeled)
            major_locator = mdates.HourLocator(interval=2)
            major_formatter = mdates.DateFormatter('%m-%d %H')
            
            # Minor ticks: every 30 minutes (no labels)
            minor_locator = mdates.HourLocator(byhour=1)
            for ax in axes:
                for t in t_dt:
                    ax.axvline(t, color="black", lw=2, alpha=0.625, zorder=0, ls='--')
                for t in X_times:
                    ax.axvline(t, color="red", lw=2, alpha=0.35, zorder=0, ls='--')
                for t in M_times:
                    ax.axvline(t, color="green", lw=2, alpha=0.35, zorder=0, ls='--')
                for t in C_times:
                    ax.axvline(t, color="blue", lw=2, alpha=0.25, zorder=0, ls='--')
                
                ax.set_xlim(xmin, xmax)
                ax.xaxis.set_major_locator(major_locator)
                ax.xaxis.set_major_formatter(major_formatter)
                ax.xaxis.set_minor_locator(minor_locator)
                ax.tick_params(axis='x', which='minor', length=4)
                #ax.title.set_fontsize(16)
                #ax.xaxis.label.set_size(15)
                #ax.yaxis.label.set_size(15)
                
            # -----------------------------
            # Titles & labels
            # -----------------------------
            ax_kur.set_title(f"{klab} Correlations")
            ax_skw.set_title(f"{slab} Correlations")
            ax_ent.set_title(f"{elab} Correlations")
            if sharp_region != None:
                ax_dL.set_title("dL (running mean ±2σ)")
                ax_dH.set_title("dH (running mean ±2σ)")
            
            ax_kur.set_ylabel("Correlation")
            ax_skw.set_ylabel("Correlation")
            ax_ent.set_ylabel("Correlation")
            if sharp_region != None:
                ax_dL.set_ylabel("Signal")
                ax_dH.set_ylabel("Signal")
                ax_dH.set_xlabel("Time")
            else:
                ax_ent.set_xlabel("Time")
            
            # -----------------------------
            # Save
            # -----------------------------
            plt.tight_layout(rect=[0, 0, 1, 0.96])
            if sharp_region != None:
                fname = outdir / f"{fname}_dLdH_{start_str}_to_{end_str}.pdf"
            else:
                fname = outdir / f"{fname}_{start_str}_to_{end_str}.pdf"
            fig.savefig(fname, format="pdf")
            print(f"Saved: {fname}")
            plt.show()
    if sharp_region != None:
        if metadata_files:
            return [dfs,correlations,time,dL,dH,deltaL,deltaH]
        else:
            return [dfs,correlations,time,dL,dH,root['dLFlux'],root['dHFlux']]
    else:
        return [dfs,correlations]
