# almanac.py
"""
almanac.py

Python translation of IDL almanac.pro

Dependencies:
    - numpy, scipy, astropy, sunpy, tqdm
    - helper_functions.py (functions: mgn, region_size_filter, pad_3d, unpad_3d, one2n, etc.)
    - get_aia_synoptic_data.py (function: get_aia_synoptic_data)

Primary entrypoint:
    almanac_wave(wave, end_time, ..., test=True, usrdir="...", savdir="...", movies=True, movdir="...")

Notes:
    - This version intentionally processes a single wavelength at a time.
    - Saves results with Python pickle (compressed).
"""

import os
import sys
import copy
from contextlib import redirect_stderr
import glob
import pickle
from pathlib import Path
from typing import Optional, Sequence, List, Tuple, Dict, Any
import imageio.v2 as imageio_v2

import numpy as np
from tqdm import tqdm
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import sunpy.map
from datetime import datetime, timedelta
import scipy
from scipy.ndimage import uniform_filter1d
import gzip

# Import helper modules (assumed to be in the same directory / PYTHONPATH)
from helper_functions import (
    mgn,
    region_size_filter,
    pad_3d,
    unpad_3d,
    one2n,
    gaussian_function,
    normalize_to_uint8,
    update_kurtosis_file
)
from get_aia_synoptic_data import get_aia_synoptic_data
import warnings
from sunpy.util.exceptions import SunpyUserWarning

warnings.filterwarnings(
    "ignore",
    message=".*off-disk coordinates.*",
    category=SunpyUserWarning
)

# Default thresholds copied/adapted from IDL script
R_THRESH = {
    #193: 2.0,
    #211: 2.3,
    193: 2.0,
    211: 2.0,
    304: 2.0,
}
PX_THRESH = {
    #193: int(15e2),
    #211: int(15e2),
    193: int(1e3),
    211: int(1e3),
    304: int(1e3),
}
MIN_VOX_VOL = 6_000  # min voxels for candidate region
MIN_TIME = 18.0  # minutes (as float)
SMO_TIME = 30.0  # smoothing time in minutes




def _ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)




def _find_fits_files(directory: str) -> List[str]:
    """
    Find FITS files (any extension) in directory and return sorted list.
    """
    p = Path(directory)
    if not p.exists():
        return []
    files = sorted([str(x) for x in p.glob("*.fits")])
    return files




def _read_frames_to_cube(filelist):
    """
    Read list of AIA FITS files into a 3D cube (nx, ny, nt)
    and collect their headers/WCS info.

    - Uses SunPy (preferred) or FITS fallback.
    - No rotation or flipping applied.
    - Keeps WCS metadata consistent with native AIA FITS.
    """

    imgs = []
    hdr_list = []

    for f in tqdm(filelist, desc="Reading FITS", unit=" file"):
        try:
            # --- Preferred SunPy read ---
            m = sunpy.map.Map(f)
            data = m.data.astype(np.float32)

            hdr = copy.deepcopy(m.fits_header)
            nax1 = hdr.get("NAXIS1", data.shape[1])
            nax2 = hdr.get("NAXIS2", data.shape[0])
            crpix1 = hdr.get("CRPIX1", nax1 / 2)
            crpix2 = hdr.get("CRPIX2", nax2 / 2)
            cdelt1 = hdr.get("CDELT1", 1.0)
            cdelt2 = hdr.get("CDELT2", 1.0)

            # Create a clean header copy (ensuring numeric consistency)
            hdr_new = hdr.copy()
            hdr_new["NAXIS1"], hdr_new["NAXIS2"] = nax1, nax2
            hdr_new["CRPIX1"], hdr_new["CRPIX2"] = crpix1, crpix2
            hdr_new["CDELT1"], hdr_new["CDELT2"] = cdelt1, cdelt2

            imgs.append(data)
            hdr_list.append(hdr_new)

        except Exception as e:
            # --- Fallback: manual FITS open ---
            try:
                with fits.open(f, memmap=True) as hdul:
                    data = None
                    hdr = None
                    for hdu in hdul:
                        if hdu.data is not None:
                            data = hdu.data.astype(np.float32)
                            hdr = hdu.header
                            break
                    if data is None:
                        print(f"Warning: {f} has no image data, skipping.")
                        continue

                    nax1 = hdr.get("NAXIS1", data.shape[1])
                    nax2 = hdr.get("NAXIS2", data.shape[0])
                    crpix1 = hdr.get("CRPIX1", nax1 / 2)
                    crpix2 = hdr.get("CRPIX2", nax2 / 2)
                    cdelt1 = hdr.get("CDELT1", 1.0)
                    cdelt2 = hdr.get("CDELT2", 1.0)

                    hdr_new = hdr.copy()
                    hdr_new["NAXIS1"], hdr_new["NAXIS2"] = nax1, nax2
                    hdr_new["CRPIX1"], hdr_new["CRPIX2"] = crpix1, crpix2
                    hdr_new["CDELT1"], hdr_new["CDELT2"] = cdelt1, cdelt2

                    imgs.append(data)
                    hdr_list.append(hdr_new)

            except Exception as ee:
                print(f"Error reading {f}: {ee}")
                continue

    if len(imgs) == 0:
        print("No valid image data found.")
        return np.array([]), []

    # --- Stack into (nx, ny, nt) ---
    cube = np.stack(imgs, axis=-1)  # (nx, ny, nt)
    return cube, hdr_list




def _rebin_cube(cube: np.ndarray, new_nx: int, new_ny: int) -> np.ndarray:
    """
    Rebin cube (ny,nx,nt) to (new_ny,new_nx,nt) using local averaging.
    Simple implementation: reshape & mean when integer factors; otherwise use slicing/resizing with interpolation.
    For simplicity, use basic block-mean rebin when possible; else use scipy zoom (if available).
    """
    nx, ny, nt = cube.shape
    if nx == new_nx and ny == new_ny:
        return cube

    # integer factor case
    if nx % new_nx == 0 and ny % new_ny == 0:
        fx = nx // new_nx
        fy = ny // new_ny
        # block mean
        cube_rb = cube.reshape(new_nx, xy, new_ny, fy, nt).mean(axis=(1, 3))
        return cube_rb
        



def _build_ratio_images(im: np.ndarray, range_k: int):
    """
    Compute temporal running-average ratio images for flare detection.

    Parameters
    ----------
    im : np.ndarray
        Input datacube of shape (nx, ny, nt), where nt is the number of time steps.
    range_k : int
        Temporal half-width of the smoothing kernel. The full kernel size is
        nker = 2*range_k + 1, corresponding to ±range_k frames.

    Returns
    -------
    diff : np.ndarray
        Absolute deviation cube |im - avg| of shape (nx, ny, nt),
        where 'avg' is the temporal running mean of each pixel.
    medim : np.ndarray
        2D median absolute deviation across the time dimension:
        median(|im - avg|, axis=2). Used as a baseline threshold image.

    Notes
    -----
    - Uses reflective temporal padding to avoid zero-edge artefacts that can
      cause false detections at the start and end of time series.
    - This implementation is fully vectorized and typically 50–100× faster
      than pixelwise convolution.
    - Ideal for AIA CoM/flare detection preprocessing where temporal smoothing
      is needed per pixel.
    """
    nker = 2 * range_k + 1

    # --- Compute temporal running mean (reflect padding avoids edge artefacts) ---
    avg = uniform_filter1d(im, size=nker, axis=2, mode="reflect")

    # --- Absolute difference from the running mean ---
    diff = np.abs(im - avg)

    # --- Median deviation across time ---
    medim = np.median(diff, axis=2)

    return diff, medim




def _apply_region_masking(ratio: np.ndarray, r_thresh: float, px_thresh: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Given ratio cube shape (nx,ny,nt), produce mask (nx,ny,nt), masksmo smoothed, and labelmask.
    Uses helper region_size_filter and gaussian smoothing similarly to IDL.
    """
    nx, ny, nt = ratio.shape

    # create mask by applying region_size_filter per frame
    mask = np.zeros_like(ratio, dtype=bool)
    for it in range(nt):
        # threshold to boolean via region_size_filter
        mask[:, :, it] = region_size_filter(ratio[:, :, it], r_thresh, int(0.7 * px_thresh), return_mask=True)

    # smoothing in x then y
    kermask = gaussian_function(1.0, normalize=True)
    # convolution across axes (separable)
    from scipy.ndimage import convolve1d

    masksmo = mask.astype(float)
    masksmo = convolve1d(masksmo, kermask, axis=0, mode="mirror")
    masksmo = convolve1d(masksmo, kermask, axis=1, mode="mirror")

    # time smoothing narrow kernel
    kert = gaussian_function(0.7, normalize=True)
    # pad kert to length (1,1,nk)
    masksmo = convolve1d(masksmo, kert, axis=2, mode="mirror")

    median_val = np.median(masksmo)
    masksmo_bool = masksmo > median_val

    # apply region_size_filter on masksmo per frame again (with more strict npix)
    mask2 = np.zeros_like(mask)
    for it in range(nt):
        mask2[:, :, it] = region_size_filter(masksmo_bool[:, :, it].astype(float), 0.8, px_thresh, return_mask=True)

    # labelmask = unpad3d(label_region(pad_3d(mask2,1), /ulong),1)*mask   in IDL
    # we will label the padded mask, then unpad and multiply by original mask
    pad = pad_3d(mask2.astype(np.uint8), 1)
    from scipy.ndimage import label as nd_label

    lbl_padded, n_lbl = nd_label(pad)
    lbl_unpad = unpad_3d(lbl_padded, 1)
    labelmask = lbl_unpad * mask

    return mask, masksmo_bool, labelmask




def _histogram_reverse_indices(labelmask: np.ndarray, minvoxvol: int) -> Tuple[np.ndarray, np.ndarray, List[int]]:
    """
    Compute histogram of labelmask (values >0), return hist (counts), unique labels above threshold,
    and a reverse-indices-like structure: returns flattened indices for each label as dictionary.
    This function is a conceptual replacement for IDL histogram(...,reverse_ind=revind).
    """
    lab = labelmask.astype(np.int64)
    lab_flat = lab.ravel()
    # bin counts for positive labels only
    labels, counts = np.unique(lab_flat[lab_flat > 0], return_counts=True)
    # pick labels with count > minvoxvol
    good_idx = np.where(counts > minvoxvol)[0]
    good_labels = labels[good_idx]
    good_counts = counts[good_idx]
    # Build reverse index map: dict label->indices list
    rev_map = {}
    for labv in good_labels:
        inds = np.where(lab_flat == labv)[0]
        rev_map[int(labv)] = inds
    return labels, counts, rev_map




def _save_event_pickle(outdir: str, timestamp_str: str, wave: int, iregion: int, payload: dict):
    """
    Save payload dict to a compressed .pkl.gz file under:
        outdir/{year}/{month}/almanac_{timestamp}_wvlnth_{wave}_{region}.pkl.gz

    Compression uses gzip (fast and widely supported).
    """
    year = timestamp_str[0:4]
    month = timestamp_str[4:6]
    dirsav = os.path.join(outdir, year, month)
    _ensure_dir(dirsav)

    fname = os.path.join(dirsav, f"almanac_{timestamp_str}_wvlnth_{wave}_{iregion}.pkl.gz")
    with gzip.open(fname, "wb") as fh:
        pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return fname




def almanac_wave(
    wave: int,
    end_time: str,
    datadir_root: str = ".",
    test: bool = False,
    usrdir: Optional[str] = None,
    flist: Optional[List[str]] = None,
    savdir: str = "./almanac_out",
    show: bool = False,
    cadence: int = 360,
    movies: bool = False,
    movdir: Optional[str] = None,
    arg: Optional[int] = None,
    nrt: bool = False,
):
    
    """
    Process a single wavelength similar to IDL almanac, but in Python.

    Parameters
    ----------
    wave : int
        Wavelength (e.g. 94, 131, 171, 193, ...)
    end_time : str
        ISO-like end time string 'YYYY/MM/DD HH:MM' (used for selecting directory/date)
    datadir_root : str
        Base data directory (if not using get_aia_synoptic_data to fetch)
    test : bool
        If True, expect files already present locally and do not attempt downloads.
    usrdir : str
        If test=True, directory of pre-downloaded files
    flist : list[str]
        Optional explicit list of files (full paths)
    savdir : str
        Directory to save outputs (pickles). Defaults to cwd.
    show : bool
        If True, print progress / show interactive (not implemented in headless).
    cadence : int
        Desired cadence (seconds) to use when downloading.
    movies : bool
        If True, create movies via moviemaker
    movdir : str
        Directory to write movies (if different than savdir).
    arg : int
        Node/processor number equivalent (keeps filename compatibility)
    nrt : bool
        Use nrt data if True (delegated to get_aia_synoptic_data)
    """
    # --- setup defaults & thresholds ---
    savdir = savdir
    movdir = movdir or savdir
    datadir_root = datadir_root or os.getcwd()
    step = cadence
    flag_synoptic = not nrt

    if wave in R_THRESH:
        r_thresh = R_THRESH[wave]
        px_thresh = PX_THRESH[wave]
    else:
        r_thresh = R_THRESH[193]
        px_thresh = PX_THRESH[193]

    # If not test, attempt download using get_aia_synoptic_data
    if not test:
        # We need time window: derive time0/time1 around end_time (IDL used same time for time0/time1)
        # We convert end_time to two strings time0 and time1 using the timestamp helper
        # For simplicity, we set time window to that day (00:00 - 23:59) around end_time
        try:
            et = end_time.strip()
            if "/" in et:
                et_parse = datetime.strptime(et, "%Y/%m/%d %H:%M")
            else:
                et_parse = datetime.strptime(et, "%Y-%m-%d %H:%M")
        except Exception:
            # try ISO fallback
            et_parse = datetime.fromisoformat(end_time)

        day0 = et_parse.replace(hour=0, minute=0, second=0)
        day1 = et_parse.replace(hour=23, minute=59, second=0)

        time0 = day0.strftime("%Y/%m/%d %H:%M")
        time1 = day1.strftime("%Y/%m/%d %H:%M")

        # Call your provided downloader to ensure files exist locally
        attempted_files, status = get_aia_synoptic_data(time0, time1, wlusr=[wave], increment=step, topsavedir=datadir_root, verbose=False, nrt=nrt)
        # The downloader saves under datadir_root/<wl3>/* filenames
        # set datadir for reading
        wl3 = f"{wave:03d}"
        datadir = os.path.join(datadir_root, wl3)
    else:
        # test mode: read from usrdir or flist
        if flist:
            read_files = list(flist)
            datadir = os.path.dirname(read_files[0]) if len(read_files) > 0 else usrdir or datadir_root
        else:
            datadir = usrdir if usrdir else os.path.join(datadir_root, f"{wave:03d}")

    # --- gather files ---
    if flist:
        read_files = list(flist)
    else:
        read_files = _find_fits_files(datadir)

    if len(read_files) == 0:
        print(f"{wave}: No FITS files found in {datadir}. Exiting.")
        return [],[],[],[]

    # Reconstruct filenames like IDL expects and read into datacube
    # read all frames into cube
    im_cube, hdrs = _read_frames_to_cube(read_files)

    if im_cube.size == 0:
        print(f"{wave}: No image data read. Exiting.")
        return [],[],[], hdrs

    # --- IDL: im=float(im>0); set negatives to zero ---
    im_cube = (im_cube > 0).astype(float) * im_cube  # set negatives to 0
    nx, ny, nt = im_cube.shape

    # convert to exposure-normalized (IDL loop dividing by exptime in headers)
    for i in range(nt):
        try:
            exptime = hdrs[i].get("EXPTIME", hdrs[i].get("exptime", 1.0))
        except Exception:
            exptime = 1.0
        if exptime is None or exptime == 0:
            exptime = 1.0
        im_cube[:, :, i] = im_cube[:, :, i] / exptime

    # Rebin to 1024x1024 if larger
    if nx > 1024 or ny > 1024:
        im_cube = _rebin_cube(im_cube, 1024, 1024)
        nx, ny, nt = im_cube.shape
        # adjust headers crudely (not perfect)
        for it in range(nt):
            hdrs[it]["NAXIS1"] = nx
            hdrs[it]["NAXIS2"] = ny
            # scale CRPIX/CDelt if present
            if "CRPIX1" in hdrs[it]:
                hdrs[it]["CRPIX1"] = hdrs[it]["CRPIX1"] * (nx / hdrs[it].get("NAXIS1", nx))
            if "CRPIX2" in hdrs[it]:
                hdrs[it]["CRPIX2"] = hdrs[it]["CRPIX2"] * (ny / hdrs[it].get("NAXIS2", ny))

    # --- remove frames that are all zeros or contain NaNs ---
    valid_inds = []
    for i in range(nt):
        frame = im_cube[:, :, i]
        if np.isnan(np.median(frame)):
            continue        # skip frames fully NaN
        if np.median(frame) == 0:
            continue        # skip frames fully zero
        valid_inds.append(i)
    
    if len(valid_inds) == 0:
        print(f"{wave}: No valid frames after removing zeros/NaNs. Exiting.")
        return [],im_cube,[],hdrs
    
    # keep only valid frames
    im_cube = im_cube[:, :, valid_inds]
    hdrs = [hdrs[i] for i in valid_inds]

    nx, ny, nt = im_cube.shape

    # require minimum frames
    # compute cadence from first few filenames (IDl did this using substrings)
    # We will estimate cadence in minutes by timestamps in headers if available
    def _frame_time_from_hdr(hdr):
        # attempt to read DATE-OBS or T_OBS etc.
        for k in ("DATE-OBS", "DATE_OBS", "DATE", "T_OBS", "T_OBS"):
            if k in hdr:
                try:
                    return datetime.fromisoformat(hdr[k])
                except Exception:
                    try:
                        return datetime.strptime(hdr[k], "%Y-%m-%dT%H:%M:%S.%f")
                    except Exception:
                        pass
        # fallback: None
        return None

    # compute estimated cadence td in minutes from header times if available
    times = [_frame_time_from_hdr(h) for h in hdrs[:4]]
    valid_times = [t for t in times if t is not None]
    if len(valid_times) >= 2:
        td_seconds = (valid_times[1] - valid_times[0]).total_seconds()
        td = td_seconds / 60.0
    else:
        td = cadence / 60.0  # use provided cadence
    min_frames = int(round(MIN_TIME / td)) if td > 0 else int(round(MIN_TIME))

    if min_frames < 1:
        print(f"{wave}: Insufficient cadence; adjust cadence. Exiting.")
        return [], im_cube, [], hdrs

    range_k = int(round(SMO_TIME / td))

    # Downsample / crop to central region (IDL trimmed edges)
    ndx = int(np.round((nx - nx * 0.65) / 2.0))
    ndy = int(np.round((nx - nx * 0.85) / 2.0))
    # ensure valid
    x0 = max(0, ndx)
    x1 = nx - x0
    y0 = max(0, ndy)
    y1 = nx - ndy if (nx - ndy) > 0 else ny
    im_cube = im_cube[x0:x1, y0:y1, :]
    nx, ny, nt = im_cube.shape
    
    # adjust headers crpix
    for it in range(nt):
        if "CRPIX1" in hdrs[it]:
            hdrs[it]["CRPIX1"] -= x0
        if "CRPIX2" in hdrs[it]:
            hdrs[it]["CRPIX2"] -= y0

    # standardize intensity (IDL: median scaling to 150)
    im_cube = np.clip(im_cube, 0, 2e3)
    val = np.median(im_cube)
    if val == 0:
        val = 1.0
    im_cube = im_cube / val
    im_cube = im_cube * 150.0

    # compute temporal average/diff and ratio
    diff, medim = _build_ratio_images(im_cube, range_k)

    # ratio = diff / rebin(medim, nx, ny, nt)  (IDL used rebin to duplicate medim across time)
    medim_rep = np.repeat(medim[:, :, np.newaxis], nt, axis=2)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = diff / medim_rep
    ratio = np.clip(ratio, 0, 30)

    # region detection & masking
    mask, masksmo_bool, labelmask = _apply_region_masking(ratio, r_thresh, px_thresh)

    # filter out tiny labelled volumes
    # compute histogram-like reverse indices
    labels, counts, rev_map = _histogram_reverse_indices(labelmask, MIN_VOX_VOL)
    if len(rev_map) == 0:
        print(f"{wave}: No candidate regions above minvoxvol threshold.")
        return [], im_cube, ratio, hdrs
      
    # Loop through candidate regions
    results = []
    # preserve original nt before region-specific cropping
    orig_nt = nt

    for idx_region, labval in enumerate(sorted(rev_map.keys())):
        inds_flat = rev_map[labval]
        # Convert flattened indices to x,y,t using one2n (helper has similar behavior)
        ix_arr, iy_arr, it_arr = one2n(inds_flat, dimensions=[nx, ny, nt])
        # one2n returns tuple arrays per dimension; ensure order matches IDL (x,y,t)
        if isinstance(ix_arr, (list, tuple)) and len(ix_arr) == 3:
            iy = ix_arr[0]
            ix = ix_arr[1]
            it = ix_arr[2]
        else:
            # fallback if one2n returns arrays directly
            ix, iy, it = ix_arr, iy_arr, it_arr

        trange = int(it.max() - it.min() + 1)
        if trange < min_frames:
            # skip too-short events
            continue

        # compute centroid weighted by ratio within this region per frame
        unique_times = np.unique(it)
        hgs_lon = []
        hgs_lat = []
        hgc_lon = []
        hgc_lat = []
        coords_list = []
        xavg_list = []
        yavg_list = []

        for tnow in unique_times:
            mask_time = (it == tnow)
            ixnow = ix[mask_time]  # horizontal (x)
            iynow = iy[mask_time]  # vertical (y)
            flat_inds_time = inds_flat[mask_time]
        
            # weight by ratio values
            weights = ratio.reshape(-1)[flat_inds_time]
            if weights.sum() == 0:
                weights = np.ones_like(weights)
            ixav = np.sum(ixnow * weights) / np.sum(weights)
            iyav = np.sum(iynow * weights) / np.sum(weights)
            xavg_list.append(ixav)
            yavg_list.append(iyav)
        
            # WCS conversion
            try:
                hdr_now = hdrs[int(tnow)]
                mp = sunpy.map.Map(im_cube[:, :, int(tnow)], hdr_now)
                coord = mp.pixel_to_world(ixav * u.pixel, iyav * u.pixel)
                hgs = coord.transform_to(sunpy.coordinates.HeliographicStonyhurst)
                hgc = coord.transform_to(sunpy.coordinates.HeliographicCarrington)
                # Due to a weird transpose issue with x and y, lon and lat are mixed-up
                # here to correct the issue for the coordinates.
                lat_hgs, lon_hgs = hgs.lon.deg, hgs.lat.deg
                lat_hgc, lon_hgc = hgc.lon.deg, hgc.lat.deg
            except Exception:
                lon_hgs = lat_hgs = lon_hgc = lat_hgc = None
                coord = None
        
            hgs_lon.append(lon_hgs)
            hgs_lat.append(lat_hgs)
            hgc_lon.append(lon_hgc)
            hgc_lat.append(lat_hgc)
            coords_list.append(coord)

        # Build output arrays like IDL's 'com' structure
        # com shape (trange, 4, 2) in IDL; we'll use dict
        com = {
            "stonyhurst_lon": hgs_lon,
            "stonyhurst_lat": hgs_lat,
            "carrington_lon": hgc_lon,
            "carrington_lat": hgc_lat,
            "coords": coords_list,
            "pixel_avg": list(zip(xavg_list, yavg_list)),
        }

        # Extract masks and images for event window
        tmin = int(it.min())
        tmax = int(it.max())
        # clamp within existing cube
        tmin = max(0, tmin)
        tmax = min(orig_nt - 1, tmax)
        cmemask = (labelmask == labval).astype(np.uint8)
        cmemask_segment = cmemask[:, :, tmin:tmax + 1]
        cmeim_segment = im_cube[:, :, tmin:tmax + 1]
        cmerat_segment = ratio[:, :, tmin:tmax + 1]

        # convert to 8bit integers between 0-255
        cmeim_segment = normalize_to_uint8(cmeim_segment)
        cmemask_segment = normalize_to_uint8(cmemask_segment)
        # invert so that black is a detection, white is background
        cmemask_segment = 255 - cmemask_segment
        
        # Prepare saved payload
        timestamp_hdr = hdrs[tmin].get("DATE-OBS", hdrs[tmin].get("DATE", "unknown"))
        # generate safe timestamp string YYYYMMDD_HHMMSS (fallback if not present)
        try:
            dt_obj = datetime.fromisoformat(timestamp_hdr)
            ts_str = dt_obj.strftime("%Y%m%d_%H%M")
        except Exception:
            ts_str = datetime.utcnow().strftime("%Y%m%d_%H%M")

        # save images as sunpy maps
        aiaim_maps = []
        maskim_maps = []
        ratio_maps = []

        for i in range(cmeim_segment.shape[2]):
            aiaim_maps.append(sunpy.map.Map(cmeim_segment[:,:,i], hdrs[tmin+i]))
            maskim_maps.append(sunpy.map.Map(cmemask_segment[:,:,i], hdrs[tmin+i]))
            ratio_maps.append(sunpy.map.Map(cmerat_segment[:,:,i], hdrs[tmin+i]))
        
        info = {
            "nx": nx,
            "ny": ny,
            "nt": tmax - tmin + 1,
            "hdr": hdrs[tmin:tmax + 1],
            "index": inds_flat,
            "values": None,  # we could place intensity values if desired
            "ratio": ratio_maps,
            "xtrim": ndx,
            "ytrim": ndy,
            "wave": str(wave),
            "CoM": com,
            "totals": None,
            "aiaim": aiaim_maps,
            "maskim": maskim_maps,
            "comboim": np.hstack((cmeim_segment, cmemask_segment)),
            "im_cube": im_cube
        }

        # Save event to pickle
        save_fname = _save_event_pickle(savdir, ts_str, wave, idx_region, info)

        # Optionally make movie for this region
        if movies:
            import matplotlib.pyplot as plt
            # use moviemaker to create mp4 from comboim (or use AIA images)
            # moviemaker expects a directory of frames or images; we will create a temporary directory
            tmp_movie_dir = os.path.join(movdir, f"almanac_{ts_str}_wvlnth_{wave}_{idx_region}")
            os.makedirs(tmp_movie_dir, exist_ok=True)

            # create temporary FITS files for frames or PNGs; but moviemaker accepts fits_dir in our earlier version
            # Here: write frames as PNGs using helper mgn if requested.
            frames_files = []
            for i_frame in range(info["nt"]):
                # Flip/prepare frame image
                frame_img = info["comboim"][:, :, i_frame]  # AIA image
                processed = mgn(frame_img, h=0.925)
                dmin, dmax = 0.15, 0.9
                proc_clipped = np.clip(processed, dmin, dmax)
                proc_byte = ((proc_clipped - dmin) * 255.0 / (dmax - dmin)).astype(np.uint8)
                
                # Get observation time from header
                hdr_now = hdrs[tmin + i_frame]
                obs_time = hdr_now.get("DATE-OBS", hdr_now.get("DATE", "unknown"))
                # Optional: format nicely
                try:
                    obs_dt = datetime.fromisoformat(obs_time)
                    obs_str = obs_dt.strftime("%Y-%m-%d %H:%M")
                except Exception:
                    obs_str = str(obs_time)
               
                # Plot with matplotlib and overlay text
                fig, ax = plt.subplots(figsize=(16, 8))
                ax.imshow(proc_byte, cmap="gray", origin="lower")
                ax.axis("off")
                ax.text(5, 5, obs_str, color="yellow", fontsize=12,
                        bbox=dict(facecolor='black', alpha=0.5, pad=2))
               
                # save PNG
                out_png = os.path.join(tmp_movie_dir, f"frame_{i_frame:04d}.png")
                fig.savefig(out_png, bbox_inches="tight", pad_inches=0)
                plt.close(fig)
                frames_files.append(out_png)

            # call moviemaker with frames directory: our moviemaker function earlier expects fits_dir, but we changed it to accept directory of fits.
            # Here we'll call moviemaker with fits_dir pointing to the PNG directory and it will attempt to read *.fits; simplest to call ffmpeg separately or adapt moviemaker.
            # For simplicity, call moviemaker using images=combo mp4 via imageio directly:im array (not implemented in that moviemaker). Instead produce a simple
            movie_path = os.path.join(movdir, f"almanac_{ts_str}_wvlnth_{wave}_{idx_region}.mp4")
            
            # suppress ffmpeg warnings
            with open(os.devnull, "w") as fnull:
                with redirect_stderr(fnull):
                    with imageio_v2.get_writer(movie_path, fps=10) as writer:
                        for fpng in frames_files:
                            img = imageio_v2.imread(fpng)
                            writer.append_data(img)

        results.append({"region": idx_region, "info": info, "pickle": save_fname})

    return results, im_cube, ratio, hdrs




def almanac_instance(args):
    end_time, fits_dir, wave, wvl_files = args
    """Run a single almanac_wave() call and return timing."""
    files = [os.path.join(f"{fits_dir}", f"{wave:03d}", fname) for fname in wvl_files[wave]]
    out = almanac_wave(wave,
                       end_time,
                       test=True,
                       movies=False,
                       flist=files,
                       movdir='./movies'
                      )
    results, im_cube, ratio, hdrs = out
    # before finding potential CME events, determine the kurtosis of im_cube and
    # ratio
    update_kurtosis_file(
        wave=wave,
        im_cube=im_cube,
        ratio_cube=ratio,
        hdrs=hdrs,
        outdir="./kurtosis"
    )
    
    return
