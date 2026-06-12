import re
import pickle
import gzip
import os
import shutil
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
from dateutil import parser
from collections import defaultdict
from scipy import ndimage
from scipy.ndimage import label, generate_binary_structure
from typing import List, Tuple, Dict, Any, Callable
import warnings
from sunpy.util.exceptions import SunpyUserWarning

warnings.filterwarnings(
    "ignore",
    message=".*off-disk coordinates.*",
    category=SunpyUserWarning
)




def build_binned_com_cube(file_list: List[str],
                          bin_size: int = 50,
                          cadence_minutes: int = 6
                          ) -> Tuple[np.ndarray, List[datetime], Dict[Tuple[int,int,int], List[Tuple[str,int]]]]:

    cadence_seconds = cadence_minutes * 60
    all_times = []

    first_valid_file = None
    first_data = None

    # --------------------------------------------
    # Single-pass: load each file once
    # --------------------------------------------
    file_data = []
    for f in file_list:
        try:
            if f.endswith(".gz"):
                with gzip.open(f, "rb") as fh:
                    data = pickle.load(fh)
            else:
                with open(f, "rb") as fh:
                    data = pickle.load(fh)
        except Exception as e:
            print(f"WARNING: failed to load '{f}': {e}")
            continue

        if first_valid_file is None:
            first_valid_file = f
            first_data = data

        nt = data.get("nt", 0)
        hdrs = data.get("hdr", [])
        pixel_avg = data.get("CoM", {}).get("pixel_avg", [])

        # collect valid times for global time axis
        times = []
        for it in range(nt):
            try:
                t = datetime.strptime(hdrs[it]["T_OBS"][:16], "%Y-%m-%dT%H:%M")
                times.append(t)
            except Exception:
                continue
        all_times.extend(times)

        file_data.append((f, nt, hdrs, pixel_avg, times))

    if not all_times:
        print("WARNING: No valid observation times found. Returning empty outputs.")
        return [], [], {}

    # --------------------------------------------
    # Build global time axis
    # --------------------------------------------
    t_min, t_max = min(all_times), max(all_times)
    nt_total = int((t_max - t_min).total_seconds() / cadence_seconds) + 1
    time_axis = [t_min + timedelta(seconds=i*cadence_seconds) for i in range(nt_total)]

    # --------------------------------------------
    # Spatial bins
    # --------------------------------------------
    nx = first_data.get("nx")
    ny = first_data.get("ny")
    if nx is None or ny is None:
        raise ValueError("First loaded file missing 'nx'/'ny' keys.")

    nx_bins = int(nx // bin_size)
    ny_bins = int(ny // bin_size)
    cube_binned = np.zeros((nx_bins, ny_bins, nt_total), dtype=np.uint16)

    provenance: Dict[Tuple[int,int,int], List[Tuple[str,int]]] = defaultdict(list)

    # --------------------------------------------
    # Process each file vectorized
    # --------------------------------------------
    for f, nt, hdrs, pixel_avg, times in file_data:
        if not times or not pixel_avg:
            continue

        # truncate pixel_avg to available frames
        n = min(len(pixel_avg), len(times))
        px_py_arr = np.array([pixel_avg[it] for it in range(n)])
        times_arr = np.array([times[it] for it in range(n)])

        # filter out invalid pixels
        valid_mask = (~np.isnan(px_py_arr[:,0])) & (~np.isnan(px_py_arr[:,1]))
        if not valid_mask.any():
            continue

        px_py_arr = px_py_arr[valid_mask]
        times_arr = times_arr[valid_mask]
        frame_indices = np.where(valid_mask)[0]

        # compute t_idx once
        t_idx_arr = ((np.array([(t - t_min).total_seconds() for t in times_arr]) / cadence_seconds).astype(int))

        # compute spatial bins
        bx_arr = (px_py_arr[:,0] // bin_size).astype(int)
        by_arr = (px_py_arr[:,1] // bin_size).astype(int)

        # filter bins in range
        in_range = (bx_arr >= 0) & (bx_arr < nx_bins) & (by_arr >= 0) & (by_arr < ny_bins) & (t_idx_arr >= 0) & (t_idx_arr < nt_total)
        bx_arr = bx_arr[in_range]
        by_arr = by_arr[in_range]
        t_idx_arr = t_idx_arr[in_range]
        frame_indices = frame_indices[in_range]

        # accumulate and update provenance
        for bx, by, t_idx, fi in zip(bx_arr, by_arr, t_idx_arr, frame_indices):
            cube_binned[bx, by, t_idx] = np.uint16(cube_binned[bx, by, t_idx] + 10)
            provenance[(bx, by, t_idx)].append((f, fi))

    return cube_binned, time_axis, provenance
    



def process_clusters(labeled, cube, time_axis, provenance, min_cluster_size=3):
    """
    Process all clusters in a binned CoM cube using global file and Stonyhurst caches.
    
    Returns:
        detections : list of (onset_time, mean_lon, mean_lat)
        assoc_files_list : list of sets of valid pixel files
    """

    num_labels = labeled.max()
    detections = []
    assoc_files_list = []

    # --- Global file + Stonyhurst cache ---
    file_cache = {}   # fname -> data dict
    stony_cache = {}  # fname -> (lon_arr, lat_arr)

    # Extract all unique files from provenance
    unique_files = {fname for lst in provenance.values() for fname, _ in lst}

    for fname in unique_files:
        try:
            with gzip.open(fname, "rb") as fh:
                data = pickle.load(fh)
                file_cache[fname] = data
                lon = np.array(data.get("CoM", {}).get("stonyhurst_lon", []))
                lat = np.array(data.get("CoM", {}).get("stonyhurst_lat", []))
                stony_cache[fname] = (lon, lat)
        except Exception:
            continue

    # --- Loop over clusters ---
    for lbl in range(1, num_labels + 1):
        coords = np.argwhere(labeled == lbl)
        voxel_count = coords.shape[0]
        if voxel_count < min_cluster_size:
            continue

        # --- Collect contributing files + frames ---
        contrib = []
        cluster_files = set()
        for x, y, z in coords:
            key = (x, y, z)
            if key in provenance:
                for fname, frame in provenance[key]:
                    contrib.append((fname, frame))
                    cluster_files.add(fname)

        if len(cluster_files) < min_cluster_size:
            continue

        # --- Pixel-space validation ---
        pxs, pys = [], []
        valid_pixel_files = []

        for fname, frame in contrib:
            data = file_cache.get(fname)
            if data is None:
                continue
            pixel_avg = data.get("CoM", {}).get("pixel_avg", [])
            if frame >= len(pixel_avg):
                continue
            try:
                px, py = pixel_avg[frame]
                px, py = float(px), float(py)
            except Exception:
                continue
            if np.isnan(px) or np.isnan(py):
                continue
            pxs.append(px)
            pys.append(py)
            valid_pixel_files.append(fname)

        if len(valid_pixel_files) < min_cluster_size:
            continue

        # --- CME onset time ---
        t_indices = coords[:, 2]
        onset_idx = t_indices.min()
        onset_time = time_axis[onset_idx]

        # --- Collect heliographic coordinates ---
        longs, lats = [], []
        for fname, frame in contrib:
            lonlat = stony_cache.get(fname)
            if lonlat is None:
                continue
            lon_arr, lat_arr = lonlat
            if frame >= len(lon_arr) or frame >= len(lat_arr):
                continue
            lon_val, lat_val = lon_arr[frame], lat_arr[frame]
            if np.isnan(lon_val) or np.isnan(lat_val):
                continue
            longs.append(lon_val)
            lats.append(lat_val)

        if len(longs) < min_cluster_size:
            continue

        # --- Robust mean without outliers ---
        def mean_no_outliers(arr):
            arr = np.asarray(arr)
            q1, q3 = np.percentile(arr, [25, 75])
            iqr = q3 - q1
            mask = (arr >= q1 - 1.5*iqr) & (arr <= q3 + 1.5*iqr)
            return float(np.mean(arr[mask]))

        mean_lon = mean_no_outliers(longs)
        mean_lat = mean_no_outliers(lats)

        detections.append((
            onset_time.strftime("%Y-%m-%d %H:%M"),
            round(mean_lon, 2),
            round(mean_lat, 2)
        ))

        assoc_files_list.append(set(valid_pixel_files))

    return detections, assoc_files_list




def detect_cme_origins_all_clusters(
    end_time: str,
    min_cluster_size: int = 3,
    bin_size: int = 80,
    time_range: int = 6
):
    """
    Detect ALL CME clusters using spatial clustering over all ALMANAC files.

    No temporal grouping. No selection of largest cluster.
    Every cluster with voxel-count >= min_cluster_size is processed.
    """

    # Load *all* ALMANAC pkl.gz files for that month
    dt_end = datetime.strptime(end_time, "%Y-%m-%d %H:%M")
    dt_start = dt_end - timedelta(hours=time_range)

    month_dir = Path("almanac_out") / f"{dt_end.year:04d}" / f"{dt_end.month:02d}"

    file_list = []
    for p in month_dir.rglob("*.pkl.gz"):
        name = p.name
        try:
            # extract YYYYMMDD_HHMM
            ts = name.split("_")[1] + "_" + name.split("_")[2]
            file_dt = datetime.strptime(ts, "%Y%m%d_%H%M")
        except Exception:
            continue  # skip malformed filenames

        if dt_start <= file_dt <= dt_end:
            file_list.append(str(p))

    file_list = sorted(file_list)
    if not file_list:
        print("No ALMANAC files found.")
        return [], [], []


    cube, time_axis, provenance = build_binned_com_cube(
        file_list, bin_size=bin_size
    )

    # label all 3D clusters
    structure = np.ones((3,3,3), dtype=bool)
    labeled, num_labels = label(cube > 0, structure=structure)

    # Process clusters efficiently
    detections, assoc_files_list = process_clusters(
        labeled, cube, time_axis, provenance,
        min_cluster_size=min_cluster_size
    )

    return detections, assoc_files_list, file_list

