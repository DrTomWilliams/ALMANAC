# moviemaker.py
import os
import subprocess
import gc
import gzip
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from matplotlib.colors import Normalize
from pathlib import Path
from glob import glob
from astropy.io import fits
from concurrent.futures import ProcessPoolExecutor
import imageio.v2 as imageio
import sunpy.map
from dateutil import parser
from datetime import timedelta
from copy import deepcopy
from multiprocessing import Pool, cpu_count
from helper_functions import *



def ensure_divisible_by_16(img):
    """Resize/pad image so both dimensions are divisible by 16 (for FFmpeg compatibility)."""
    h, w = img.shape[:2]
    new_h = int(np.ceil(h / 16) * 16)
    new_w = int(np.ceil(w / 16) * 16)

    if (new_h, new_w) != (h, w):
        # Handle grayscale, RGB, and RGBA images
        if img.ndim == 2:  # grayscale
            padded = np.zeros((new_h, new_w), dtype=img.dtype)
            padded[:h, :w] = img
        else:  # color image (3 or 4 channels)
            c = img.shape[2]
            padded = np.zeros((new_h, new_w, c), dtype=img.dtype)
            padded[:h, :w, :] = img
        return padded

    return img
    
    
    

def round_to_cadence(dt, cadence_minutes=6):
    """
    Round datetime dt to nearest multiple of `cadence_minutes`.
    """
    total_minutes = dt.hour*60 + dt.minute
    remainder = total_minutes % cadence_minutes
    if remainder < cadence_minutes / 2:
        delta = -remainder
    else:
        delta = cadence_minutes - remainder
    rounded_dt = dt + timedelta(minutes=delta, seconds=-dt.second, microseconds=-dt.microsecond)
    return rounded_dt




def load_maps_for_time(info, frame_idx):
    """
    Safely retrieve AIA and mask maps for a specific frame index.

    If frame_idx is None (no data for that time):
        - AIA map is black
        - Mask map is white (bright, not dark)
    """
    ref_map = info["aiaim"][0]  # reference map for shape & metadata
    shape = ref_map.data.shape
    meta = ref_map.meta

    # --- Handle missing frame ---
    if frame_idx is None:
        # Black AIA map
        aia_map = sunpy.map.Map(np.zeros(shape, dtype=np.float32), meta)
        # White mask map
        mask_map = sunpy.map.Map(np.full(shape, 255, dtype=np.float32), meta)
        return aia_map, mask_map

    # --- Use actual maps if available ---
    aia_map = info["aiaim"][frame_idx]
    mask_map = info["maskim"][frame_idx]

    return aia_map, mask_map




def get_frame_index(info, target_time):
    """
    Find frame index corresponding to target_time (datetime object).
    Returns None if the time is outside the range.
    """
    times = info["times"]
    if target_time < times[0] or target_time > times[-1]:
        return None
    # Find closest time
    for i, t in enumerate(times):
        if abs((t - target_time).total_seconds()) < 60:  # allow up to 60s tolerance
            return i
    return None




def generate_frame(args):
    t, files, assoc_info = args

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    import matplotlib.patheffects as PathEffects
    import gc

    fig_width = assoc_info["_fig_width"]
    fig_height = assoc_info["_fig_height"]
    pairs_per_row = assoc_info["_pairs_per_row"]

    nrows = int(np.ceil(len(files) / pairs_per_row))
    pair_width = 1.0 / pairs_per_row
    half_width = pair_width / 2.
    pair_height = 1.0 / nrows

    fig = plt.figure(figsize=(fig_width, fig_height))

    for idx, f in enumerate(files):
        row = idx // pairs_per_row
        col_in_row = idx % pairs_per_row

        left_aia = col_in_row * pair_width
        left_mask = left_aia + half_width
        bottom = 1.0 - (row + 1) * pair_height

        info = assoc_info["info_by_file"][f]
        wave = info["wave"]
        frame_idx = get_frame_index(info, t)
        aia_map, mask_map = load_maps_for_time(info, frame_idx)

        aia_arr = mgn(aia_map.data)
        lo, hi = np.nanpercentile(aia_arr, (20, 99.5))
        norm = Normalize(vmin=lo, vmax=hi)

        ax_aia = fig.add_axes([left_aia, bottom, half_width, pair_height])
        ax_aia.imshow(aia_arr, cmap=plt.get_cmap(f"sdoaia{wave}"),
                      norm=norm, origin="lower")
        ax_aia.set_xticks([]); ax_aia.set_yticks([])
        for spine in ax_aia.spines.values():
            spine.set_linewidth(0)

        txt = ax_aia.text(0.125, 0.08, f"AIA {wave}",
                          transform=ax_aia.transAxes,
                          ha='center', va='top',
                          fontsize=10, color='white')
        txt.set_path_effects(
            [PathEffects.withStroke(linewidth=2.5, foreground='black')]
        )

        mask_arr = mask_map.data.astype(float)
        ax_mask = fig.add_axes([left_mask, bottom, half_width, pair_height])
        ax_mask.imshow(mask_arr, cmap="gray", vmin=0, vmax=1,
                       origin="lower")
        ax_mask.set_xticks([]); ax_mask.set_yticks([])
        for spine in ax_mask.spines.values():
            spine.set_linewidth(0)

        if idx == len(files) - 1:
            ax_mask.text(0.25, 0.08, t.strftime('%Y-%m-%d %H:%M'),
                         transform=ax_mask.transAxes,
                         ha='center', va='top',
                         fontsize=10, color='red')

    # --- Extract RGB array (Matplotlib ≥3.8 safe) ---
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba())
    frame_array = buf[..., :3].copy()

    plt.close(fig)
    gc.collect()

    return frame_array
    
    
    

def build_assoc_info(file_list, cadence_minutes=6):
    info_by_file = {}
    global_times = set()
    file_order = []

    for f in file_list:
        with gzip.open(f, "rb") as fp:
            info = pickle.load(fp)

        hdrs = info["hdr"]
        raw_times = [parser.parse(h["T_OBS"]) for h in hdrs]
        times = [round_to_cadence(t, cadence_minutes) for t in raw_times]

        info_by_file[f] = {
            "wave": info["wave"],
            "aiaim": info["aiaim"],       # keep SunPy Maps intact
            "maskim": info["maskim"],
            "times": times,
            "nx": info.get("nx", 1024),
            "ny": info.get("ny", 1024),
        }

        file_order.append(f)
        global_times.update(times)

    return {
        "file_list": file_order,
        "frame_times": sorted(global_times),
        "info_by_file": info_by_file,
    }




def plot_cme_multiwave_parallel(
    assoc_info,
    detection,
    movies_dir="./movies",
    pairs_per_row=2
):
    import re
    import imageio.v2 as imageio
    from multiprocessing import Pool, cpu_count
    from pathlib import Path
    import numpy as np

    name = (
        f'{detection[0].replace("-","_").replace(" ","_").replace(":","")}'
        f'_{("N"+str(round(detection[2]))).replace("N-","S")}'
        f'{("W"+str(round(detection[1]))).replace("W-","E")}'
    )

    files = sorted(
        assoc_info["file_list"],
        key=lambda f: tuple(map(int, re.findall(r'_wvlnth_(\d+)_(\d+)', f)[0]))
    )

    nwaves = len(files)

    movies_dir = Path(movies_dir)
    movies_dir.mkdir(exist_ok=True)
    movie_path = movies_dir / f"{name}_channels_{nwaves}.mp4"

    # --- FIGURE SIZE SCALING ---
    nrows = int(np.ceil(nwaves / pairs_per_row))
    ref_figsize = (12, 4.625)
    fig_width = ref_figsize[0]
    fig_height = ref_figsize[1] * (nrows / 2)

    assoc_info["_fig_width"] = fig_width
    assoc_info["_fig_height"] = fig_height
    assoc_info["_pairs_per_row"] = pairs_per_row

    tasks = [(t, files, assoc_info) for t in assoc_info["frame_times"]]

    print("Starting parallel frame generation + streaming to video...")

    with Pool(min(4, cpu_count())) as pool, \
         imageio.get_writer(
             movie_path,
             format="FFMPEG",     # force ffmpeg
             mode="I",
             fps=30,
             codec="libx264",
             pixelformat="yuv420p"
         ) as writer:

        for frame in pool.imap_unordered(generate_frame, tasks, chunksize=2):

            # Ensure divisible by 16 (ffmpeg requirement)
            h, w = frame.shape[:2]
            new_h = int(np.ceil(h / 16) * 16)
            new_w = int(np.ceil(w / 16) * 16)

            if (new_h, new_w) != (h, w):
                padded = np.zeros((new_h, new_w, 3), dtype=np.uint8)
                padded[:h, :w] = frame
                frame = padded

            writer.append_data(frame)

    print(f"Video saved to: {movie_path}\n")
