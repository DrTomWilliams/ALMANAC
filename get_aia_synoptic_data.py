"""
get_aia_synoptic_data.py

Translated from provided IDL script to Python.

Dependencies:
    - requests

# ---------- Example usage ----------
if __name__ == "__main__":
    attempted_files, status = get_aia_synoptic_data(
        "2011-03-04 00:00", "2011-03-11 00:00",
        wlusr=[94, 131, 171, 193, 211, 335, 304],
        increment=10 * 60,
        verbose=False,
        nproc=8  # set to number of CPUs you want to use.
    )
    print(f"\nDownloaded {sum(status)} of {len(status)} successfully.")
"""

from datetime import datetime, timedelta
from urllib.parse import urlparse
from pathlib import Path
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import os
import time
import requests
from typing import List, Tuple, Optional




# ---------- Utility helpers ----------
def _zero_pad_wavelengths(wl: List[int], width: int) -> List[str]:
    fmt = "{:0" + str(width) + "d}"
    return [fmt.format(int(w)) for w in wl]




def _make_time_grid(tstart: datetime, tend: datetime, seconds: int) -> List[datetime]:
    times = []
    t = tstart
    while t <= tend:
        times.append(t)
        t += timedelta(seconds=seconds)
    return times




def _format_req_strings(dt: datetime, nrt: bool) -> Tuple[str, str, str]:
    reqdate11 = dt.strftime("%Y/%m/%d")
    reqdate08 = dt.strftime("%Y%m%d")
    reqtime08 = dt.strftime("%H%M%S") if nrt else dt.strftime("%H%M")
    return reqdate11, reqdate08, reqtime08




# ---------- Worker function ----------
def _download_task(args):
    """Worker: download a single FITS file"""
    url, localpath, verbose = args
    try:
        localdir = os.path.dirname(localpath)
        os.makedirs(localdir, exist_ok=True)
        if os.path.exists(localpath):
            if verbose:
                print(f"Exists: {localpath}")
            return (os.path.basename(localpath), True)

        with requests.get(url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(localpath, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
        if verbose:
            print(f"Downloaded: {localpath}")
        time.sleep(0.2)  # courtesy pause for JSOC
        return (os.path.basename(localpath), True)
    except Exception as e:
        if verbose:
            print(f"Failed: {url}\n   Error: {e}")
        return (os.path.basename(localpath), False)




# ---------- Main downloader ----------
def get_aia_synoptic_data(
    tstartusr: str,
    tendusr: str,
    wlusr: Optional[List[int]] = None,
    test: bool = False,
    increment: int = 60,
    topsavedir: Optional[str] = None,
    verbose: bool = False,
    nrt: bool = False,
    nproc: int = None,
) -> Tuple[List[str], List[bool]]:
    """
    Sequential-by-wavelength, parallel-by-time AIA synoptic downloader.
    """

    # Defaults
    if wlusr is None or len(wlusr) == 0:
        wlusr = [94, 131, 171, 193, 211, 335, 304]

    # Helper: parse dates
    def parse_datetime(s: str, end=False) -> datetime:
        s = s.strip().replace("/", "-")
        if len(s) < 11:
            s += " 23:59" if end else " 00:00"
        for fmt in ("%Y-%m-%d %H:%M", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d"):
            try:
                return datetime.strptime(s, fmt)
            except ValueError:
                continue
        raise ValueError(f"Unrecognized date format: {s}")

    tstart = parse_datetime(tstartusr)
    tend = parse_datetime(tendusr, end=True)
    
    # specify the top directory for data
    if topsavedir is None:
        topsavedir = os.path.join('Data',tend.strftime("%Y_%m_%d"))
    
    # Time grid
    req_times = _make_time_grid(tstart, tend, increment)
    reqdate11, reqdate08, reqtime08, reqhour = [], [], [], []
    for dt in req_times:
        r11, r08, rtime08 = _format_req_strings(dt, nrt)
        reqdate11.append(r11)
        reqdate08.append(r08)
        reqtime08.append(rtime08)
        reqhour.append("H" + rtime08[:2] + "00")

    savedir = topsavedir or os.getenv("AIA_SYNOPTIC") or os.getcwd()
    if nrt:
        savedir = os.path.join(savedir, "nrt")

    wlstr3 = _zero_pad_wavelengths(wlusr, 3)
    wlstr4 = _zero_pad_wavelengths(wlusr, 4)
    remotetop = "http://jsoc2.stanford.edu/data/aia/synoptic"
    if nrt:
        remotetop += "/nrt"

    attempted_files = []
    status_list = []

    # Number of processes
    nproc = nproc or max(1, cpu_count() - 2)

    for iwl, w in enumerate(wlusr):
        localdir = os.path.join(savedir, wlstr3[iwl])

        # Build all tasks for this wavelength
        tasks = []
        for i, dt in enumerate(req_times):
            fname = f"AIA{reqdate08[i]}_{reqtime08[i]}_{wlstr4[iwl]}.fits"
            url = f"{remotetop}/{reqdate11[i]}/{reqhour[i]}/{fname}"
            localpath = os.path.join(localdir, fname)
            attempted_files.append(fname)

            if test:
                print(localpath)
                print(url)
                status_list.append(False)
                continue

            tasks.append((url, localpath, verbose))

        if test:
            continue

        # Parallel download over time samples for this wavelength
        with Pool(processes=nproc) as pool:
           # add tqdm for progress tracking
            results = list(
                tqdm(
                    pool.imap(_download_task, tasks),
                    total=len(tasks),
                    desc=f"Downloading {wlusr[iwl]} Å",
                    unit=" file"
                )
            )

        # Collect results
        for fname, ok in results:
            status_list.append(ok)

    return attempted_files, status_list
