def run_almanac(
    end_time="2013-05-15 07:12",
    waves=[94, 131, 171, 193, 211, 304, 335],
    cadence=6,
    tRange=6,
    minCluster=4,
    bin_size=16,
    kurtosis=False
):
    """
    Run the full ALMANAC workflow: download AIA files, process ALMANAC, detect CMEs, 
    and generate multiwave movies.

    Parameters
    ----------
    end_time : str
        Target time in format "YYYY-MM-DD HH:MM"
    waves : list of int
        Wavelengths to process, e.g. [131], [94,131,193]
    cadence : int
        Cadence in minutes
    tRange : int
        Time range in hours for the ALMANAC window
    minCluster : int
        Minimum number of channels for a detection
    """
    import warnings
    from datetime import datetime
    import os, time, pickle
    from concurrent.futures import ProcessPoolExecutor, as_completed

    # Custom routines
    from get_aia_synoptic_data import get_aia_synoptic_data
    from helper_functions import split_by_wavelength, compute_time_strings, snap_to_cadence
    from almanac import almanac_instance
    from SpaceTimeContinuum import detect_cme_origins_all_clusters
    from moviemaker import build_assoc_info, plot_cme_multiwave_parallel
    from sunpy.util.exceptions import SunpyUserWarning

    warnings.filterwarnings(
        "ignore",
        category=SunpyUserWarning,
        message="The conversion of these 2D helioprojective coordinates to 3D is all NaNs.*"
    )

    time1 = time.perf_counter()

    # Ensure end_time is compatible with cadence
    end_time_dt = datetime.strptime(end_time, "%Y-%m-%d %H:%M")
    end_time_snapped = snap_to_cadence(end_time_dt, cadence_minutes=cadence)
    end_time = end_time_snapped.strftime("%Y-%m-%d %H:%M")
    
    # Download files for desired wavelength(s)
    start_time, fits = compute_time_strings(end_time, offset=tRange)
    fits = os.path.join("Data",fits)
    attempted_files, status = get_aia_synoptic_data(
        start_time,
        end_time,
        wlusr=waves,
        increment=cadence*60,
        nproc=min(8, os.cpu_count())
    )
    
    # dict of successful files per wavelength
    wvl_files = split_by_wavelength(attempted_files, status)
    time2 = time.perf_counter()
    
    # Assign cores up to number of AIA channels and process ALMANAC
    max_workers = min(len(waves), os.cpu_count())
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(almanac_instance, [end_time, fits, wave, wvl_files]): wave for wave in waves}

        for future in as_completed(futures):
            wave = futures[future]
            try:
                future.result()  # just ensure it finished
                print(f"Channel {wave} completed")
            except Exception as e:
                print(f"Channel {wave} failed: {e}")
    
    time3 = time.perf_counter()
    
    # Cluster the detections, remove non-events, and find origin of all possible CMEs
    print('')
    detections, assoc_files,file_list = detect_cme_origins_all_clusters(end_time, 
                                        min_cluster_size=minCluster,
                                        bin_size=bin_size)
    
    # clean-up files not associated with a detection
    for f in file_list:
        if all(f not in file_set for file_set in assoc_files):
            print(f'Deleted: {f} not in any cluster > {minCluster} files')
            try:
                os.remove(f)
            except:
                continue
    print('')
    
    time4 = time.perf_counter()
    
    # Turn clusters into movies and output the detections
    det_dir = "detections"
    os.makedirs(det_dir, exist_ok=True)
    
    for i, det in enumerate(detections):
        try:
            # associated files with each detection
            file_list = assoc_files[i]
            
            # Detection time assumed to be first element
            det_time = det[0]
        
            # Make filename safe
            if hasattr(det_time, "strftime"):
                time_str = det_time.strftime("%Y%m%d_%H%M%S")
            else:
                time_str = str(det_time).replace(":", "").replace(" ", "_")

            fname = os.path.join(det_dir, f"detection_{time_str}.pkl")
            
            with open(fname, "wb") as f:
                pickle.dump(
                    {
                        "detection": det,
                        "associated_files": file_list,
                        "end_time": end_time,
                        "waves": waves,
                        "cadence": cadence,
                        "tRange": tRange,
                        "minCluster": minCluster,
                        "bin_size": bin_size
                    },
                    f,
                    protocol=pickle.HIGHEST_PROTOCOL
                )
            print(f"Saved detection {i} -> {fname}")

            assoc_info = build_assoc_info(file_list)
            plot_cme_multiwave_parallel(assoc_info, det)
        except Exception as e:
            print(f"Detection {i} failed: {e}")
            continue
    time5 = time.perf_counter()

    # Print timings
    print('\n===========================================================')
    print(f'Total Download Time: {time2 - time1:.2f} seconds')
    print(f'ALMANAC Processing Time: {time3 - time2:.2f} seconds')
    print(f'CME Detection/Grouping Time: {time4 - time3:.2f} seconds')
    print(f'Movie(s) Generation Time: {time5 - time4:.2f} seconds')
    print(f'Total Processing Time: {time5 - time1:.2f} seconds')
    
    return detections, assoc_files
