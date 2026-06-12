# run_single_cme_worker.py
import json
import sys
import os
from run_almanac import run_almanac
from datetime import datetime

if __name__ == "__main__":
    cme_time = sys.argv[1]

    # create logs directory
    log_dir = "logs"
    os.makedirs(log_dir, exist_ok=True)
    
    log_file = os.path.join(log_dir, f"worker_{cme_time.replace(':','').replace(' ','_')}.log")
    
    def log(msg):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(log_file, "a") as f:
            f.write(f"[{timestamp}] {msg}\n")

    try:
        log(f"Starting CME worker for {cme_time}")

        detect, assoc = run_almanac(
            end_time=cme_time,
            cadence=6,
            tRange=6,
            minCluster=4
        )

        assoc_json = [list(s) for s in assoc]

        log(f"Finished CME worker for {cme_time} successfully")

        print(json.dumps([detect, assoc_json]))

    except Exception as e:
        error = {"error": str(e), "type": type(e).__name__}
        print(json.dumps([error, None]))
        log(f"ERROR: {type(e).__name__}: {str(e)}")
        sys.exit(1)

