# ALMANAC
## Automated Detection of CoronaL MAss Ejecta origiNs for Space Weather AppliCations

[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)]()
[![License](https://img.shields.io/badge/license-MIT-green.svg)]()

ALMANAC is an automated framework for detecting, tracking, and characterising the early signatures and source regions of coronal mass ejections (CMEs) using multi-wavelength solar observations.

The pipeline combines solar imaging observations, automated feature detection, spatio-temporal clustering, and statistical analysis to identify CME-related disturbances and provide quantitative event properties suitable for space-weather applications.

---

# Publications and Usage

If you use ALMANAC in scientific research, please cite:

> **T. Williams & H. Morgan (2022) Space Weather, 20 (11), e2022SW003253**
> **T. Williams, C. B. Prior, D. MacTaggart, H. Morgan (2026), ApJ, In Review**

Additional publications using or extending ALMANAC:

> **O. P. M. Aslam, D. MacTaggart, T. Williams, L. Fletcher, P. Romano (2024) MNRAS, 534, 444**
> **T. Williams, C. B. Prior, D. MacTaggart (2025) ApJ, 980, 102**


---

# Overview

ALMANAC processes solar observations through the following stages:

1. Retrieval of solar observations
2. Image preprocessing and normalisation
3. Multi-frame ratio image generation
4. Automated CME signature detection
5. Spatio-temporal clustering
6. CME source region identification
7. Event association and characterisation
8. Output generation and visualisation

The pipeline is designed for both:

- retrospective scientific analysis
- near-real-time space-weather monitoring

---

# Pipeline Structure

The repository is organised as follows:

```text
ALMANAC/
│
├── run_almanac.py
│   Main ALMANAC execution script and pipeline interface
│
├── almanac.py
│   Core ALMANAC processing routines:
│   - image preparation
│   - ratio image generation
│   - CME feature extraction
│   - kurtosis time-series generation
│
├── SpaceTimeContinuum.py
│   CME detection and event association:
│   - spatio-temporal clustering
│   - CME origin identification
│   - cluster processing
│   - event association products
│
├── moviemaker.py
│   Visualisation utilities:
│   - frame generation
│   - CME detection movies
│
├── get_aia_synoptic_data.py
│   Solar data acquisition:
│   - AIA data retrieval
│   - parallel download utilities
│
├── helper_functions.py
│   Shared numerical and analysis utilities
│
├── notebooks/
│   │
│   ├── ALMANAC Pipeline.ipynb
│   │   Main workflow examples:
│   │   - running ALMANAC
│   │   - batch CME processing
│   │   - output inspection
│   │
│   └── kurtosis.ipynb
│       Analysis of kurtosis-based CME signatures
│
├── data/
│   Input solar observations
│
├── output/
│   Generated ALMANAC products:
│   - detections
│   - event associations
│   - diagnostic files
│
├── logs/
│   Runtime logs from batch processing
│
└── README.md
```


---

# Installation

## Requirements

ALMANAC requires:

- Python >= 3.10
- SunPy
- Astropy
- NumPy
- SciPy
- Pandas
- Matplotlib
- Scikit-image
- ImageIO
- Zarr
- tqdm
- Multiprocessing support

---

## Create environment

Using conda:
**bash**
```
conda create -n almanac python=3.10
conda activate almanac

pip install numpy scipy pandas matplotlib
pip install astropy sunpy
pip install scikit-image imageio zarr tqdm
```

# Running ALMANAC

## Single Event Processing

A single ALMANAC run can be executed using:

```python
from run_almanac import run_almanac

lasco_time = '2025-11-11 10:30'

detection, files = run_almanac(
    end_time=lasco_time,
    cadence=6,
    tRange=6,
    minCluster=4
)
```

where:

- `end_time` is the CME observation time
- `cadence` is the temporal sampling interval (minutes)
- `tRange` defines the analysis window (hours)
- `minCluster` defines the minimum number of detections required for clustering

The function returns:

- `detection`: detected CME candidate information
- `files`: associated output products

---

## Batch Processing

For processing multiple CME events, a wrapper script using subprocess execution and JSON logging is recommended.

The batch workflow:

- isolates individual CME runs
- stores stdout/stderr logs
- records failed events
- enables large-scale catalogue generation

An example implementation is provided in:

```
notebooks/ALMANAC Pipeline.ipynb
```

---

## Kurtosis Analysis

ALMANAC generates kurtosis-based time-series products describing statistical changes in solar image sequences.

Examples of analysing these products are provided in:

```
notebooks/kurtosis.ipynb
```
