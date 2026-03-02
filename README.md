# Fixed Income Assignment 1


## Quick Start

### 1. Create Environment
```powershell
python -m venv venv
.\venv\Scripts\Activate.ps1
```

### 2. Install Dependencies
```powershell
pip install -r requirements.txt
```

### 3. Run Notebook
- **Part 1**: Open `Assignment_1.ipynb` and run all 

## Codebase Structure

```
src/fi_pricing/
├── curves/              # Yield curve models
│   ├── base.py          # Abstract base class
│   ├── nss.py           # Nelson-Siegel-Svensson
│   ├── calibrator.py    # NSS calibration
│   └── zcy_extractor.py # Zero-coupon yield extraction
├── models/              # Interest rate models
│   ├── affine.py        # Affine model base
│   ├── one_factor.py    # Vasicek, CIR, Hull-White
│   └── twoFG.py         # Two-Factor Gaussian (G2++)

src/part2/
└── analysis.py          # CIR/EKF estimation & weekly yields

src/helpers/
└── two_factor_gaussian_calibration.py  # Cap pricing utilities

Data/                    # Market data file
```
