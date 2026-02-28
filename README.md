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

### 3. Run Main Notebook
Open `Assignment_1.ipynb` in Jupyter or VS Code and run all cells.

## Codebase Structure

```
src/fi_pricing/
├── curves/              # Yield curve models
│   ├── base.py          # Abstract base class for curves
│   ├── nss.py           # Nelson-Siegel-Svensson implementation
│   ├── calibrator.py    # NSS parameter calibration
│   └── zcy_extractor.py # Extract zero-coupon yields from bonds
│
├── models/              # Interest rate models
│   ├── affine.py        # Abstract base for affine models
│   ├── one_factor.py    # Vasicek, CIR, Hull-White (1-factor)
│   └── twoFG.py         # Two-Factor Gaussian (G2++)
│
scripts/
└── part2_run.py         # Batch runner for Part 2 exercises

Data/                    # Market data (bonds, swaps, caps)
part2_outputs/           # Generated plots and results
```

### Key Features
- **NSS Calibration**: Fit yield curves to market bond prices
- **Affine Models**: Closed-form pricing for bonds and options (Jamshidian decomposition)
- **G2++ Model**: Two-factor Gaussian model with cap/swaption calibration
- **Vectorized**: Numpy-based batch operations for performance