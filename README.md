# SSMM — Simplified Soil Moisture Model

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)

**SSMM** is a one-dimensional, multi-layer soil water balance model that simulates soil
moisture dynamics using Darcy's Law and Richards' Equation. It was built for reproducibility:
give it daily precipitation, potential evapotranspiration (PET), and a handful of soil
properties, and it returns daily volumetric soil moisture for a 5-layer soil column.

A small sample dataset is included so you can run the model immediately — see
[Quickstart](#quickstart) below.

---

## Model overview

| | |
|---|---|
| Dimensions | 1D vertical soil column |
| Layers | 5 |
| Layer thicknesses | `300, 350, 400, 450, 500` mm (2000 mm total depth) |
| Time step | 1 hour internally (24 sub-steps per forcing day), leapfrog integration |
| Forcing | Daily precipitation and PET |
| ET extraction | First layer only |
| Initial soil moisture | Field capacity |

**Units**

| Quantity | Unit |
|---|---|
| Precipitation | mm/day |
| PET | mm/day |
| Volumetric water content | fraction (0–1) per layer |
| Hydraulic conductivity (Ks) | mm/hour |
| Matric potential (Ψ) | mm |

**Lower boundary conditions**

| Option | Description |
|---|---|
| `gravitational` (default) | Free drainage — gradient set by gravity only |
| `ground_water` | Hydraulic head gradient of 0 at the lower boundary |
| `no_flow` | Sealed bottom — zero flux at the lower boundary |

---

## Installation

```bash
git clone https://github.com/Li-N-K/SSMModel.git
cd SSMModel
python3 -m venv .venv
source .venv/bin/activate      # on Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Requires Python 3.8+, `numpy`, and `pandas`.

> **Using Anaconda?** Install into a fresh virtual environment as shown above rather than
> the `base` environment or an existing conda env. A `base` environment with a stale/broken
> `pyarrow` install can make `import pandas` fail with a NumPy 1.x/2.x ABI error that has
> nothing to do with this repo — a clean venv sidesteps it entirely.

---

## Quickstart

Run the model on the included sample data (two years of daily precipitation and PET for
Birmingham, AL):

```bash
python SSMModel.py \
  --precip sample_data/precip.csv \
  --pet sample_data/pet.csv \
  --soil sample_data/soil.json \
  --output soil_moisture_output.csv
```

This writes `soil_moisture_output.csv` with daily volumetric soil moisture for each of the
5 layers plus total column soil moisture.

---

## Running on your own data

### 1. Precipitation CSV

Two columns: a date and daily precipitation in mm/day.

```csv
Date,Precip
2020-01-01,4.2
2020-01-02,0.0
```

### 2. PET CSV

Two columns: a date and daily potential evapotranspiration in mm/day.

```csv
Date,PET
2020-01-01,3.1
2020-01-02,4.5
```

Both files need a parseable date column and should cover the same date range. If you only
have daily temperature, you'll need to derive PET yourself first (e.g. Hargreaves-Samani,
Penman-Monteith) — SSMM takes PET as an input, it does not compute it.

### 3. Soil properties JSON

```json
{
  "Qs": 0.45,
  "Field Capacity": 0.30,
  "Wilting Point": 0.15,
  "Ks(mm/h)": 10.0,
  "Psi(mm)": 100.0,
  "b": 4.9
}
```

| Key | Meaning |
|---|---|
| `Qs` | Saturated volumetric water content |
| `Field Capacity` | Volumetric water content at field capacity |
| `Wilting Point` | Volumetric water content at wilting point |
| `Ks(mm/h)` | Saturated hydraulic conductivity |
| `Psi(mm)` | Matric potential at saturation |
| `b` | Pore size distribution index (Clapp-Hornberger `b`) |

Then run:

```bash
python SSMModel.py \
  --precip precip.csv \
  --pet pet.csv \
  --soil soil.json \
  --lower_boundary gravitational \
  --output model_output.csv
```

---

## Output

A CSV with daily volumetric soil moisture per layer and total column soil moisture in mm:

```csv
Date,layer1,layer2,layer3,layer4,layer5,Total_Soil_Moisture
2020-01-01,0.30,0.30,0.30,0.30,0.30,600.00
2020-01-02,0.28,0.29,0.30,0.30,0.30,595.35
```

---

## Repository structure

```
SSMModel/
├── SSMModel.py       # model + CLI entry point
├── sample_data/       # example precip/PET/soil inputs (Birmingham, AL, 2021–2022)
├── requirements.txt
└── LICENSE
```

---

## Reproducibility notes

- An internal dummy row is appended so the final real day is simulated correctly, then dropped
  from the output.
- The model is deterministic and scriptable for batch workflows.
- Designed to be modular enough to slot into larger hydrological or climate pipelines.

### Fixes since the original thesis version

- **Pandas compatibility**: forcing series were indexed positionally (`series[j]`) relying on
  pandas' old integer-fallback behavior on a `DatetimeIndex`. This was removed in pandas 2.x and
  raised a `KeyError`. Switched to explicit `.iloc[j]`.
- **Leapfrog bootstrap bug**: on the very first sub-step of the simulation, the central-difference
  update referenced `theta[i - 2]`, which at `i = 1` wraps around (via NumPy's negative indexing)
  to the last, still-uninitialized row of the array instead of a real prior state. This made every
  run collapse to the wilting point within the first day. The first step now falls back to a
  forward-difference update using the actual initial condition; all subsequent steps are
  unchanged.

## License

MIT — see [LICENSE](LICENSE).
