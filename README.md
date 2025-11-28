# SLM - Stochastic Landscape Model

**SLM** is a stochastic landscape-scale hydrological-biogeochemical model designed to simulate discharge and concentration time series in river networks for different landscape and hydroclimatic scenarios.

## 📁 Project Structure

```text
variland/
├── 📂 images/                    # example GIF
├── 📂 river_network/            
│   ├── OCN_example.RData         # Example river network 
├── 📂 scripts/                   # Core simulation scripts
│   ├── 🐍 python/               
│   │   ├── functions.py          # Core Python functions
│   └── 📊 r/                    
│       ├── functions.R           # Core R functions
│       └── run_SLM.R             # simulation runner
├── 📂 simulations/               # Configuration files
│   ├── config_example.json       # Example configuration template
├── 📂 renv/                      # R environment management
├── 🐳 environment.yml            # Python environment specification
├── 🔧 run_SLM.sh                 # execution script
├── 📄 LICENSE                    # License information
└── 📋 README.md                  # This documentation
```

## Quick Start

### Prerequisites

- **R** (≥4.2) with `renv` package
- **Python** (≥3.11) with Conda

### 1. Environment Setup

#### 🐍 Python Environment

Create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate SLM
```

#### 📊 R Environment

Restore the R dependencies:

```r
renv::restore()
```

### 2. Configuration Files

The `simulations/` folder contains `.json` config files. These config files specify system paths, model parameters and simulation settings. An example configuration file is provided in `simulations/config_example.json`.

### 3. Running Simulations

Execute the main simulation:

```bash
./run_SLM.sh
```

**Output**: Results are saved as `.parquet` files in a results folder defined in the `config_example.json`.

## Example Results

![Variland Simulation Animation](images/example_animation.gif)

Example animation showing spatiotemporal water quality dynamics in an artificial river network.

### Complete Animation Gallery

Explore **32 exemplary animations** showcasing different scenarios:

**[View Animation Gallery →](https://syncandshare.desy.de/index.php/s/oSM5ynsJgq4dH3x)**

## 📄 License & Copyright

### Copyright Notice

Copyright © 2025 **UFZ-Helmholtz-Centre for Environmental Research**

### License

- This software is licensed under the **GNU General Public License v3.0 (GPL-3.0)**. Please see the [LICENSE](LICENSE) file for the complete license text and legal terms.
- The documentation and example data are licensed under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license.

### Citation

This software is archived on Zenodo: [![DOI](https://zenodo.org/badge/1103193816.svg)](https://doi.org/10.5281/zenodo.17749061)

A publication using this software is in review.

### Disclaimer

This software is provided "as is", without warranty of any kind, express or implied. The authors and copyright holders shall not be liable for any claim, damages, or other liability arising from the use of this software.

### Contact

For questions, collaboration inquiries, or reporting issues:

- **Author**: Linus Schauer
- **Institution**: UFZ - Helmholtz-Centre for Environmental Research
- **GitHub Issues**: [Report a bug or request a feature](https://github.com/linusschauer/SLM/issues)