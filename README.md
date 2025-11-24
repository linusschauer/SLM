# 🌊 SLM - Stochastic Landscape Model

**SLM** is a stochastic landscape-scale hydrological-biogeochemical model designed to simulate discharge and concentration time series in river networks. The explorative model was developed to analyze the **space-time variance of water quality** in river networks in different landscape and hydroclimatic scenarios.

## 📁 Project Structure

```text
variland/
├── 📂 images/                    # Project assets and visualizations
├── 📂 river_network/            
│   ├── OCN_example.RData         # Example river network 
├── 📂 scripts/                   # Core simulation scripts
│   ├── 🐍 python/               # Python routing algorithms
│   │   ├── functions.py          # Core Python functions
│   └── 📊 r/                     # R simulation framework
│       ├── functions.R           # Core R functions
│       └── run_variland.R        # Main simulation runner
├── 📂 simulations/               # Configuration files for different scenarios
│   ├── config_example.json      # Example configuration template
├── 📂 renv/                      # R environment management
├── 🐳 environment.yml            # Conda environment specification
├── 🔧 run_variland.sh            # Main execution script
├── 📄 LICENSE.md                # License information
└── 📋 README.md                  # This documentation
```

## 🚀 Quick Start

### Prerequisites

- **R** (≥4.2) with `renv` package
- **Python** (≥3.11) with Conda

### 1. Environment Setup

#### 🐍 Python Environment

Create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate stocha
```

> The Python environment handles computationally intensive discharge and load routing through the river network.

#### 📊 R Environment

Restore the R dependencies:

```r
renv::restore()
```

> The R environment manages model simulations and river network topology operations.


### 2. Configuration Files

The `simulations/` folder contains pre-configured `.json` config files. These config files specify system paths, model parameters and simulation settings. An example configuration file is provided in `simulations/config_example.json`.

### 3. Running Simulations

Execute the main simulation:

```bash
./run_SLM.sh
```

**Output**: Results are saved as `.parquet` files in a results folder defined in the `config_example.json`.

## 📊 Example Results

![Variland Simulation Animation](images/example_animation.gif)

Example animation showing spatiotemporal water quality dynamics in an artificial river network.

### 🎬 Complete Animation Gallery

Explore **32 exemplary animations** showcasing different scenarios:

**[View Animation Gallery →](https://syncandshare.desy.de/index.php/s/oSM5ynsJgq4dH3x)**

## 📄 License

Please see the file [LICENSE.md](LICENSE.md) for further information about how the content is licensed.