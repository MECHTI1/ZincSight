ZincSight
=========
**Interpretable prediction of zinc ion locations in proteins**

ZincSight is a tool for predicting zinc ion binding sites in protein structures. It is designed for **high-throughput analysis**, while keeping predictions **accurate** and **interpretable**.  

Although trained on zinc-binding sites, ZincSight also detects potential **transition-metal binding sites**, including sites for **iron (Fe), manganese (Mn), cobalt (Co), copper (Cu), and nickel (Ni)**, in experimental structures and predicted models such as AlphaFold (see the [paper](https://doi.org/10.1002/pro.70350)).

## Features

### Flexible input options
1. Protein IDs (e.g., PDB, AlphaFold, or ESM Metagenomic Atlas) 
2. Uploaded structure files in PDB or MMCIF format

### Comprehensive output:
  - Detailed CSV files including (among others):
    * predicted zinc-ion x,y,z coordinates
    * predicted binding residues
    * probability for being TP site
    * references for experimentally solved zinc-binding templates used in prediction
  - PyMOL sessions for easy visualization

For guidance on interpreting the output, see [Interpreting predictions](docs/INTERPRETING_RESULTS.md).

## Quick Start with Google Colab
The fastest way to get started with ZincSight is through our interactive Google Colab notebooks:

- **ZincSight Colab:** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/MECHTI1/ZincSight/blob/master/ZincSight.ipynb)
- **ZincSight BatchDrive Colab:** Built for massively scalable screening across millions of structures, with native Google Drive integration for large input/output workflows  
  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/MECHTI1/ZincSight/blob/master/ZincSight_MaxDrive.ipynb)

  
## Local Installation

### Downloading the repository
```bash
git clone https://github.com/MECHTI1/ZincSight.git
cd ZincSight
```

### Dependencies installation
```bash
pip install -r requirements.txt
```

### ***Optional***: adding a new zinc-binding motif template to the training set

In order to add a new zinc-binding motif from a PDB structure that was not available at the time of training, use the `add_template.py` script with the following arguments:

| Argument | Description                             |
|----------|-----------------------------------------|
| PDB_ID   | PDB identifier                          |
| RESIDUES | Comma-separated pairs of CHAIN_RESIDUEs |

Example execution:
```shell
python ./add_templates 4f3w B_89,B_56,B_92
```

### ZincSight execution

After installing the environment, you can run ZincSight using the `zincsight.py` script with the following arguments:

| Argument                                          | Description                                  |
|---------------------------------------------------|----------------------------------------------|
| `-r`, `--rotamers`                                | Include histidine rotamers (default: no)     |
| `-i <IDENTIFIERS>`, `--identifiers <IDENTIFIERS>` | Comma-separated identifiers                  |
| `-f <FILE>`, `--file <FILE>`                      | Path of a text file with identifier per line |
| `-c <NUMBER>`, `--cores <NUMBER>`                 | Number of cores (default: 2)                 |


At least one of `IDENTIFIERS` or `FILE` is required.

Script execution with the included example input file:
```shell
python ./zincsight.py -f example_input.txt
```

## Published benchmark (zinc-binding sites)

Results reported for the original ZincSight method and comparison tools on 27 non-redundant AlphaFold2 models with 44 experimentally verified zinc-binding sites ([Mechtinger et al., *Protein Science*, 2025, Table 1](https://doi.org/10.1002/pro.70350)). Reference zinc positions were mapped from the corresponding experimental structures; predictions were evaluated at 2 Å and 5 Å distance thresholds.

| Method | PR-AUC (2 Å) | PR-AUC (5 Å) |
|---|---:|---:|
| **ZincSight (with histidine rotamer sampling)** | **0.85** | **0.98** |
| ZincSight (without histidine rotamer sampling) | 0.83 | 0.90 |
| Metal3D (3+) | 0.82 | 0.93 |
| Master of Metals | 0.68 | 0.85 |
| PinMyMetal | 0.68 | 0.75 |
| AlphaFold 3 (1 Zn) | 0.39 | 0.77 |
| AlphaFold 3 (2 Zn) | 0.58 | 0.76 |
| AlphaFold 3 (3 Zn) | 0.59 | 0.72 |
| AlphaFold 3 (4 Zn) | 0.64 | 0.69 |

On this dataset (mean protein length: 485 amino acids), the study reported mean runtimes of 0.87 s per structure for ZincSight with histidine rotamer sampling and 0.47 s without, using two Intel i7 CPU cores. Reported comparison-tool means were 2.56 s for Master of Metals, 17.4 s for PinMyMetal and 39.8 s for Metal3D (Section 2.7). These are study measurements, not Tamarind turnaround estimates.

## Large-scale application: N-Zinclopedia

ZincSight was used to screen **1,334,703 protein and domain structure representatives** for the [N-Zinclopedia database](https://bentallab.github.io/nzinclopedia/), as reported in the [2026 N-Zinclopedia preprint](https://doi.org/10.64898/2026.10.18.759936). This total counts structures screened, rather than proteins confirmed to bind a metal. The database combines ZincSight predictions with MAHOMES-II classifications of candidate catalytic and non-catalytic sites. Experiments are needed to confirm metal binding and catalytic activity at individual sites.

The demonstrated screening scale, CPU inference and experimental coordination templates make ZincSight useful for large protein batches and for inspecting the basis of site predictions.

## Citation

If you use ZincSight or its results, please cite:

> Mechtinger G, Axel G, Kolodny R, Ben-Tal N. Interpretable prediction of zinc ion location in proteins with ZincSight. *Protein Science.* 2025;34(11):e70350. https://doi.org/10.1002/pro.70350
