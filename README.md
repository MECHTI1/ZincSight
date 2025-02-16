ZincSight
=========

A rapid and precise tool for large-scale zinc ion location prediction in proteins. The tool is optimized for high-throughput analyses while maintaining accuracy.

## Features

- Fast processing: Approximately 0.5–1 second per protein query on Google Colab—achieved by utilizing two threads on an Intel Xeon® CPU @ 2.20GHz, or alternatively, using two local 11th Gen Intel® Core™ i7-1165G7 CPU cores.
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

## 1. Quick Start with Google Colab
The fastest way to get started with ZincSight is through our interactive Google Colab notebook:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/MECHTI1/ZincSight/blob/master/ZincSight.ipynb)

## 2. Local Installation

```bash
git clone https://github.com/MECHTI1/ZincSight.git
cd ZincSight
conda env create -f environment.yml
conda activate zincsight_env
python -c "print('ZincSight Environment OK')"
