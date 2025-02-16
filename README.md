ZincSight
=========

**ZincSight** is a rapid and precise tool for large-scale zinc ion location prediction in proteins. It can be installed locally on Ubuntu using a Conda environment and used interactively via Google Colab.

ZincSight predicts zinc-binding sites in proteins. It accepts input as:

- **Protein IDs** (e.g., PDB, AlphaFold, or ESM Metagenomic Atlas)
- **Uploaded structure files** in PDB or MMCIF format

The tool is optimized for high-throughput analyses while maintaining accuracy.

## Installation & Usage Guide 
Google Colab
-----------------------------
For an interactive session, use our dedicated Google Colab notebook:

**[ZincSight Colab Notebook](https://colab.research.google.com/github/MECHTI1/ZincSight/blob/master/ZincSight.ipynb)** 

(Ubuntu + Conda) 
-----------------------------

1. **Clone the Repository**

   ```bash
   git clone https://github.com/your-username/ZincSight.git
   cd ZincSight
To set up the required environment for ZincSight, use the following commands:

```bash
conda env create -f environment.yml
conda activate zincsight_env
```

## Contact
- **Email:** mechtinger1@mail.tau.ac.il  
