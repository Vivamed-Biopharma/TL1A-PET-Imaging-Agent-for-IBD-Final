# TL1A In Silico — Environment & Installation

This package is designed to be lightweight and runnable on a laptop. Choose one of the two setups below. Option A is minimal (no RDKit). Option B includes RDKit via conda-forge.

## Option A — Minimal (no RDKit)

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install --upgrade pip
pip install pandas biopython numpy scipy matplotlib seaborn SALib
# Optional but small:
pip install anarci
pip install jupyter ipykernel
python -m ipykernel install --user --name tl1a-min --display-name "Python (tl1a-min)"
```

Notes:
- This option runs all notebooks except the RDKit-specific cells in 11_rdkit_chelator.ipynb.
- If ANARCI fails to install on your platform, you can skip Notebook 02.

## Option B — With RDKit (recommended)

```bash
conda create -n tl1a python=3.10 -y
conda activate tl1a
conda install -c conda-forge rdkit prody pdb2pqr -y
pip install pandas biopython numpy scipy matplotlib seaborn SALib anarci jupyter ipykernel
python -m ipykernel install --user --name tl1a --display-name "Python (tl1a)"
```

Notes:
- RDKit is used in 11_rdkit_chelator.ipynb for computing small-molecule descriptors on NOTA and p-SCN-Bn-NOTA.
- prody and pdb2pqr are optional extras for future structure-based work.

## Verifying the environment

```bash
python -c "import pandas, Bio, numpy, scipy, matplotlib, SALib; print('OK')"
```

If you plan to run the optional visualization:

```bash
pip install py3Dmol
```

## Running the notebooks

From the repository root:

```bash
jupyter notebook
```

Open the notebooks under `tl1a_in_silico/` in order. All notebooks are self-contained and print tables or summary text.
