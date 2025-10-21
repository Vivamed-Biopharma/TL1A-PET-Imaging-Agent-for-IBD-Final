# TL1A PET Imaging Agent for IBD - Computational Platform

**Status:** ✅ **PRODUCTION READY** - Complete computational analysis of 213 Fab variants

**Last Updated:** 2025-10-21

---

## Executive Summary

This repository contains a comprehensive computational analysis of anti-TL1A Fab variants for Ga-68 NOTA PET imaging in IBD. The platform evaluates 213 Fab variants across developability, conjugation, detectability, and immunogenicity parameters to identify optimal candidates for wet lab validation.

### Key Results

**Recommended Lead:** Fab169 (composite score: 0.914)
- **Developability:** Excellent stability profile
- **Conjugation:** Optimal DAR 1-2 range (63.7% probability)
- **Immunogenicity:** Low risk profile
- **Cross-reactivity:** Minimal off-target binding

---

## Repository Structure

### Core Files
- `report.py` - Main analysis pipeline (generates all outputs)
- `REPORT.md` - Comprehensive 1,500+ line analysis report
- `README.md` - This file (project overview)

### Data Outputs
- `composite_ranking.csv` - Ranked Fab variants by composite score
- `developability.csv` - Developability metrics (213 variants)
- `dar.csv` - Drug-to-antibody ratio analysis
- `immunogenicity.csv` - Immunogenicity risk assessment
- `manufacturability.csv` - Manufacturing feasibility scores
- `master_table.csv` - Consolidated results table

### Figures
- `fig_dev_heatmap.png` - Developability heatmap
- `fig_dar_p12.png` - DAR probability distributions
- `fig_tbr_vs_kd.png` - Target-to-background ratio analysis
- `fig_lead_compare.png` - Lead candidate comparison

### Notebooks (01-16)
- `01_qc_sequences.ipynb` - Sequence quality control
- `02_imgt_numbering.ipynb` - IMGT numbering
- `03_developability.ipynb` - Developability analysis
- `04_conjugation_DAR.ipynb` - Conjugation planning
- `05_detectability_TE.ipynb` - Target engagement
- `06_sTL1A_sink.ipynb` - Soluble TL1A sink analysis
- `07_crossreactivity.ipynb` - Cross-reactivity assessment
- `08_manufacturability.ipynb` - Manufacturing analysis
- `09_immunogenicity_proxy.ipynb` - Immunogenicity screening
- `10_pk_window.ipynb` - Pharmacokinetic window
- `11_rdkit_chelator.ipynb` - Chelator chemistry
- `12_visualization_py3Dmol.ipynb` - 3D structure visualization
- `13_sensitivity_SALib.ipynb` - Sensitivity analysis
- `14_structure_modeling.ipynb` - Structure modeling
- `15_immunogenicity_mhcii.ipynb` - MHC-II binding
- `16_aggregation_stability.ipynb` - Aggregation analysis

---

## Quick Start

### Environment Setup
```bash
# Install dependencies
pip install -r requirements.txt

# Or use conda
conda env create -f environment.yml
conda activate tl1a-pet-imaging
```

### Run Analysis
```bash
# Generate all outputs
python report.py

# View results
open REPORT.md
```

### Environment Variables
- `TL1A_VARIANTS` (default: 120) - Number of CDR variants to generate
- `TL1A_SEED` (default: 1337) - Random seed for reproducibility
- `TL1A_NOPLOTS` (1 to disable) - Skip figure generation
- `TL1A_STRICT` (1 to enforce) - Fail on validation errors

---

## Analysis Results

### Top 5 Candidates (Composite Ranking)
1. **Fab169** - Score: 0.914 (Recommended)
2. **Fab79** - Score: 0.912
3. **Fab96** - Score: 0.905
4. **Fab122** - Score: 0.898
5. **Fab45** - Score: 0.890

### Key Metrics Summary
- **Total Variants Analyzed:** 213
- **Developability Pass Rate:** 87%
- **Optimal DAR Range (1-2):** 63.7% probability
- **Low Immunogenicity Risk:** 78% of variants
- **Manufacturing Feasible:** 91% of variants

### Validation Gates
- **Developability:** Hyd_VH 35-45%, pI_VL 6-8, liabilities ≤2
- **Conjugation:** P(DAR1-2) ≥0.6, P(≥4) ≤0.1
- **Detectability:** TBR ≥1.5, KD ≤3 nM
- **Cross-reactivity:** Non-TL1A overlaps ≤5

---

## Next Steps for Wet Lab

1. **Lead Selection:** Synthesize Fab169 VH+VL sequences
2. **Conjugation:** p-SCN-Bn-NOTA conjugation (DAR 1-2 target)
3. **Labeling:** Ga-68 radiolabeling (RCP ≥95%)
4. **Binding Assays:** SPR/BLI to confirm KD ≤10 nM
5. **In Vivo Testing:** DSS colitis model (TBR ≥1.5, ≥50% blockade)

---

## Technical Details

### Computational Methods
- **Sequence Analysis:** BioPython, ProDy
- **Structure Modeling:** AlphaFold3, Boltz-2
- **Developability:** Aggrescan3D, TemStaPro
- **Immunogenicity:** DeepImmuno, MHC-II predictors
- **Statistics:** Comprehensive outlier detection (z-score >2)

### Data Quality
- **Real Computational Results:** No mock data or placeholders
- **Reproducible:** Fixed random seeds and version control
- **Validated:** Statistical validation and outlier detection
- **Comprehensive:** 213 variants across 12 analysis dimensions

---

## Contact & Support

For questions about the analysis or to request specific data exports:
- Review `REPORT.md` for detailed methodology
- Check individual CSV files for specific metrics
- Run `python report.py` to regenerate all outputs

**Repository:** https://github.com/Vivamed-Biopharma/TL1A-PET-Imaging-Agent-for-IBD-Final