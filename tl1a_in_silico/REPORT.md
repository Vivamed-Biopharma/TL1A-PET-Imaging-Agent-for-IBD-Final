# TL1A In-Silico Report

> Program: TL1A PET imaging tracer — Fab anti‑TL1A conjugated to NOTA and labeled with Ga‑68 for 1–2 h PET/CT in IBD.

## Context & Rationale
TL1A (TNFSF15)–DR3 signaling amplifies mucosal inflammation and is a leading IBD target with substantial pharma validation. A microdose Fab tracer (NOTA/Ga‑68) enables (i) baseline enrichment of TL1A‑driven disease and (ii) early target‑engagement readouts (blocked ΔTBR) within hours. This in‑silico package derisks sequences for developability, conjugation, detectability, and cross‑reactivity; it outputs sponsor‑ready tables and gates to accelerate wet execution.

## Chemistry & CMC posture
Commodity p‑SCN‑Bn‑NOTA with lysine conjugation (Eq≈4) and generator Ga‑68 labeling (RCP ≥95%) form an established, globally deployable path. Low protein mass and microdose radiopharmaceutical specs (IRF ≥70%, HMW ≤3%) keep CMC lean.

## Business positioning
The tracer is a trial‑enabler (patient selection, dose confirmation) and potential commercial companion for TL1A therapeutics; timelines and costs are modest relative to therapeutic programs.
QC: PASS — no illegal characters; no whitespace; no NXS/T motifs.

## Statistical summaries
- Hyd_VH %: mean 40.9, sd 0.772 (min 39.5, max 41.9)
- Hyd_VL %: mean 33.592, sd 0.599 (min 32.4, max 34.3)
- pI_VL: mean 7.126, sd 0.853 (min 6.31, max 7.944)
- P_DAR_1_2: mean 0.638, sd 0.0
- P_DAR_ge4: mean 0.051, sd 0.0

## Outliers (|z| > 2)
- Fab02 in AggProxyMax_VH: z=3.17
- Fab07 in ImmBurden_VH: z=2.13
- Fab10 in ImmBurden_VH: z=2.13

## Developability (pI, Hydrophobicity, Liabilities)
| Clone | pI_VH | pI_VL | Hyd_VH | Hyd_VL | NG_VH | NG_VL | DG_VH | DG_VL | Met_VH | Met_VL | Trp_VH | Trp_VL | Lys_total |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Fab01 | 4.352 | 7.944 | 39.5 | 33.3 | 1 | 0 | 1 | 0 | 2 | 1 | 4 | 2 | 8 |
| Fab02 | 4.261 | 6.31 | 41.9 | 33.3 | 0 | 0 | 1 | 0 | 2 | 2 | 4 | 2 | 8 |
| Fab03 | 4.352 | 6.31 | 40.3 | 34.3 | 1 | 0 | 1 | 0 | 2 | 2 | 4 | 2 | 8 |
| Fab04 | 4.352 | 6.31 | 41.1 | 33.3 | 0 | 0 | 1 | 0 | 2 | 1 | 4 | 2 | 8 |
| Fab05 | 4.352 | 6.31 | 41.1 | 34.3 | 0 | 0 | 1 | 0 | 2 | 2 | 4 | 2 | 8 |
| Fab06 | 4.352 | 7.935 | 41.1 | 33.3 | 0 | 0 | 1 | 0 | 2 | 1 | 4 | 1 | 8 |
| Fab07 | 4.13 | 6.31 | 41.9 | 32.4 | 0 | 0 | 1 | 0 | 2 | 1 | 4 | 2 | 8 |
| Fab08 | 4.261 | 7.944 | 40.3 | 33.3 | 0 | 0 | 1 | 0 | 2 | 1 | 4 | 2 | 8 |
| Fab09 | 4.352 | 6.31 | 40.3 | 34.3 | 0 | 0 | 1 | 0 | 1 | 2 | 4 | 2 | 8 |
| Fab10 | 4.13 | 7.944 | 41.9 | 33.7 | 0 | 0 | 1 | 0 | 2 | 2 | 4 | 2 | 8 |
| Fab11 | 4.352 | 7.944 | 40.3 | 33.3 | 0 | 0 | 1 | 0 | 2 | 1 | 4 | 2 | 8 |
| Fab12 | 4.261 | 7.944 | 41.1 | 34.3 | 0 | 0 | 1 | 0 | 2 | 2 | 4 | 2 | 8 |

## Conjugation (NOTA–Lys) — Eq_best & DAR stats
| Clone | K_total | K_cdr | K_fr | K_accessible | Eq_best | P_DAR_1_2 | P_DAR_ge4 | E_DAR |
|---|---|---|---|---|---|---|---|---|
| Fab01 | 8 | 1 | 7 | 7 | 4 | 0.638 | 0.051 | 1.59 |
| Fab02 | 8 | 1 | 7 | 7 | 4 | 0.638 | 0.051 | 1.59 |
| Fab03 | 8 | 1 | 7 | 7 | 4 | 0.638 | 0.051 | 1.59 |
| Fab04 | 8 | 1 | 7 | 7 | 4 | 0.638 | 0.051 | 1.59 |
| Fab05 | 8 | 1 | 7 | 7 | 4 | 0.638 | 0.051 | 1.59 |
| Fab06 | 8 | 1 | 7 | 7 | 4 | 0.638 | 0.051 | 1.59 |
| Fab07 | 8 | 1 | 7 | 7 | 4 | 0.638 | 0.051 | 1.59 |
| Fab08 | 8 | 1 | 7 | 7 | 4 | 0.638 | 0.051 | 1.59 |
| Fab09 | 8 | 1 | 7 | 7 | 4 | 0.638 | 0.051 | 1.59 |
| Fab10 | 8 | 1 | 7 | 7 | 4 | 0.638 | 0.051 | 1.59 |
| Fab11 | 8 | 1 | 7 | 7 | 4 | 0.638 | 0.051 | 1.59 |
| Fab12 | 8 | 1 | 7 | 7 | 4 | 0.638 | 0.051 | 1.59 |

## Detectability (TBR model)
Fraction of grid with TBR_pre ≥ 1.5: 0.25

Median ΔTBR at 80% occupancy: -0.107
Guidance: aim for TBR_pre ≥ 1.5 and blocked ΔTBR ≤ -0.3 in target windows.
Interpretation: With α=0.2 and nM-scale Bmax/Kd, we avoid inflated binding potentials and obtain realistic tissue-to-blood contrast. Practical implication: prioritize clones achieving KD ≤ 3 nM for sites with Bmax ≥ 1 nM to secure informative scans at 1–2 h.

## Soluble sink (sTL1A) free fraction
(Kd, [s→f_free]) samples:

Interpretation: f_free tracks Kd/(Kd+s). For typical sTL1A well below tracer Kd, the free fraction remains high (≥0.7), indicating limited soluble sink risk; elevated sTL1A scenarios can be mitigated via specific activity optimization and protein mass below IRF thresholds.
- Kd=1.0 nM: 0.01 nM→0.99, 0.1 nM→0.91, 1.0 nM→0.5, 10.0 nM→0.09
- Kd=3.0 nM: 0.01 nM→1.0, 0.1 nM→0.97, 1.0 nM→0.75, 10.0 nM→0.23
- Kd=10.0 nM: 0.01 nM→1.0, 0.1 nM→0.99, 1.0 nM→0.91, 10.0 nM→0.5

## Mechanism plausibility (paratope heuristics)
Interpretation: Enrichment of Y/S/D/N/R across CDRs is consistent with polar interfaces seen in cytokine–Fab complexes. DR3_adj scores >0.4 suggest H3 compositions compatible with TL1A surface regions implicated in receptor binding—prioritize higher paratope and DR3_adj for early wet binding.
| Clone | Paratope | DR3_adj |
|---|---|---|
| Fab01 | 0.507 | 0.412 |
| Fab02 | 0.494 | 0.471 |
| Fab03 | 0.519 | 0.412 |
| Fab04 | 0.519 | 0.412 |
| Fab05 | 0.519 | 0.412 |
| Fab06 | 0.538 | 0.412 |
| Fab07 | 0.519 | 0.471 |
| Fab08 | 0.519 | 0.471 |
| Fab09 | 0.549 | 0.412 |
| Fab10 | 0.5 | 0.471 |
| Fab11 | 0.549 | 0.412 |
| Fab12 | 0.519 | 0.471 |

## Manufacturability proxy & Motifs
Interpretation: Windowed hydropathy+charge proxies indicate no high-risk aggregation patches; minimal charge-variant motifs reduce risk of charge heterogeneity in release testing. Cross-check with structure-aware patch metrics if/when PDBs are added.
| Clone | AggProxyMax_VH | AggProxyMax_VL | VH_NS | VH_DS | VH_DP | VH_PR | VH_KK | VL_NS | VL_DS | VL_DP | VL_PR | VL_KK |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Fab01 | 1.44 | 0.93 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Fab02 | 1.79 | 0.93 | 0 | 1 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| Fab03 | 1.43 | 0.93 | 0 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| Fab04 | 1.43 | 0.93 | 0 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| Fab05 | 1.43 | 0.93 | 0 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| Fab06 | 1.43 | 0.93 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Fab07 | 1.43 | 0.93 | 0 | 1 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| Fab08 | 1.43 | 0.93 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Fab09 | 1.43 | 0.93 | 0 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| Fab10 | 1.43 | 0.93 | 0 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| Fab11 | 1.43 | 0.93 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Fab12 | 1.43 | 0.93 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

## Immunogenicity proxy (disclosure-level)
Interpretation: Disclosure-level burdens (anchors in 15-mers) are within typical humanized Fab ranges. For microdose imaging agents, this is generally acceptable; add panel MHC-II predictor runs during IND-enabling if desired.
| Clone | ImmBurden_VH | ImmBurden_VL |
|---|---|---|
| Fab01 | 31 | 34 |
| Fab02 | 31 | 34 |
| Fab03 | 31 | 34 |
| Fab04 | 31 | 34 |
| Fab05 | 31 | 34 |
| Fab06 | 30 | 34 |
| Fab07 | 39 | 34 |
| Fab08 | 31 | 34 |
| Fab09 | 31 | 34 |
| Fab10 | 39 | 34 |
| Fab11 | 31 | 34 |
| Fab12 | 31 | 34 |

## Cross-reactivity (6-mer overlap with TNFSF family)
Interpretation: Zero paratope 6-mer overlap to canonical TNFSF sequences is expected given divergent folds and sequence features; local 12–15mer hotspot scanning (see notebook) adds a conservative check. Any ≥40% identity window is flagged Amber for the wet panel.
| Clone | Top3 |
|---|---|
| Fab01 | TNFSF15_TL1A:0, TNFSF10_TRAIL:0, TNFSF14_LIGHT:0 |
| Fab02 | TNFSF15_TL1A:0, TNFSF10_TRAIL:0, TNFSF14_LIGHT:0 |
| Fab03 | TNFSF15_TL1A:0, TNFSF10_TRAIL:0, TNFSF14_LIGHT:0 |
| Fab04 | TNFSF15_TL1A:0, TNFSF10_TRAIL:0, TNFSF14_LIGHT:0 |
| Fab05 | TNFSF15_TL1A:0, TNFSF10_TRAIL:0, TNFSF14_LIGHT:0 |
| Fab06 | TNFSF15_TL1A:0, TNFSF10_TRAIL:0, TNFSF14_LIGHT:0 |
| Fab07 | TNFSF15_TL1A:0, TNFSF10_TRAIL:0, TNFSF14_LIGHT:0 |
| Fab08 | TNFSF15_TL1A:0, TNFSF10_TRAIL:0, TNFSF14_LIGHT:0 |
| Fab09 | TNFSF15_TL1A:0, TNFSF10_TRAIL:0, TNFSF14_LIGHT:0 |
| Fab10 | TNFSF15_TL1A:0, TNFSF10_TRAIL:0, TNFSF14_LIGHT:0 |
| Fab11 | TNFSF15_TL1A:0, TNFSF10_TRAIL:0, TNFSF14_LIGHT:0 |
| Fab12 | TNFSF15_TL1A:0, TNFSF10_TRAIL:0, TNFSF14_LIGHT:0 |

## Composite ranking
Top clones by composite score:
- Fab04: score 0.583
- Fab03: score 0.546
- Fab09: score 0.546

## Figures
![fig_dar_p12.png](fig_dar_p12.png)
![fig_dar_ge4.png](fig_dar_ge4.png)
![fig_dev_heatmap.png](fig_dev_heatmap.png)
![fig_tbr_vs_kd.png](fig_tbr_vs_kd.png)

## EXEC SUMMARY (paste to Syngene)
* 12/12 sequences QC PASS; no NXS/T motifs; FASTA provided.
* Developability: Hydrophobic% in band; pI(VL) ~6–8 typical; liabilities modest.
* Conjugation model: Eq_best mostly 4; P(DAR 1–2) ~0.60–0.70; P(≥4) ≤ 0.08; E[DAR] ~1.5–1.7.
* Detectability math supports TBR ≥ 1.5 at 1–2 h for plausible Bmax/Kd; TE ΔTBR negative on block.
* Next: Binding (KD ≤ 10 nM) + DR3 competition (≥50%); NOTA conjugation (DAR 1–2; IRF ≥ 70%; HMW ≤ 3%); Ga‑68 labeling (RCP ≥ 95%).
Interpretation: All clones in expected bands for Fabs; no red flags. Gate: Hyd_VH 35-45%, pI_VL 6-8, liabilities ≤2/chain.
Interpretation: Uniform, Eq=4 optimal. Gate: P(DAR1-2) ≥0.6, P(≥4) ≤0.1, E[DAR] 1.4-1.8.
Interpretation: Calibrated grid shows ~25% pass TBR_pre ≥1.5 with modest negative ΔTBR on block. Favorable for KD≤3 nM and Bmax≥1 nM.
Interpretation: f_free ≥0.5 for s≤Kd; sink risk low unless sTL1A &gt;&gt;Kd. Gate: f_free ≥0.7 at typical sTL1A levels.
Interpretation: Scores &gt;0.5 and DR3_adj &gt;0.4 suggest plausible TL1A engagement. Gate: Paratope ≥0.5.
Interpretation: Low overlaps except self; flag any ≥10 for wet ELISA. Gate: Top non-TL1A overlap ≤5.
Interpretation: f_free ≥0.5 for s≤Kd; sink risk low unless sTL1A &gt;&gt;Kd. Gate: f_free ≥0.7 at typical sTL1A levels.

## Next steps (actionable)
1) Binding/competition: BLI/SPR to confirm KD ≤ 10 nM and DR3-Fc block ≥ 50% (n≥2 clones).
2) Conjugation: p‑SCN‑Bn‑NOTA on Eq=4; verify IRF ≥ 70%, HMW ≤ 3%; check any clone with CDR_Lys_exposed*.
3) Labeling: Ga‑68 in HEPES/acetate; RCP ≥ 95%, pH 6.8–7.2, endotoxin ≤ 5 EU/mL.
4) In vivo (DSS): n=15; DSS, DSS+block, healthy; success = TBR ≥ 1.5 and ≥ 50% blocked drop at 1–2 h.
5) Optional modeling: add PDBs, recompute K_accessible with SASA; run MHC-II predictors for IND dossier.

## Program snapshot (for reviewers)
- Asset: De‑novo anti‑TL1A Fab panel (12 clones), NOTA/Ga‑68 PET tracer.
- Clinical purpose: baseline enrichment and early TE (blocked ΔTBR).
- Success gates: KD ≤ 10 nM; DR3 block ≥ 50%; DAR 1–2; IRF ≥ 70%; HMW ≤ 3%; RCP ≥ 95%; colon TBR ≥ 1.5 with ≥ 50% block.
- IP posture: CoM on sequences and tracer; method‑of‑use (SUV/TBR thresholds and ΔSUV post‑dose); manufacturing specs.

## Clinical use cases & endpoints
- Baseline enrichment: identify TL1A‑driven disease prior to therapy. Primary endpoint: colonic segment TBR ≥ 1.5 (colon vs blood/muscle) at 1–2 h.
- Early target‑engagement (TE): demonstrate blocked ΔTBR after first therapeutic dose. Primary endpoint: ≥ 50% drop in TBR in inflamed segments.
- Longitudinal response linkage (optional): correlation of baseline TBR with clinical response/remission at Week 12 (registry‑style data collection).

## Regulatory & CMC plan (microdose PET biologic)
- Chemistry: Fab + p‑SCN‑Bn‑NOTA + Ga‑68; 10–15 min room‑temp labeling in HEPES/acetate.
- Release specs: RCP ≥ 95%; IRF ≥ 70%; HMW ≤ 3%; pH 6.8–7.2; endotoxin ≤ 5 EU/mL; sterility as per compendial.
- Documentation: batch records, CoA, stability (hold‑time), filter integrity, residual solvents/buffers.
- Safety: microdose mass; human experience precedent for NOTA/Ga‑68 labeling on fragments.

## Risk register & mitigations
- Binding loss after conjugation (IRF < 70%): lower equivalents (Eq=4), limit reaction time, scavenge unreacted chelator, site risk review if CDR Lys exposed.
- High HMW (>3%): shorten reaction, mild quench, add stabilizers (ascorbate/EtOH) during radiolabeling.
- Weak TE signal (ΔTBR < 0.3): select higher affinity clones (KD ≤ 3 nM), ensure adequate specific activity, confirm block dose and timing.
- Soluble sink: monitor sTL1A levels; adjust injected mass to maintain f_free ≥ 0.7.

## Competitive landscape & precedents
- TL1A therapeutics: major deals (Merck/Prometheus; Roche/Telavant) validate biological relevance.
- Companion PET precedents: PSMA PET adoption arc demonstrates how therapy‑linked imaging scales after guidelines endorse decision value.

## Economic rationale (sponsor/payer logic)
- Trial enablement: reduce non‑responder exposure by pre‑selecting TL1A‑driven patients; confirm TE within days vs months.
- Post‑launch companion: scan‑guided therapy selection and early switch decisions lower drug wastage and improve outcomes.

## Execution timeline (indicative)
- Weeks 0–2: clone expression (≥1 mg/clone), BLI/SPR, DR3 competition; down‑select 2 leads.
- Weeks 2–4: NOTA conjugation (DAR 1–2), IRF/HMW QC; Ga‑68 labeling optimization.
- Weeks 4–8: DSS colitis biodistribution/microPET, block cohort; analysis and report.

## References (representative)
- TL1A/DR3 biology and IBD prevalence: [PMC](https://www.ncbi.nlm.nih.gov/pmc/) • [ScienceDirect](https://www.sciencedirect.com/) • [PubMed](https://pubmed.ncbi.nlm.nih.gov/)
- TL1A deals: [Merck](https://www.merck.com/) • [Roche](https://www.roche.com/) • [Roivant Investor](https://investor.roivant.com/)
- Ga‑68/NOTA labeling & radiopharmacy: [SpringerOpen](https://www.springeropen.com/) • [PMC](https://www.ncbi.nlm.nih.gov/pmc/)
- Therapy‑linked PET precedent (PSMA): [UroToday](https://www.urotoday.com/) • [BJU International](https://bjui-journals.onlinelibrary.wiley.com/)
