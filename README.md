# Drug-Resistance Prediction from Genotype — *Staphylococcus aureus*

A genomic-AMR (antimicrobial resistance) pipeline that tries to predict real, lab-measured drug
resistance (**methicillin**, **erythromycin**) directly from a bacterial genome's curated
AMR-gene content — real genome, gene-annotation, and phenotype-test data fetched **live from the
BV-BRC (Bacterial and Viral Bioinformatics Resource Center) Data API**, joined and modeled with an
MLST-aware (clonal-lineage-aware) train/test split. Implemented in **Python**
(scikit-learn, primary) and **R** (`glm`/`randomForest`, twin).

## Aim

Given a bacterial genome's curated set of known AMR genes, can a classical-ML model predict whether
that isolate will test resistant to a given antibiotic — and if it can't, why not?

## Objective

Fetch real genome, curated AMR-gene, and lab-measured phenotype data for a real bacterial species
from BV-BRC, join them into a genome × gene feature matrix, train classical ML models under an
MLST-aware (clonal-lineage-leakage-aware) train/test split, and honestly report whether gene
presence/absence predicts measured resistance — including tracing *why* it doesn't, if it doesn't,
back to the real data itself rather than blaming the model.

## Data fetch

Two real data pivots happened before landing on this design — first species (*Klebsiella
pneumoniae* → *S. aureus*, after directly measuring that only ~7-18 K. pneumoniae genomes had both
a real phenotype label and a curated gene call), then feature design (hand-picked mecA/erm columns
→ the full curated AMR-gene set, after finding even S. aureus's mecA/erm-specific overlap was only
8-2 genomes). Full diagnostic detail lives in the project's private teaching files (not pushed
here).

Final real, live-fetched data, from BV-BRC (taxon_id `1280`, *Staphylococcus aureus*), for two
drugs: methicillin (classic driver: the *mecA*/*mecC* cassette) and erythromycin (classic driver:
the *erm* gene family).

- **`genome_amr`** (real lab phenotype-test results, broth microdilution/disk diffusion AST):
  25,000 raw rows → 5,697 after filtering to these two drugs → 221 methicillin-labeled isolates
  (51.6% resistant), 441 erythromycin-labeled isolates (39.2% resistant).
- **`sp_gene`** (curated AMR-gene calls from genome annotation): 25,000 raw rows → 5,597 rows with
  `property == "Antibiotic Resistance"` → 189 distinct curated genes.
- **`genome`** metadata (incl. MLST sequence type, used for the leakage-aware split): 25,000 rows,
  94.1% with non-null MLST.

## Data describe

`sp_gene` pivoted into a genome × gene presence/absence matrix, kept to genes seen in ≥5 genomes
(119 genes total — a broad, unfiltered feature set, not hand-picked to mecA/erm; see Results).
LEFT-joined to the phenotype labels (a genome with no gene hits is kept as a valid all-zero row,
not dropped — an earlier inner-join version silently deleted susceptible isolates). Final modeling
frames: **meth_frame (219 × 122)**, **ery_frame (439 × 122)**.

## Methods / Workflow — what we did

```
BV-BRC Data API (live RQL query, taxon_id=1280)
   │ genome_amr: filter to methicillin/erythromycin, evidence="Laboratory Method"
   │ sp_gene: filter to property="Antibiotic Resistance", pivot wide (genome x gene)
   │ genome: MLST sequence type per genome
   ▼
   LEFT join labels -> gene features, fillna(0)
   ▼
   meth_frame (219 x 122), ery_frame (439 x 122)
   ▼
   MLST-aware split (no clonal group crosses train/test) + random split, both with class-balance retry
   ▼
   Logistic Regression + Random Forest (Python: scikit-learn; R: glm + randomForest)
   ▼
   evaluate: sensitivity, specificity, AUROC  +  BV-BRC's own AdaBoost baseline (Step 4b)
   ▼
   feature importance -- does mecA/erm actually rise to the top on their own?
   ▼
   direct verification: mecA/erm presence-fraction by class -- the real, final finding
```

1. Fetch real phenotype labels, curated AMR-gene calls, and MLST metadata live from BV-BRC.
2. LEFT-join into genome × gene feature matrices (methicillin and erythromycin separately),
   deliberately keeping gene-negative genomes rather than dropping them.
3. Split each frame two ways: an MLST-aware split (no clonal lineage crosses train/test) and a
   plain random split, both with a class-balance retry.
4. Train Logistic Regression and Random Forest (Python: scikit-learn; R: `glm` + `randomForest`),
   compared against BV-BRC's own AdaBoost baseline.
5. Evaluate with sensitivity, specificity, AUROC.
6. Inspect feature importance — do mecA/erm actually rise to the top on their own?
7. Directly verify the textbook-expected causal genes' presence-fraction by class as the final,
   decisive check (not inferred from metrics alone).

## Results

| Drug | Split | Model | Sensitivity | Specificity | AUROC |
|---|---|---|---|---|---|
| methicillin | mlst-aware | LogReg / RF | 1.00 | 0.00 | 0.50 |
| erythromycin | mlst-aware | LogReg / RF | 0.00 | ~0.97-0.99 | ~0.47-0.49 |
| — | R (`glm`/`randomForest`), all 8 drug x model x split combos | | | | 0.42-0.58 |

Every combination, in both languages, landed at or near 0.5 (chance) — with the models degenerating
to "always predict one class" rather than discriminating. Feature importance in both languages was
dominated by universal housekeeping/target genes (`gyrA`, `gyrB`, `rpoB`, `rpoC`, ribosomal protein
genes) rather than the textbook-expected `mecA`/`erm` genes.

**The real finding, confirmed by direct measurement, not inferred from the poor metrics:**
`meth_frame.loc[meth_frame['label']==1, 'mecA'].mean()` (and every other mec-family/erm-family
column, both classes) reads **exactly 0.000**. Not weak signal — zero genomes anywhere in either
219-row methicillin or 439-row erythromycin modeling frame carry any mecA-cassette or erm-family
gene, resistant or not.

R's logistic regression additionally surfaced a distinct statistical artifact worth naming
separately: **quasi-complete separation** — several unrelated genes' coefficients converged to the
exact same arbitrary value (e.g. `17.63586` shared by 8 methicillin genes), a known small-sample
`glm()` fitting-procedure artifact, confirmed distinct from the missing-gene finding because
`randomForest`'s importances on the identical data did not show the same pattern.

## Biology interpretation of results

The causal resistance genes are essentially absent from this specific BV-BRC join of
phenotype-tested genomes, because phenotype testing (wet-lab AST) and curated gene annotation
(computational, from WGS) are two independently-run pipelines that only partially overlap in which
genomes they cover. This is a **data-availability finding**, not a modeling failure or a claim that
gene-presence can't predict AMR in general — confirmed independently in both Python and R, ruling
out a language-specific bug as the explanation. Biologically, mecA-cassette carriage and erm-family
methylase genes are the textbook-established causal drivers of methicillin and erythromycin
resistance in *S. aureus* respectively; their total absence from both modeling frames means the
models were never given the actual causal signal to learn from, which is exactly why they
degenerated to predicting the majority class — a model cannot discover a relationship whose
predictor variable never varies in the data it's shown, no matter how sophisticated the algorithm.
The models correctly learned that unrelated housekeeping/target genes (`gyrA`, `gyrB`, `rpoB`,
`rpoC`) carried whatever weak, likely-spurious signal was present in this particular join, rather
than any biologically causal AMR mechanism.

## Learning through project

A negative or null result, run down to its root cause by direct measurement rather than left as an
unexplained "the model didn't work," is a real and useful finding — here, tracing a chance-level
AUROC all the way back to a **zero-variance causal feature** (mecA/erm presence exactly 0.000 in
both classes) turns an apparent modeling failure into a specific, actionable finding about this
data join's coverage gap between two independently-run BV-BRC pipelines (wet-lab phenotype testing
vs. computational gene annotation). Two dataset pivots (species, then feature design) before landing
on a workable design is itself a normal and expected part of real bioinformatics work, not a sign of
a flawed plan — the earlier attempts directly informed why this final design was chosen. Confirming
a surprising negative result independently in a second language (R glm/randomForest vs. Python
scikit-learn) is what makes it trustworthy — it rules out a single-language implementation bug as
the explanation, and R additionally surfaced a genuinely distinct statistical artifact
(quasi-complete separation) worth naming separately rather than conflating with the missing-gene
finding.

## Limitations

This is a negative result specific to this exact join, snapshot, and species/drug pair — not a
general claim about genomic AMR prediction (which has real, published successes on datasets
purpose-built to jointly have both phenotype and genotype data, e.g. CRyPTIC for TB). `sp_gene`
presence/absence reflects BV-BRC's own curation pipeline and reference-database version, and a
genome absent from `sp_gene` is treated identically to a genome confirmed gene-negative (BV-BRC's
schema doesn't distinguish the two at the row level). MLST is a coarse clonal-lineage proxy and
reduces but doesn't eliminate leakage between close genetic relatives.

## Reproduce

**Python:** run `amr_saureus_prediction.ipynb` top to bottom (installs its own dependencies via
`%pip install`). Requires internet access (live BV-BRC API queries) and produces `data_py/` and
`results_py/`.

**R:** run the Python notebook first (produces `data_py/meth_frame.csv` and `data_py/ery_frame.csv`),
then run `amr_saureus_prediction.R`. Requires `randomForest` and `pROC`:

```r
install.packages(c("randomForest", "pROC"))
```

## Tech

`Python` (scikit-learn) · `R` (`glm`, `randomForest`) · BV-BRC Data API · MLST-aware leakage-safe
splitting · Logistic Regression · Random Forest

## Files

```
amr_saureus_prediction.ipynb   # Python pipeline: fetch -> label -> feature matrix -> split -> model -> evaluate -> interpret
amr_saureus_prediction.R       # R twin: same pipeline on the same joined bridge CSVs, glm + randomForest
results_py/                    # model_performance.csv, feature_importance.csv (Python)
results_R/                     # model_performance_R.csv, feature_importance_R.csv (R)
```

`data_py/` (the fetched/joined intermediate CSVs the R script reads) is not tracked — it's fully
re-derivable by running the Python notebook, which fetches live from the BV-BRC API.

## License

All rights reserved — see `LICENSE`. This repository is public for portfolio/demonstration purposes
only; no permission is granted to copy, modify, or reuse any part of it.
