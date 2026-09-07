# CTC_Cluster_43: Reproducible Analysis

Reviewer-facing code and frozen cohort manifest for the manuscript:

**A circulating tumor cell cluster-associated transcriptional program shows reproducible donor-level enrichment but limited external prognostic evidence.**

## Scope

This repository contains only the final analysis workflow used for the manuscript. The discovery analysis treats biological CTC clusters as the experimental unit and uses one paired observation per donor. The external analysis evaluates the frozen 43-gene programme in a prespecified 35-patient TCGA primary-tumour cohort using progression-free interval (PFI).

## Main analysis checkpoints

- GSE109761 expression matrix: 19,955 genes × 357 expression samples
- Pure CTC-cluster cells: 106
- Pure CTC-single cells: 131
- Biological cluster units: 93
- Paired donors: 10
- Genome-wide biological genes tested: 19,870
- Exact sign configurations: 1,024
- Genes with FDR < 0.05: 0
- Frozen programme: `CTC_Cluster_43` (43 genes)
- Discovery programme effect: 1.151135
- Discovery exact two-sided P: 0.001953125
- TCGA patients: 35
- PFI events: 14
- HR per 1-SD score: 1.327032
- 95% CI: 0.706724–2.491799
- P: 0.378766
- C-index: 0.552469
- Median-split log-rank P: 0.660190

## Repository contents

```text
.
├── CTC_Cluster_43_reproducible_analysis.ipynb
├── CTC_Cluster_43_reproducible_analysis.py
├── requirements.txt
├── data_manifest/
│   └── tcga_frozen_manifest.csv
└── results/
    └── README.md
```

The large public expression files are **not** stored in Git. They are downloaded from the public GEO/GDC resources by the analysis workflow.

## Reproduction

### Google Colab

1. Open `CTC_Cluster_43_reproducible_analysis.ipynb` in Google Colab.
2. Install the packages from `requirements.txt`.
3. Run the notebook from top to bottom on CPU.
4. The notebook writes intermediate and final tables/checkpoints to `results/`.

### Local Python

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python CTC_Cluster_43_reproducible_analysis.py
```

Windows PowerShell:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python CTC_Cluster_43_reproducible_analysis.py
```

## Frozen external cohort

`data_manifest/tcga_frozen_manifest.csv` records the 35 selected TCGA primary-tumour samples and their exact GDC file UUIDs. This prevents a changing GDC catalogue from silently changing the validation cohort.

The external score is calculated as:

1. `log2(TPM + 1)`
2. gene-wise standardisation across the 35 TCGA samples (`ddof=1`)
3. equal-weight mean across the 43 frozen genes

Clinical PFI is merged at the patient level from the TCGA-CDR resource.

## Discovery design

Cells are not treated as independent biological replicates. CTC-cluster cells are grouped by donor and cluster ID; records with missing cluster IDs are retained as singleton biological units. Cluster-unit expression is averaged within donor, single-CTC expression is averaged within donor, and the paired cluster-minus-single difference is the donor-level observation.

The primary genome-wide test enumerates all 1,024 possible sign configurations of the 10 paired donor differences. Multiple testing is controlled with Benjamini–Hochberg FDR.

The 43-gene programme was frozen before external testing. No outcome-driven gene selection is performed in the TCGA cohort.

## Interpretation

The discovery result supports a reproducible donor-level transcriptional programme associated with CTC clustering. The external TCGA PFI association is not statistically significant and is not presented as validation of a prognostic biomarker. The TCGA analysis is an external test with limited sample size (35 patients; 14 PFI events).

The programme also contains platelet-associated biology. The final 43 genes have zero overlap with the predefined eight-gene canonical platelet panel used for the sensitivity analysis; this does not establish that platelet-related signal is absent.

## Data availability

The workflow uses public resources:

- GEO accession: `GSE109761`
- TCGA/GDC primary-tumour RNA-seq data
- TCGA-CDR clinical follow-up/PFI data

No private patient data are included in this repository.

## License

Code is released under the MIT License. Public datasets remain subject to their original database terms.
