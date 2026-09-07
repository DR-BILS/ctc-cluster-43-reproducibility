# CTC_Cluster_43 corrected reproduction script
# Generated from CTC_Cluster_43_reproducible_analysis.ipynb
# Run top-to-bottom. The notebook is the primary reviewer interface.


# ======================================================================
# # CTC_Cluster_43 — corrected reproducible analysis
# 
# This notebook is the **reviewer-facing reproduction workflow for the corrected manuscript**:
# 
# > *A circulating tumor cell cluster-associated transcriptional program shows reproducible donor-level enrichment but limited external prognostic evidence.*
# 
# The notebook supersedes the previous cell-level/UVM–PAAD–KIRC workflow. It implements the corrected donor-level discovery design, freezes the 43-gene programme before external testing, and reproduces the 35-patient TCGA external PFI analysis.
# 
# **Important:** run from top to bottom. Do not reuse the old 50/42-gene or cell-level analyses.
# ======================================================================


# Cell 2
# 1. Environment and dependencies
import sys, os, platform, subprocess, importlib, gzip, io, tarfile, tempfile, hashlib, itertools, time, json, math
required = {
    "pandas":"pandas", "numpy":"numpy", "scipy":"scipy",
    "statsmodels":"statsmodels", "matplotlib":"matplotlib",
    "lifelines":"lifelines", "requests":"requests", "GEOparse":"GEOparse",
    "openpyxl":"openpyxl"
}
missing=[]
for pip_name, mod_name in required.items():
    try:
        importlib.import_module(mod_name)
    except ImportError:
        missing.append(pip_name)
if missing:
    subprocess.run([sys.executable, "-m", "pip", "install", "-q"] + missing, check=True)

import pandas as pd, numpy as np, scipy
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import requests
from lifelines import CoxPHFitter, KaplanMeierFitter
import GEOparse

SEED = 42
np.random.seed(SEED)
DATA_DIR = "data"
RESULTS_DIR = "results"
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

print("Python:", sys.version.split()[0])
print("CPU-only workflow; GPU not required.")



# ======================================================================
# ## 2. Frozen CTC_Cluster_43 definition
# 
# The 43 genes below are frozen **before** the TCGA external test. Equal weights are used.
# 
# The original robust set contained 44 genes. `NABP1` was removed because it formed a singleton heterogeneous module. No TCGA outcome was used to alter the final gene list.
# ======================================================================


# Cell 4
SIGNATURE_GENES = ['RSU1', 'UBA7', 'MOB3C', 'ENDOD1', 'TMEM40', 'KCNK6', 'ZYX', 'GNB5', 'C10orf54', 'IL8', 'MIR22HG', 'ERN1', 'TST', 'ACTN1', 'TPST2', 'MAFG', 'RAB32', 'CD97', 'GNA13', 'MFN2', 'FBXO18', 'PHTF2', 'LCN2', 'PRKCD', 'STX16', 'TANGO2', 'BNIP2', 'CD68', 'RAP2B', 'SLC6A6', 'EIF4G3', 'SEPT5-GP1BB', 'HIPK2', 'KIAA0513', 'AMPD2', 'SERPINB1', 'TAB2', 'IL17RA', 'CENPU', 'SMOX', 'STK4', 'POU2F1', 'TLN1']
GENCODE_V36 = {'RSU1': 'ENSG00000148484', 'UBA7': 'ENSG00000182179', 'MOB3C': 'ENSG00000142961', 'ENDOD1': 'ENSG00000149218', 'TMEM40': 'ENSG00000088726', 'KCNK6': 'ENSG00000099337', 'ZYX': 'ENSG00000159840', 'GNB5': 'ENSG00000069966', 'C10orf54': 'ENSG00000108219', 'IL8': 'ENSG00000169429', 'MIR22HG': 'ENSG00000186594', 'ERN1': 'ENSG00000178607', 'TST': 'ENSG00000128311', 'ACTN1': 'ENSG00000072110', 'TPST2': 'ENSG00000128294', 'MAFG': 'ENSG00000197063', 'RAB32': 'ENSG00000118508', 'CD97': 'ENSG00000123146', 'GNA13': 'ENSG00000120063', 'MFN2': 'ENSG00000116688', 'FBXO18': 'ENSG00000134452', 'PHTF2': 'ENSG00000006576', 'LCN2': 'ENSG00000148346', 'PRKCD': 'ENSG00000163932', 'STX16': 'ENSG00000124222', 'TANGO2': 'ENSG00000183597', 'BNIP2': 'ENSG00000140299', 'CD68': 'ENSG00000129226', 'RAP2B': 'ENSG00000181467', 'SLC6A6': 'ENSG00000131389', 'EIF4G3': 'ENSG00000075151', 'SEPT5-GP1BB': 'ENSG00000284874', 'HIPK2': 'ENSG00000064393', 'KIAA0513': 'ENSG00000135709', 'AMPD2': 'ENSG00000116337', 'SERPINB1': 'ENSG00000021355', 'TAB2': 'ENSG00000055208', 'IL17RA': 'ENSG00000177663', 'CENPU': 'ENSG00000151725', 'SMOX': 'ENSG00000088826', 'STK4': 'ENSG00000101109', 'POU2F1': 'ENSG00000143190', 'TLN1': 'ENSG00000137076'}
ALIASES = {'C10orf54': 'TSPAN14', 'IL8': 'CXCL8', 'CD97': 'ADGRE5', 'FBXO18': 'FBH1', 'SEPT5-GP1BB': 'AC000093.1'}
PLATELET_PANEL = ["ITGA2B","PPBP","PF4","ITGB3","GP9","TUBB1","GP1BA","SELP"]

assert len(SIGNATURE_GENES) == 43
assert set(SIGNATURE_GENES) == set(GENCODE_V36)
print("Frozen signature:", len(SIGNATURE_GENES), "genes")
print(SIGNATURE_GENES)



# ======================================================================
# ## 3. GSE109761 discovery data
# 
# The processed normalized matrix contains 19,955 genes and 357 expression columns.
# 
# Primary discovery inclusion:
# 
# - `CTC-cluster`: 106 cells
# - `CTC-single`: 131 cells
# - `CTC-cluster-WBC` and `CTC-single-WBC`: excluded from the core comparison
# - cluster cells sharing donor + cluster ID: one biological cluster unit
# - missing cluster ID: singleton biological unit
# - donor-level cluster mean versus donor-level single-CTC mean
# - primary paired cohort: 10 donors
# 
# This explicitly prevents treating cells from the same biological cluster or donor as independent replicates.
# ======================================================================


# Cell 6
# Download processed GEO matrix and metadata.
FTP_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109761/suppl/GSE109761_processed_normalized_matrix_hs.txt.gz"
gz_path = os.path.join(DATA_DIR, "GSE109761_processed_normalized_matrix_hs.txt.gz")
matrix_path = os.path.join(DATA_DIR, "GSE109761_processed_normalized_matrix_hs.txt")

if not os.path.exists(matrix_path):
    r = requests.get(FTP_URL, timeout=120, headers={"User-Agent":"Mozilla/5.0 reviewer-reproduction"})
    r.raise_for_status()
    assert r.content[:2] == b"\x1f\x8b", "Expected gzip data from NCBI FTP mirror."
    open(gz_path,"wb").write(r.content)
    with gzip.open(gz_path,"rb") as fin, open(matrix_path,"wb") as fout:
        fout.write(fin.read())

expr = pd.read_csv(matrix_path, sep="\t", index_col=0)
print("Expression matrix:", expr.shape)
assert expr.shape == (19955, 357)

gse = GEOparse.get_GEO(geo="GSE109761", destdir=DATA_DIR, silent=True)
meta = gse.phenotype_data.copy()
print("Metadata records:", len(meta))



# Cell 7
# Classify samples using the validated `sample type` characteristic.
sample_type_col = [c for c in meta.columns if "sample type" in c.lower()][0]
donor_col_candidates = [c for c in meta.columns if "donor" in c.lower()]
cluster_id_candidates = [c for c in meta.columns if "cluster" in c.lower() and "id" in c.lower()]
print("sample type field:", sample_type_col)
print("donor fields:", donor_col_candidates)
print("cluster ID fields:", cluster_id_candidates)

# The title field is the expression-column key in this dataset.
meta2 = meta.copy()
meta2["sample_type"] = meta2[sample_type_col].astype(str).str.strip()
assert "title" in meta2.columns

# Identify donor and cluster-ID fields by exact semantic content where possible.
def choose_field(candidates, preferred_terms):
    for term in preferred_terms:
        for c in candidates:
            if term in c.lower():
                return c
    return candidates[0] if candidates else None

donor_col = choose_field(donor_col_candidates, ["donor"])
cluster_id_col = choose_field(cluster_id_candidates, ["ctc-cluster id", "cluster id"])

meta2 = meta2.set_index("title")
matched = meta2.index.intersection(expr.columns)
meta2 = meta2.loc[matched].copy()
print("Matched metadata/expression columns:", len(matched))

cluster_meta = meta2[meta2["sample_type"] == "CTC-cluster"].copy()
single_meta = meta2[meta2["sample_type"] == "CTC-single"].copy()

print("CTC-cluster cells:", len(cluster_meta))
print("CTC-single cells:", len(single_meta))
assert len(cluster_meta) == 106
assert len(single_meta) == 131



# Cell 8
# Construct biological cluster units.
cluster_meta["donor"] = cluster_meta[donor_col].astype(str)
single_meta["donor"] = single_meta[donor_col].astype(str)

def clean_cluster_id(x):
    if pd.isna(x) or str(x).strip() in {"", "nan", "None", "NA"}:
        return None
    return str(x).strip()

cluster_meta["cluster_id_clean"] = cluster_meta[cluster_id_col].map(clean_cluster_id)

# Known ID: donor + cluster ID = one biological unit.
# Missing ID: singleton unit unique to that cell.
cluster_meta["_cluster_unit"] = [
    f"{d}|{cid}" if cid is not None else f"{d}|MISSING|{i}"
    for i, (d, cid) in enumerate(zip(cluster_meta["donor"], cluster_meta["cluster_id_clean"]))
]

n_known = cluster_meta["cluster_id_clean"].notna().sum()
n_missing = cluster_meta["cluster_id_clean"].isna().sum()
n_units = cluster_meta["_cluster_unit"].nunique()

print("Known cluster-ID cells:", n_known)
print("Missing cluster-ID cells:", n_missing)
print("Biological cluster units:", n_units)

assert n_known == 84
assert n_missing == 22
assert n_units == 93



# Cell 9
# Build donor-level expression matrices.
# Cluster: average cells within each biological unit, then average units within donor.
# Single CTC: average single cells within donor.
cluster_unit_expr = expr[cluster_meta.index].T.groupby(cluster_meta["_cluster_unit"]).mean()
cluster_unit_donor = cluster_meta.groupby("_cluster_unit")["donor"].first()
cluster_unit_expr["donor"] = cluster_unit_donor
donor_cluster = cluster_unit_expr.groupby("donor").mean(numeric_only=True)

donor_single = expr[single_meta.index].T.groupby(single_meta["donor"]).mean()

paired_donors = sorted(set(donor_cluster.index) & set(donor_single.index))
print("Paired donors:", paired_donors)
print("Number paired:", len(paired_donors))
assert len(paired_donors) == 10

donor_cluster = donor_cluster.loc[paired_donors]
donor_single = donor_single.loc[paired_donors]

# Remove ERCC features.
biological_genes = [g for g in donor_cluster.columns if not str(g).upper().startswith("ERCC")]
assert len(biological_genes) == 19870

donor_diff = donor_cluster[biological_genes] - donor_single[biological_genes]
print("Genome-wide donor-difference matrix:", donor_diff.shape)
assert donor_diff.shape == (10, 19870)



# ======================================================================
# ## 4. Exact genome-wide donor-level inference
# 
# For each gene, all 1,024 sign configurations are enumerated. The exact two-sided P value is the fraction of configurations with absolute mean difference at least as large as observed.
# 
# Because only 10 donors are paired, the smallest non-zero exact two-sided P value is 0.001953125.
# ======================================================================


# Cell 11
# Exact sign-flip test for every gene.
signs = np.array(list(itertools.product([-1, 1], repeat=len(paired_donors))), dtype=float)
D = donor_diff.to_numpy(dtype=float)  # genes x donors
obs = D.mean(axis=1)
null_means = (D @ signs.T) / len(paired_donors)
p_exact = np.mean(np.abs(null_means) >= (np.abs(obs) - 1e-12), axis=1)

genome = pd.DataFrame({
    "Gene": donor_diff.columns,
    "Effect": obs,
    "Exact_P": p_exact,
})
genome["FDR_BH"] = multipletests(genome["Exact_P"], method="fdr_bh")[1]
genome["Direction"] = np.where(genome["Effect"] >= 0, "Positive", "Negative")
genome["abs_effect"] = genome["Effect"].abs()
genome = genome.sort_values(["Effect"], ascending=False).reset_index(drop=True)

print("Genes tested:", len(genome))
print("Exact P < 0.05:", (genome["Exact_P"] < 0.05).sum())
print("Exact P < 0.01:", (genome["Exact_P"] < 0.01).sum())
print("FDR < 0.05:", (genome["FDR_BH"] < 0.05).sum())
print("Minimum exact P:", genome["Exact_P"].min())
assert len(genome) == 19870
assert (genome["FDR_BH"] < 0.05).sum() == 0
assert np.isclose(genome["Exact_P"].min(), 0.001953125)
genome.to_csv(os.path.join(RESULTS_DIR,"Table_S1_Full_Genomewide_Results.csv"), index=False)



# ======================================================================
# ## 5. Reproducibility-based programme construction
# 
# Four analyses are used:
# 
# **A.** Primary 10-donor paired analysis.
# 
# **B.** Exclude donor CD.
# 
# **C.** Restrict cluster cells to known cluster identifiers only.
# 
# **D.** Exclude the predefined eight-gene platelet panel.
# 
# For each analysis, genes must be positive, rank ≤100, and have median Cohen's d_z ≥1 across A–D. This produced 44 genes. `NABP1` was then removed because it formed a singleton heterogeneous module, giving the frozen 43-gene programme.
# ======================================================================


# Cell 13
def build_donor_diff(cluster_meta_in, single_meta_in, expr_in, donors_keep=None, known_only=False, exclude_genes=None):
    cm = cluster_meta_in.copy()
    sm = single_meta_in.copy()
    if donors_keep is not None:
        cm = cm[cm["donor"].isin(donors_keep)]
        sm = sm[sm["donor"].isin(donors_keep)]
    if known_only:
        cm = cm[cm["cluster_id_clean"].notna()].copy()
    cm["_unit"] = [
        f"{d}|{cid}" if cid is not None else f"{d}|MISSING|{i}"
        for i, (d, cid) in enumerate(zip(cm["donor"], cm["cluster_id_clean"]))
    ]
    cu = expr_in[cm.index].T.groupby(cm["_unit"]).mean()
    unit_donor = cm.groupby("_unit")["donor"].first()
    dc = cu.groupby(unit_donor).mean(numeric_only=True)
    ds = expr_in[sm.index].T.groupby(sm["donor"]).mean()
    paired = sorted(set(dc.index) & set(ds.index))
    dc, ds = dc.loc[paired], ds.loc[paired]
    genes = [g for g in dc.columns if not str(g).upper().startswith("ERCC")]
    if exclude_genes:
        genes = [g for g in genes if g not in set(exclude_genes)]
    return dc[genes] - ds[genes]

def exact_gene_table(diff):
    s = np.array(list(itertools.product([-1,1], repeat=diff.shape[0])), dtype=float)
    x = diff.to_numpy(dtype=float)
    effect = x.mean(axis=0)
    null = (s @ x.T) / diff.shape[1]
    p = np.mean(np.abs(null) >= (np.abs(effect) - 1e-12), axis=0)
    dz = effect / x.std(axis=0, ddof=1)
    out = pd.DataFrame({"Gene":diff.columns,"Effect":effect,"Exact_P":p,"Cohen_dz":dz})
    out["Rank"] = out["Effect"].rank(method="min", ascending=False).astype(int)
    return out

diff_A = build_donor_diff(cluster_meta, single_meta, expr)
diff_B = build_donor_diff(cluster_meta, single_meta, expr, donors_keep=[d for d in paired_donors if d != "CD"])
diff_C = build_donor_diff(cluster_meta, single_meta, expr, known_only=True)
diff_D = build_donor_diff(cluster_meta, single_meta, expr, exclude_genes=PLATELET_PANEL)

res = {"A": exact_gene_table(diff_A), "B": exact_gene_table(diff_B), "C": exact_gene_table(diff_C), "D": exact_gene_table(diff_D)}
for k,v in res.items():
    print(k, "donors =", len({"A":diff_A,"B":diff_B,"C":diff_C,"D":diff_D}[k].index), "genes =", len(v))



# Cell 14
# Merge analyses and apply the prespecified robust-set rule.
common = set.intersection(*[set(v["Gene"]) for v in res.values()])
rank_cols = pd.DataFrame({"Gene": sorted(common)})
for k in ["A","B","C","D"]:
    t = res[k].set_index("Gene")
    rank_cols[k+"_rank"] = t.loc[rank_cols["Gene"], "Rank"].to_numpy()
    rank_cols[k+"_effect"] = t.loc[rank_cols["Gene"], "Effect"].to_numpy()
    rank_cols[k+"_dz"] = t.loc[rank_cols["Gene"], "Cohen_dz"].to_numpy()

rank_cols["median_rank"] = rank_cols[[f"{k}_rank" for k in "ABCD"]].median(axis=1)
rank_cols["median_dz"] = rank_cols[[f"{k}_dz" for k in "ABCD"]].median(axis=1)
rank_cols["positive_all"] = (rank_cols[[f"{k}_effect" for k in "ABCD"]] > 0).all(axis=1)
rank_cols["rank_le_100_all"] = (rank_cols[[f"{k}_rank" for k in "ABCD"]] <= 100).all(axis=1)

robust44 = rank_cols[
    rank_cols["positive_all"] &
    rank_cols["rank_le_100_all"] &
    (rank_cols["median_dz"] >= 1)
].copy().sort_values("median_rank")

print("Robust set:", len(robust44))
assert len(robust44) == 44

robust44.to_csv(os.path.join(RESULTS_DIR,"robust_44_genes.csv"), index=False)

frozen_genes = [g for g in robust44["Gene"] if g != "NABP1"]
assert len(frozen_genes) == 43
assert frozen_genes == SIGNATURE_GENES
print("Frozen CTC_Cluster_43 reproduced exactly.")
pd.DataFrame({"Order":range(1,44),"Gene":frozen_genes,"Alias":[ALIASES.get(g,"") for g in frozen_genes]}).to_csv(
    os.path.join(RESULTS_DIR,"Table_2_Frozen_CTC_Cluster_43.csv"), index=False)



# ======================================================================
# ## 6. Discovery programme score
# 
# For discovery, gene-wise z-scores are calculated across the 20 donor-condition observations using the **population standard deviation (`ddof=0`)**. The programme score is the equal-weight mean of the 43 z-scores.
# 
# The paired programme effect is cluster score minus single-CTC score.
# ======================================================================


# Cell 16
# 20 donor-condition observations, frozen 43 genes.
conditions = pd.concat([
    donor_single.loc[paired_donors, frozen_genes].assign(Condition="Single CTC"),
    donor_cluster.loc[paired_donors, frozen_genes].assign(Condition="CTC cluster")
], axis=0)
# Rebuild with explicit unique row index.
conditions = pd.concat([
    donor_single.loc[paired_donors, frozen_genes].assign(Donor=paired_donors, Condition="Single CTC"),
    donor_cluster.loc[paired_donors, frozen_genes].assign(Donor=paired_donors, Condition="CTC cluster")
], ignore_index=True)

Z = (conditions[frozen_genes] - conditions[frozen_genes].mean(axis=0)) / conditions[frozen_genes].std(axis=0, ddof=0)
programme_scores = Z.mean(axis=1)
conditions["CTC_Cluster_43"] = programme_scores

single_scores = conditions.loc[conditions["Condition"]=="Single CTC"].set_index("Donor")["CTC_Cluster_43"]
cluster_scores = conditions.loc[conditions["Condition"]=="CTC cluster"].set_index("Donor")["CTC_Cluster_43"]
programme_diff = cluster_scores.loc[paired_donors] - single_scores.loc[paired_donors]

obs_score = programme_diff.mean()
exact_programme_p = np.mean(
    np.abs((programme_diff.to_numpy()[None,:] * signs).mean(axis=1)) >= abs(obs_score)-1e-12
)
dz = programme_diff.mean() / programme_diff.std(ddof=1)

print("Programme mean effect:", obs_score)
print("Median:", programme_diff.median())
print("SD:", programme_diff.std(ddof=1))
print("Cohen dz:", dz)
print("Positive donors:", (programme_diff>0).sum(), "/", len(programme_diff))
print("Exact sign-flip P:", exact_programme_p)

assert np.isclose(obs_score, 1.151134569989992, atol=1e-6)
assert np.isclose(dz, 1.774589, atol=1e-5)
assert exact_programme_p == 0.001953125

pd.DataFrame({"Donor":paired_donors,"Programme_difference":programme_diff.values}).to_csv(
    os.path.join(RESULTS_DIR,"Table_S2_Programme_donor_differences.csv"), index=False)



# Cell 17
# Sensitivity analyses on the frozen 43-gene programme.
def score_effect_for_diff(diff, genes=frozen_genes):
    # Score from the same 20-condition discovery observations after applying the requested donor/gene restriction.
    donors = diff.index.tolist()
    c = donor_cluster.loc[donors, genes]
    s = donor_single.loc[donors, genes]
    both = pd.concat([s.assign(Condition="Single"), c.assign(Condition="Cluster")])
    z = (both[genes] - both[genes].mean()) / both[genes].std(ddof=0)
    scores = z.mean(axis=1)
    n = len(donors)
    ds = scores.iloc[:n].to_numpy()
    dc = scores.iloc[n:].to_numpy()
    d = dc-ds
    ss = np.array(list(itertools.product([-1,1], repeat=n)))
    p = np.mean(np.abs((d[None,:]*ss).mean(axis=1)) >= abs(d.mean())-1e-12)
    return pd.Series(d,index=donors), float(d.mean()), float(d.median()), float(d.std(ddof=1)), float(d.mean()/d.std(ddof=1)), float(p)

# A
dA = programme_diff
# B
b_diff = diff_B[frozen_genes]
dB, meanB, medB, sdB, dzB, pB = score_effect_for_diff(b_diff)
# C
c_diff = diff_C[frozen_genes]
dC, meanC, medC, sdC, dzC, pC = score_effect_for_diff(c_diff)
# D is identical because the frozen 43 has zero overlap with the platelet panel.
dD, meanD, medD, sdD, dzD, pD = score_effect_for_diff(diff_D[frozen_genes])

sens = pd.DataFrame([
    ["Primary",len(dA),dA.mean(),dA.median(),dA.std(ddof=1),(dA>0).sum(),(dA<0).sum(),exact_programme_p,dz],
    ["Exclude CD",len(dB),meanB,medB,sdB,(dB>0).sum(),(dB<0).sum(),pB,dzB],
    ["Known cluster IDs",len(dC),meanC,medC,sdC,(dC>0).sum(),(dC<0).sum(),pC,dzC],
    ["Platelet sensitivity",len(dD),meanD,medD,sdD,(dD>0).sum(),(dD<0).sum(),pD,dzD],
], columns=["Analysis","Donors","Mean","Median","SD","Positive","Negative","Exact_P","Cohen_dz"])
print(sens.to_string(index=False))
sens.to_csv(os.path.join(RESULTS_DIR,"Table_3_Sensitivity_Analyses.csv"),index=False)

assert np.isclose(meanB,1.180084,atol=1e-5)
assert np.isclose(dzB,1.933329,atol=1e-5)
assert pB == 0.00390625



# ======================================================================
# ## 7. Frozen TCGA external cohort
# 
# The external cohort is the frozen 35-patient GDC selection used in the manuscript. The manifest records the exact patient, primary-tumour sample barcode and GDC STAR-Counts file UUID.
# 
# This avoids re-discovering a different cohort when the GDC catalogue changes.
# ======================================================================


# Cell 19
# Frozen 35-sample GDC manifest.
manifest = pd.DataFrame([('TCGA-06-0211', 'TCGA-06-0211-01A', '54bde4b3-b158-4836-a593-68b414076832', 'TCGA-GBM', '9aaeecd'), ('TCGA-44-2656', 'TCGA-44-2656-01A', '77819b24-ea90-427b-8a99-c9a90417c903', 'TCGA-LUAD', 'eecca6e2e475aab335bd7c365025ca0d0daa144b'), ('TCGA-44-3917', 'TCGA-44-3917-01A', 'b0b5b05d-f782-449d-94cf-8fd74a4bd30c', 'TCGA-LUAD', 'eecca6e2e475aab335bd7c365025ca0d0daa144b'), ('TCGA-44-3918', 'TCGA-44-3918-01A', 'cc1174cd-f902-46f0-9e69-b285b2b1317c', 'TCGA-LUAD', 'eecca6e2e475aab335bd7c365025ca0d0daa144b'), ('TCGA-44-4112', 'TCGA-44-4112-01A', '9af8b2eb-3dbb-4688-a5fd-15abc0460f57', 'TCGA-LUAD', '5a02c12930eb23436322a7494047b907b97570c5'), ('TCGA-44-6146', 'TCGA-44-6146-01A', '074b0792-df3c-4b59-9f50-793bc14bcb81', 'TCGA-LUAD', 'eecca6e2e475aab335bd7c365025ca0d0daa144b'), ('TCGA-44-6147', 'TCGA-44-6147-01A', '9a64a124-9787-4378-ba01-b79aa183f25f', 'TCGA-LUAD', 'eecca6e2e475aab335bd7c365025ca0d0daa144b'), ('TCGA-A6-2672', 'TCGA-A6-2672-01A', '9a6cef0c-1db2-4ec8-9918-b261ce6f5a86', 'TCGA-COAD', '470801dfbc085958058a5b50518e4f3a5e4721c3'), ('TCGA-A6-2674', 'TCGA-A6-2674-01A', 'd84aa7f4-7210-4121-9e4a-e6000d30ff90', 'TCGA-COAD', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-A6-2677', 'TCGA-A6-2677-01A', '78e565d9-4c8c-44b1-8d58-473381e92c52', 'TCGA-COAD', '5a02c12930eb23436322a7494047b907b97570c5'), ('TCGA-A6-2684', 'TCGA-A6-2684-01A', '97dbbc5d-67e7-43a0-bc8f-d3425e0f3ce4', 'TCGA-COAD', '470801dfbc085958058a5b50518e4f3a5e4721c3'), ('TCGA-A6-3809', 'TCGA-A6-3809-01A', 'c741f3fd-133d-470c-8021-cd3432c61ad0', 'TCGA-COAD', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-A6-3810', 'TCGA-A6-3810-01A', 'cd10bd5e-b75b-4d47-a103-9c259dcf1a19', 'TCGA-COAD', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-A6-5656', 'TCGA-A6-5656-01A', '08872f38-d20f-4610-898c-aff8fe15bfe2', 'TCGA-COAD', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-A6-5659', 'TCGA-A6-5659-01A', 'f01575cb-c8fa-4a21-9e61-768d25cabc58', 'TCGA-COAD', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-A6-5661', 'TCGA-A6-5661-01A', '53494104-fea8-48d7-881e-03c760435509', 'TCGA-COAD', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-A6-5665', 'TCGA-A6-5665-01A', '824e0db3-9c7a-4cc7-ae84-6319a4e53e86', 'TCGA-COAD', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-A6-6650', 'TCGA-A6-6650-01A', 'cee77d7f-4428-4df1-8532-7782481459a1', 'TCGA-COAD', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-A6-6780', 'TCGA-A6-6780-01A', 'e065cd2d-aecb-4c38-ae7e-70ec3da7430f', 'TCGA-COAD', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-A6-6781', 'TCGA-A6-6781-01A', 'd2bb361b-c78c-482f-b4d4-bd77ac7c0f37', 'TCGA-COAD', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-A7-A0DC', 'TCGA-A7-A0DC-01A', '0e1999e0-9c02-4804-88b6-60b46d100739', 'TCGA-BRCA', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-A7-A13E', 'TCGA-A7-A13E-01A', '0a511373-8fd3-433a-a5e1-877530d6a239', 'TCGA-BRCA', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-B2-3923', 'TCGA-B2-3923-01A', '4134fa64-6365-4a43-b990-48342be44405', 'TCGA-KIRC', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-B2-5633', 'TCGA-B2-5633-01A', 'bdb74506-8766-42e1-922e-e03da4637a70', 'TCGA-KIRC', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-BK-A0CA', 'TCGA-BK-A0CA-01A', '220a61db-a613-42ef-8c4b-7c68bd310d36', 'TCGA-UCEC', '470801dfbc085958058a5b50518e4f3a5e4721c3'), ('TCGA-BK-A0CC', 'TCGA-BK-A0CC-01A', '63d0b3c6-f57d-4bbb-be7a-162afa1fc7a2', 'TCGA-UCEC', '470801dfbc085958058a5b50518e4f3a5e4721c3'), ('TCGA-BK-A139', 'TCGA-BK-A139-01A', 'a2fb0c94-9c06-4d3c-a392-6f3e1e1517ec', 'TCGA-UCEC', '470801dfbc085958058a5b50518e4f3a5e4721c3'), ('TCGA-BK-A26L', 'TCGA-BK-A26L-01A', 'a8569b16-beb6-4d8c-ba39-ffbfd1769663', 'TCGA-UCEC', '470801dfbc085958058a5b50518e4f3a5e4721c3'), ('TCGA-BL-A0C8', 'TCGA-BL-A0C8-01A', '49a6de74-9360-4347-a594-d64d5aa5ce9e', 'TCGA-BLCA', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-BL-A13I', 'TCGA-BL-A13I-01A', 'f36b3cb3-b632-43c6-a6e3-c325c9e7026f', 'TCGA-BLCA', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-BL-A13J', 'TCGA-BL-A13J-01A', '0d50ca90-4dfc-4dc8-b353-0d2af65649d6', 'TCGA-BLCA', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-HC-7740', 'TCGA-HC-7740-01A', 'a01b7ba6-acd9-40cb-8359-c6565714cdab', 'TCGA-PRAD', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-HC-8258', 'TCGA-HC-8258-01A', '107b0ba8-f479-42dc-a916-385b845a162c', 'TCGA-PRAD', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-HC-8261', 'TCGA-HC-8261-01A', '2b5963dc-6bbd-4541-a6ed-19fcc046fc10', 'TCGA-PRAD', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f'), ('TCGA-HC-8265', 'TCGA-HC-8265-01A', 'cc7216dc-6b4a-4e8e-a69e-fa32c570c3a3', 'TCGA-PRAD', '5d8c131bbff59fb0c969217fc1d44e6d1503cd1f')], columns=["patient","sample_barcode","file_id","project_id","workflow_version"])
manifest.to_csv(os.path.join(RESULTS_DIR,"tcga_frozen_manifest.csv"),index=False)
print(manifest.to_string(index=False))
assert len(manifest)==35 and manifest["patient"].nunique()==35 and manifest["sample_barcode"].nunique()==35



# Cell 20
# Download and parse the 35 frozen GDC STAR-Counts files.
GDC_DATA = "https://api.gdc.cancer.gov/data/"

def download_gdc_file(file_id, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return
    r = requests.get(GDC_DATA + file_id, timeout=180)
    r.raise_for_status()
    open(dest,"wb").write(r.content)

def parse_star_counts(path):
    # STAR Counts files are tab-delimited; comment lines begin with ##.
    df = pd.read_csv(path, sep="\t", comment="#")
    cols = {c.lower():c for c in df.columns}
    tpm = next((cols[c] for c in cols if c == "tpm_unstranded"), None)
    gid = next((c for c in df.columns if c.lower() == "gene_id"), df.columns[0])
    gname = next((c for c in df.columns if c.lower() == "gene_name"), None)
    if tpm is None:
        raise ValueError(f"No tpm_unstranded column in {path}: {list(df.columns)}")
    out = df[[gid] + ([gname] if gname else []) + [tpm]].copy()
    out = out.rename(columns={gid:"gene_id", tpm:"tpm_unstranded"})
    if gname: out = out.rename(columns={gname:"gene_name"})
    return out

sample_gene = {}
for _, r in manifest.iterrows():
    fp = os.path.join(DATA_DIR, f"{r.file_id}.tsv")
    download_gdc_file(r.file_id, fp)
    tab = parse_star_counts(fp)
    sample_gene[r.sample_barcode] = tab

print("Downloaded/parsed:", len(sample_gene), "GDC files")
assert len(sample_gene)==35



# Cell 21
# Resolve the frozen 43 genes using GENCODE v36 Ensembl IDs.
def gene_value(tab, gene, ensembl_id):
    ids = tab["gene_id"].astype(str).str.replace(r"\.\d+$","",regex=True)
    hit = tab.loc[ids == ensembl_id]
    if len(hit) == 1:
        return float(hit["tpm_unstranded"].iloc[0])
    # Fallback to gene_name for documented aliases only.
    if "gene_name" in tab.columns:
        candidates = [gene] + ([ALIASES[gene]] if gene in ALIASES else [])
        hit = tab[tab["gene_name"].astype(str).isin(candidates)]
        if len(hit) >= 1:
            return float(hit["tpm_unstranded"].iloc[0])
    raise KeyError(f"{gene} / {ensembl_id} not resolved")

tcga_expr = pd.DataFrame(index=manifest["sample_barcode"], columns=frozen_genes, dtype=float)
for _, r in manifest.iterrows():
    tab = sample_gene[r.sample_barcode]
    for g in frozen_genes:
        tcga_expr.loc[r.sample_barcode,g] = gene_value(tab,g,GENCODE_V36[g])

print("TCGA expression:", tcga_expr.shape)
print("Missing values:", int(tcga_expr.isna().sum().sum()))
assert tcga_expr.shape == (35,43)
assert tcga_expr.isna().sum().sum() == 0



# Cell 22
# External programme score:
# log2(TPM+1) -> gene-wise z-score across 35 samples using ddof=1 -> equal-weight mean.
log_expr = np.log2(tcga_expr + 1.0)
tcga_z = (log_expr - log_expr.mean(axis=0)) / log_expr.std(axis=0, ddof=1)
tcga_score = tcga_z.mean(axis=1)
tcga_score.name = "CTC_Cluster_43"

print("TCGA score mean:", tcga_score.mean())
print("TCGA score median:", tcga_score.median())
print("TCGA score SD:", tcga_score.std(ddof=1))
print("TCGA score range:", tcga_score.min(), tcga_score.max())
assert np.isclose(tcga_score.mean(),0,atol=1e-12)
assert np.isclose(tcga_score.median(),0.1689783578,atol=1e-6)
assert np.isclose(tcga_score.std(ddof=1),0.5355102148,atol=1e-6)

tcga_score.to_csv(os.path.join(RESULTS_DIR,"tcga_CTC_Cluster_43_scores.csv"),header=True)



# ======================================================================
# ## 8. TCGA-CDR PFI and survival analysis
# 
# The clinical source is TCGA-CDR. The primary analysis is continuous Cox regression. The reported hazard ratio is per 1-s.d. increase in the TCGA programme score.
# 
# Median split is secondary/descriptive and uses the prespecified cohort median.
# ======================================================================


# Cell 24
# TCGA-CDR workbook.
TCGA_CDR_UUID = "1b5f413e-a8d1-4d10-92eb-7c4ae739ed81"
cdr_path = os.path.join(DATA_DIR,"TCGA-CDR-SupplementalTableS1.xlsx")
if not os.path.exists(cdr_path):
    r = requests.get(GDC_DATA + TCGA_CDR_UUID, timeout=180)
    r.raise_for_status()
    open(cdr_path,"wb").write(r.content)

cdr = pd.read_excel(cdr_path)
print("TCGA-CDR shape:", cdr.shape)
print("Columns:", list(cdr.columns))
assert "PFI" in cdr.columns and "PFI.time" in cdr.columns



# Cell 25
# Match patient IDs and fit the primary continuous Cox model.
tcga = manifest[["patient","sample_barcode"]].copy()
tcga["CTC_Cluster_43"] = tcga["sample_barcode"].map(tcga_score)

# TCGA-CDR patient identifier is usually `bcr_patient_barcode`.
patient_col = "bcr_patient_barcode" if "bcr_patient_barcode" in cdr.columns else next(
    c for c in cdr.columns if "patient" in c.lower() and "barcode" in c.lower()
)
clinical = cdr[[patient_col,"PFI","PFI.time"]].copy().rename(columns={patient_col:"patient"})
clinical["PFI"] = pd.to_numeric(clinical["PFI"], errors="coerce")
clinical["PFI.time"] = pd.to_numeric(clinical["PFI.time"], errors="coerce")

tcga_pfi = tcga.merge(clinical,on="patient",how="inner").dropna(subset=["CTC_Cluster_43","PFI","PFI.time"])
print("Matched TCGA patients:", len(tcga_pfi))
print("PFI events:", int(tcga_pfi["PFI"].sum()))
assert len(tcga_pfi)==35
assert int(tcga_pfi["PFI"].sum())==14

# HR per raw score unit, then convert to HR per 1-SD score.
cph = CoxPHFitter()
cph.fit(tcga_pfi[["PFI.time","PFI","CTC_Cluster_43"]], duration_col="PFI.time", event_col="PFI")
beta = cph.params_["CTC_Cluster_43"]
hr_raw = np.exp(beta)
ci_raw = np.exp(cph.confidence_intervals_.loc["CTC_Cluster_43"].to_numpy())
p = cph.summary.loc["CTC_Cluster_43","p"]
score_sd = tcga_pfi["CTC_Cluster_43"].std(ddof=1)
hr_1sd = np.exp(beta * score_sd)
ci_1sd = np.exp(cph.confidence_intervals_.loc["CTC_Cluster_43"].to_numpy() * score_sd)
cindex = cph.concordance_index_

print("HR raw score:", hr_raw)
print("HR per 1 SD:", hr_1sd)
print("95% CI per 1 SD:", ci_1sd)
print("P:", p)
print("C-index:", cindex)

assert np.isclose(hr_1sd,1.327032,atol=1e-5)
assert np.isclose(p,0.378766,atol=1e-5)
assert np.isclose(cindex,0.552469,atol=1e-5)



# Cell 26
# Secondary descriptive median split.
cutoff = tcga_pfi["CTC_Cluster_43"].median()
tcga_pfi["group"] = np.where(tcga_pfi["CTC_Cluster_43"] >= cutoff, "High","Low")
high = tcga_pfi[tcga_pfi["group"]=="High"]
low = tcga_pfi[tcga_pfi["group"]=="Low"]

km_high = KaplanMeierFitter().fit(high["PFI.time"], high["PFI"], label="High")
km_low = KaplanMeierFitter().fit(low["PFI.time"], low["PFI"], label="Low")

# lifelines logrank test
from lifelines.statistics import logrank_test
lr = logrank_test(high["PFI.time"], low["PFI.time"], event_observed_A=high["PFI"], event_observed_B=low["PFI"])

print("Median cutoff:", cutoff)
print("High:", len(high), "events", int(high.PFI.sum()))
print("Low:", len(low), "events", int(low.PFI.sum()))
print("Log-rank P:", lr.p_value)

assert np.isclose(cutoff,0.1689783578,atol=1e-6)
assert len(high)==18 and int(high.PFI.sum())==8
assert len(low)==17 and int(low.PFI.sum())==6
assert np.isclose(lr.p_value,0.660190,atol=1e-5)

tcga_pfi[["patient","CTC_Cluster_43","PFI","PFI.time"]].to_csv(
    os.path.join(RESULTS_DIR,"Table_S4_TCGA_Patient_Level_PFI.csv"),index=False)



# ======================================================================
# ## 9. Reproducibility checkpoint
# 
# The following assertions are the minimum reviewer checkpoints. If any fail, the run is not equivalent to the manuscript analysis.
# ======================================================================


# Cell 28
checkpoints = {
    "GSE genes": expr.shape[0] == 19955,
    "GSE expression samples": expr.shape[1] == 357,
    "CTC-cluster cells": len(cluster_meta) == 106,
    "CTC-single cells": len(single_meta) == 131,
    "Biological cluster units": cluster_meta["_cluster_unit"].nunique() == 93,
    "Paired donors": len(paired_donors) == 10,
    "Genome-wide genes": len(genome) == 19870,
    "Exact configurations": len(signs) == 1024,
    "No genome-wide FDR<0.05 genes": int((genome.FDR_BH < 0.05).sum()) == 0,
    "Frozen signature": len(frozen_genes) == 43,
    "TCGA patients": len(tcga_pfi) == 35,
    "TCGA events": int(tcga_pfi.PFI.sum()) == 14,
    "TCGA genes complete": tcga_expr.shape == (35,43),
    "Primary programme P": np.isclose(exact_programme_p,0.001953125),
    "Primary programme effect": np.isclose(obs_score,1.151134569989992,atol=1e-6),
    "TCGA HR/SD": np.isclose(hr_1sd,1.327032,atol=1e-5),
    "TCGA P": np.isclose(p,0.378766,atol=1e-5),
    "TCGA log-rank P": np.isclose(lr.p_value,0.660190,atol=1e-5),
}
for k,v in checkpoints.items():
    print(f"{'PASS' if v else 'FAIL'}: {k}")
assert all(checkpoints.values())
print("\nALL MANUSCRIPT REPRODUCIBILITY CHECKPOINTS PASSED.")



# ======================================================================
# ## 10. Outputs
# 
# The notebook writes reproducible CSV outputs to `results/`, including:
# 
# - `Table_S1_Full_Genomewide_Results.csv`
# - `Table_2_Frozen_CTC_Cluster_43.csv`
# - `Table_3_Sensitivity_Analyses.csv`
# - `Table_S4_TCGA_Patient_Level_PFI.csv`
# - `tcga_frozen_manifest.csv`
# - `tcga_CTC_Cluster_43_scores.csv`
# 
# The supplementary enrichment/module characterisation is descriptive and is not used as independent differential-expression confirmation. The central inferential results are the donor-level exact analysis, the frozen programme, and the prespecified external TCGA PFI test.
# ======================================================================

