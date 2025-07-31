# run_bb_gsea_CC.py

import pandas as pd
import mygene
import gseapy as gp
import os

# -------- Load and clean expression data --------
bb_expression = pd.read_csv("/scratch/users/k24057496/GSEA/merged_cellular_BB.txt", sep='\t', quotechar='"')
bb_expression.columns = bb_expression.columns.str.strip()
bb_expression.set_index('gene', inplace=True)
bb_expression = bb_expression.rename_axis(None).transpose()
bb_expression = bb_expression.loc[:, bb_expression.sum(axis=0) >= 10]

# -------- Normalize expression --------
def deseq2_norm(df):
    size_factors = df.sum(axis=1) / df.sum(axis=1).median()
    norm_df = df.div(size_factors, axis=0)
    return norm_df, size_factors

bb_expression_normalised, _ = deseq2_norm(bb_expression)
bb_expression_normalised = bb_expression_normalised.transpose()

# -------- Load metadata and prepare CLS file --------
bb_metadata = pd.read_csv('/scratch/users/k24057496/GSEA/phenotypeStarting_BB.txt', sep='\t')
bb_metadata.columns = bb_metadata.columns.str.strip()
bb_metadata['Status'] = bb_metadata['Status'].astype(str)
bb_metadata['label'] = bb_metadata['Status'] + "_" + bb_metadata['Individual'].astype(str)
bb_metadata = bb_metadata.sort_values(by='label')

bb_samples_order = bb_metadata['Individual'].tolist()
indiv_2_label = dict(zip(bb_metadata['Individual'], bb_metadata['label']))

# Ensure matching samples
missing_samples = [s for s in bb_samples_order if s not in bb_expression_normalised.columns]
if missing_samples:
    raise ValueError(f"Missing samples in expression: {missing_samples}")

bb_expression_normalised = bb_expression_normalised[bb_samples_order]
bb_expression_normalised.rename(columns=indiv_2_label, inplace=True)

# Save CLS file
cls_path = "/scratch/users/k24057496/GSEA/bb_labels.cls"
with open(cls_path, "w") as f:
    f.write(f"{len(bb_metadata)} {len(set(bb_metadata['Status']))} 1\n")
    f.write(f"# {' '.join(sorted(set(bb_metadata['Status'])))}\n")
    f.write(' '.join(bb_metadata['Status']) + "\n")

# -------- Convert Ensembl IDs to gene symbols --------
mg = mygene.MyGeneInfo()
query = mg.querymany(bb_expression_normalised.index.tolist(), scopes='ensembl.gene', fields='symbol', species='human')
id_map = {entry['query']: entry.get('symbol') for entry in query if 'symbol' in entry}
bb_expression_normalised.rename(index=id_map, inplace=True)
bb_expression_normalised = bb_expression_normalised[~bb_expression_normalised.index.duplicated(keep='first')]
bb_expression_normalised.dropna(inplace=True)

# -------- Run GSEA for GO Cellular Component --------
outdir = "/scratch/users/k24057496/GSEA/3_GSEA_results_CC/"
os.makedirs(outdir, exist_ok=True)

bb_gsea_res = gp.gsea(
    data=bb_expression_normalised,
    gene_sets=["GO_Cellular_Component_2025"],
    cls=cls_path,
    sample_norm_method='rank',
    permutation_type='phenotype',
    permutation_num=1000,
    outdir=outdir,
    threads=8,
    min_size=5,
    max_size=1000,
    seed=0
)

# -------- Save results --------
results_df = bb_gsea_res.res2d
results_df['Term'] = results_df['Term'].str.replace("GO_Cellular_Component_2025__", "", regex=False)
results_df.to_csv("/scratch/users/k24057496/GSEA/gsea_results_CC_summary.csv", index=False)
