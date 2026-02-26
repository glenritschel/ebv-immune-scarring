from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pandas as pd

# counts: genes x samples DataFrame (raw integer counts)
# metadata: samples x covariates DataFrame with 'cluster' column

dds = DeseqDataSet(
    counts=counts_df,        # rows=genes, cols=samples
    metadata=metadata_df,    # rows=samples
    design_factors="cluster"
)
dds.deseq2()

stats = DeseqStats(dds, contrast=["cluster", "0", "1"])
stats.summary()

results = stats.results_df  # padj, log2FoldChange, etc.
results.sort_values("padj").head(20)
