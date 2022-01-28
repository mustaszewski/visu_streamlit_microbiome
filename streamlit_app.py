from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
import streamlit as st

# load abundance data
# 
df_16s_abu = pd.read_excel("data/original/GSE113690_Autism_16S_rRNA_OTU_assignment_and_abundance.xls", index_col=None)

df_16s_abu.rename(columns={"OTU":"OTU_ID"}, inplace=True)

# define labels of taxonomic ranks
rank_labels = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]

# split string in taxonomy columns into multiple columns, each representing one taxonomic level
taxonomic_ranks = df_16s_abu["taxonomy"].str.split(";", expand=True)
taxonomic_ranks.columns = rank_labels # relabel taxonomy level columns
# cleanup columns by stripping the _k__, _c__ etc. suffixes
taxonomic_ranks = taxonomic_ranks.replace("^_?[^_]__", "", regex=True)

taxonomy = pd.concat([df_16s_abu.iloc[:, 0:2], taxonomic_ranks], axis=1)

metadata_colnames = ["sample_id", "stage", "gender", "age", "has_16S_RNA_sequencing", "has_metagenomic_sequencing", "has_metabonomic_analysis", "constipation"]

# load metadata from Excel file into DataFrame
df_meta = pd.read_excel("data/original/Table_S1_Sample_information.xlsx", names=metadata_colnames)

# cleanup some string values by replacement
df_meta = df_meta.replace("Autism", "ASD")
df_meta = df_meta.replace("Yes", "yes")
df_meta = df_meta.replace("No", "no")
df_meta.drop(["has_16S_RNA_sequencing", "has_metagenomic_sequencing", "has_metabonomic_analysis"], axis=1, inplace=True)


# merge abundance data with metadata
ids_otu = df_16s_abu["OTU_ID"]
ids_samples = df_16s_abu.columns[2:]

df_16s_abu_transposed = df_16s_abu.drop(["taxonomy"], axis=1).T
df_16s_abu_transposed.columns = ids_otu
df_16s_abu_transposed = df_16s_abu_transposed.iloc[1: , :] # drop first line

df_16s_abu_transposed.index.rename("Sample_ID", inplace=True)
df_16s_abu_transposed.columns.rename("OTU_ID", inplace=True)

df_16s_abu_transposed = df_16s_abu_transposed.apply(pd.to_numeric) # convert all columns to integer

df_16s_abu_transposed

df_full = df_meta.merge(df_16s_abu_transposed, how="inner", left_on="sample_id", right_on=df_16s_abu_transposed.index)

# convert into long format

# melt wide format into long format
df_long = df_full[["sample_id"] +  ids_otu.to_list()].melt('sample_id', var_name='OTU_ID', value_name='OTU_Abundance')

# perform left-join to add  sample metadata
df_long_full = df_long.merge(df_full.iloc[:, 0:5], how="left", on="sample_id")
# perform left join to add OTU taxonomy information
df_long_full = df_long_full.merge(taxonomy, how="left", on="OTU_ID")


# phylum abundance data

phylum_abundance_per_sample = df_long_full.groupby(["sample_id", "phylum"])["OTU_Abundance"].sum()
phylum_abundance_per_sample = pd.DataFrame(phylum_abundance_per_sample).reset_index()
phylum_abundance_per_sample = phylum_abundance_per_sample.merge(df_meta, on="sample_id", how="left")
phylum_abundance_per_sample


# compute pca

pca_3d_model = PCA(n_components=3)
pca_3d = pca_3d_model.fit_transform(df_full[ids_otu])

# append PCA coordinates to data frame as separate columns
df_full[["pca1", "pca2", "pca3"]] = pca_3d
df_full.head()

# interactive compound pca with 2 histograms
def draw_pca():
    interval = alt.selection_interval()
    
    scatter_base = alt.Chart(df_full).mark_point().encode(
        x=alt.X('pca1', title='PCA #1'),
        color=alt.condition(interval, 'stage:N', alt.value('lightgray')),
        tooltip=['sample_id']
        ).add_selection(
            interval
            )

    scatter_pca = scatter_base.encode(y=alt.Y('pca2', title='PCA #2')) | scatter_base.encode(y=alt.Y('pca3', title='PCA #3'))
    
    histo =alt.Chart(phylum_abundance_per_sample).mark_bar().encode(
        y='stage:N',
        x=alt.X('sum(OTU_Abundance)', stack='normalize'),#'mean(OTU_Abundance):Q', #    x=alt.X('mean(OTU_Abundance)', stack='zero'),

        color='phylum:N',
        tooltip=['phylum:N', alt.Tooltip('sum(OTU_Abundance):Q', format=',.1f')],
        #tooltip=['mean(OTU_Abundance):Q', 'phylum:N'],
    ).transform_lookup(
        lookup='sample_id',
        from_=alt.LookupData(data=df_full, key='sample_id',
        fields=['pca1', 'pca2'])
    ).transform_filter(interval)

    compound_pcas_histo = alt.vconcat(
        scatter_pca, histo,
        center=True
    ).resolve_scale(
        color='independent'
    )

    return compound_pcas_histo



### run app code

st.write(draw_pca())