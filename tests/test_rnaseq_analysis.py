
from biorag.rnaseq_analysis import RNASeqAnalysis


x = RNASeqAnalysis("tests/human_gene_v2.2.h5")

#adata = x.create_anndata_from_samples(["GSM1000981, GSM1002540","GSM1009637", "GSM1158460"])

element_list = ["GSM1000981, GSM1002540","GSM1009637", "GSM1158460"]
geneset = ["IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6"]
res = x.perform_enrichment_on_samples_batched(element_list, geneset)

print(res)
print(res.columns)

# metadata = a4.data.series("tests/human_gene_v2.2.h5", "GSE64016")
# print(metadata)

def test_rnaseq_analysis_adata_samples():
    adata = x.create_anndata_from_samples(["GSM1000981, GSM1002540","GSM1009637", "GSM1158460"]) # some samples have no counts
    assert adata.obs.shape[0] == 2

def test_rnaseq_analysis_adata_samples():
    adata = x.create_anndata_from_series("GSE64016")
    assert adata.obs.shape[0] == 460

def test_list_to_dc_geneset():
    df = x.list_to_dc_geneset(["IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6"])
    assert df.shape == (6, 2)

def test_perform_enrichment_on_samples_batched():
    element_list = ["GSM1000981, GSM1002540","GSM1009637", "GSM1158460"]
    geneset = ["IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6"]
    res = x.perform_enrichment_on_samples_batched(element_list, geneset)
    assert res.shape == (2, 8)