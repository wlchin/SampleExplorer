
from biorag.rnaseq_analysis import RNASeqAnalysis
import pytest

x = RNASeqAnalysis("tests/human_gene_v2.2.h5")
#res = x.perform_enrichment_on_samples_batched(element_list, geneset)
#adata = x.create_anndata_from_samples(["GSM1000981, GSM1002540","GSM1009637", "GSM1158460"])

element_list = ["GSM1000981, GSM1002540","GSM1009637", "GSM1158460"]
geneset = ["IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6"]


# print(res)
# print(res.columns)

# metadata = a4.data.series("tests/human_gene_v2.2.h5", "GSE64016")
# print(metadata)

@pytest.mark.xfail
def test_rnaseq_analysis_adata_samples():
    adata = x.create_anndata_from_samples(["GSM1000981, GSM1002540","GSM1009637", "GSM1158460"]) # some samples have no counts
    assert adata.obs.shape[0] == 2

@pytest.mark.xfail
def test_rnaseq_analysis_adata_samples():
    adata = x.create_anndata_from_series("GSE64016")
    assert adata.obs.shape[0] == 460

@pytest.mark.xfail
def test_rnaseq_analysis_test_list_to_dc_geneset():
    df = x.list_to_dc_geneset(["IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6"])
    assert df.shape == (6, 2)

@pytest.mark.xfail
def test_rnaseq_analysis_test_perform_enrichment_on_samples_batched():
    element_list = ["GSM1000981, GSM1002540","GSM1009637", "GSM1158460"]
    geneset = ["IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6"]
    res = x.perform_enrichment_on_samples_batched(element_list, geneset)
    assert res.shape == (2, 9)