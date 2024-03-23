
from biorag.transcriptome_embedding import Transcriptome_embedding
import anndata as ad

x = ad.read_h5ad("tests/test_transcriptome_db.h5ad", backed = "r")

x1 = Transcriptome_embedding(x.obs, x.obsm["embedding"])

def test_get_closest_transcriptional_studies():
    res = x1.get_closest_transcriptional_studies("GSE47774,GSE47792", k = 5)
    assert res.shape[0] == 5


def test_ebedding_retreival():
    res = x1.search_relevant_embeddings_by_series(["GSE47774,GSE47792"])
    assert res.shape[0] == 1413

single_sample = ["GSM1000984"]

def test_embedding_single_sample_retreival():
    res = x1.get_closest_transcriptional_studies_by_sample(single_sample, k = 5)
    assert res.shape[0] == 5


