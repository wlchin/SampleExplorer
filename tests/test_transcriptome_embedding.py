
from biorag.transcriptome_embedding import Transcriptome_embedding
import anndata as ad

ad.read_h5ad("tests/test_transcriptome_db.h5ad", backed = "r")

