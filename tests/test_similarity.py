import anndata as ad
import numpy as np

from sample_explorer.similarity import cosine_similarity
from tests.common import SEMANTIC_DB


def test_matches_formula():
    rng = np.random.default_rng(0)
    a = rng.normal(size=(5, 8))
    b = rng.normal(size=(3, 8))
    expected = (a / np.linalg.norm(a, axis=1, keepdims=True)) @ (
        b / np.linalg.norm(b, axis=1, keepdims=True)
    ).T
    result = cosine_similarity(a, b)
    assert result.shape == (5, 3)
    np.testing.assert_allclose(result, expected, atol=1e-12, rtol=0)


def test_zero_vector_gives_zero():
    assert (cosine_similarity(np.zeros((1, 4)), np.ones((2, 4))) == 0).all()
    assert (cosine_similarity(np.ones((2, 4)), np.zeros((1, 4))) == 0).all()


def test_accepts_h5py_dataset():
    backed = ad.read_h5ad(SEMANTIC_DB, backed="r")
    try:
        dataset = backed.X
        assert not isinstance(dataset, np.ndarray)
        in_memory = np.asarray(dataset[:])
        query = dataset[[0, 1],]
        result = cosine_similarity(dataset, query)
        np.testing.assert_allclose(result, cosine_similarity(in_memory, np.asarray(query)))
    finally:
        backed.file.close()


def test_float16_is_computed_in_float64():
    v = np.full((1, 1000), 3, dtype=np.float16)
    result = cosine_similarity(v, v)
    assert result.dtype == np.float64
    np.testing.assert_allclose(result, [[1.0]], atol=1e-12, rtol=0)


def test_float32_stays_float32():
    v = np.ones((2, 4), dtype=np.float32)
    assert cosine_similarity(v, v).dtype == np.float32
