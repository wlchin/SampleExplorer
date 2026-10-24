"""MiniLMEncoder reproduces the sentence-transformers paraphrase-MiniLM-L6-v2 vectors."""
from pathlib import Path

import anndata as ad
import numpy as np

from sample_explorer.embedding import MiniLMEncoder

SEMANTIC_DB = Path(__file__).parent / "test_semantic_db.h5ad"

LONG_TEXT = " ".join(f"gene{i} expression" for i in range(200))


def test_matches_stored_semantic_vectors():
    sem = ad.read_h5ad(SEMANTIC_DB)
    emb = MiniLMEncoder().encode(sem.obs["rag_text"].tolist())
    assert np.abs(emb - np.asarray(sem.X)).max() < 1e-5


def test_shapes_and_dtype():
    enc = MiniLMEncoder()
    assert enc.encode("x").shape == (384,)
    assert enc.encode(["a", "b"]).shape == (2, 384)
    assert enc.encode([]).shape == (0, 384)
    assert enc.encode("").shape == (384,)
    for out in (enc.encode("x"), enc.encode(["a", "b"]), enc.encode([])):
        assert out.dtype == np.float32


def test_padding_invariance():
    enc = MiniLMEncoder()
    together = enc.encode(["short", LONG_TEXT])
    assert np.abs(together[0] - enc.encode("short")).max() < 1e-6
    assert np.abs(together[1] - enc.encode(LONG_TEXT)).max() < 1e-6


def test_truncates_at_128_tokens():
    enc = MiniLMEncoder()
    assert np.array_equal(enc.encode(LONG_TEXT), enc.encode(LONG_TEXT + " extra words"))


def test_download_model_fetches_pinned_files(monkeypatch):
    from sample_explorer import embedding
    calls = []

    def fake_download(repo, filename, revision):
        calls.append((repo, filename, revision))
        return f"/cache/{filename}"

    monkeypatch.setattr(embedding, "hf_hub_download", fake_download)
    assert embedding.download_model() == ["/cache/tokenizer.json", "/cache/onnx/model.onnx"]
    assert calls == [(embedding.MODEL_REPO, "tokenizer.json", embedding.MODEL_REVISION),
                     (embedding.MODEL_REPO, "onnx/model.onnx", embedding.MODEL_REVISION)]
