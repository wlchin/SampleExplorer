"""Torch-free sentence embeddings: paraphrase-MiniLM-L6-v2 run with ONNX Runtime."""
from collections.abc import Sequence

import numpy as np
import onnxruntime as ort
from huggingface_hub import hf_hub_download
from tokenizers import Tokenizer

MODEL_REPO = "sentence-transformers/paraphrase-MiniLM-L6-v2"
MODEL_REVISION = "c9a2bfebc254878aee8c3aca9e6844d5bbb102d1"
MAX_SEQ_LENGTH = 128
EMBEDDING_DIM = 384


class MiniLMEncoder:
    """Encodes text to mean-pooled, unnormalised MiniLM vectors, as sentence-transformers does."""

    def __init__(self, repo: str = MODEL_REPO, revision: str = MODEL_REVISION):
        self.tokenizer = Tokenizer.from_file(
            hf_hub_download(repo, "tokenizer.json", revision=revision)
        )
        self.tokenizer.enable_truncation(max_length=MAX_SEQ_LENGTH)
        self.tokenizer.enable_padding()
        self.session = ort.InferenceSession(
            hf_hub_download(repo, "onnx/model.onnx", revision=revision),
            providers=["CPUExecutionProvider"],
        )

    def encode(self, texts: str | Sequence[str], batch_size: int = 32) -> np.ndarray:
        """Return float32 vectors: shape (384,) for a str, (n, 384) for a sequence."""
        single = isinstance(texts, str)
        items = [texts] if single else list(texts)
        batches = [
            self._encode_batch(items[i : i + batch_size])
            for i in range(0, len(items), batch_size)
        ]
        out = (
            np.vstack(batches)
            if batches
            else np.empty((0, EMBEDDING_DIM), dtype=np.float32)
        )
        return out[0] if single else out

    def _encode_batch(self, batch: list[str]) -> np.ndarray:
        enc = self.tokenizer.encode_batch(batch)
        mask = np.array([e.attention_mask for e in enc], dtype=np.int64)
        feeds = {
            "input_ids": np.array([e.ids for e in enc], dtype=np.int64),
            "attention_mask": mask,
            "token_type_ids": np.array([e.type_ids for e in enc], dtype=np.int64),
        }
        hidden = self.session.run(None, feeds)[0]
        m = mask[..., None].astype(np.float32)
        pooled = (hidden * m).sum(axis=1) / np.clip(m.sum(axis=1), 1e-9, None)
        return pooled.astype(np.float32)
