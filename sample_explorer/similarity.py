import numpy as np


def _as_float_2d(x, dtype):
    x = np.asarray(x)
    if x.ndim == 1:
        x = x.reshape(1, -1)
    return x.astype(dtype, copy=False)


def _normalise_rows(x):
    norms = np.linalg.norm(x, axis=1, keepdims=True)
    norms[norms == 0] = 1.0  # zero rows stay zero, so their similarity is 0
    return x / norms


def cosine_similarity(a, b):
    """
    Cosine similarity between every row of ``a`` and every row of ``b``.

    Accepts array-likes, including h5py datasets. Rows with zero norm have
    similarity 0, as in scikit-learn. Like scikit-learn, the computation is in
    float32 only when both inputs are float32; every other dtype (float16,
    integers, mixed) is computed in float64.

    Returns:
        np.ndarray of shape (len(a), len(b)).
    """
    a = np.asarray(a)
    b = np.asarray(b)
    dtype = np.float32 if a.dtype == np.float32 and b.dtype == np.float32 else np.float64
    a = _as_float_2d(a, dtype)
    b = _as_float_2d(b, dtype)
    return _normalise_rows(a) @ _normalise_rows(b).T
