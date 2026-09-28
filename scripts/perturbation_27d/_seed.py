"""Per-dataset bootstrap seeds: rng_for(key) is a generator seeded from 20250908 and the key (dataset + test), so a
bootstrap SE does not depend on which other datasets were processed before it."""
import hashlib
import numpy as np
def rng_for(key):
    return np.random.default_rng([20250908, int(hashlib.md5(key.encode()).hexdigest()[:8], 16)])
