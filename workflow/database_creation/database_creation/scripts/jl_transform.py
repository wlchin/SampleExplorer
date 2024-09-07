import numpy as np
from sklearn import random_projection
import pickle

memmap_filename = "results/study_counts.dat"

rows = 287553
cols = 67186

# Reload the memmap file in read-only mode
reloaded_memmap_matrix = np.memmap(memmap_filename, dtype='float16', mode='r', shape=(rows, cols))

# Now you can work with the reloaded memmap matrix
print("Reloaded Memmap Matrix:")
print(reloaded_memmap_matrix)

transformer = random_projection.GaussianRandomProjection(n_components=1000)
X_new = transformer.fit_transform(reloaded_memmap_matrix)
X_new.shape

pickle_file_path = "results/matrix_v2.pkl"

# Save the NumPy matrix using pickle
with open(pickle_file_path, 'wb') as file:
    pickle.dump(X_new, file)