import numpy as np
import pandas as pd
import pygalmesh

# Load grid coordinates from R
grid_coords = pd.read_csv("grid_coords.csv").values
d = grid_coords.shape[1]

# Compute bounding box (domain) from grid
mins = grid_coords.min(axis=0)
maxs = grid_coords.max(axis=0)

# Define signed distance function for bounding box
class BoxDomain(pygalmesh.DomainBase):
    def __init__(self, mins, maxs):
        super().__init__(len(mins))
        self.mins = mins
        self.maxs = maxs
        self.center = (mins + maxs) / 2
        self.half_size = (maxs - mins) / 2

    def eval(self, x):
        return max(abs(xi - ci) - hi for xi, ci, hi in zip(x, self.center, self.half_size))

# Build the domain
domain = BoxDomain(mins, maxs)

# Generate the mesh
mesh = pygalmesh.generate_mesh(
    domain,
    max_cell_circumradius=0.2,  # adjust resolution
    num_lloyd_steps=5
)

# mesh.points: mesh vertices (n_points, d)
# mesh.cells[0].data: simplices (n_cells, d+1)

# Save for use in R or Python
np.savez("pygalmesh_mesh.npz", points=mesh.points, simplices=mesh.cells[0].data)
