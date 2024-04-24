import numpy as np
from goo import goo
from importlib import reload

reload(goo)
goo.setup_world(seed=1)

def random_point_inside_sphere(radius):
    """Generate a random point inside a sphere."""
    r = (radius - 1) * np.cbrt(np.random.rand())
    theta = np.random.uniform(0, 2*np.pi)
    phi = np.random.uniform(0, np.pi)
    x = r * np.sin(phi) * np.cos(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(phi)
    return x, y, z

def check_min_distance(new_point, points, min_distance):
    """Check if the new point meets the minimum distance requirement."""
    for point in points:
        distance = np.linalg.norm(np.array(new_point) - np.array(point))
        if distance < min_distance:
            return False
    return True

# Create cells within the sphere with minimum distance constraint
num_cells = 6  # Number of cells
radius = 4  # Sphere radius
min_distance = 2.5  # Minimum distance between cells

cells = []

while len(cells) < num_cells:
    new_point = random_point_inside_sphere(radius)
    if check_min_distance(new_point, cells, min_distance):
        cells.append(new_point)
        cell_name = f"cell_A{len(cells)}"
        goo.make_cell(cell_name, loc=new_point, type="cellsA")

# Add forces and motions
homoA = -10000
motionA = -2500

for i in range(num_cells):
    cell_name = f"cell_A{i+1}"
    goo.add_homo_adhesion(cell_name, homoA)
    goo.add_motion(cell_name, motionA)

goo.add_sphere_boundaries(loc=(0, 0, 0), radius=radius)

# Simulation setup
handlers = goo.handler_class()
handlers.launch_simulation(start=1,
                           end=300,
                           filepath="/Users/antoine/Harvard/MegasonLab/GPU_backup/AntoineRuzette/goo/data/motion_diffusion/20240422_cell_mixing/img2_",
                           adhesion=True,
                           data=True,
                           growth='tissue',
                           motility=True,
                           division=False,
                           target_volume=50)
