from goo import goo
import sys
from importlib import reload
import numpy as np


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


def main(seed, motion, adhesion):
    # reload(goo)
    goo.setup_world(seed=seed)

    # Create cells within the sphere with minimum distance constraint
    num_cells = 10  # Number of cells
    radius = 13  # Sphere radius
    min_distance = 2.5  # Minimum distance between cells

    cells = []

    while len(cells) < num_cells:
        new_point = random_point_inside_sphere(radius)
        if check_min_distance(new_point, cells, min_distance):
            cells.append(new_point)
            cell_name = f"cell_A{len(cells)}"
            goo.make_cell(cell_name, loc=new_point, type="cellsA")

    # Add forces and motions
    homoA = -motion
    motionA = -adhesion

    for i in range(num_cells):
        cell_name = f"cell_A{i+1}"
        goo.add_homo_adhesion(cell_name, homoA)
        goo.add_motion(cell_name, motionA)

    goo.add_sphere_boundaries(loc=(0, 0, 0), radius=radius)

    # Simulation setup
    handlers = goo.handler_class()
    handlers.launch_simulation(start=1,
                            end=300,
                            filepath="/Users/antoine/Harvard/MegasonLab/GPU_backup/AntoineRuzette/goo/data/motion_diffusion/20240422_cell_mixing/img2_seed{s}_adh{adhesion}_motion{motion}",
                            adhesion=True,
                            data=True,
                            growth='tissue',
                            motility=True,
                            division=False,
                            target_volume=50, 
                            growth_rate=2
                            headless=True)


if __name__ == "__main__":
    # Check if the "--seed" argument is provided
    if "--seed" in sys.argv:
        seed_idx = sys.argv.index("--seed") + 1
        adh_idx = sys.argv.index("--adhesion") + 1
        motion_idx = sys.argv.index("--motion") + 1

        seed = int(sys.argv[seed_idx])
        adhesion = int(sys.argv[adh_idx])
        motion = int(sys.argv[motion_idx])

        # Call the main function with the provided seed
        main(seed=seed, adhesion=adhesion, motion=motion)
    else:
        print("Error: Seed argument is missing. Please provide the --seed argument.")
