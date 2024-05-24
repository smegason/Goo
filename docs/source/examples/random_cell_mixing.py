from importlib import reload
import goo
from goo import *

reload(goo)
reset_modules()
reset_scene()


def random_point_inside_sphere(radius):
    """Generate a random point inside a sphere."""
    r = (radius - 1) * np.cbrt(np.random.rand())
    theta = np.random.uniform(0, 2 * np.pi)
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

create_boundary((0, 0, 0), size=radius * 1.2)

cellsA = SimpleType("A")
cellsA.homo_adhesion_strength = 2500
cells = []

while len(cells) < num_cells:
    new_point = random_point_inside_sphere(radius)
    if check_min_distance(new_point, cells, min_distance):
        cells.append(new_point)
        cell_name = f"cell_A{len(cells)}"
        color = tuple(np.random.random_sample(3))
        cellsA.create_cell(cell_name, new_point, color=color)

sim = Simulator([cellsA])
sim.setup_world()
sim.toggle_gravity(True)
sim.add_handlers(
    [
        GrowthPIDHandler(target_volume=30),
        AdhesionLocationHandler(),
        RandomMotionHandler(distribution=ForceDist.CONSTANT, max_strength=2500),
        DataExporter(
            path="/tmp/out.json", options=DataFlag.TIMES | DataFlag.CONTACT_AREAS
        ),
        SceneExtensionHandler(300),
    ]
)
sim.run(300)
