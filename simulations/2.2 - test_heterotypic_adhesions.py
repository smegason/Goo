from importlib import reload
import goo
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()

locs_A = [
    (-1, 0, 0),
    (-3, 0, 0),
    (-1, 2, 0),
    (-3, 2, 0),
]

locs_B = [
    (1, 0, 0),
    (3, 0, 0),
    (1, 2, 0),
    (3, 2, 0),
]

celltypeA = goo.SimpleType("A")
celltypeA.homo_adhesion_strength = -1000

celltypeB = goo.SimpleType("B")
celltypeB.homo_adhesion_strength = -1000
celltypeA.set_hetero_adhesion(celltypeB, 3500)

for i, loc in enumerate(locs_A):
    cell = celltypeA.create_cell("cellA" + str(i), loc)
    cell.stiffness = 15
    cell.pressure = 5

for i, loc in enumerate(locs_B):
    cell = celltypeB.create_cell("cellB" + str(i), loc)
    cell.stiffness = 15
    cell.pressure = 5

sim = goo.Simulator([celltypeA, celltypeB])
sim.setup_world()
sim.add_handler(AdhesionLocationHandler())
