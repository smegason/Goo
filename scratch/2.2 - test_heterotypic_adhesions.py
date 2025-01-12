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

celltypeA = goo.CellType("A", pattern="simple", homo_adhesion_strength=100, size=1)

celltypeB = goo.CellType("B", pattern="simple", homo_adhesion_strength=100, size=1)
celltypeA.set_hetero_adhesion_strength(celltypeB, 350)

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
sim.add_handler(RecenterHandler())
