from importlib import reload
import goo
from goo.handler import *
from goo.cell import rsetattr

reload(goo)

goo.reset_modules()
goo.reset_scene()

locs_A = [
    (-1, 0, 0),
    (-3, 0, 0),
    (-1, 2, 0),
    (-3, 2, 0),
]
celltypeA = goo.create_celltype("A", physics_on=True)
celltypeA.homo_adhesion_strength = 1000

locs_B = [
    (1, 0, 0),
    (3, 0, 0),
    (1, 2, 0),
    (3, 2, 0),
]
celltypeB = goo.create_celltype("B", physics_on=True)
celltypeB.homo_adhesion_strength = 1000

celltypeA.set_hetero_adhesion(celltypeB, 3500)

for i, loc in enumerate(locs_A):
    cell = celltypeA.create_cell("cellA" + str(i), loc)
    cell.stiffness = 5
for i, loc in enumerate(locs_B):
    celltypeB.create_cell("cellB" + str(i), loc)
    cell.stiffness = 5

sim = goo.Simulator(celltypes=[celltypeA, celltypeB])
sim.add_handler(ForceUpdateHandler())
sim.toggle_gravity(False)
sim.run_simulation(start=1, end=100)
