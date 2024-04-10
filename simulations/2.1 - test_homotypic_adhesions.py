from importlib import reload
import goo
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()

celltype = goo.OpaqueType("default", physics_on=True)
celltype.set_homo_adhesion(5000)

locs = [
    (-1.36, -1, 0.2),
    (-0.2, -1, -1.36),
    (1.36, -1, -0.2),
    (0.2, -1, 1.36),
    (-0.73, 1, 1.1),
    (1.2, 1, 0.84),
    (0.92, 1, -1.08),
    (-1, 1, -0.82),
]

for i, loc in enumerate(locs):
    cell = celltype.create_cell("cell" + str(i), loc)
    cell.stiffness = 15

sim = goo.Simulator(celltypes=[celltype])
sim.setup_world()
sim.add_handler(AdhesionLocationHandler())
