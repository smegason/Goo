from importlib import reload
import goo
from goo.handler import *

reload(goo)

goo.reset_modules()
goo.reset_scene()

celltype = goo.create_celltype("default", physics_on=True)
celltype.set_homo_adhesion(2000)

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
    celltype.create_cell("cell" + str(i), loc)

sim = goo.Simulator(celltypes=[celltype])
sim.add_handler(ForceUpdateHandler())
sim.toggle_gravity(False)
sim.run_simulation(start=1, end=250)
