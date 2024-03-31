from importlib import reload
import goo
from goo.cell import *
from goo.reloader import *
from goo.division import *
from goo.simulator import *
from goo.force import *

reload(goo)

reset_modules()
reset_scene()

celltype = CellType("default", physics_on=True)
locs = [
    (0, 0, 0),
    (0, 2, 0),
    (2, 0, 0),
    (2, 2, 0),
    (1, -0.5, 2),
    (2.5, 1, 2),
    (1, 2.5, 2),
    (-0.5, 1, 2),
]

for i, loc in enumerate(locs):
    celltype.create_cell("cell" + str(i), loc)
celltype.homo_adhesion_strength = 1000

sim = Simulator()
sim.toggle_gravity(False)
