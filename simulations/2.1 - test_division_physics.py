from importlib import reload
import goo
from goo.division import *
from goo.handler import ForceUpdateHandler

reload(goo)

goo.reset_modules()
goo.reset_scene()

goo.make_force("Force", (-1, 1, 0), 500)
goo.make_force("Force", (1, -1, 0), 500)
celltype = goo.CellType("default", physics_on=True)
cell = celltype.create_cell("cell", (1, 1, 0), size=2)

sim = goo.Simulator(celltypes=[celltype])
division_handler = TimeDivisionPhysicsHandler(BisectDivisionLogic, mu=20)
# TODO: integrate division and origin change handling
sim.add_handler(division_handler)

sim.toggle_gravity(False)
sim.run_simulation(start=1, end=49)
