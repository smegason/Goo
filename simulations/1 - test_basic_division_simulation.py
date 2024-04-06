from importlib import reload
import goo
from goo.division import *

reload(goo)

goo.reset_modules()
goo.reset_scene()

celltype = goo.create_celltype("default", physics_on=False)
cell = celltype.create_cell("cell", (0, 0, 0), size=5, subdivisions=5)

sim = goo.Simulator(celltypes=[celltype])
division_handler = TimeDivisionHandler(BisectDivisionLogic, mu=10)
sim.add_handler(division_handler)

sim.run_simulation(start=1, end=80)
