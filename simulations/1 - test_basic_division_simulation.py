from importlib import reload
import goo
from goo.cell import *
from goo.reloader import *
from goo.division import *
from goo.simulator import *

reload(goo)

reset_modules()
reset_scene()

celltype = CellType("default", physics_on=False)
cell = celltype.create_cell("cell", (0, 0, 0), size=5)

sim = Simulator(celltypes=[celltype])
division_handler = TimeDivisionHandler(BisectDivisionLogic)
sim.add_handler(division_handler)

sim.run_simulation(start=1, end=80)
