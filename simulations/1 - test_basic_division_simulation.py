from importlib import reload
import goo
from goo.handler import *
from goo.division import *

reload(goo)
goo.reset_modules()
goo.reset_scene()


celltype = goo.CellType("cells", physics_enabled=False)
cell = celltype.create_cell("cell", (0, 0, 0), size=5, subdivisions=5)

sim = goo.Simulator(celltypes=[celltype])
division_handler = TimeDivisionHandler(BisectDivisionLogic, mu=10)
sim.add_handler(division_handler)
