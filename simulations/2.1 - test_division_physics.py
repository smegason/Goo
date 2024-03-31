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

make_force("Force", (-1, 1, 0), 500)
make_force("Force", (1, -1, 0), 500)
celltype = CellType("default", physics_on=True)
# celltype.homo_adhesion_strength = 0
cell = celltype.create_cell("cell", (1, 1, 0), size=2)

sim = Simulator(celltypes=[celltype])
division_handler = TimeDivisionPhysicsHandler(BisectDivisionLogic, mu=20)
sim.add_handler(division_handler)

sim.toggle_gravity(False)
sim.run_simulation(start=1, end=49)
