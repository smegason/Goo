from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)

goo.reset_modules()
goo.reset_scene()

celltype = goo.create_celltype("cellA")
celltype.homo_adhesion_strength = 2000

cell = celltype.create_cell("cell", (1, 1, 0), scale=(1, 1, 1))

sim = goo.Simulator(celltypes=[celltype])
sim.toggle_gravity(False)

sim.add_handler(TimeDivisionHandler(BisectDivisionLogic, mu=51))
sim.add_handler(AdhesionLocationHandler())
sim.add_handler(RemeshHandler(freq=5))

sim.run_simulation(start=1, end=250)
