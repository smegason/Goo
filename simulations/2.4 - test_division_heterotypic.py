from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)

goo.reset_modules()
goo.reset_scene()

celltype = goo.create_celltype("A")
celltype.homo_adhesion_strength = -5000

celltype2 = goo.create_celltype("B")
celltype2.homo_adhesion_strength = -5000
celltype.set_hetero_adhesion(celltype2, 5000)

celltype.create_cell("cellA", (0, 0, 0))
celltype2.create_cell("cellB", (1, 1, 0))

sim = goo.Simulator(celltypes=[celltype, celltype2])
sim.toggle_gravity(False)

sim.add_handlers(
    [
        TimeDivisionHandler(BisectDivisionLogic, mu=50),
        AdhesionLocationHandler(),
        RemeshHandler(freq=10),
    ]
)

sim.run_simulation(start=1, end=250)
