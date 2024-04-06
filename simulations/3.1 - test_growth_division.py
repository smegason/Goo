from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)

goo.reset_modules()
goo.reset_scene()

celltype = goo.create_celltype("A")
cell = celltype.create_cell("cellA", (0, 0, 0))
cell.stiffness = 1

sim = goo.Simulator(celltypes=[celltype])
sim.toggle_gravity(False)

sim.add_handlers(
    [
        GrowthPIDHandler(target_volume=30),
        SizeDivisionHandler(BisectDivisionLogic, threshold=29.5),
        RemeshHandler(),
        AdhesionLocationHandler(),
    ]
)

sim.run_simulation(start=1, end=250)
