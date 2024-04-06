from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)

goo.reset_modules()
goo.reset_scene()

sim = goo.Simulator()

celltype = goo.create_celltype("A")
cell = celltype.create_cell("cellA", (0, 0, 0))
cell.stiffness = 1

sim.setup_world()
sim.add_celltype(celltype)
sim.add_handlers(
    [
        GrowthPIDHandler(target_volume=30),
        SizeDivisionHandler(BisectDivisionLogic, threshold=29.5),
        RemeshHandler(),
        AdhesionLocationHandler(),
    ]
)

sim.run_simulation(start=1, end=250)
