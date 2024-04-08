from importlib import reload
import goo
from goo import *

reload(goo)

reset_modules()
reset_scene()

sim = Simulator()

celltype = create_celltype("A")
cell = celltype.create_cell("cellA", (0, 0, 0))
cell.stiffness = 15

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

sim.render(start=1, end=10)
