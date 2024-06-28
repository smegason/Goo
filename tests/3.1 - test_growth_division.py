from importlib import reload
import goo
from goo import *

reload(goo)

reset_modules()
reset_scene()

sim = Simulator()

celltype = CellType("A")
celltype.homo_adhesion_strength = 5000
cell = celltype.create_cell("cellA", (0, 0, 0))
cell.stiffness = 15

sim.setup_world()
sim.add_celltype(celltype)
sim.add_handlers(
    [
        SizeDivisionHandler(BisectDivisionLogic, threshold=29.5),
        GrowthPIDHandler(target_volume=30),
        RemeshHandler(),
        AdhesionLocationHandler(),
        # ColorizeHandler(Colorizer.PRESSURE),
        # DataExporter(
        #     path="", options=DataFlag.TIMES | DataFlag.VOLUMES | DataFlag.PRESSURES
        # ),
    ]
)

# sim.render(start=1, end=10)
