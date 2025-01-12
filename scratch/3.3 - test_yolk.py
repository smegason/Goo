from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)

goo.reset_modules()
goo.reset_scene()

yolks = goo.YolkType("yolks")
yolks.create_cell("yolkA", (0, 0, 0))

cells = goo.CellType("A")
cells.create_cell("cellA", (10, 0, 0))

sim = goo.Simulator([cells, yolks])
sim.setup_world()
sim.add_handlers(
    [
        RemeshHandler(voxel_size=0),
        RecenterHandler(),
    ]
)
sim.add_handlers(
    [
        GrowthPIDHandler(target_volume=30),
        SizeDivisionHandler(BisectDivisionLogic, threshold=29.5),
    ],
    celltypes=[cells],
)
