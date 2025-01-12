from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()

celltype = goo.CellType("A", pattern="standard")
cell = celltype.create_cell("cellA", (0, 0, 0), target_volume=50)
cell.pressure = 20

sim = goo.Simulator([celltype])
sim.setup_world()

sim.add_handlers(
    [
        GrowthPIDHandler(),
        RemeshHandler(),
        # DataExporter(path=None, options=DataFlag.VOLUMES),
    ]
)

for i in range(25):
    bpy.context.scene.frame_set(i + 1)
    print(cell.volume(), cell.growth_controller.next_volume, cell.pressure)
