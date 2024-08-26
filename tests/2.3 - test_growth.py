from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()

celltype = goo.CellType("A", pattern="simple")
cell = celltype.create_cell("cellA", (0, 0, 0))

sim = goo.Simulator([celltype])
sim.setup_world()

sim.add_handlers(
    [
        GrowthPIDHandler(),
        # RemeshHandler(freq=1, voxel_size=0.25),
        # DataExporter(path=None, options=DataFlag.VOLUMES),
    ]
)

for i in range(50):
    bpy.context.scene.frame_set(i + 1)
    print(cell.volume())
