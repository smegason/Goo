from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()

celltype = goo.SimpleType("A")
cell = celltype.create_cell("cellA", (0, 0, 0))

sim = goo.Simulator([celltype])
sim.setup_world()

sim.add_handlers(
    [
        GrowthPIDHandler(
            target_volume=30,
            growth_type=Growth.LINEAR,
            growth_rate=1,
            initial_pressure=0.01,
            Kp=0.05,
        ),
        RemeshHandler(freq=1, voxel_size=0.25),
        # DataExporter(path=None, options=DataFlag.VOLUMES),
    ]
)

for i in range(50):
    bpy.context.scene.frame_set(i + 1)
    print(cell.volume())
