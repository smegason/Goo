from importlib import reload
import goo
from goo.handler import *
from mathutils import Euler

reload(goo)
goo.reset_modules()
goo.reset_scene()


celltype = goo.SimpleType("cellsA")
cell = celltype.create_cell("cell", (0, 0, 0), size=2)
cell.stiffness = 15
cell.pressure = 1

sim = goo.Simulator([celltype])
sim.setup_world()
sim.add_handlers(
    [
        GrowthPIDHandler(),
        RecenterHandler(),
        RandomMotionHandler(ForceDist.CONSTANT, max_strength=8000),
        RemeshHandler(),
        DataExporter(
            path="/tmp/out.json",
            options=DataFlag.TIMES | DataFlag.MOTION_PATH | DataFlag.FORCE_PATH,
        ),
    ]
)

# for i in range(200):
#     bpy.context.scene.frame_set(i + 1)
