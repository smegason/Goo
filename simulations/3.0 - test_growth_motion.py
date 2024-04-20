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
cell.motion_force.strength = 8000

sim = goo.Simulator([celltype])
sim.setup_world()
sim.add_handlers(
    [
        GrowthPIDHandler(),
        AdhesionLocationHandler(),
        MotionHandler(),
        RemeshHandler(),
    ]
)
