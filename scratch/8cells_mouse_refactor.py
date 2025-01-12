from importlib import reload
import goo
from goo import *

reload(goo)
reset_modules()
reset_scene()

sim = Simulator()

celltype = SimpleType("A")
celltype.homo_adhesion_strength = 10000
celltype.motion_strength = 500

locs = [
    (-1.36, -1, 0.2),
    (-0.2, -1, -1.36),
    (1.36, -1, -0.2),
    (0.2, -1, 1.36),
    (-0.73, 1, 1.1),
    (1.2, 1, 0.84),
    (0.92, 1, -1.08),
    (-1, 1, -0.82),
]

for i, loc in enumerate(locs):
    celltype.create_cell(f"cell_A{i}", loc)

sim.setup_world()
sim.add_celltype(celltype)
sim.add_handlers(
    [
        GrowthPIDHandler(),
        RecenterHandler(),
        RandomMotionHandler(),
        SceneExtensionHandler(end=500),
    ]
)

# sim.render(end=500)
