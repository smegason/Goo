from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()

celltype = goo.SimpleType("A")
celltype.homo_adhesion_strength = -5000

celltype2 = goo.SimpleType("B")
celltype2.homo_adhesion_strength = -5000
celltype.set_hetero_adhesion(celltype2, 10000)

cell = celltype.create_cell("cellA", (0, 0, 0))
cell2 = celltype2.create_cell("cellB", (1.4, 1.4, 0))

sim = goo.Simulator(celltypes=[celltype, celltype2])
sim.setup_world()

sim.add_handlers(
    [
        # TimeDivisionHandler(BisectDivisionLogic, mu=50),
        GrowthPIDHandler(target_volume=5),
        RecenterHandler(),
        RemeshHandler(sphere_factor=0.25),
    ]
)
