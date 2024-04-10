from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()

celltype = goo.OpaqueType("A")
celltype.homo_adhesion_strength = -5000

celltype2 = goo.OpaqueType("B")
celltype2.homo_adhesion_strength = -5000
celltype.set_hetero_adhesion(celltype2, 5000)

cell = celltype.create_cell("cellA", (0, 0, 0))
cell.stiffness = 15
cell2 = celltype2.create_cell("cellB", (1, 1, 0))
cell2.stiffness = 15

sim = goo.Simulator(celltypes=[celltype, celltype2])
sim.setup_world()

sim.add_handlers(
    [
        TimeDivisionHandler(BisectDivisionLogic, mu=50),
        AdhesionLocationHandler(),
        RemeshHandler(freq=10),
    ]
)
