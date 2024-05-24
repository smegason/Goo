from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()

force = goo.create_force("Force", (0, 0, 0), 5000)
force.obj.field.shape = "PLANE"

celltype = goo.SimpleType("cellA")
celltype.homo_adhesion_strength = 5000
cell = celltype.create_cell("cell", (0, 0, 3))
cell.stiffness = 15
cell.pressure = 5

sim = goo.Simulator(celltypes=[celltype])
sim.setup_world()
sim.add_handlers(
    [
        TimeDivisionHandler(BisectDivisionLogic, mu=20),
        AdhesionLocationHandler(),
        RemeshHandler(),
    ]
)
