from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)

goo.reset_modules()
goo.reset_scene()

force = goo.create_force("Force", (0, 0, 0), 1000)

celltype = goo.OpaqueType("cellA")
celltype.homo_adhesion_strength = 5000
cell = celltype.create_cell("cell", (1, 1, 0), mesh_kwargs={"scale": (1, 1, 1)})
cell.stiffness = 15

sim = goo.Simulator(celltypes=[celltype])
sim.setup_world()
sim.add_handlers(
    [
        TimeDivisionHandler(BisectDivisionLogic, mu=51),
        AdhesionLocationHandler(),
        # RemeshHandler(freq=5),
    ]
)
