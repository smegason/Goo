from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)

goo.reset_modules()
goo.reset_scene()

celltype = goo.CellType("A")
cell = celltype.create_cell("cellA", (0, 0, 0))
goo.create_turbulence("turbulence", (2, 0, 0), -3000)

sim = goo.Simulator(celltypes=[celltype])
sim.toggle_gravity(False)

sim.add_handlers(
    [
        GrowthPIDHandler(target_volume=30),
        RemeshHandler(voxel_size=0),
        AdhesionLocationHandler(),
    ]
)
