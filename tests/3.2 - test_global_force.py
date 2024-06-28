from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)

goo.reset_modules()
goo.reset_scene()

goo.create_turbulence("turbulence", (2, 0, 0), -3000)

celltype = goo.CellType("A")
cell = celltype.create_cell("cellA", (0, 0, 0))
# cell.stiffness = 15
# cell.pressure = 5

sim = goo.Simulator(celltypes=[celltype])
sim.toggle_gravity(False)

sim.add_handlers(
    [
        GrowthPIDHandler(target_volume=5),
        # RemeshHandler(),
        AdhesionLocationHandler(),
    ]
)
