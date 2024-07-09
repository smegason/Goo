from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()

# Defining cells
celltype = goo.SimpleType("cellA")
celltype.homo_adhesion_strength = 1000
cell1 = celltype.create_cell(name="cell", loc=(0, -1.5, 0), color=(0, 1, 1))
cell2 = celltype.create_cell(name="cell", loc=(0, 1.5, 0), color=(0, 1, 1))
cell1.stiffness = 2
cell1.pressure = 5
cell2.stiffness = cell1.stiffness
cell2.pressure = cell1.pressure

molA = goo.Molecule("molA", conc=1, D=1)
diffusionsystem = goo.DiffusionSystem(molecules=[molA])

sim = goo.Simulator(celltypes=[celltype], diffsystems=[diffusionsystem])
sim.setup_world()
sim.add_handlers(
    
    [
        GrowthPIDHandler(target_volume=50),
        AdhesionLocationHandler(),
        RandomMotionHandler(distribution=ForceDist.CONSTANT, max_strength=5000),
        DiffusionHandler()
    ]
)