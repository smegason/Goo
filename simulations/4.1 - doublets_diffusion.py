from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()

celltype = goo.SimpleType("cellA")
celltype.homo_adhesion_strength = 5000
cell1 = celltype.create_cell(name="cell", loc=(0, -1.5, 0), color=(0, 1, 1))
cell2 = celltype.create_cell(name="cell", loc=(0, 1.5, 0), color=(0, 1, 1))
cell1.stiffness = 5
cell1.pressure = 5
cell2.stiffness = cell1.stiffness
cell2.pressure = cell1.pressure

molA = goo.Molecule("molA", conc=1, D=1)
molB = goo.Molecule("molB", conc=0.5, D=2)

diffusionsystem = goo.DiffusionSystem(molecules=[molA, molB])

sim = goo.Simulator(celltypes=[celltype])
sim.setup_world()
sim.add_handlers(
    
    [
        GrowthPIDHandler(target_volume=50),
        # SizeDivisionHandler(BisectDivisionLogic, mu=60, sigma=10),
        AdhesionLocationHandler(),
        RemeshHandler(),
        # RandomMotionHandler(distribution=ForceDist.CONSTANT, max_strength=0),
    ]
)

# Run simulation headless
output_path = "/Users/antoine/Harvard/MegasonLab/GPU_backup/AntoineRuzette/goo/data/adhesion/20240626_adhesion_surfaces/20240626_adhesion_POINTS"
sim.render(start=1, end=50, path=output_path)