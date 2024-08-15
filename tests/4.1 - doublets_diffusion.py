from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()

celltype = goo.SimpleType("cellA")
celltype.homo_adhesion_strength = 0
cell1 = celltype.create_cell(name="cell", loc=(0, 5, 0), color=(0, 1, 1))
# cell2 = celltype.create_cell(name="cell1", loc=(3, 0, 0), color=(0, 1, 1))
cell1.stiffness = 1
cell1.pressure = 5

molA = goo.Molecule("molA", conc=1, D=1, gradient="linear")
molB = goo.Molecule("molB", conc=0.5, D=2, gradient="random")

diffusionsystem = goo.DiffusionSystem(molecules=[molA, molB])

sim = goo.Simulator(celltypes=[celltype], diffsystem=[diffusionsystem], time=500)
sim.setup_world()
sim.add_handlers(
    [
        GrowthPIDHandler(target_volume=50),
        # SizeDivisionHandler(BisectDivisionLogic, mu=60, sigma=10),
        RecenterHandler(),
        # RemeshHandler(),
        # RandomMotionHandler(distribution=ForceDist.CONSTANT, max_strength=2000),
        DiffusionHandler(diffusionSystem=diffusionsystem),
        # MolecularSensingHandler(diffusionSystem=diffusionsystem),
        DataExporter(
            path="/Users/antoine/Harvard/MegasonLab/GPU_backup/AntoineRuzette/goo/data/diffusion/20240701_molecules_secretion/linear_gradient.json",
            options=DataFlag.CELL_CONCENTRATIONS | DataFlag.MOTION_PATH,  # or ALL flag
        ),
    ]
)
