from importlib import reload
import goo
from goo.gene import *
from goo.division import *
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()

celltype = goo.CellType("cellA", pattern="simple", target_volume=70)
celltype.homo_adhesion_strength = 500
cell = celltype.create_cell(name="cell1", loc=(0, 0, 0), color=(0, 0, 0))
cell.stiffness = 5

# Defining genes
x = Gene("x")
y = Gene("y")
z = Gene("z")

network1 = GeneRegulatoryNetwork()
network1.load_circuits(
    DegFirstOrder(x, 0.1),
    DegFirstOrder(y, 0.1),
    DegFirstOrder(z, 0.1),
    ProdRepression(y, x, kcat=0.4, n=3),
    ProdRepression(z, y, kcat=0.4, n=3),
    ProdRepression(x, z, kcat=0.4, n=3),
)
cell.grn = network1
cell.metabolites = {x: 2, y: 0.1, z: 0.1}

sim = goo.Simulator(celltypes=[celltype], time=200, physics_dt=1)
sim.setup_world(seed=2025)
sim.add_handlers(
    [
        goo.GrowthPIDHandler(),
        goo.SizeDivisionHandler(goo.BisectDivisionLogic, mu=60, sigma=2),
        goo.RecenterHandler(),
        goo.RemeshHandler(),
        goo.NetworkHandler(),
        goo.ColorizeHandler(goo.Colorizer.GENE, x, range=(1, 2)),
    ]
)