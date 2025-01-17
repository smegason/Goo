from importlib import reload
import goo
from goo.gene import *
from goo.division import *
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()

# Defining cells
x = Gene("x")
y = Gene("y")
z = Gene("z")

celltype = goo.CellType("cellA", pattern="simple", target_volume=70)
celltype.homo_adhesion_strength = 500
celltype.motion_strength = 1000
cell = celltype.create_cell(name="cell1", loc=(0, 0, 0), color=(0, 0, 0))
cell.stiffness = 5

network1 = GeneRegulatoryNetwork()
network1.load_circuits(
    DegFirstOrder(x, 0.1),
    DegFirstOrder(y, 0.1),
    DegFirstOrder(z, 0.1),
    ProdRepression(y, x, kcat=0.2, n=3),
    ProdRepression(z, y, kcat=0.2, n=3),
    ProdRepression(x, z, kcat=0.2, n=3)
    )
cell.grn = network1
cell.metabolites = {x: 2, y: 0.1, z: 0.1}
cell.link_gene_to_property(gene=x, property="motion_strength")

sim = goo.Simulator(celltypes=[celltype], time=500, physics_dt=1)
sim.setup_world(seed=2025)
sim.add_handlers(
    [
        goo.GrowthPIDHandler(),
#        goo.SizeDivisionHandler(goo.BisectDivisionLogic, mu=60, sigma=2),
        goo.RecenterHandler(),
#        goo.RemeshHandler(),
        goo.NetworkHandler(),
        goo.ColorizeHandler(goo.Colorizer.GENE, x, range=(1, 2)),
        goo.RandomMotionHandler(goo.ForceDist.GAUSSIAN)
    ]
)