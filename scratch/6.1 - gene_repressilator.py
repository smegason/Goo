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


celltype = goo.CellType("cellA", pattern="simple")

cell1 = celltype.create_cell(
    name="cell", loc=(0, -1.5, 0), color=(0, 0, 0), physics_enabled=False
)
network1 = GeneRegulatoryNetwork()
network1.load_circuits(
    DegFirstOrder(x, 0.1),
    DegFirstOrder(y, 0.1),
    DegFirstOrder(z, 0.1),
    ProdRepression(y, x, kcat=0.4, n=3),
    ProdRepression(z, y, kcat=0.4, n=3),
    ProdRepression(x, z, kcat=0.4, n=3),
)
cell1.grn = network1
cell1.metabolites = {x: 2, y: 0.1, z: 0.1}

cell2 = celltype.create_cell(
    name="cell", loc=(0, 1.5, 0), color=(0, 0, 0), physics_enabled=False
)
network2 = GeneRegulatoryNetwork()
cell2.grn = network2
cell2.metabolites = {x: 2}

sim = goo.Simulator(celltypes=[celltype])
sim.setup_world()
sim.add_handlers([NetworkHandler(), ColorizeHandler(Colorizer.GENE, x, range=(0, 2))])

print("done loading")
