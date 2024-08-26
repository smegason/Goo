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


celltype = goo.SimpleType("cellA", physics_enabled=False)

cell1 = celltype.create_cell(name="cell", loc=(0, -1.5, 0), color=(0, 0, 0))
network1 = GeneRegulatoryNetwork()
network1.load_circuits(DegFirstOrder(x, 0.08))
cell1.grn = network1
cell1.metabolite_concs = {x: 2}

cell2 = celltype.create_cell(name="cell", loc=(0, 1.5, 0), color=(0, 0, 0))
network2 = GeneRegulatoryNetwork()
cell2.grn = network2
cell2.metabolite_concs = {x: 2}

sim = goo.Simulator(celltypes=[celltype])
sim.setup_world()
sim.add_handlers([NetworkHandler(), ColorizeHandler(Colorizer.GENE, x, range=(0, 2))])

print("done loading")
