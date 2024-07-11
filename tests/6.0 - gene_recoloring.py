from importlib import reload
import goo
from goo.circuits import *
from goo.division import *
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()

# Defining cells
x = Gene("x")
grn = Circuit(deg_first_order(x, 0.1))
celltype = goo.SimpleType("cellA", physics_enabled=False, grn=grn)

cell1 = celltype.create_cell(name="cell", loc=(0, -1.5, 0))
cell2 = celltype.create_cell(name="cell", loc=(0, 1.5, 0))

cell1.genes_conc['x'] = 3.0

print(cell1.circuit_tellurium)

sim = goo.Simulator(celltypes=[celltype])
sim.setup_world()
sim.add_handlers([GeneHandler()])
