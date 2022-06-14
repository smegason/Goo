# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo

goo.setup_world()

# make cells
cell1 = Goo.Cell("cell_", loc=(1, 0, 0))
cell2 = Goo.Cell("cell_", loc=(-1, 0, 0))
