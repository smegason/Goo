# cleavage.py - simulates the first 12 cell cycles of zebrafish development

import bpy, treelib as tr
from custom import Goo

Goo.setup_world()

cell = Goo.Cell("cell_", loc = (0, 0, 0))
cell_tree = tr.Tree()
cell_tree.create_node("cell_", "cell_", data = cell.data)