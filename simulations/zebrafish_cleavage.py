# zebrafish_cleavage.py - simulates the first 12 cell cycles of zebrafish development

import Goo from Goo

Goo.setup_world()

#make the blastomere
cell = Goo.Cell("cell_", loc = (0, 0, 0))

#make the yolk
yolk = Goo.Cell("yolk", loc = (0,0, -1))