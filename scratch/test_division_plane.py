from importlib import reload
import goo
from goo import *

reload(goo)
reset_modules()
reset_scene()

cellsA = SimpleType("A")

cells = [
    cellsA.create_cell("cell_A1", (0, +0, 0), scale=(1, 0.8, 0.8), rotation=(0, 45, 0)),
    cellsA.create_cell("cell_A2", (0, +4, 0), scale=(0.8, 0.8, 1), rotation=(45, 0, 0)),
    cellsA.create_cell("cell_A3", (0, -4, 0), scale=(0.8, 1, 0.8), rotation=(0, 0, 45)),
]

for cell in cells:
    logic = BisectDivisionLogic(margin=0.025)
    cell.divide(logic)
    logic.flush()
