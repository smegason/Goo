# Goal: test the ability to import and reload modules

from importlib import reload
import goo
from goo import *

reload(goo)
reset_modules()

reset_scene()
cell = create_cell("cell", (0, 0, 0))
