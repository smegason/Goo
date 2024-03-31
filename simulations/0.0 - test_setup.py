# Goal: test the ability to import and reload modules

import bmesh
from importlib import reload
import goo
from goo import goo as g
from goo.cell import *
from goo.reloader import *

reload(goo)
reset_modules()
reset_scene()

cell = create_cell("cell", (0, 0, 0), physics_on=False)
