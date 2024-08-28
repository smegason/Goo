import sys
import scipy
import numpy
import bpy
import pytest

import goo
from goo.cell import CellType, YolkType, SimpleType
from goo.reloader import *
from goo.simulator import *
from goo.force import *
from goo.division import *
from goo.handler import *
from goo.molecule import * 
from goo.boundary import *


def test_libraries_loaded():
    assert sys.modules.get('bpy') is not None, "bpy is not loaded"
    assert sys.modules.get('goo') is not None, "goo is not loaded"
    assert sys.modules.get('scipy') is not None, "scipy is not loaded"
    assert sys.modules.get('numpy') is not None, "numpy is not loaded"
    assert sys.modules.get('pytest') is not None, "pytest is not loaded"
