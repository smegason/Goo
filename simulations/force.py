import sys
from re import T

import bpy
import numpy as np


class Force():
    """
    A class for representing forces between cells in Goo. Currently, only adhesion forces are supported.  

    :param force_name: Name of force
    :param cell_name: Name of cell
    :param strength: Strength of force
    :param falloff_power: Power of the falloff of force
    """
    def __init__(self, force_name, cell_name, strength, falloff_power):
        self.name = force_name
        self.strength = strength
        self.associated_cell = cell_name
        self.falloff_power = falloff_power
        self.falloff_type = 'SPHERE'
        #self.guide_clump_amount = clumping

def make_force(force):
        """
        Makes a Blender force from the Goo force. 

        :param force: a Goo force

        :return: None
        """
        # Add a force object
        cell = force.associated_cell
        bpy.ops.object.effector_add(type='FORCE',
                                    enter_editmode=False,
                                    align='WORLD',
                                    location=bpy.data.objects[cell].location,
                                    scale=(1, 1, 1))

        # Add force parameters
        bpy.context.object.field.strength = force.strength
        bpy.context.object.name = force.name
        bpy.context.object.field.falloff_power = force.falloff_power
        #bpy.context.object.field.guide_clump_amount = force.guide_clump_amount





'''def initialize_cell_sheet():
    bpy.ops.mesh.primitive_grid_add(size=2,
                                    enter_editmode=False,
                                    align='WORLD',
                                    location=(0, 0, 0),
                                    scale=(8.28558, 8.28558, 8.28558))
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.transform.shear(value=-0.5,
                            orient_axis='Z',
                            orient_axis_ortho='Y',
                            orient_type='GLOBAL',
                            orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                            orient_matrix_type='GLOBAL',
                            mirror=True,
                            use_proportional_edit=False,
                            proportional_edit_falloff='SMOOTH',
                            proportional_size=1,
                            use_proportional_connected=False,
                            use_proportional_projected=False,
                            release_confirm=True)
    bpy.ops.object.mode_set(mode='OBJECT')
    
    bpy.ops.transform.resize(value=(8.28558, 8.28558, 8.28558),
                             orient_type='GLOBAL',
                             orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                             orient_matrix_type='GLOBAL',
                             mirror=False,
                             use_proportional_edit=False,
                             proportional_edit_falloff='SMOOTH',
                             proportional_size=1,
                             use_proportional_connected=False,
                             use_proportional_projected=False)
    
    c = Cell("cell", loc=(0, 0, 0))
    make_cell(c)
    bpy.data.objects["Grid"].select_set(False)
    bpy.data.objects["cell"].select_set(False)
    bpy.data.objects["cell"].select_set(True)
    bpy.data.objects["Grid"].select_set(True)
    bpy.ops.object.parent_set(type='OBJECT')
    bpy.context.object.instance_type = 'VERTS'
    bpy.ops.object.duplicates_make_real()
    bpy.ops.outliner.item_activate(deselect_all=True)
    bpy.data.objects["Grid"].select_set(True)
    bpy.ops.object.delete(use_global=False)'''