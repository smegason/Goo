# Divides 2 types of cells

import bpy
from goo import goo
from importlib import reload
reload(goo)

goo.setup_world()

master_coll = bpy.context.view_layer.layer_collection

collection = bpy.context.blend_data.collections.new(name='type1')
bpy.context.collection.children.link(collection)
layer_collection = bpy.context.view_layer.layer_collection.children[collection.name]
bpy.context.view_layer.active_layer_collection = layer_collection
bpy.app.handlers.frame_change_post.clear()
cell = goo.Cell(name_string="cell1", loc=(0, 0, 0))
goo.make_cell(cell)

bpy.context.view_layer.active_layer_collection = master_coll

collection = bpy.context.blend_data.collections.new(name='type2')
bpy.context.collection.children.link(collection)
layer_collection = bpy.context.view_layer.layer_collection.children[collection.name]
bpy.context.view_layer.active_layer_collection = layer_collection
cell = goo.Cell("cell_(1, 3, 3)", loc=(1, 3, 3))
goo.make_cell(cell)

handlers = goo.handler_class()
handlers.active_cell_types = ["type1", "type2"]
handlers.set_division_rate("type1", 2)
handlers.set_division_rate("type2", 3)

bpy.app.handlers.frame_change_post.append(handlers.div_handler)
