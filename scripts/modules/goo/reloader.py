import sys
import bpy


def reset_modules():
    to_delete = []
    for modname, _ in sys.modules.items():
        if modname.startswith("goo"):
            to_delete.append(modname)
    for modname in to_delete:
        del sys.modules[modname]


def reset_scene():
    bpy.app.handlers.frame_change_pre.clear()
    bpy.app.handlers.frame_change_post.clear()
    bpy.context.scene.frame_set(1)
    bpy.context.scene.frame_start = 1
    bpy.context.scene.frame_end = 250

    if bpy.context.active_object and bpy.context.active_object.mode != "OBJECT":
        bpy.ops.object.mode_set(mode="OBJECT")

    for obj in bpy.context.scene.objects:
        if obj.type in ["CAMERA", "LIGHT"]:
            continue
        bpy.data.objects.remove(obj, do_unlink=True)

    for mesh in bpy.data.meshes:
        bpy.data.meshes.remove(mesh, do_unlink=True)

    for col in bpy.data.collections:
        bpy.data.collections.remove(col)
