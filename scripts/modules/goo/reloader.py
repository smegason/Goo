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
    bpy.app.handlers.frame_change_post.clear()
    bpy.context.scene.frame_set(1)
    try:
        bpy.ops.object.mode_set(mode="OBJECT")
    except:
        pass
    for obj in bpy.context.scene.objects:
        bpy.data.objects.remove(obj, do_unlink=True)

    for mesh in bpy.data.meshes:
        bpy.data.meshes.remove(mesh, do_unlink=True)
