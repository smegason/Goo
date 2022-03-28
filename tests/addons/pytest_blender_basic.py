import bpy


bl_info = {
    "name": "Move X Axis",
    "category": "Object",
}


class PytestBlenderObjectMoveX(bpy.types.Operator):
    """My Object Moving Script"""

    bl_idname = "object.move_x"
    bl_label = "Move X by One"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        scene = context.scene
        for obj in scene.objects:
            obj.location.x += 1.0

        return {"FINISHED"}


def register():
    bpy.utils.register_class(PytestBlenderObjectMoveX)


def unregister():
    bpy.utils.unregister_class(PytestBlenderObjectMoveX)


if __name__ == "__main__":
    register()