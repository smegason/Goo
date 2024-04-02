import bpy


class BlenderObject:
    def __init__(self, obj):
        self.obj = obj

    @property
    def name(self):
        return self.obj.name

    @name.setter
    def name(self, name):
        self.obj.name = name

    @property
    def loc(self):
        return self.obj.location

    @loc.setter
    def loc(self, loc):
        self.obj.location = loc


class Force(BlenderObject):
    def __init__(self, obj):
        super(Force, self).__init__(obj)
        self.type = "force"

    def enable(self):
        if not self.obj.field:
            with bpy.context.temp_override(active_object=self.obj):
                bpy.ops.object.forcefield_toggle()
        else:
            self.obj.field.type = "FORCE"

    def disable(self):
        if self.obj.field:
            self.obj.field.type = "NONE"

    def enabled(self):
        return self.obj.field and self.obj.field.type == "FORCE"

    @property
    def strength(self):
        return self.obj.field.strength

    @strength.setter
    def strength(self, strength):
        self.obj.field.strength = strength

    @property
    def falloff(self):
        return self.obj.field.falloff

    @falloff.setter
    def falloff(self, falloff):
        self.obj.field.falloff_power = falloff

    @property
    def min_dist(self):
        if self.obj.field.use_min_distance:
            return self.obj.field.min_dist
        return None

    @min_dist.setter
    def min_dist(self, min_dist):
        if min_dist is None:
            self.obj.field.use_min_distance = False
        else:
            self.obj.field.use_min_distance = True
            self.obj.field.distance_min = min_dist

    @property
    def max_dist(self):
        if self.obj.field.use_max_distance:
            return self.obj.field.max_dist
        return None

    @max_dist.setter
    def max_dist(self, max_dist):
        if max_dist is None:
            self.obj.field.use_max_distance = False
        else:
            self.obj.field.use_max_distance = True
            self.obj.field.distance_max = max_dist


class AdhesionForce(Force):
    def __init__(self, obj):
        super(AdhesionForce, self).__init__(obj)
        self.type = "adhesion"

    @property
    def strength(self):
        return -self.obj.field.strength

    @strength.setter
    def strength(self, strength):
        self.obj.field.strength = -strength


def make_force(name, loc, strength, falloff=0, min_dist=None, max_dist=None):
    obj = bpy.data.objects.new(name, None)
    obj.location = loc

    force = Force(obj)
    force.enable()
    force.strength = strength
    force.falloff = falloff
    force.min_dist = min_dist
    force.max_dist = max_dist

    bpy.context.collection.objects.link(obj)
    return force


def make_homo_adhesion(cell, strength=2000):
    homo_force = AdhesionForce(cell)
    homo_force.enable()

    homo_force.strength = strength
    homo_force.min_dist = 0.6
    homo_force.max_dist = 1.4
    homo_force.falloff = 0.5
    return homo_force
