import bpy

def add_material(mat):
    bpy.data.materials.new(name = mat.name)
    bpy.data.materials[mat.name].node_tree.nodes["Mix Shader"].inputs[0].default_value = 0.079
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[0].default_value = (0.00749999, 0.020955, 0.3, 1)
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[1].default_value = 0
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[2].default_value[0] = 1.1
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[2].default_value[1] = 0.2
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[2].default_value[2] = 0.2
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[3].default_value = (0.8, 0.8, 0.8, 1)
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[4].default_value = 0.136
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[5].default_value = 0.5
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[6].default_value = 0.555
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[7].default_value = 0.318
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[8].default_value = 0.041
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[9].default_value = 0.048
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[10].default_value = 0.052
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[11].default_value = 0.03
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[12].default_value = 0.114
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[13].default_value = 0.123
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[14].default_value = 1.45
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[15].default_value = 0.882
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[16].default_value = 0
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[17].default_value = (0, 0, 0, 1)
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[18].default_value = 1
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[19].default_value = 0.414
    bpy.data.materials[mat.name].node_tree.nodes["Hue Saturation Value"].inputs[0].default_value = 0.8
    bpy.data.materials[mat.name].node_tree.nodes["Hue Saturation Value"].inputs[1].default_value = 2
    bpy.data.materials[mat.name].node_tree.nodes["Hue Saturation Value"].inputs[2].default_value = 2
    bpy.data.materials[mat.name].node_tree.nodes["Hue Saturation Value"].inputs[3].default_value = 1
    bpy.data.materials[mat.name].node_tree.nodes["Hue Saturation Value"].inputs[4].default_value = (0.8, 0.8, 0.8, 1)
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[1].default_value = 0
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[2].default_value[0] = 1
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[2].default_value[1] = 0.2
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[2].default_value[2] = 0
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[3].default_value = (0.8, 0.8, 0.8, 1)
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[4].default_value = 0
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[5].default_value = 0.5
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[6].default_value = 0
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[7].default_value = 0.482
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[8].default_value = 0
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[9].default_value = 0
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[10].default_value = 0
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[11].default_value = 0
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[12].default_value = 0
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[13].default_value = 0
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[14].default_value = 1.45
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[15].default_value = 1
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[16].default_value = 0
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[17].default_value = (0, 0, 0, 1)
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[18].default_value = 1
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[19].default_value = 0.555
    bpy.context.object.active_material.blend_method = 'BLEND'
    bpy.context.object.active_material.shadow_method = 'NONE'
    bpy.context.object.active_material.refraction_depth = 0
    bpy.context.object.active_material.pass_index = 0
    bpy.context.object.active_material.diffuse_color = (0.8, 0.8, 0.8, 1)
    bpy.context.object.active_material.metallic = 0
    bpy.context.object.active_material.roughness = 0.4

class Material():
    def __init__(self):
        self.name = "Cell Material"
        self.BDSF_1_fac = 0.079
        self.BDSF_1_color = (0.00749999, 0.020955, 0.3, 1)
        self.BDSF_1_subsurface = 0
        self.BDSF_1_subsurf_radius = (1.1, 0.2, 0.2)
        self.BDSF_1_subsurf_color = (0.8, 0.8, 0.8, 1)
        self.BDSF_1_metallic = 0.136
        self.BDSF_1_specular = 0.5
        self.BDSF_1_specular_tint = 0.555
        self.BDSF_1_roughness = 0.318
        self.BDSF_1_anisotropic = 0.041
        self.BDSF_1_anisotropic_rot = 0.048
        self.BDSF_1_sheen = 0.052
        self.BDSF_1_sheen_tint = 0.03
        self.BDSF_1_clearcoat = 0.114
        self.BDSF_1_clear_rough = 0.123
        self.BDSF_1_IOR = 1.45
        self.BDSF_1_transmission = 0.882
        self.BDSF_1_transmission_rough = 0
        self.BDSF_1_emission_color = (0, 0, 0, 1)
        self.BDSF_1_emission_strength = 1
        self.BDSF_1_alpha = 0.414

        self.BDSF_2_hue = 0.8
        self.BDSF_2_saturation = 2
        self.BDSF_2_value = 2
        self.BDSF_2_fac = 1
        self.BDSF_2_color = (0.8, 0.8, 0.8, 1)
        self.BDSF_2_subsurface = 0
        self.BDSF_2_subsurf_radius = (1.1, 0.2, 0.2)
        self.BDSF_2_subsurf_color = (0.8, 0.8, 0.8, 1)
        self.BDSF_2_metallic = 0
        self.BDSF_2_specular = 0.5
        self.BDSF_2_specular_tint = 0
        self.BDSF_2_roughness = 0.482
        self.BDSF_2_anisotropic = 0
        self.BDSF_2_anisotropic_rot = 0
        self.BDSF_2_sheen = 0
        self.BDSF_2_sheen_tint = 0
        self.BDSF_2_clearcoat = 0
        self.BDSF_2_clear_rough = 0
        self.BDSF_2_IOR = 1.45
        self.BDSF_2_transmission = 1
        self.BDSF_2_transmission_rough = 0
        self.BDSF_2_emission_color = (0, 0, 0, 1)
        self.BDSF_2_emission_strength = 1
        self.BDSF_2_alpha = 0.555

        self.blend_method = 'BLEND'
        self.shadow_method = 'NONE'
        self.refraction_depth = 0
        self.pass_index = 0
        self.diffuse_color = (0.8, 0.8, 0.8, 1)
        self.metallic = 0
        self.roughness = 0.4

def setup_world():
    

