"""Extract the watertight table body from table.fbx as an ASCII STL (run with Blender).

    blender --background --python extract_surface.py

The FBX holds several objects — the solid table body ("Cube", a watertight manifold),
a thin non-watertight tabletop overlay, and some decorative/floor planes. Only the
watertight body can be voxelised reliably, so we keep that one and write table_solid.stl
next to this script. fbx_table_stress.py consumes that STL.
"""

import os
import bpy

HERE = os.path.dirname(os.path.abspath(bpy.data.filepath or __file__))
FBX = os.path.join(HERE, "table.fbx")
OUT = os.path.join(HERE, "table_solid.stl")
SOLID_BODY = "Cube"   # the watertight object in this asset

bpy.ops.wm.read_factory_settings(use_empty=True)
bpy.ops.import_scene.fbx(filepath=FBX)

for obj in list(bpy.data.objects):
    if obj.type != "MESH" or obj.name != SOLID_BODY:
        bpy.data.objects.remove(obj, do_unlink=True)

for obj in bpy.data.objects:
    obj.select_set(True)
bpy.context.view_layer.objects.active = bpy.data.objects[SOLID_BODY]

try:                                   # Blender >= 4.1
    bpy.ops.wm.stl_export(filepath=OUT, ascii_format=True, export_selected_objects=True)
except AttributeError:                 # Blender < 4.1
    bpy.ops.export_mesh.stl(filepath=OUT, ascii=True, use_selection=True)

print(f"Wrote {OUT}")
