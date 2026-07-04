# CalculiX++ examples

Runnable examples that exercise the library through its Python bindings. Build the
Python module first, then point `PYTHONPATH` at it:

```bash
cmake --build build -j
export PYTHONPATH=$PWD/build/python
```

| Example | What it shows |
|---|---|
| [`table_stress/`](table_stress/) | Linear-static stress of a 3D table under a 100 kg central object, with von Mises visualisation (matplotlib + ParaView `.vtk`). Two variants: a parametric mesh, and the real `table.fbx` art asset voxelised into a solid FE mesh. |
| [`Antenna_1/`](Antenna_1/) | Stress of a tapered cell-tower antenna mast (`Antenna.obj`) under its own weight + a 300 kg top payload — a hollow thin-walled tube (surface-raster fill → eroded core), mass-calibrated, direct solver. |

[`femkit.py`](femkit.py) is the shared toolkit behind these: mesh IO (STL/OBJ),
voxelisation (winding-number inside test, or surface-fill for thin shells), the C3D8
hex-mesh builder + deck assembler, and the von Mises / VTK / render helpers.
