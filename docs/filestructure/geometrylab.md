# geometrylab

``` { .sh .no-copy }
.
├─ geometrylab/
│  ├─ geometry/
│  │  ├─ meshpy.py
│  │  ├─ meshprimitives.py
│  │  ├─ meshutilities.py
│  │  ├─ frame.py
│  │  ├─ polyline.py
│  │  ├─ circle.py
│  │  └─bspline.py
│  ├─ gui/
│  │  ├─ geolabcomponent.py
│  │  ├─ geolabgui.py
│  │  ├─ geolabmesh.py
│  │  ├─ geolabpoints.py
│  │  ├─ geolabscene.py
│  │  ├─ handler.py
│  │  ├─ multiscenemanager.py
│  │  ├─ scenemanager.py
│  │  └─ tools.py
│  ├─ optimization/
│  │  ├─ gridshell.py
│  │  ├─ guidedprojection.py
│  │  ├─ guidedprojectionbase.py
│  │  └─ combnormals.py
│  ├─ vtkplot/
│  │  ├─ scenemanager.py
│  │  ├─ plotmanager.py
│  │  ├─ meshplotmanager.py
│  │  ├─ plotutilities.py
│  │  ├─ glyphs.py
│  │  ├─ glyphsource.py
│  │  ├─ pointsource.py
│  │  ├─ pointsplotmanager.py
│  │  ├─ edgesource.py
│  │  ├─ polylinesource.py
│  │  ├─ bsplineplotmanager.py
│  │  ├─ vectorsource.py
│  │  ├─ vector3dsource.py
│  │  ├─ meshvectorsource.py
│  │  ├─ facesource.py
│  │  ├─ toolbar.py
│  │  ├─ selector.py
│  │  ├─ viewer.py
│  │  └─ check.py
│  ├─ utilities/
│  │  └─ utilities.py
│  ├─ fitting/
│  │  ├─ jetfitting.py
│  │  ├─ linearregression.py
│  │  ├─ cluster.py
│  │  └─ bspline.py
│  └─ test/
│     └─ paneling.py
└─ archgeolab/
```

## geometry

|  Files              | Functions                               |
| ------------------- | --------------------------------------- |
| `meshpy.py`         | very basic mesh connection by half-edges|
| `meshprimitives.py` | plane, cylinder, sphere, torus, arrows  |
| `meshutilities.py`  | clean boundary, add noise               |
| `frame.py`          | frame vectors $(origin,e_1,e_2,e_3)$    |
| `polyline.py`       | construction of polyline                |
| `circle.py`         | construction of circles                 |
| `bspline.py`        | construction of B-spline                |


## gui

|  Files                 | Functions                                             |
| ---------------------- | ----------------------------------------------------- |
| `geolabcomponent.py`   | GeolabComponent                                       |
| `geolabgui.py`         | Toolbar button settings                               |
| `geolabpoints.py`      | points plotting settings                              |
| `geolabscene.py`       | basic Mayavi GUI viewer setting                       |
| `handler.py`           | Handler                                               |
| `multiscenemanager.py` | multi-scenes settings                                 |
| `scenemanager.py`      | callbacks of selection, hide, adding, removing        |
| `tools.py`             | loads, remesh, corner tolerance, save mesh in Toolbar |


## optimization

|  Files                    | Functions                                  |
| ------------------------- | ------------------------------------------ |
| `gridshell.py`            | basic geometric constraints settings       |
| `guidedprojection.py`     | Guided Projection algorithm settings       |
| `guidedprojectionbase.py` | basice algorithm settings                  |
| `combnormals.py`          | constraints for equilibrium shellstructure |

## vtkplot

|  Files                  | Functions                                        |
| ----------------------- | ------------------------------------------------ |
│ `scenemanager.py`       | scene manager class                              |
│ `plotmanager.py`        | plot manager class                               |
│ `meshplotmanager.py`    | mesh plottint manager class                      |
│ `plotutilities.py`      | plotting attributes                              |  
│ `glyphs.py`             | plot arrows, spheres, pipes, discs, cones, rings |
│ `glyphsource.py`        | Glyph3D                                          |
│ `pointsource.py`        | Points                                           |
│ `pointsplotmanager.py`  | plot points                                      |
│ `edgesource.py`         | plot edge plot                                   |
│ `polylinesource.py`     | plot polylines                                   |
│ `bsplineplotmanager.py` | plot B-spline curves                             |
│ `vectorsource.py`       | vector settings                                  |
│ `vector3dsource.py`     | plot vectors                                     |
│ `meshvectorsource.py`   | plot vectors as meshes                           |
│ `facesource.py`         | plot faces                                       |
│ `toolbar.py`            | Toolbar settings                                 |
│ `selector.py`           | interactive selections                           |
│ `viewer.py`             | GUI viewer layout                                |
│ `check.py`              | check the changes of mesh connectivity           |

## utilities

|  Files          | Functions            |
| --------------- | -------------------- |
| `utilities.py`  | array data structure |

## fitting

|  Files                | Functions                                     |
| --------------------- | --------------------------------------------- |
| `jetfitting.py`       | quadratic fitting to get principal curvatures |
| `linearregression.py` | linear regression, sphere fitting             |
| `cluster.py`          | K-means cluster of given points               |
| `bspline.py`          | construction of B-spline curves               |

## test

|  Files          | Functions                                     |
| --------------- | --------------------------------------------- |
| `paneling.py`   | basic frame integrated with algorithm and GUI |
