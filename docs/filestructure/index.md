# File Structure


``` { .sh .no-copy }
.
├─ geometrylab/
│  ├─ geometry/
│  │  ├─ meshpy.py
│  │  ├─ meshprimitives.py
│  │  ├─ meshutilities.py
│  │  ├─ frame.py
│  │  ├─ points.py
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
├─ archgeolab/
│  ├─ readfile_orthonet.py
│  ├─ guidedprojection_orthonet.py
│  ├─ opt_gui_orthonet.py
│  ├─ archgeometry/
│  │  ├─ quadrings.py
│  │  ├─ gui_basis.py
│  │  ├─ gridshell_new.py
│  │  ├─ orient.py
│  │  ├─ orthogonalVectors.py
│  │  ├─ curves.py
│  │  └─ conicSection.py
│  └─ constraints/
│  │  ├─ constraints_basis.py
│  │  ├─ constraints_net.py
│  │  ├─ constraints_fairness.py
│  │  ├─ constraints_glide.py
      └─ constraints_equilibrium.py
```

<!-- ## geometrylab

## archgeolab
### readfile_orthonet
### opt_gui_orthonet
### guidedprojection_orthonet

### archgeometry
#### quadrings
#### orient
#### curves
#### conicSection
#### orthogonalVectors
#### gridshell_new
#### gui_basic
### constraints
#### constraints_basic
#### constraints_net
#### constraints_fairness
#### constraints_glide
#### constraints_equilibrium -->