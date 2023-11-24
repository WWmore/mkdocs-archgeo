# Plot Mesh

#### Plot mesh
``` py
from geometrylab.vtkplot.edgesource import Edges
from geometrylab.vtkplot.facesource import Faces
def plot_mesh(self, mesh, name):
    showe = Edges(mesh,color ='black',tube_radius=0.5*self.meshmanager.r,name=name+'e')
    showf = Faces(mesh,color = (77,77,77),opacity=0.1,name=name+'f')
    self.meshmanager.add([showe, showf])
```

## Quad Planarity


## Gaussian Curvature


## Mean Curvature


## Curve network


## Isolines



## Support Structures
