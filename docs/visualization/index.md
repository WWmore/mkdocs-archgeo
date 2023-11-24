# Visualization

ArchGeo provides an interactive design GUI environment, which is based on the Mayavi library. 

![File](../assets/mayavi.png)


## Mayavi

[Mayavi](https://docs.enthought.com/mayavi/mayavi/) is an interactive scientific data visualization and 3D plotting in Python. 
Its development is base on [TraitsUI](https://docs.enthought.com/traitsui/) and [vtk](https://vtk.org/).
Mayavi library provides Python code and example gallery of visulazation of points, curves, surfaces, meshes, vector field and animations.


## Read data

#### Open an obj file
``` py
## Instantiate the sample component
from opt_gui_orthonet import OrthoNet
component = OrthoNet()

## Instantiate the main geolab application
from archgeolab.archgeometry.gui_basic import GeolabGUI
GUI = GeolabGUI()

## Add the component to geolab
GUI.add_component(component)

## Open an obj file e.g. file_path=path+r'\heart.obj'
GUI.open_obj_file(file_path)

## Open another obj file
#GUI.open_obj_file(reffile)

## Start geolab main loop
GUI.start()
```


## Plot data

#### Plot mesh
``` py
from geometrylab.vtkplot.edgesource import Edges
from geometrylab.vtkplot.facesource import Faces
def plot_mesh(self, mesh, name):
    showe = Edges(mesh,color ='black',tube_radius=0.5*self.meshmanager.r,name=name+'e')
    showf = Faces(mesh,color = (77,77,77),opacity=0.1,name=name+'f')
    self.meshmanager.add([showe, showf])
```

#### Plot meshedges
``` py
def plot_mesh_edges(self):
    self.meshmanager.plot_edges(color=(157,157,157),tube_radius=0.4*self.meshmanager.r)
```

#### Hide meshedges
``` py
def hide_mesh_edges(self):
    self.meshmanager.hide_edges()
```

#### Plot colored edges
``` py
def plot_colored_edges(self, poly, data, name):
    val = np.max(data)*1.2
    self.meshmanager.plot_polyline(polyline=poly,edge_data=data,color='Blues',lut_range=[0,val],tube_radius=1.2*self.meshmanager.r,name=name)       
```

#### Plot meshfaces
``` py
def plot_mesh_faces(self):
    self.meshmanager.plot_faces(color='white',glossy=1,opacity=1)
```

#### Hide meshfaces
``` py
def hide_mesh_faces(self):
    self.meshmanager.hide_faces()
```

#### Plot colored faces
``` py
def plot_colored_faces(self, i_faces, name):
    ## chosen indices of mesh faces are plotted differently from the left
    data = np.zeros(self.mesh.F)
    data[i_faces] = 1
    self.meshmanager.plot_faces(face_data=data,color=[(8,45,130), (198,110,203)],opacity=[0,0.5],name=name)       
```

#### Plot colored points
``` py
def plot_points(self, Points, name):
    self.meshmanager.plot_glyph(points=Points,color='brg',lut_range='-:0:+',radius=2*self.meshmanager.r,name=name)   
```

#### Plot polyline
``` py
from geometrylab.geometry import Polyline
def plot_polyline(self, Points, name):
    poly = Polyline(Points,closed=False)  
    self.meshmanager.plot_polyline(poly,color=(138,43,226),glossy=1,tube_radius=0.5*self.meshmanager.r,name=name)
```

#### Plot vectors
```py
def plot_vectors(self, an, vn, name):
    self.meshmanager.plot_vectors(anchor=an,vectors=vn,position='tail',color = (255,0,255),name=name)  
```



## Save data

``` py
def save_obj(new_mesh, new_mesh_name, save_path):
    name = ('{}').format(new_mesh_name)  
    completeName = os.path.join(save_path, name)   
    new_mesh.make_obj_file(completeName)
    print('\n\n NOTE: <'+new_mesh_name+'> mesh has been saved in <'+completeName+'>\n')
```

```py 
def save_csv(vertices, name):
    ## save x, y, z list of the vertices list
    import csv
    x, y, z = np.transpose(vertices)
    with open(name+'_x.csv', 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(x)
    with open(name+'_y.csv', 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(y)
    with open(name+'_z.csv', 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(z)
    print('NOTE: .csv files has been saved!')
```