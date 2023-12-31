# Home

Oirginal source code visits [DOS](https://github.com/WWmore/DOS). 
It contains the implementation associated with the paper ["Discrete Orthogonal Structures"](https://doi.org/10.1016/j.cag.2023.05.024). Please cite the paper if you use this code in your project.

ArchGeo is a Python library for processing manifold and orientable meshes based on half-edge data structures.
It provides not only a GUI to visualize the 3D models, including point cloud, curves, surfaces, meshes and vector field, but also an optimization framework to solve geometry processing problems. 
ArchGeo is a powerful tool to do interactive design.

ArchGeo offers four core features:

* Visualizing any point cloud, curves, surfaces, meshes and vector field;
* Geometrically processing any manifold and orientable mesh;
* Solving geometric optimization problems;
* Fast interactive design.

## Triangular meshes
Triangular meshes are the most commonly used meshes in geometry processing and computer graphics, since they can explicitly represent the model shape.
Other good properties include representing flat surfaces efficiently, allocating more faces to the areas with fine detail, easily attaching  data (e.g. RGB colors, texture coordinates, normal vectors, etc.) on vertices and interpolating over the whole surface, and simplier way to store triangular faces in data structure compared with other polygon faces. 

**Regular** triangular mesh has all the vertices of valence 6.
Three angles of 60 degrees corresponds to equilateral triangles, which are the most visually preferred but very limited ones in freeform shape. 
Skinny triangles with more extreme angles are not ideal in both geometry processing and industrial applications.

## Why quad meshes?
Quad meshes are extensively used in animation, architectural designs, and industrial analysis.
They are also majorly discussed objects with fruitful discrete theories in Discrete Differential Geometry (DDG) [1].

Each vertex of a **regular** quad mesh has valence 4.

Compared with triangular meshes, quad meshes have the following good features:

* (usually) small number of valence at vertex star, which leads to less supporting weight in structural application
* can have torsion-free structure
* have nice geometry theory in discretization of (smooth) differential geometry

## Half-edge data structure
Our core processing object is (quad) mesh. We use the popular [half-edge data structure](https://cs184.eecs.berkeley.edu/sp20/article/17/an-introduction-to-half-edge-dat) to represent the geometric relations between vertices.

## Architectural Geometry

The research area of **Architectural Geometry** traces back to the research paper ["Geometric modeling with conical meshes and developable surfaces"](https://doi.org/10.1145/1141911.1141941) [2] in 2006 by Prof. Helmut Pottmann and colleagues. 

In 2007, the book ["Architectural geometry"](http://www.architecturalgeometry.at/) [3] was published, and since then, numerous high-quality research papers have been presented in top conferences and journals every year. 

In 2008, the first [Advances in Architectural Gometry](https://www.architecturalgeometry.org/aag08/) (AAG) conference was launched in Vienna, followed by many follow-up application papers in architecture, structure, CAD and CAM.

In 2015, a survey paper titled ["Architectural geometry"](https://doi.org/10.1016/j.cag.2014.11.002) [4] presented the fruitful achievements of this research area and also listed some promising research directions. 

This area, which combines theory and application, is still active with many researchers dedicated to advancing its development.


-----------------------------------------------------------
<span style="font-family:Papyrus; font-size:0.8em;">[1] Alexander Bobenko, Suris Yuri. 2008. Discrete differential geometry: Integrable structure. Vol. 98. American Mathematical Soc.</span>

<span style="font-family:Papyrus; font-size:0.8em;">[2] Yang Liu, Helmut Pottmann, Johannes Wallner, Yongliang Yang, Wenping Wang. 2006. Geometric modeling with conical meshes and developable surfaces. ACM Trans. Graphics 25, 3, 681--689.</span>

<span style="font-family:Papyrus; font-size:0.8em;">[3] Helmut Pottmann, Andreas Asperl, Axel Kililan. 2007. Architectural geometry. Bentley Institute Press.</span>

<span style="font-family:Papyrus; font-size:0.8em;">[4] Helmut Pottmann, Michael Eigensatz, Amir Vaxman, Johannes Wallner. 2015. Architectural geometry. Computers & Graphics. 47, 145--164.</span>

<!-- [1] Yang Liu, Helmut Pottmann, Johannes Wallner, Yongliang Yang, Wenping Wang. 2006. Geometric modeling with conical meshes and developable surfaces. ACM Trans. Graphics 25, 3, 681--689.

[2] Helmut Pottmann, Andreas Asperl, Axel Kililan. 2007. Architectural geometry. Bentley Institute Press.

[3] Helmut Pottmann, Michael Eigensatz, Amir Vaxman, Johannes Wallner. 2015. Architectural geometry. Computers & Graphics. 47, 145--164. -->


