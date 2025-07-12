# Home

ArchGeo is a Python-based geometry processing library equipped with integrated visualization and optimization capabilities. 
It is designed to facilitate efficient design and construction by leveraging the knowledge of Discrete Differential Geometry (DDG). 
Utilizing a half-edge data structure, ArchGeo adeptly handles a variety of meshes. The library not only provides visualization and rendering for geometric objects such as point clouds and curves but also excels in mesh optimization design. 
For visualization tasks, it employs the powerful Mayavi3D library, while its optimization processes, particularly for quadrilateral meshes, are driven by a Gauss-Newton algorithm coupled with an efficient solver.

To gain a clear understanding of how ArchGeo operates, we recommend exploring the [DOS](https://github.com/WWmore/DOS) project, which is closely related to the published research paper ["Discrete Orthogonal Structures"](https://doi.org/10.1016/j.cag.2023.05.024) [1]. This project delves into the comprehensive implementations and showcases the capabilities of ArchGeo. It covers the relevant theory within DDG and optimized implementations, including discrete developable surfaces, discrete minimal surfaces, discrete constant mean curvature surfaces, principal meshes, and principal stress meshes, all based on the analysis and construction of discrete orthogonal quad meshes.

Additional features of ArchGeo will be open-sourced in the near future. If you utilize this library in your projects, please cite the paper [1].


**ArchGeo offers four core features**:

1. **Comprehensive Visualization**: ArchGeo provides detailed and intuitive visualizations of various geometric objects, including point clouds, curves, surfaces, meshes, and vector fields.

2. **Advanced Geometric Processing**: ArchGeo processes manifold and orientable meshes, with specialized tools for analyzing and manipulating quad meshes.

3. **Powerful Optimization Solutions**: ArchGeo leverages Discrete Differential Geometry (DDG) to solve geometric optimization problems for quad meshes, delivering efficient and effective solutions to complex geometric challenges.

4. **Rapid Interactive Design**: ArchGeo supports real-time editing and interactive design, allowing for immediate adjustments and instant visual feedback to enhance the creative process.

<!-- * Visualizing any point cloud, curves, surfaces, meshes and vector field;
* Geometrically processing any manifold and orientable mesh;
* Solving geometric optimization problems;
* Fast interactive edit and design. -->

## Triangular meshes
Triangular meshes are the most commonly used meshes in geometry processing and computer graphics, since they can explicitly represent the model shape.
Other good properties include representing flat surfaces efficiently, allocating more faces to the areas with fine detail, easily attaching data (e.g. RGB colors, texture coordinates, normal vectors, etc.) on vertices and interpolating over the whole surface, and simplier way to store triangular faces in data structure compared with other polygon faces. 

**Regular triangular mesh** has all the vertices of valence 6.
Three angles of 60 degrees corresponds to equilateral triangles, which are the most visually preferred but very limited ones in freeform shape. 
Skinny triangles with more extreme angles are not ideal in both geometry processing and industrial applications.

## Why quad meshes?
Quad meshes are extensively used in animation, architectural designs, and industrial analysis.
TThey are also a major focus in Discrete Differential Geometry (DDG), where they have given rise to a wealth of discrete theories [2].

Each vertex of a **regular quad mesh** has valence 4.

Compared with triangular meshes, quad meshes have the following good features:

* (usually) small number of valence at vertex star, which leads to less supporting weight in structural application
* can have torsion-free structure
* have nice geometry theory in discretization of (smooth) differential geometry

## Half-edge data structure
ArchGeo processes manifold and orientable meshes, with a particular focus on quad meshes, utilizing the [half-edge data structure](https://cs184.eecs.berkeley.edu/sp20/article/17/an-introduction-to-half-edge-dat) to represent the geometric relationships between vertices. This structure enables efficient traversal and manipulation of mesh elements, facilitating a wide range of geometry processing tasks.

## Architectural Geometry

The research area of **Architectural Geometry** originated from the 2006 paper ["Geometric modeling with conical meshes and developable surfaces"](https://doi.org/10.1145/1141911.1141941) [3] by Prof. Helmut Pottmann and colleagues. 

In 2007, the book ["Architectural geometry"](http://www.architecturalgeometry.at/) [4] was published, and since then, numerous high-quality research papers have been presented in top conferences such as ACM SIGGRAPH (Asia) and AAD, as well as in leading journals like ACM TOG, CAD, and CAGD.

In 2008, the first [Advances in Architectural Gometry](https://www.architecturalgeometry.org/aag08/) (AAG) conference was launched in Vienna, followed by many follow-up application papers in architecture, structure, CAD and CAM.

In 2015, a survey paper titled ["Architectural geometry"](https://doi.org/10.1016/j.cag.2014.11.002) [5] presented the fruitful achievements of this research area and also listed some promising research directions. 

This area, which combines theory and application, remains active with many researchers dedicated to advancing its development.


-----------------------------------------------------------


<span style="font-family:Papyrus; font-size:0.8em;">[1] Felix Dellinger, Xinye Li, Hui Wang. 2023. Discrete orthogonal structures. Computers & Graphics. 114, 126-137.</span>

<span style="font-family:Papyrus; font-size:0.8em;">[2] Alexander Bobenko, Suris Yuri. 2008. Discrete differential geometry: Integrable structure. Vol. 98. American Mathematical Soc.</span>

<span style="font-family:Papyrus; font-size:0.8em;">[3] Yang Liu, Helmut Pottmann, Johannes Wallner, Yongliang Yang, Wenping Wang. 2006. Geometric modeling with conical meshes and developable surfaces. ACM Trans. Graphics 25, 3, 681-689.</span>

<span style="font-family:Papyrus; font-size:0.8em;">[4] Helmut Pottmann, Andreas Asperl, Axel Kililan. 2007. Architectural geometry. Bentley Institute Press.</span>

<span style="font-family:Papyrus; font-size:0.8em;">[5] Helmut Pottmann, Michael Eigensatz, Amir Vaxman, Johannes Wallner. 2015. Architectural geometry. Computers & Graphics. 47, 145-164.</span>

<!-- [1] Yang Liu, Helmut Pottmann, Johannes Wallner, Yongliang Yang, Wenping Wang. 2006. Geometric modeling with conical meshes and developable surfaces. ACM Trans. Graphics 25, 3, 681--689.

[2] Helmut Pottmann, Andreas Asperl, Axel Kililan. 2007. Architectural geometry. Bentley Institute Press.

[3] Helmut Pottmann, Michael Eigensatz, Amir Vaxman, Johannes Wallner. 2015. Architectural geometry. Computers & Graphics. 47, 145--164. -->


