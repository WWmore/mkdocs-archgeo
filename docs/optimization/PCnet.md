# Principal curvature net

Chinese explaination refers to the Section 3.5 in the [PhD thesis](https://www.huiwang.me/assets/pdf/hui-phd-thesis.pdf).

## Definition

**Principal curvature net** is a curve network formed by principal curvature lines, whose directions at each point follow the maximum and minimum curvature of the surface.
Any net that is conjugate and orthogonal is a principal curvature net [1]. 
The corresponding discretization is orthogonal planar quad mesh, i.e. orthogonal PQ mesh.

Other discretizations include circular mesh [1] and conical mesh [2].

## Constraint of PQ mesh

The representation of planar quad faces is similar to the planar vertex stars, but adding additional quad face normals $f_n$ as auxiliary variables.
$f_n$ are unit normals orthogonal to vectors $v_{i+1}-v_i (i=1,2,3)$:

$$
f_n ^2 = 1, f_n \cdot (v_{i+1}-v_i) = 0, i=1,2,3.
$$

The number of all variables is $|X| = 3|V| + 3|F|$ and the number of hard constraints is $N = |F| + 4|F|$.

| Variable     | Symbol       | Number           |
| ------------ | ------------ | ---------------- |
| `vertices`   | $v  \in R^3$ | $3\vert V \vert$ |
| `normals`    | $f_n\in R^3$ | $3\vert F \vert$ |

The function for PQ mesh is `DOS/archgeolab/constraints/constraints_basic.py/con_planarity_constraints()`.

<!-- [![PQ](../assets/pq.png)](https://www.youtube.com/embed/m-CFC0XZ488) -->

-----------------------------------------------------------
<span style="font-family:Papyrus; font-size:0.8em;">[1] Alexander Bobenko, Suris Yuri. 2008. Discrete differential geometry: Integrable structure. Vol. 98. American Mathematical Soc.</span>

<span style="font-family:Papyrus; font-size:0.8em;">[2] Yang Liu, Helmut Pottmann, Johannes Wallner, Yongliang Yang, Wenping Wang. 2006. Geometric modeling with conical meshes and developable surfaces. ACM Trans. Graphics 25, 3, 681--689.</span>