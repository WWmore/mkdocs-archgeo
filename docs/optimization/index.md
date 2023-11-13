# Optimization

ArchGeo provides an optimization framework to solve geometry processing problems.

It utilizes a Gauss-Newton algorithm, so-called Guided Projection Algorithm as dissused in the paper ["Form-finding with polyhedral meshes made simple"](https://doi.org/10.1145/2601097.2601213), to produce constrained quadmesh models.


<details>
<summary><span style="font-weight: bold;">Abstract of the paper 'Form-finding with Polyhedral Meshes Made Simple'</span></summary>

  We solve the form-finding problem for polyhedral meshes in a way which combines form, function and fabrication; taking care of user-specified constraints like boundary interpolation, planarity of faces, statics, panel size and shape, enclosed volume, and last, but not least, cost. Our main application is the interactive modeling of meshes for architectural and industrial design. Our approach can be described as guided exploration of the constraint space whose algebraic structure is simplified by introducing auxiliary variables and ensuring that constraints are at most quadratic. Computationally, we perform a projection onto the constraint space which is biased towards low values of an energy which expresses desirable "soft" properties like fairness. We have created a tool which elegantly handles difficult tasks, such as taking boundary-alignment of polyhedral meshes into account, planarization, fairing under planarity side conditions, handling hybrid meshes, and extending the treatment of static equilibrium to shapes which possess overhanging parts.

</details>
<br>


For the Chinese explaination, please refer to the Section 3.7 in the [PhD thesis](https://www.huiwang.me/assets/pdf/hui-phd-thesis.pdf).

## Guided Projection Algorithm

Guided Projection Algorithm solves the geometry optimization problem by representing all the constraints to be equations with no more than quadratic degrees.
If there are $N$ constraints, then there exist symmetric matrices $A_i$, vectors $b_i$ and constants $c_i$ such that all the constraints can be represented in the following form

$$
\varphi_i(X) = \frac{1}{2}X^T A_i X + b_i^T X +c_i = 0, i=1,\cdots,N,
$$

where $X$ is a vector including all variables.

$X$ can be extended once more geometry constraints are added. 
Additional variables as auxiliary variables may be added into $X$ when lowering higher order (more than quadratic) equations to at-most quadratic equations, which requires geometric understanding.
The whole optimization is an iterative process until a satisfied solver $X$ is obtained. 

Suppose the solver in the last iteration is $X_n$, above equations can be linearized using Taylor expansion

$$
\varphi_i(X) \thickapprox \varphi_i(X_n) + \nabla \varphi_i(X_n)^T(X-X_n) = 0, i=1,\cdots,N,
$$

which can be written as $H \cdot X  = r$, where 

\[ H =
\left[\begin{array}{cc}
\nabla\varphi_1(X_{n})^T\\
\nabla\varphi_2(X_{n})^T\\
\vdots \\
\nabla\varphi_N(X_{n})^T
\end{array}\right]
 =
\left[\begin{array}{cc}
(A_1\cdot X+b_1)^{T}\\
(A_2\cdot X+b_2)^{T}\\
\vdots \\
(A_N\cdot X+b_N)^{T}
\end{array}\right],
\]

\[
r =
\left[\begin{array}{cc}
-\varphi_1(X_n) + \nabla\varphi_1(X_n)^TX_n\\
-\varphi_2(X_n) + \nabla\varphi_2(X_n)^TX_n\\
\vdots \\
-\varphi_N(X_n) + \nabla\varphi_N(X_n)^TX_n
\end{array}\right]
=
\left[\begin{array}{cc}
\frac{1}{2}\cdot X^T\cdot A_1\cdot X - c_1\\
\frac{1}{2}\cdot X^T\cdot A_2\cdot X - c_2\\
\vdots \\
\frac{1}{2}\cdot X^T\cdot A_{N}\cdot X -c_N
\end{array}\right].
\]

We do not solve $H \cdot X  = r$ directly, since this linear system is typically underdetermined.
Usually there is enough solution space, and we add fairness term (explained later) and 
a controlled solver distance from the previous value $X_n$ as regularizers.
Then we solve 

\[
\|HX - r\|^2 + \|KX - s\|^2 + \epsilon^2\|X - X_n\|^2 \to min ,
\]

where $\|KX - s\|$ and $\|X - X_n\|$ are regularizers, and $\epsilon=0.001$ for almost all the optimization cases.

Furthermore, we only solve the linear system

\[
(H^T H +  K^T K + \epsilon^2 I)X = H^T r +  K^T s + \epsilon^2 X_N,
\]

which can be solved fast by the [SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.spsolve.html#scipy.sparse.linalg.spsolve) sparse matrix solver.
In this codebase, we use [PyPardiso](https://pypi.org/project/pypardiso/) to increase the computing.



### Hard constraints


List of variables.


### Fairness term

The fairness term $\|KX - s\|$ is a simple and efficient soft constraint that plays a crucial role in smoothing out polylines in visual appearance.
Soft constraints typically have smaller and controllable weights compared to hard constraints. 

While the smoothness of vertices in triangular mesh is achieved using Laplacian fairness, the focus here is on the smoothness of polylines.

For each vertex $v$ of valence 4 in a quad mesh, we consider 4 neighbouring connected vertices be $v_i(i=1,\cdots,4)$ in clockwise order and use vanishing second order differences to represent the smoothness of any three consecutive vertices

$$
v_1 -2 v + v_3=0, v_2 -2 v + v_4=0.
$$

### Sparse matrix construction

