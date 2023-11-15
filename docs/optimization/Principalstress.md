# Principal stress net

## Definition

A mesh is in **equilibrium** if the sum of the forces at each vertex is zero, which means no bending forces appear in physical realization.
If an orthogonal quad mesh is in equilibrium, it is a discrete version of the **net of principal stress lines** in the surface, i.e. **Principal stress net** [1].

## Constraint

If a vertical load $p_i$ is applied in an unsupported vertex $v_i$, the equilibrium condition in $v_i$ reads

$$
\sum_{j=1}^4 w_{ij}(v_i-v_{ij}) - (0,0,p_i) = (0,0,0),
$$

where the sum is over 4 neighbouring connected vertices $v_{ij}$ of $v_i$, and $w_{ij}$ denotes the force density in the edge from $v_i$ to $v_j$. 

The force densities $w_{ij}$ are introduced as auxiliary variables and meet $w_{ij} = -w_{ji}$. 
The vertical load $p_i$ can depend on the mesh in which case we update it after every iteration. 

The equilibrium function refers to `DOS/archgeolab/constraints/constraints_equilibrium.py/equilibrium_constraints()`.

<!-- [![Funicular](../assets/funicular.png)](https://www.youtube.com/embed/sOzjRHIrR-s) -->

-----------------------------------------------------------
<span style="font-family:Papyrus; font-size:0.8em;">[1] Martin Kilian, Davide Pellis, Johannes Wallner, Helmut Pottmann. Material-minimizing forms and structures. ACM Trans. Graphics 36,  6 (2017): 1-12.</span>

