# Optimization

ArchGeo provides an optimization framework to solve geometry processing problems.

It utilizes a Gauss-Newton algorithm, so-called Guided Projection Algorithm as dissused in the paper ["Form-finding with polyhedral meshes made simple"](https://doi.org/10.1145/2601097.2601213), to produce constrained quadmesh models.


<details>
<summary><span style="font-weight: bold;">Abstract of the paper 'Form-finding with Polyhedral Meshes Made Simple'</span></summary>

  We solve the form-finding problem for polyhedral meshes in a way which combines form, function and fabrication; taking care of user-specified constraints like boundary interpolation, planarity of faces, statics, panel size and shape, enclosed volume, and last, but not least, cost. Our main application is the interactive modeling of meshes for architectural and industrial design. Our approach can be described as guided exploration of the constraint space whose algebraic structure is simplified by introducing auxiliary variables and ensuring that constraints are at most quadratic. Computationally, we perform a projection onto the constraint space which is biased towards low values of an energy which expresses desirable "soft" properties like fairness. We have created a tool which elegantly handles difficult tasks, such as taking boundary-alignment of polyhedral meshes into account, planarization, fairing under planarity side conditions, handling hybrid meshes, and extending the treatment of static equilibrium to shapes which possess overhanging parts.

</details>
<br>


## Guided Projection Algorithm

