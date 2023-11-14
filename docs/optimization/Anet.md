# A-net



## Definition


## Constraint

| Variable     | Symbol          | Number     |
| ------------ | --------------- | ---------- |
| `vertices`   | $v_i \in R^3$   | $3|V|$     |



| $H \cdot X  = r$  | Representation                                                                               |
| ----------------- | -------------------------------------------------------------------------------------------- |
| `H: shape`        | $(N,3|V|)$                                                                                   |
| `H: row`          | np.tile(np.arange($N$),12)                                                                   |
| `H: col`          | $[col_1,col_2,col_3,col_4]$                                                                  |
| `H: data`         | $2[X_n[col_1]-X_n[col_3],X_n[col_4]-X_n[col_2],X_n[col_3]-X_n[col_1],X_n[col_2]-X_n[col_4]]$ |
| `r`               | $(X_n[col_1] - X_n[col_3])^2 - (X_n[col_2] - X_n[col_4])^2$                                  |


The function is `DOS/archgeolab/constraints/constraints_basic.py/con_equal_length()`.


For the Chinese explaination, please refer to the Section 3.6 in the [PhD thesis](https://www.huiwang.me/assets/pdf/hui-phd-thesis.pdf).


## Minimal net



[![Anet](../assets/anet.png)](https://www.youtube.com/embed/KQbJ2e_Ow7M)






