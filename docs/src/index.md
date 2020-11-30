# MultiObjectiveProblems.jl Documentation 

This package provides some test problems for multiobjective optimization.
All problems are considered to be of the form
```math
\min_{\mathbf x\in \mathcal{X}}
\begin{bmatrix}
f_1(\mathbf x)\\ \vdots \\f_k(\mathbf x)
\end{bmatrix}
```
For now only unconstrained and finitely box constrained problems are provided, 
i.e. ``\mathcal{X} = ℝ^n`` or 
``\mathcal{X} = \{\mathbf x \in ℝ^n : \mathbf l \le \mathbf x \le \mathbf u \}``.

Each problem is a subtype of `MOP`. 
```@docs 
MOP
```

Consider for example the canonical 2-Parabola example
```math 
\min_{\mathbf x\in ℝ^2} 
\begin{bmatrix}
(x_1 - 1)^2 + (x_2 - 1)^2 \\
(x_1 + 1)^2 + (x_2 + 1)^2
\end{bmatrix}.
```
It is provided as `TwoParabolasUnconstrained`.
Simply initialize an instance via
```
mop = TwoParabolasUnconstrained();
```

This problem has 2 objectives
```
julia> num_objectives( mop ) 
2
```
We can retrieve them individually as 
```
f1, f2 = get_objectives(mop);
```
They take vector valued arguments of type `Vector{R} where R<:Real`.