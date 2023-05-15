# NonLinearProg


NonLinearProg is a Julia package for Non Linear Programming. It is primarily a light wrapper (akin to the [fmincon](https://uk.mathworks.com/help/optim/ug/fmincon.html) function in Matlab) to set up and solve nonlinear optimization problems featuring a non linear objective function, linear and nonlinear equality and inequality constraints, and bound constraints. The package features other useful functions for nonlinear programming. 

The package contains the following functions:
- `fmincon`: a wrapper for [MathOptInterface.jl](https://github.com/jump-dev/MathOptInterface.jl) to solve non linear optimization problems of the form

$$
\begin{align*}
\min_x f(x)&\\
s.t.\\
lb \leq x \leq ub\\
A\times x \leq b\\
A_{eq}\times x = b_{eq}\\
lb_h \leq h(x) \leq ub_h
\end{align*}
$$

where $f$ and $h$ are nonlinear functions.

- `derivative`: a function to numerically compute the gradient of scalar-valued functions or the jacobian of vector-valued functions.

- `derivative_par`: a parallelized version of the `derivative` function.