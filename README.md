# NonLinearProg


A Julia package for Non Linear Programming. The package contains the following functions:
- `fmincon`: a wrapper for [MathOptInterface.jl](https://github.com/jump-dev/MathOptInterface.jl) to solve non linear optimization problems of the form

$$
\begin{align}
\min_x f(x)&\\
s.t.\\
lb \leq x \leq ub\\
A\times x \leq b\\
A_{eq}\times x = b_{eq}\\
lb_h \leq h(x) \leq ub_h
\end{align}
$$

- `derivative`: a function to numerically compute the gradient of scalar-valued functions or the jacobian of vector-valued functions.

