# This is a MathOptInterface wrapper to solve non-linear optimization problem of the form
# min f(x) 
# s.t 
# lb <= x <= ub
# Aeq*x = beq
# A*x  <= b
# cons_lb <= h(x) <= cons_ub

module NonLinearProg

import MathOptInterface
import Ipopt
const MOI = MathOptInterface

include("fmincon.jl")


end