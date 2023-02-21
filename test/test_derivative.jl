using NonLinearProg
using LinearAlgebra
using ForwardDiff
using TickTock

fun = function(x)
    Ix = I(length(x))
    fval =  x'*Ix*x
    return fval[1]
end


x0 = rand(5)

tick()
NonLinearProg.derivative(fun,x0)
tock()

tick()
ForwardDiff.gradient(fun, x0)
tock()


