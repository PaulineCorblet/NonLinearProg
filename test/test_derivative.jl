using NonLinearProg


### Tests
# R -> R
fun1 = function(x)
    return x
end
NonLinearProg.derivative(fun1, 1.0)
NonLinearProg.derivative(fun1, [1.0])


# R -> Rn
fun2 = function(x)
    return [x; 2*x]
end
NonLinearProg.derivative(fun2, 1.0)
NonLinearProg.derivative(fun2, [1.0])

# Rm -> R
fun3 = function(x)
    return sum(x)
end
NonLinearProg.derivative(fun3, [1.0; 2.0])

# Rm -> R
fun4 = function(x)
    return [sum(x) ; sum(2*x)]
end
NonLinearProg.derivative(fun4, [1.0; 2.0])


### Parallel version
addprocs(2)

@everywhere function fun(x)
    return sum(x.*x)
end

NonLinearProg.derivative_par(fun,rand(10))

function fun(x)
    return sum(x.*x)
end
