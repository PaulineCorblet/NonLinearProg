using NonLinearProg
using LinearAlgebra
using ForwardDiff
using TickTock

fun = function(x)
    Ix = I(length(x))
    fval =  x'*Ix*x
    return fval[1]
end


x0 = rand(50)

tick()
NonLinearProg.derivative(fun,x0)
tock()

tick()
ForwardDiff.jacobian(fun, x0)
tock()

setA = [1 4]
fun = function(x)

    alpha = x[setA]

    fval = [exp(alpha[1])*2 + alpha[2]*3 ; x[3]]
   

    return fval
end



fun = function(x)
    Ix = I(length(x))
    fval =  x'*Ix*x
    return fval[1]
end
x0 = rand(5)



function insert(x,xk,k)
    new_x = copy(x)
    new_x[k] = xk
    return new_x
end


function derivative_par(f,x)

    for k=1:length(x)
    r = NonLinearProg.derivative(xk -> f(insert(x,xk,k)),x[k])
    println(r)
    end

end

NonLinearProg.derivative(fun2,1.0)

derivative_par(fun,x0)

    
addprocs(2)
pmap(x -> [2*x; x], 1:4)


function fun(x)
return x
end
