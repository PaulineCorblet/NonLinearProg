using NonLinearProg

# Test with Rosenbrock function
function fun(x)
    return (1-x[1])^2+100*(x[2]-x[1]^2)^2
end

function gfun(x)
    return [-2*(1-x[1])-200*x[1]*(x[2]-x[1]^2); 200*(x[2]-x[1]^2)] 
end

x0 = [0.9266578756641997; 0.4065842040265859]
xopt, objopt, termstat = NonLinearProg.fmincon(fun,x0,g=gfun)
xopt, objopt, termstat = NonLinearProg.fmincon(fun,x0,g=gfun,tol=1e-4, max_iter=10)


# Add constraint: on the unit disk
function h(x)
    return [x[1]^2+x[2]^2]
end

function J(x)
    return [2*x[1] 2*x[2]]
end

cons_ub = fill(1., 1)
lb = fill(.9, 2)
xopt, objopt, termstat = NonLinearProg.fmincon(fun, x0; g=gfun, h=h, J=J, nlcon_ub=cons_ub)
xopt, objopt, termstat = NonLinearProg.fmincon(fun,x0; g=gfun, h=h, J=J, nlcon_ub=cons_ub, lb=lb)
