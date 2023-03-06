# This is a MathOptInterface wrapper to solve non-linear optimization problem of the form
# min f(x) 
# s.t 
# lb <= x <= ub
# Aeq*x = beq
# A*x  <= b
# nlcon_lb <= h(x) <= nlcon_ub

mutable struct OptProblemFmincon{f,g,h,J,J_struct} <: MOI.AbstractNLPEvaluator
    obj::f
    obj_grad::g
    cons::h
    cons_jac::J
    cons_jac_struct::J_struct
end

function MOI.initialize(prob::OptProblemFmincon, requested_features::Vector{Symbol})
    for feat in requested_features
        if !(feat in MOI.features_available(prob))
            error("Unsupported feature $feat")
        end
    end
end

function MOI.features_available(prob::OptProblemFmincon)
    if isnothing(prob.obj_grad) & isnothing(prob.cons_jac)
        return
    elseif !isnothing(prob.obj_grad) & isnothing(prob.cons_jac)
        return [:Grad]
    elseif isnothing(prob.obj_grad) & !isnothing(prob.cons_jac)
        return [:Jac]
    else
        return [:Grad, :Jac]
    end
end

# Objective
MOI.eval_objective(prob::OptProblemFmincon, x) = prob.obj(x)

# Objective gradient
function MOI.eval_objective_gradient(prob::OptProblemFmincon, grad_f, x) 
    grad_temp = prob.obj_grad(x)
    for j=1:length(x)
        grad_f[j] = grad_temp[j]
    end
end

# Non linear constraints
function MOI.eval_constraint(prob::OptProblemFmincon, cons_h, x)
    cons_temp = prob.cons(x)
    cons_h .= cons_temp
end

# Non linear constraints jacobian
function MOI.eval_constraint_jacobian(prob::OptProblemFmincon, jac_J, x)
    cons_jac_temp = prob.cons_jac(x)
    jac_J .= vec(cons_jac_temp)
end

# Non linear constraints jacobian structure
MOI.jacobian_structure(prob::OptProblemFmincon) = prob.cons_jac_struct


# Wrapper
function fmincon(f, x0, optimizer; A = nothing, b = nothing, Aeq = nothing, beq = nothing, 
                                  g = nothing, h = nothing, J = nothing, 
                                  lb = nothing, ub = nothing, nlcon_lb = nothing, nlcon_ub = nothing, tol = nothing)

    # Initialize
    MOI.empty!(optimizer)
    println("Problem summary:")

    # Add variables
    n = length(x0)
    x = MOI.add_variables(optimizer, n)
    println(string("Number of variables:                     ", n,"."))

    # Add starting value
    for i=1:n
        MOI.set(optimizer, MOI.VariablePrimalStart(), x[i], x0[i])
    end

    # Add linear inequality constraints
    if ~isnothing(A) & ~isnothing(b)
        for i=1:size(A,1)
            MOI.add_constraint(optimizer, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(A[i,:], x), 0.0), MOI.LessThan(b[i]))
        end
        println(string("Number of linear inequality constraints: ",size(A,1),"."))
    else
        println(string("Number of linear inequality constraints: ",0,"."))
    end

    # Add linear inequality constraints
    if ~isnothing(Aeq) & ~isnothing(beq)
        for i=1:size(Aeq,1)
            MOI.add_constraint(optimizer, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(Aeq[i,:], x), 0.0), MOI.EqualTo(beq[i]))
        end
        println(string("Number of linear equality constraints:   ",size(Aeq,1),"."))
    else
        println(string("Number of linear equality constraints:   ",0,"."))
    end

    # Add variable bounds
    if !isnothing(lb)
        for i in 1:n
            if lb[i] > -Inf
                MOI.add_constraint(optimizer, x[i], MOI.GreaterThan(lb[i]))
            end
        end
    end
    if !isnothing(ub)
        for i in 1:n
            if ub[i] < Inf
                MOI.add_constraint(optimizer, x[i], MOI.LessThan(ub[i]))
            end
        end
    end

    # Add nonlinear constraints
    if ~isnothing(h)

        if ~isnothing(nlcon_lb) & isnothing(nlcon_ub)
            nlcon_ub = fill(Inf, length(nlcon_lb))
        end
        if  isnothing(nlcon_lb) & ~isnothing(nlcon_ub)
            nlcon_lb = fill(-Inf, length(nlcon_ub))
        end
        if isnothing(nlcon_lb) & isnothing(nlcon_ub)
            K       = length(h(x0)) #or perhaps we should throw an error
            nlcon_lb = fill(-Inf, K)
            nlcon_ub = fill(Inf, K)
        end

        K        = length(nlcon_lb)
        J_struct = Matrix{Tuple{Int64,Int64}}(undef,K,n)
        for i=1:K
            for j=1:n
                J_struct[i,j] = (i,j)
            end
        end
        J_struct = J_struct[:]
        println(string("Number of nonlinear constraints:         ",K,"."))
    else
        J_struct = nothing
        nlcon_lb  = Float64[]
        nlcon_ub  = Float64[]

        println(string("Number of nonlinear constraints:         ",0,"."))
    end

    # Pass non-linear constraints & objective function
    prob       = OptProblemFmincon(f,g,h,J,J_struct) 
    block_data = MOI.NLPBlockData(MOI.NLPBoundsPair.(nlcon_lb, nlcon_ub), prob, true)
    MOI.set(optimizer, MOI.NLPBlock(), block_data)
  
    # Optimize
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    if ~isnothing(tol)
        MOI.set(optimizer, MOI.RawOptimizerAttribute("tol"), tol)
    end

    MOI.optimize!(optimizer)
    xval     = MOI.get(optimizer, MOI.VariablePrimal(), x)
    objval   = MOI.get(optimizer, MOI.ObjectiveValue())
    termstat = MOI.get(optimizer, MOI.TerminationStatus())

    return xval, objval, termstat
end

