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


mutable struct OptProblem{f,g,h,J,J_struct} <: MOI.AbstractNLPEvaluator
    obj::f
    obj_grad::g
    cons::h
    cons_jac::J
    cons_jac_struct::J_struct
end

function MOI.initialize(prob::OptProblem, requested_features::Vector{Symbol})
    for feat in requested_features
        if !(feat in MOI.features_available(prob))
            error("Unsupported feature $feat")
        end
    end
end

function MOI.features_available(prob::OptProblem)
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
MOI.eval_objective(prob::OptProblem, x) = prob.obj(x)

# Objective gradient
function MOI.eval_objective_gradient(prob::OptProblem, grad_f, x) 
    grad_temp = prob.obj_grad(x)
    for j=1:length(x)
        grad_f[j] = grad_temp[j]
    end
end

# Non linear constraints
function MOI.eval_constraint(prob::OptProblem, cons_h, x)
    cons_temp = prob.cons(x)
    cons_h .= cons_temp
end

# Non linear constraints jacobian
function MOI.eval_constraint_jacobian(prob::OptProblem, jac_J, x)
    cons_jac_temp = prob.cons_jac(x)
    jac_J .= vec(cons_jac_temp)
end

# Non linear constraints jacobian structure
MOI.jacobian_structure(prob::OptProblem) = prob.cons_jac_struct


# Wrapper
function fmincon(f, x0; A = nothing, b = nothing, Aeq = nothing, beq = nothing, 
                                  g = nothing, h = nothing, J = nothing, 
                                  lb = nothing, ub = nothing, cons_lb = nothing, cons_ub = nothing,
                                  optimizer = Ipopt.Optimizer())

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

        if ~isnothing(cons_lb) & isnothing(cons_ub)
            cons_ub = fill(Inf, length(cons_lb))
        end
        if  isnothing(cons_lb) & ~isnothing(cons_ub)
            cons_lb = fill(-Inf, length(cons_ub))
        end
        if isnothing(cons_lb) & isnothing(cons_ub)
            K       = length(h(x0)) #or perhaps we should throw an error
            cons_lb = fill(-Inf, K)
            cons_ub = fill(Inf, K)
        end

        K        = length(cons_lb)
        J_struct = Matrix{Tuple{Int64,Int64}}(undef,K,n)
        for i=1:K
            for j=1:n
                J_struct[i,j] = (i,j)
            end
        end
        J_struct = J_struct[:]
        println(string("Number of nonlinear constraints:         ",K,"."))
        prob       = OptProblem(f,g,h,J,J_struct) 
        block_data = MOI.NLPBlockData(MOI.NLPBoundsPair.(cons_lb, cons_ub), prob, true)
        MOI.set(optimizer, MOI.NLPBlock(), block_data)
    else
        # J_struct = nothing
        # cons_lb  = Float64[]
        # cons_ub  = Float64[]
        println(string("Number of nonlinear constraints:         ",0,"."))
    end
  
    # Optimize
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.optimize!(optimizer)
    xval     = MOI.get(optimizer, MOI.VariablePrimal(), x)
    objval   = MOI.get(optimizer, MOI.ObjectiveValue())
    termstat = MOI.get(optimizer, MOI.TerminationStatus())

    return xval, objval, termstat
end

end