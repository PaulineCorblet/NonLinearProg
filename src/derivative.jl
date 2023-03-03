### A function to compute the gradient of scalar-valued functions or the jacobian of vector-valued functions.

function derivative(fun::Function, x0; epsilon = nothing, method = "forward", print_level = 0)

    # Dimensions
    f0 = fun(x0)
    nI = length(f0)
    nJ = length(x0)
    J = zeros(nI,nJ)

    # Step size (based on Calculus.jl and Numerical recipes)
    scale = maximum([ones(nJ) abs.(x0)], dims=2)[:]

    # Compute numerical derivatives
    if method == "forward"
        if isnothing(epsilon)
            epsilon = sqrt.(eps.(eltype.(x0))) .* scale
        else
            epsilon = fill(epsilon, nJ)
        end

        if nJ == 1
            x1      = x0 .+ epsilon[1]
            f1      = fun(x1)
            J[:,1] .= (f1 - f0)/epsilon[1]
        else
            for j=1:nJ
                if print_level == 1
                    print(string("Computing derivatives using ", method," method, ",j,"/",nJ,". epsilon = ",epsilon[j],".","\r"))
                end

            x1      = copy(x0)
            x1[j]   = x1[j] + epsilon[j]
            f1      = fun(x1)
            J[:,j] .= (f1 - f0)/epsilon[j]
            end
        end

    elseif method == "central"
        if isnothing(epsilon)
            epsilon = cbrt.(eps.(eltype.(x0))) .* scale
        else
            epsilon = fill(epsilon, nJ)
        end

        if nJ == 1
            x1      = x1 .- epsilon[1]
            f1      = fun(x1)
            x2      = x2 .+ epsilon[1]
            f2      = fun(x2)
            J[:,1] .= (f2 - f1)/(2*epsilon[1])
        else
            for j=1:nJ
                if print_level == 1
                    print(string("Computing derivatives using ", method," method, ",j,"/",nJ,". epsilon = ",epsilon[j],".","\r"))
                end

            x1      = copy(x0)
            x1[j]   = x1[j] - epsilon[j]
            f1      = fun(x1)
            x2      = copy(x0)
            x2[j]   = x2[j] + epsilon[j]
            f2      = fun(x2)
            J[:,j] .= (f2 - f1)/(2*epsilon[j])
            end 
        end
    end

    if nI == 1
        J = J[:]
    end
    if nJ>1 & print_level == 1
        println()
        println("Done.")
    end
    return J
end



function derivative_par(f::Function, x; epsilon = nothing, method = "forward", print_level = 0)

    # Insert function
    function insert(x,xk,k)
        new_x = copy(x)
        new_x[k] = xk
        return new_x
    end


    # Parallel computing
    J = Distributed.pmap(k -> NonLinearProg.derivative(xk -> f(insert(x,xk,k)), x[k], epsilon = epsilon, method = method, print_level = print_level), collect(1:length(x)))

    # Return
    return reduce(vcat, J)
end