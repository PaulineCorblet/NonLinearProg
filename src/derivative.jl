### A function to compute the gradient of scalar-valued functions or the jacobian of vector-valued functions.

function derivative(fun::Function, x0; epsilon = 1e-6, method = "forward")

    # Dimensions
    f0 = fun(x0)
    nI = length(f0)
    nJ = length(x0)
    J = zeros(nI,nJ)

    # Compute numerical derivatives
    if method == "forward"

        if nJ == 1
            x1      = x0 .+ epsilon
            f1      = fun(x1)
            J[:,1] .= (f1 - f0)/epsilon
        else
            for j=1:nJ
            print(string("Computing derivatives using ", method," method, ",j,"/",nJ,".","\r"))
            x1      = copy(x0)
            x1[j]   = x1[j] + epsilon
            f1      = fun(x1)
            J[:,j] .= (f1 - f0)/epsilon
            end
        end

    elseif method == "central"

        if nJ == 1
            x1      = x1 .- epsilon
            f1      = fun(x1)
            x2      = x2 .+ epsilon
            f2      = fun(x2)
            J[:,1] .= (f2 - f1)/(2*epsilon)
        else
            for j=1:nJ
            print(string("Computing derivatives using ", method," method, ",j,"/",nJ,".","\r"))
            x1      = copy(x0)
            x1[j]   = x1[j] - epsilon
            f1      = fun(x1)
            x2      = copy(x0)
            x2[j]   = x2[j] + epsilon
            f2      = fun(x2)
            J[:,j] .= (f2 - f1)/(2*epsilon)
            end 
        end
    end

    if nI == 1
        J = J[:]
    end
    if nJ>1
        println()
        println("Done.")
    end
    return J
end



function derivative_par(f,x)

    # Insert function
    function insert(x,xk,k)
        new_x = copy(x)
        new_x[k] = xk
        return new_x
    end


    # Parallel computing
    J = pmap(k -> NonLinearProg.derivative(xk -> f(insert(x,xk,k)), x[k]), collect(1:length(x)))

    # Return
    return reduce(vcat, J)
end