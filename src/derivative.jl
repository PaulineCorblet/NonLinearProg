### A function to compute the gradient of scalar-valued functions or the jacobian of vector-valued functions.

function derivative(fun::Function, x0; epsilon = 1e-6, method = "forward")

    # Dimensions
    f0 = fun(x0)
    nI = length(f0)
    nJ = length(x0)
    J = zeros(nI,nJ)

    # Compute numerical derivatives
    if method == "forward"
        for j=1:nJ
            sleep(1)
            print(string("Computing derivatives using ", method," method, ",j,"/",nJ,".","\r"))
            x1     = copy(x0)
            x1[j]  = x1[j] + epsilon
            f1     = fun(x1)
            J[:,j] .= (f1 - f0)/epsilon
        end
    elseif method == "central"
        for j=1:nJ
            x1     = copy(x0)
            x1[j]  = x1[j] - epsilon
            f1     = fun(x1)
            x2     = copy(x0)
            x2[j]  = x2[j] + epsilon
            f2     = fun(x2)
            J[:,j] .= (f2 - f1)/(2*epsilon)
        end 
    end

    if nI == 1
        J = J[:]
    end
    println()
    println("Done.")

    return J
end
