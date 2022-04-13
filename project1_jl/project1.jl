#=
        AA 222 Project 1, Unconstrained Optimization
        File Name: project1.jl
        File Description: TODO
        Author: Betty Wan, boyuwan@stanford.edu
        Date: Apr 11st, 2022

        References: 
        [1] Textbook, TODO
        [2] AA 222 Project Tutorial, https://www.youtube.com/watch?v=ZPZ9xknXGcQ
=#

#=
    If you want to use packages, please do so up here.
    Note that you may use any packages in the julia standard library
    (i.e. ones that ship with the julia language) as well as Statistics
    (since we use it in the backend already anyway)
=#

# Julia Packages:
using LinearAlgebra
using Statistics

# File Includes:
# include("plotting.jl")


"""
    optimize(f, g, x0, n, prob)

Arguments:
    - `f`: Function to be optimized
    - `g`: Gradient function for `f`
    - `x0`: (Vector) Initial position to start from
    - `n`: (Int) Number of evaluations allowed. Remember `g` costs twice of `f`
    - `prob`: (String) Name of the problem. So you can use a different strategy for each problem. E.g. "simple1", "secret2", etc.

Returns:
    - The location of the minimum
"""
function optimize(f, g, x0, n, prob)
    if prob == "simple1" 
        xhistory, fhistory = hookeJeeves(f, g, x0, n, 0.75, 0.001, 0.5)
    elseif prob == "simple2"
        xhistory, fhistory = hookeJeeves(f, g, x0, n, 1.5, 0.001, 0.5)
    elseif prob == "simple3"
        xhistory, fhistory = hookeJeeves(f, g, x0, n, 0.5, 0.001, 0.5)
    else 
        xhistory, fhistory = hookeJeeves(f, g, x0, n, 1, 0.001, 0.5)
    end

    x_best = xhistory[argmin(fhistory)]
    return x_best
end

basis(i, n) = [k == i ? 1.0 : 0.0 for k in 1 : n]

"""
    hookeJeeves(f, g, x, n, a, ε, γ=0.5)   

Arguments:
    - `f`: Function to be optimized
    - `g`: Gradient function for `f`
    - `x`: (Vector) Initial position to start from
    - `n`: (Int) Number of evaluations allowed.
    - `a`: The starting step size, reads as alpha
    - `ε`: Tolerance value, at which the function stops running if the step size is less than ε
    - `γ`: Decay step. The factor of reduction on step size if no improvement is seen

Returns:
    - The history of x and f explored
"""

function hookeJeeves(f, g, x, n, a, ε, γ=0.5) 

    xhistory = [x]
    fhistory = [f(x)]
    y, dim = fhistory[1], length(x)

    while a > ε

        improved = false 
        x_best, y_best = x, y

        for i in 1 : dim
            for sgn in (-1,1)
                if count(f, g) >= n
                    return xhistory, fhistory
                end
                xnew = x + sgn * a * basis(i, dim) 
                ynew = f(xnew)
                if ynew < y_best
                    push!(xhistory, xnew)
                    push!(fhistory, ynew)
                    x_best, y_best, improved = xnew, ynew, true
                end
            end 
        end
        x, y = x_best, y_best

        if !improved
            a *= γ 
        end
    
    end
    return xhistory, fhistory 

end