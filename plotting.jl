# import Pkg
# Pkg.add("Plots")
using Plots

include("project1_jl/helpers.jl")
include("project1_jl/project1.jl")
include("project1_jl/simple.jl")

# function rosenbrock(x::Vector)
#     return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
# end

function rosenbrock(x,y)
    return (1.0 - x)^2 + 100.0 * (y - x^2)^2
end

# function rosenbrock_gradient(x::Vector)
#     storage = zeros(2)
#     storage[1] = -2.0 * (1.0 - x[1]) - 400.0 * (x[2] - x[1]^2) * x[1]
#     storage[2] = 200.0 * (x[2] - x[1]^2)
#     return storage
# end

# Plotting optimizer function
xhistory, fhistory = hookeJeeves(rosenbrock, rosenbrock_gradient, [1, 1], 20, 0.75, 0.001, 0.5)  
println(xhistory)

# Contour Plot of Rosenbrock
xr = -2:0.1:2
yr = -2:0.1:2

contour(xr, yr, rosenbrock, levels = [10,25,50,100,200,250,300], colorbar = false, c = cgrad(:viridis, rev = true), legend = false, xlims = (-2, 2), ylims = (-2, 2), xlabel = "x1", ylabel = "x2", aspectratio = :equal, clim = (2, 500))
plot!([xhistory[i][1] for i = 1:length(xhistory)], [xhistory[i][2] for i = 1:length(xhistory)], color = :black)
savefig("example_contour.png")

# Convergence Plot
plot(collect(1:length(fhistory)), fhistory, xlabel = "Iteration", ylabel = "f(x)")
savefig("example_convergence.png")