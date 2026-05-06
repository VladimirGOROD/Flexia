using ForwardDiff
using GLMakie
using Flexia

include("solvers.jl")


function robertson(x, params)
    μ = params[1]
    dx = similar(x)
    dx[1] = x[2]
    dx[2] = μ*(1-x[1]^2)*x[2]-x[1] 
    return dx       
end

func(x) = robertson(x, (5.0,))
jac(x) = ForwardDiff.jacobian(func, x)

x0 = Float64[2,0.1]
t = 0:0.01:50
XX = Matrix{Float64}(undef, 2, length(t))
YY = Matrix{Float64}(undef, 2, length(t))
XX[:, 1] = x0
for i in 2:length(t)
    XX[:,i] = XX[:,i-1] + step(t) * func(XX[:,i-1])
end

cros!(YY, x0, Float64[1 0; 0 1], func, jac, step(t))

# fig = Figure(size = (800, 600))
lines( t, XX[1,:], label= "euler")
lines!(t, YY[1,:], label= "cros")
current_figure()
##
