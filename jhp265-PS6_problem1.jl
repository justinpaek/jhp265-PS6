# Code for CHEME 7770 PS6 Problem 1
# Justin Paek

#Problem 1b

# Plot the nullclines for the Collins toggle switch
function u_nullcline(alpha, u, n)
    return (alpha ./ u .- 1).^(1/n)
end

function v_nullcline(alpha, u, n)
    return alpha/(1 .+ u.^n)
end

alpha = 10
u = LinRange(0,10,50)
# v = LinRange(0,1,50)

using Plots
# For n = 1
n = 1
plt1 = plot(u, [u_nullcline.(alpha, u, n), v_nullcline.(alpha, u, n)], 
title = "Nullclines: n=1", xaxis="u",  yaxis="v", ylims=(0,10))
display(plt1)
savefig(plt1, "Nullclines_1")
# For n = 2
n = 2
plt2 = plot(u, [u_nullcline.(alpha, u, n), v_nullcline.(alpha, u, n)], 
title = "Nullclines: n=2", xaxis="u",  yaxis="v", ylims=(0,10))
display(plt2)
savefig(plt2,"Nullclines_2")

# Problem 1c
# Generate a phase portrait for the system
using Makie
using AbstractPlotting
using AbstractPlotting.MakieLayout
AbstractPlotting.inline!(true)

function toggle_switch(u, v)
    du_dt = alpha/(1 .+ v.^n) .- u
    dv_dt = alpha/(1 .+ u.^n) .- v 
    return Point(du_dt, dv_dt)
end

# Construct the streamplot
plt3 = Scene(resolution =(400,400))
streamplot!(plt3, toggle_switch, 0..4, 0..4, colormap = :plasma, 
    gridsize= (32,32), arrow_size = 0.05)

display(plt3)
save("phaseplot.png", plt3)

# Problem 1e
# Define functions for the elements of the Jacobian
alpha = 10
n=1
vs= 2.6
us = 2.6

f_u = -1.0
f_v = (-alpha*n*vs^(n-1))/((1+vs^n)^2)
g_u = (-alpha*n*us^(n-1))/((1+us^n)^2)
g_v = -1.0

J = [f_u f_v; g_u g_v]
using LinearAlgebra
eigs1 = eigvals(J)
det1 = det(J)
trace1 = tr(J)

n=2
vs= 2
us = 2

f_u = -1.0
f_v = (-alpha*n*vs^(n-1))/((1+vs^n)^2)
g_u = (-alpha*n*us^(n-1))/((1+us^n)^2)
g_v = -1.0
J = [f_u f_v; g_u g_v]

eigs2 = eigvals(J)
det2 = det(J)
trace2 = tr(J)