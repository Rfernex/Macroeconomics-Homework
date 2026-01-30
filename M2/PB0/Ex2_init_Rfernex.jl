##### PLUTO NOTEBOOK LINK####
print("Pluto notebook link : https://pluto.land/n/c9gymj2g")
##############################

# package imports 
using Pkg
Pkg.add(["DataFrames", "PlotlyJS", "Plots", "GLMakie","Distributions","Random"])
using GLMakie
using Plots
using PlotlyJS
using DataFrames
using Distributions
using Random

Random.seed!(2) # anchors randomly drawn values
filepath = "/Users/rfernex/Documents/Education/SciencesPo/Courses/M2/S1/Macro III/TD/PB0/outputs/Exercise_2/" # For exporting output graphs


# Defining parameters
beta = 0.9
r = 0.2
Ny = 50        # Number of Grid Points
Na = 100        # Number of Grid Points
b = 1


################ INITIALIZATION ##############################################################################

## Defining the Grid for the Endogenous State Variable: Capital
ymin = 5; ymax = 10;    # Bounds for Grid
grid_y = collect(range(ymin, ymax, length = Ny)); # Grid for income

amin = -5; amax = 5;    # Bounds for Grid
grid_a = collect(range(amin, amax, length = Na)); # Grid for assets

# No credit constraint 
v_opt = zeros(Ny,Ny);
index_a_opt= Array{Int64,2}(undef,Ny,Ny); # Hint! 
c_1_opt = zeros(Ny,Ny);
c_2_opt = zeros(Ny,Ny);
a_opt = zeros(Ny,Ny);

# With credit constraint
v_opt_bis = zeros(Ny,Ny);
c_1_opt_bis = zeros(Ny,Ny);
c_2_opt_bis = zeros(Ny,Ny);
a_opt_bis = zeros(Ny,Ny);


################ NO CREDIT CONSTRAINT ##############################################################################

function compute_budget(y1,y2,a)
    u = zeros(2)
    u[1] = y1-a
    u[2] = y2 + (1+r)*a
    return u
end

function ufun_TPM(c1,c2,beta)
    u = log(c1) + beta*log(c2)
    return u
end

# Fill the V array for all possible choice of a 
V = zeros(Ny,Ny,Na);
for ia in 1:Na
    for iy1 in 1:Ny
        for iy2 in 1:Ny
            c1, c2 = compute_budget(grid_y[iy1],grid_y[iy2],grid_a[ia])
            if c1>0 && c2>0 
                V[iy1,iy2,ia] = ufun_TPM(c1,c2,beta)
            else
                V[iy1,iy2,ia] = -10^9
            end
        end
    end
end


# Max and recover a_opt
for iy1 in 1:Ny
    for iy2 in 1:Ny
        v_opt[iy1,iy2], idx_a_opt = findmax(V[iy1,iy2,:])
        a_opt[iy1,iy2] = grid_a[idx_a_opt]
        c_1_opt[iy1,iy2] = grid_y[iy1]-a_opt[iy1,iy2]
        c_2_opt[iy1,iy2] = grid_y[iy2]+(1+r)*a_opt[iy1,iy2]
    end
end


y1 = rand((1:length(grid_y)), 4)

# plots a as a function of y2 with y1 fixed
plots_a = []
for i in 1:length(y1)
    P = Plots.plot(grid_y,a_opt[y1[i],:],title=join(["Optimal a as a function of y2 for y1 = ",round(grid_y[y1[i]]; digits = 2)]," "),legendfontsize=6, titlefont = font(7) ,label="No credit constraint",tickfontsize = 7)
    plots_a = push!(plots_a,P)
end
Plots.plot(plots_a..., layout = 4)
Plots.savefig("/Users/rfernex/Documents/Education/SciencesPo/Courses/M2/S1/Macro III/TD/PB0/outputs/Exercise_2/q4_optimal_a.pdf")

# plots v as a function of y2 with y1 fixed
plots_v = []
for i in 1:length(y1)
    P = Plots.plot(grid_y,v_opt[y1[i],:],title=join(["Optimal v as a function of y2 for y1 = ",round(grid_y[y1[i]]; digits = 2)]," "),legendfontsize=6, titlefont = font(7) ,label="No credit constraint",tickfontsize = 7)
    plots_v = push!(plots_v,P)
end
Plots.plot(plots_v..., layout = 4)
Plots.savefig("/Users/rfernex/Documents/Education/SciencesPo/Courses/M2/S1/Macro III/TD/PB0/outputs/Exercise_2/q4_max_utility.pdf")


waves = PlotlyJS.surface(z=v_opt, x=grid_y, y=grid_y, colorscale=colors.viridis)
fig=PlotlyJS.Plot(waves, Layout(width=600, height=400, 
                   scene= attr(aspectmode="data",xaxis_title="y1",yaxis_title="y2", zaxis_title="utility",
                               camera_eye=attr(x=2.55, y=2.55, z=1.4))))
relayout!(fig, scene_aspectmode="cube")
PlotlyJS.savefig(fig, join([filepath,"/q4_3D_graph.html"],""))

################ WITH CREDIT CONSTRAINT ##############################################################################

a_opt_bis = copy(a_opt);
for iy1 in 1:Ny
    for iy2 in 1:Ny
        if a_opt_bis[iy1,iy2] < b
            a_opt_bis[iy1,iy2] = b
        end
        c_1_opt_bis[iy1,iy2] = grid_y[iy1]-a_opt_bis[iy1,iy2]
        c_2_opt_bis[iy1,iy2] = grid_y[iy2]+(1+r)*a_opt_bis[iy1,iy2]
        v_opt_bis[iy1,iy2] = ufun_TPM(c_1_opt_bis[iy1,iy2], c_2_opt_bis[iy1,iy2], beta)
    end
end

waves = PlotlyJS.surface(z=v_opt_bis, x=grid_y, y=grid_y, colorscale=colors.viridis)
fig=PlotlyJS.Plot(waves, Layout(width=600, height=400, 
                   scene= attr(aspectmode="data",xaxis_title="y1",yaxis_title="y2", zaxis_title="utility",
                               camera_eye=attr(x=2.55, y=2.55, z=1.4))))
relayout!(fig, scene_aspectmode="cube")
PlotlyJS.savefig(fig, join([filepath,"q7_3D_graph_CC.html"],""))


################ PLOTS ##############################################################################
using Plots

y2 = rand((1:length(grid_y)), 4)

# plots a as a function of y1 with y2 fixed
plots_a = []
for i in 1:length(y2)
    P = Plots.plot(grid_y,a_opt[:,y2[i]],title=join(["Optimal a as a function of y1 for y2 = ",round(grid_y[y2[i]]; digits = 2)]," "),legendfontsize=6, titlefont = font(7) ,label="No credit constraint",tickfontsize = 7)
    Plots.plot!(grid_y,a_opt_bis[:,y2[i]],title=join(["Optimal a as a function of y1 for y2 = ",round(grid_y[y2[i]]; digits = 2)]," "), legendfontsize=6, titlefont = font(7) ,label="With credit constraint",tickfontsize = 7)
    plots_a = push!(plots_a,P)
end
Plots.plot(plots_a..., layout = 4)
Plots.savefig(join([filepath,"q7_optimal_a_CC.pdf"],""))

# plots v as a function of y1 with y2 fixed
plots_v = []
for i in 1:length(y2)
    P = Plots.plot(grid_y,v_opt[:,y2[i]],title=join(["Optimal a as a function of y1 for y2 = ",round(grid_y[y2[i]]; digits = 2)]," "),legendfontsize=6, titlefont = font(7) ,label="No credit constraint",tickfontsize = 7)
    Plots.plot!(grid_y,v_opt_bis[:,y2[i]],title=join(["Optimal a as a function of y1 for y2 = ",round(grid_y[y2[i]]; digits = 2)]," "), legendfontsize=6, titlefont = font(7) ,label="With credit constraint",tickfontsize = 7)
    plots_v = push!(plots_v,P)
end
Plots.plot(plots_v..., layout = 4)
Plots.savefig(join([filepath,"q7_max_utility_CC.pdf"],""))
