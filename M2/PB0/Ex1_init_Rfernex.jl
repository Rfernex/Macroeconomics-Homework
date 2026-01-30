##### PLUTO NOTEBOOK LINK####
print("Pluto notebook link : https://pluto.land/n/f1xtwht6")
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

Random.seed!(1) # anchors randomly drawn values
filepath = "/Users/rfernex/Documents/Education/SciencesPo/Courses/M2/S1/Macro III/TD/PB0/outputs/Exercise_1/" # For exporting output graphs

###### Utility functions

# c is consumption (cannot be negative)
# l is labor (between 0 and 1)
alpha = Observable(0.3) # sets alpha as an observable variable for reactive programming (modifying alpha at any point in the programm will change its value in the related utility functions)
beta = 0.4 
sigma = 2.0 
psi = 1.0
phi = 1.0
cbar = 0.1
mu = (sigma-1)/sigma

# 1. Write the other utility functions, with c and l as only inputs.

# Cobb–Douglas utility
function ufun_CD(c::Float64, l::Float64)
    u = @lift(c^$alpha * (1-l)^(beta))
    return u
end

function ufun_CES(c::Float64, l::Float64)
    u = (c^mu + (1-l)^(1-mu))^(1/mu)
    return u
end

function ufun_FE(c::Float64, l::Float64)
    u = c^(1-sigma)/(1-sigma) - (1/psi) * l^(1+1/phi)/(1+phi)
    return u
end

function ufun_SG(c::Float64, l::Float64)
    if c>cbar 
        u = @lift((c-cbar)^$alpha * (1-l)^beta)
    else
        u = 0
    end
    return u
end
# others... 

# 2. Try!

####### Budget constraint
y = 0.1
w = 1.0

function budgeting(l,w,y)
    budget = w*l + y
    return budget
end;

l = [0.2 0.4 1];
c = budgeting.(l,w::Float64,y::Float64);

list = zeros(length(c)*length(l),4)
iter = 1
for i in c
    for j in l
        list[iter,1] = ufun_CD(i,j)[]
        list[iter,2] = ufun_CES(i,j)
        list[iter,3] = ufun_FE(i,j)
        list[iter,4] = ufun_SG(i,j)[]
        iter = iter+1
    end
end

utilities_sample = DataFrame(list, :auto)
rename!(utilities_sample,[:CD,:CES,:FRISCH,:SG]);


# 3. Define c as a function of l through the budget constraint 
# for function version, check 2.

print("see budgeting function in 2. ")

# 4. Try! 

####### Brute-force search for optimal labor

nl = 101
lgrid = collect(range(0,1,nl)) # create 'all' possible choices of labor 
v_CD = zeros(nl) # initialize the vector of utility of all possible choices 
print("See question 5 for detailed approach")

# 5. Compute the utility for each function and for all choices 
utilities = DataFrame(l = Float64[], c = Float64[], CD = Float64[], CES = Float64[], FRISCH = Float64[], SG = Float64[])
for k in 1:nl
    l = lgrid[k]
    c = budgeting(l,w,y)
    CD = ufun_CD(c,l)[]
    CES = ufun_CES(c,l)
    FRISCH = ufun_FE(c,l)
    SG = ufun_SG(c,l)[]
    push!(utilities, (l,c,CD,CES,FRISCH,SG)) 
end
utilities

# 6. Plot it
Plots.plot(lgrid, v_CD)
Plots.plot(utilities.l,[utilities.CD,utilities.CES,utilities.FRISCH,utilities.SG], labels = ["CD" "CES" "FRISCH" "SG"], title = "Utilities as a function of labor supply",titlefont = font(8), tickfontsize = 7, legendfontsize = 7)
Plots.savefig(join([filepath,"/q6_utility_comparison.pdf"],""))


# 7. Find the maximum
utility_max, idx_utility_max = findmax([utilities.CD utilities.CES utilities.FRISCH utilities.SG])
l_opt_cross = utilities.l[idx_utility_max[1]]
Plots.plot(utilities.l,[utilities.CD,utilities.CES,utilities.FRISCH,utilities.SG], labels = ["CD" "CES" "FRISCH" "SG"], title = "Utilities as a function of labor supply",titlefont = font(8), tickfontsize = 7, legendfontsize = 7)
vline!([l_opt_cross], label="Optimal l")
Plots.savefig(join([filepath,"/q7_max_utility.pdf"],""))

######### Function to compute the optimal labor choice for a given utility specification 

# 8. Create a function to compute v, with ufun as input
y = 0.1
w = 1.0 
nl = 101

function compute_v(ufun::Function, w::Float64, y::Float64, lgrid, nl)
    v = zeros(nl)
    for k in 1:nl
        l = lgrid[k]
        c = budgeting(l,w,y)
        v[k] = ufun(c,l)[]
    end
    v_max, i_max = findmax(v)
    l_opt = lgrid[i_max]
    return l_opt, v_max, i_max
end;

function max_util(w, y, lgrid, nl)
    l_CD, v_max_CD, i_max_CD = compute_v(ufun_CD,w,y,lgrid,nl)
    l_CES, v_max_CES, i_max_CES = compute_v(ufun_CES,w,y,lgrid,nl)
    l_FE, v_max_FE, i_max_FE = compute_v(ufun_FE,w,y,lgrid,nl)
    l_SG, v_max_SG, i_max_SG = compute_v(ufun_SG,w,y,lgrid,nl)
    opt_ls = [l_CD, l_CES, l_FE, l_SG]
    u_max, idx_u_max = findmax([v_max_CD, v_max_CES, v_max_FE, v_max_SG])
    cross_opt_l = opt_ls[idx_u_max]
    return cross_opt_l, u_max
end;

######### Comparative statics over income and wages

# 9. create a grid for wealth and income

ny = 6
nw = 11
ygrid = LinRange(0.0,0.5,ny);
wgrid = LinRange(1.0,2.0,nw);
V_opt = zeros(nw,ny);
l_opt = zeros(nw,ny);
i_opt = zeros(Int64,nw,ny);

# 10. Compute the optimal labor supply for each (w,y)

df_optimal_supply = DataFrame(w = Float64[], y = Float64[], opt_l = Float64[],max_utility = Float64[])

for iy in 1:ny
    for iw in 1:nw
        w = wgrid[iw]
        y = ygrid[iy]
        l_opt[iw,iy], V_opt[iw,iy] =  max_util(w,y,lgrid,nl)
        push!(df_optimal_supply, (w,y,V_opt[iw,iy],l_opt[iw,iy]))
    end
end

df_optimal_supply

####### Plots for comparative statics

# 11. plot optimal labor supply vs. wealth for fixed wages

alpha = 0.4 # please modify this parameter and rerun all lines below to try for different values of alpha (the plot files created have the value of alpha in their names)

V_opt = zeros(nw,ny,4);
l_opt = zeros(nw,ny,4);
i_opt = zeros(Int64,nw,ny,4);
df_optimal_supply_all = DataFrame(w = Float64[], y = Float64[], opt_l_CD = Float64[], opt_l_CES = Float64[], opt_l_FE = Float64[], opt_l_SG = Float64[])
for iy in 1:ny
    for iw in 1:nw
        w = wgrid[iw]
        y = ygrid[iy]
        l_opt[iw,iy,1], V_opt[iw,iy,1], i_opt[iw,iy,1] =  compute_v(ufun_CD,w,y,lgrid,nl)
        l_opt[iw,iy,2], V_opt[iw,iy,2], i_opt[iw,iy,2] =  compute_v(ufun_CES,w,y,lgrid,nl)
        l_opt[iw,iy,3], V_opt[iw,iy,3], i_opt[iw,iy,3] =  compute_v(ufun_FE,w,y,lgrid,nl)
        l_opt[iw,iy,4], V_opt[iw,iy,4], i_opt[iw,iy,4] =  compute_v(ufun_SG,w,y,lgrid,nl)
        push!(df_optimal_supply_all, (w,y,l_opt[iw,iy,1],l_opt[iw,iy,2],l_opt[iw,iy,3],l_opt[iw,iy,4]))
    end
end

# fixed income with variable wealth
plots_w = []
w1 = rand((1:length(wgrid)), 4)
for i in 1:length(w1)
    income_idx = df_optimal_supply_all.w .== wgrid[i]; 
    df_plot = df_optimal_supply_all[income_idx,:]
    P = Plots.plot(df_plot.y, [df_plot.opt_l_CD,df_plot.opt_l_CES,df_plot.opt_l_FE,df_plot.opt_l_SG],titlefont = font(5), labels = ["CD" "CES" "FRISCH" "SG"],title = join(["Optimal labor supply as a function of wealth for income w1 = ",round(wgrid[i]; digits = 2)]," "), tickfontsize = 5, legendfontsize = 4)
    plots_w = push!(plots_w,P)
end
Plots.plot(plots_w..., layout = length(w1))
filename_string = join([filepath,"q11_plot_w_",alpha,".pdf"],"") 
Plots.savefig(filename_string)


# 12. Compare the optimal choices for each utility function!

# fixed wealth with variable income
plots_y = []
y1 = rand((1:length(ygrid)), 4) 
for i in 1:length(y1)
    wage_idx = df_optimal_supply_all.y .== ygrid[i]; 
    df_plot = df_optimal_supply_all[wage_idx,:]
    P = Plots.plot(df_plot.w, [df_plot.opt_l_CD,df_plot.opt_l_CES,df_plot.opt_l_FE,df_plot.opt_l_SG],titlefont = font(5), labels = ["CD" "CES" "FRISCH" "SG"],title = join(["Optimal labor supply as a function of income for wealth y1 = ",round(ygrid[i]; digits = 2)]," "), tickfontsize = 5, legendfontsize = 4)
    plots_y = push!(plots_y,P)
end
Plots.plot(plots_y..., layout = length(y1))
filename_string = join([filepath,"q12_plot_y_",alpha,".pdf"],"")
Plots.savefig(filename_string)


# We observe that, for the CES and the Cobb-Douglas utility functions, the optimal labor supply increases markedly as a function of wages. This implies a strong substitution effect (or a strong positive income effect) as the opportunity cost of leisure rises sharply in reponse to a rise in wages for these utility functions. 
# When wealth is strictly positive, we observe a similar substitution effect for the Stone-Great utility. 
# In comparison, the utility function with the Frisch elasticity of labor exhibits a strong negative income effect. In other words, when labor income goes up, people prefer to devote more of their time to leisurely activities, at the expense of working more to get additional income. 

# Looking at the plots for the Cobb-Douglas and the Stone Geary function for different values of alpha, we see that it does not strongly impact the trend observed above, however the base level of the optimal labor supply increases sharply in alpha.

