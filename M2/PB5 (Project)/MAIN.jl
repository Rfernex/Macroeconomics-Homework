
############################ WHAT'S NEW ? + INFO ##########################################
# - employment is now a choice : the agent chooses between two different forms of employment (self-employment vs wage-employment)
# - labor is determined endogenously and chosen by the the agent when self-employed (under wage-employment the labor supply is exogenous)
# - runtime : ~4-5 min (per set-up)
# - MAKE SURE YOU HAVE DOWNLOADED THE FUNCTIONS.jl SCRIPT AS WELL BEFORE RUNNING THIS FILE
###########################################################################################

using Pkg
Pkg.add(["PlotlyJS", "Plots","Distributions","Random","CSV", "DataFrames"])
using CSV, DataFrames, Plots

############################ INITIALIZATION ###############################################

### EXPORT + IMPORT PATH : Change this line to reflect the folder where you put both the MAIN and the FUNCTIONS scripts 
root_path = "/Users/rfernex/Documents/Education/SciencesPo/Courses/M2/S1/Macro III/TD/PB5"
###########################################

### exogenous Parameters
struct par_model
    beta::Float64 # Discount factor
    mu::Float64  # Risk aversion from CRRA parameter
    r::Float64 # Interest rate 
    b::Float64 # Low productivity gain under wage-employment
    yh_e::Float64 # High productivity gain under wage-employment 
    yh_se::Float64 # High productivity gain under self-employment 
    l_bar::Int64 # Hours worked under wage contract
    F::Int64 # Business start up cost 
    maxits::Int64 # Maximum number of iterations
    tol::Float64 # Distance cap to end VFI
    psi::Float64 # Labor disutility weight
    gamma::Float64 # Labor disutility parameter 
end

setups = [
    (name="High Persistence, Low Start-up Cost, Low revenue gap",),
    (name="Low Persistence, High Start-up Cost, Large revenue gap",)
]

println("List of set ups:")
for (i, s) in enumerate(setups)
    println("$i. $(s.name)")
end

print("Enter the number of the set up you want (type it in then press enter - DO BOTH TWICE): \n") # WARNING : need to type and enter twice for choice to be registered (problem with readline function innate to Julia)
choice = parse(Int,readline())
println("Running setup $(choice)...")

if choice == 1 
    par = par_model(0.95,1.1,0.02,2,10,12,25,50,3000,0.01,0.01,1.0);
    amax = 30;
    prob = [0.9 0.1; # Define transition matrix for states of the world  (productivity)
            0.1 0.9];
else
    par = par_model(0.95,1.1,0.02,2,8,12,25,70,3000,0.01,0.01,1.1);
    amax = 30;
    prob = [0.7 0.3; # Define transition matrix for states of the world  (productivity)
            0.3 0.7];
end 

total_start_time = time() # used to track runtime

include(joinpath(root_path, "FUNCTIONS.jl"))

### Define value grids (productivity & assets)

# Assets 
Na = 300;
amin = 0;
agrid = collect(range(amin, length = Na, stop = amax));

# Productivity 
prod_grid = [0.0    par.b;
             par.yh_se  par.yh_e] # first column for self employed | second column for wage-employed

# Labor supply
Nl = 51
lmin = 0;
lmax = 50;
lgrid = round.(collect(range(lmin, length = Nl, stop = lmax)));
idx_lbar = findall(==(par.l_bar), lgrid)

# Employment grid 
empgrid = ["Self-Employed", "Wage-Employed"]
 
### Compute the utility grid 

U = zeros(Nl,Na,Na,2,2,2)

for job = 1:2, status = 1:2, next_job = 1:2, ia = 1:Na, ia_p = 1:Na, il = 1:Nl
    c = cons(ia,ia_p,il,job,status,next_job,par)
    if c < 0 
        U[il,ia,ia_p,job,status,next_job] = -1.0e12
    elseif job == 1
        U[il,ia,ia_p,job,status,next_job] = uval(c,par) - ldis(lgrid[il],par) 
    else
        U[il,ia,ia_p,job,status,next_job] = uval(c,par) - ldis(par.l_bar,par)
    end
end

##########################################################################################


############################ BLOCK 1 : VFI ###############################################

### Get the value function and the three policy function indexes using VFI 
V0 = zeros(2,Na)
P0 = zeros(2,Na)
policy_idx_job = Array{Int64}(undef, 2, 2, Na)
policy_idx_l = Array{Int64}(undef, 2, 2, Na)
policy_idx_a = Array{Int64}(undef, 2, 2, Na)
Vopt = Array{Float64}(undef, 2, 2, Na)
policy_idx_l[2, :, :] .= idx_lbar
Vopt[1,:,:], Vopt[2,:,:], policy_idx_a[1,:,:], policy_idx_a[2,:,:],  policy_idx_job[1,:,:], policy_idx_job[2,:,:], policy_idx_l[1,:,:] = VFI(V0,P0,prob,par,U)

### Retrieve the values of the six policy function (labor, job, assets), three for each job 

se_policy_l, se_policy_a, se_policy_job = get_pol_func(policy_idx_a,policy_idx_l,policy_idx_job,1)
e_policy_l, e_policy_a, e_policy_job = get_pol_func(policy_idx_a,policy_idx_l,policy_idx_job,2)

##########################################################################################

############################ BLOCK 2 : plotting value/policy functions ###################

graph_folder = joinpath(root_path, "Graphs/Results_set_up$(choice)")
mkpath(graph_folder)

## Asset policy function 

# Plot for Employed (E)
p1 = Plots.plot(agrid, Vector(e_policy_a[1, :]),label="Low Productivity",
    color=:red,linewidth=2,
    title="Asset Policy Function: Employed (E)",
    xlabel="Current Assets",
    ylabel="Next Period Assets")
Plots.plot!(p1, agrid, Vector(e_policy_a[2, :]),label="High Productivity",
    color=:blue,linewidth=2)
savefig(p1, joinpath(graph_folder, "asset_policy_employed.png"))

# Plot for Self-Employed (SE)
p2 = Plots.plot(agrid[2:end], Vector(se_policy_a[1,2:end]), label="Low Productivity",
    color=:red,linewidth=2,
    title="Asset Policy Function: Self-Employed (SE)",
    xlabel="Current Assets",
    ylabel="Next Period Assets")
Plots.plot!(p2, agrid[2:end], Vector(se_policy_a[2,2:end]), label="High Productivity",
    color=:blue,linewidth=2)
savefig(p2, joinpath(graph_folder, "asset_policy_self_employed.png"))

## Labor policy function 

# Plot for Self-Employed (SE)
p2 = Plots.plot(agrid, Vector(se_policy_l[1, :]), label="Low Productivity",
    color=:red,linewidth=2,
    title="Labor Policy Function: Self-Employed (SE)",
    xlabel="Current Assets",
    ylabel="Labor Supply")
Plots.plot!(p2, agrid, Vector(se_policy_l[2, :]), label="High Productivity",
    color=:blue,linewidth=2)
savefig(p2, joinpath(graph_folder, "labor_policy_self_employed.png"))

## Job policy function 

# Converting job types to numeric 
function map_job_to_numeric(val)
    s = string(val) 
    if occursin("Self-Employed", s) 
        return 1.0
    else 
        return 0.0
    end
end

e_policy_job_num = map(map_job_to_numeric, Matrix(e_policy_job))
se_policy_job_num = map(map_job_to_numeric, Matrix(se_policy_job))


# Plot for Employed (E)
p1 = plot(agrid, e_policy_job_num[1, :],label="Low Productivity",
    color=:red,linewidth=2,seriestype=:steppost,  
    title="Job Policy: Employed (E)",
    xlabel="Current Assets",
    ylabel="Job Choice (for next period)",
    yticks=([0, 1], ["Wage Emp", "Self-Emp"]),
    ylim=(-0.1, 1.1),
    legend=:right)
plot!(p1, agrid, e_policy_job_num[2, :],label="High Productivity",
    color=:blue,linewidth=2,seriestype=:steppost)
savefig(p1, joinpath(graph_folder, "job_policy_employed.png"))

# Plot for Self-Employed (SE)
p2 = plot(agrid[2:end], se_policy_job_num[1, 2:end],label="Low Productivity",
    color=:red,linewidth=2,seriestype=:steppost,  
    title="Job Policy: Self-Employed (SE)",
    xlabel="Current Assets",
    ylabel="Job Choice (for next period)",
    yticks=([0, 1], ["Wage Emp", "Self-Emp"]),
    ylim=(-0.1, 1.1),
    legend=:right)
plot!(p2, agrid[2:end], se_policy_job_num[2, 2:end],label="High Productivity",
    color=:blue,linewidth=2,seriestype=:steppost)
savefig(p2, joinpath(graph_folder, "job_policy_self_employed.png"))


## Value Functions

# Value function for Employed (E)
p1 = Plots.plot(agrid[2:end], Vector(Vopt[1,1,2:end]),label="Low Productivity",
    color=:red,linewidth=2,
    title="Value Function : Employed (E)",
    xlabel="Current Assets",
    ylabel="Value")
Plots.plot!(p1, agrid[2:end], Vector(Vopt[1,2,2:end]),label="High Productivity",
    color=:blue,linewidth=2)
savefig(p1, joinpath(graph_folder, "value_function_employed.png"))

# Value function for Self-Employed (E)
p2 = Plots.plot(agrid, Vector(Vopt[2,1,:]),label="Low Productivity",
    color=:red,linewidth=2,
    title="Value Function : Self-Employed (SE)",
    xlabel="Current Assets",
    ylabel="Value")
Plots.plot!(p2, agrid, Vector(Vopt[2,2,:]),label="High Productivity",
    color=:blue,linewidth=2)
savefig(p2, joinpath(graph_folder, "value_function_self_employed.png"))

##########################################################################################

############# OPTIONAL :  Export result tables to csv for inspection #######################

table_folder = joinpath(root_path, "Tables/Results_set_up$(choice)")

mkpath(table_folder)

save_matrix_to_csv(joinpath(table_folder, "se_policy_assets.csv"), se_policy_a, table_folder)
save_matrix_to_csv(joinpath(table_folder, "se_policy_labor.csv"),  se_policy_l, table_folder)
save_matrix_to_csv(joinpath(table_folder, "se_policy_job.csv"),    se_policy_job, table_folder)

save_matrix_to_csv(joinpath(table_folder, "e_policy_assets.csv"),  e_policy_a, table_folder)
save_matrix_to_csv(joinpath(table_folder, "e_policy_labor.csv"),   e_policy_l, table_folder)
save_matrix_to_csv(joinpath(table_folder, "e_policy_job.csv"),     e_policy_job, table_folder)

save_matrix_to_csv(joinpath(table_folder, "Value_SE.csv"), Vopt[1,:,:], table_folder)
save_matrix_to_csv(joinpath(table_folder, "Value_E.csv"),  Vopt[2,:,:], table_folder)


println("\nAll files saved in folder: $table_folder")

##########################################################################################

println("Total runtime: ", time() - total_start_time, " seconds")