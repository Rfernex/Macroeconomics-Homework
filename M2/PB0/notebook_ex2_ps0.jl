### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ 060841e3-0e96-4b40-a350-39059dd60ff9
begin
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
end

# ╔═╡ 7773cb15-eb89-4b01-bcfe-edb89ebb1e87
md"""
# Problem Set 0
"""

# ╔═╡ b653ec78-916a-11f0-2ac4-45523d3e870b
md"""
## Exercise 2
"""

# ╔═╡ 7c34757d-78ca-4b69-a021-359d79a4e713
md"""
#### Importing the relevant packages
"""

# ╔═╡ 8457dd67-4a38-40ca-8d11-85a4ebb113b1
md"""
#### Defining parameters
"""

# ╔═╡ 660a0a36-d463-48df-944e-25fa4df86555
begin
	beta = 0.9
	r = 0.2
	Ny = 50        # Number of Grid Points
	Na = 100        # Number of Grid Points
	b = 1
end

# ╔═╡ ec6633c8-a5b5-4cc3-93f5-e7f13c812d96
md"""
### Initialization
"""

# ╔═╡ 970e3372-0eed-443d-b361-625dba2ecf11
begin
	# Defining the Grid for the Endogenous State Variable: Capital
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
end

# ╔═╡ c81737db-fd12-41e5-bfb9-372fc0f04dda
md"""
### No credit constraint
"""

# ╔═╡ 743fa7b6-4748-42fa-924a-ada2f7324ded
begin
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
end

# ╔═╡ ee6204e5-1110-4dd6-918c-6467becb778a
md"""
#### Plot "a" as a function of y2 with y1 fixed
"""

# ╔═╡ 1df44dd1-d83e-497f-907a-cd59f000cf9c
begin
	y1 = rand((1:length(grid_y)), 4) # define random value grid
	local plots_a = []
	for i in 1:length(y1)
	    P = Plots.plot(grid_y,a_opt[y1[i],:],title=join(["Optimal a as a function of y2 for y1 = ",round(grid_y[y1[i]]; digits = 2)]," "),legendfontsize=6, titlefont = font(7) ,label="No credit constraint",tickfontsize = 7)
	    plots_a = push!(plots_a,P)
	end
	Plots.plot(plots_a..., layout = 4)
	Plots.savefig("/Users/rfernex/Documents/Education/SciencesPo/Courses/M2/S1/Macro III/TD/PB0/outputs/Exercise_2/q4_optimal_a.pdf")
	Plots.plot!()
end

# ╔═╡ 8fbc4803-5332-4f01-94e0-20fd7aeb5c7e
md"""
#### Plot "v" as a function of y2 with y1 fixed
"""

# ╔═╡ c0afcea5-ca43-4048-8db0-79266fb98e85
begin
	local plots_v = []
	for i in 1:length(y1)
	    P = Plots.plot(grid_y,v_opt[y1[i],:],title=join(["Optimal v as a function of y2 for y1 = ",round(grid_y[y1[i]]; digits = 2)]," "),legendfontsize=6, titlefont = font(7) ,label="No credit constraint",tickfontsize = 7)
	    plots_v = push!(plots_v,P)
	end
	Plots.plot(plots_v..., layout = 4)
	Plots.savefig("/Users/rfernex/Documents/Education/SciencesPo/Courses/M2/S1/Macro III/TD/PB0/outputs/Exercise_2/q4_max_utility.pdf")
	Plots.plot!()
end

# ╔═╡ cc700ee1-0569-4f71-98ea-170fef8a24c1
begin
	# 3D graph of optimal utility (no constraint) 
	waves = PlotlyJS.surface(z=v_opt, x=grid_y, y=grid_y, colorscale=colors.viridis)
	fig=PlotlyJS.Plot(waves, Layout(width=600, height=400, 
	                   scene= attr(aspectmode="data",xaxis_title="y1",yaxis_title="y2", zaxis_title="utility",
	                               camera_eye=attr(x=2.55, y=2.55, z=1.4))))
	relayout!(fig, scene_aspectmode="cube")
	PlotlyJS.savefig(fig, join([filepath,"/q4_3D_graph.html"],""))
	fig
end

# ╔═╡ 373ed927-1ff9-46e1-ba31-07d1c3c3ab20
md"""
### With credit constraint
"""

# ╔═╡ 3c6ba866-025c-45ec-b9ca-9508facf6e6f
begin
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
end

# ╔═╡ 723c2b4f-c6e7-4480-b9f6-83b32dcd4c92
md"""
### Plots
"""

# ╔═╡ bb247608-f26b-4d0e-8c40-be4e91a44c46
begin
	# 3D graph of optimal utility (with constraint) 
	wavesCC = PlotlyJS.surface(z=v_opt_bis, x=grid_y, y=grid_y, colorscale=colors.viridis)
	fig2=PlotlyJS.Plot(wavesCC, Layout(width=600, height=400, 
	                   scene= attr(aspectmode="data",xaxis_title="y1",yaxis_title="y2", zaxis_title="utility",
	                               camera_eye=attr(x=2.55, y=2.55, z=1.4))))
	relayout!(fig2, scene_aspectmode="cube")
	PlotlyJS.savefig(fig, join([filepath,"q7_3D_graph_CC.html"],""))
	fig2
end

# ╔═╡ e806f6b3-bc00-4815-9db2-b72691132e1c
md"""
#### Plot "a" as a function of y1 with y2 fixed (with Constraint)
"""

# ╔═╡ 6b56bd70-c84e-4c7e-86f1-6b69e3c6d65c
begin
	y2 = rand((1:length(grid_y)), 4) # define random value grid
	local plots_aCC = []
	for i in 1:length(y2)
	    P = Plots.plot(grid_y,a_opt[:,y2[i]],title=join(["Optimal a as a function of y1 for y2 = ",round(grid_y[y2[i]]; digits = 2)]," "),legendfontsize=6, titlefont = font(7) ,label="No credit constraint",tickfontsize = 7)
	    Plots.plot!(grid_y,a_opt_bis[:,y2[i]],title=join(["Optimal a as a function of y1 for y2 = ",round(grid_y[y2[i]]; digits = 2)]," "), legendfontsize=6, titlefont = font(7) ,label="With credit constraint",tickfontsize = 7)
	    plots_aCC = push!(plots_aCC,P)
	end
	Plots.plot(plots_aCC..., layout = 4)
	Plots.savefig(join([filepath,"q7_optimal_a_CC.pdf"],""))
	Plots.plot!()
end

# ╔═╡ e205f0d4-e6cb-420e-ba29-a67f5c2575d9
md"""
#### Plot "v" as a function of y1 with y2 fixed (with Constraint)
"""

# ╔═╡ 654e9565-92fe-4022-80b0-2638a06aef8f
begin
	local plots_vCC = []
	for i in 1:length(y2)
	    P = Plots.plot(grid_y,v_opt[:,y2[i]],title=join(["Optimal a as a function of y1 for y2 = ",round(grid_y[y2[i]]; digits = 2)]," "),legendfontsize=6, titlefont = font(7) ,label="No credit constraint",tickfontsize = 7)
	    Plots.plot!(grid_y,v_opt_bis[:,y2[i]],title=join(["Optimal a as a function of y1 for y2 = ",round(grid_y[y2[i]]; digits = 2)]," "), legendfontsize=6, titlefont = font(7) ,label="With credit constraint",tickfontsize = 7)
	    plots_vCC = push!(plots_vCC,P)
	end
	Plots.plot(plots_vCC..., layout = 4)
	Plots.savefig(join([filepath,"q7_max_utility_CC.pdf"],""))
	Plots.plot!()
end

# ╔═╡ Cell order:
# ╠═7773cb15-eb89-4b01-bcfe-edb89ebb1e87
# ╠═b653ec78-916a-11f0-2ac4-45523d3e870b
# ╠═7c34757d-78ca-4b69-a021-359d79a4e713
# ╠═060841e3-0e96-4b40-a350-39059dd60ff9
# ╠═8457dd67-4a38-40ca-8d11-85a4ebb113b1
# ╠═660a0a36-d463-48df-944e-25fa4df86555
# ╠═ec6633c8-a5b5-4cc3-93f5-e7f13c812d96
# ╠═970e3372-0eed-443d-b361-625dba2ecf11
# ╠═c81737db-fd12-41e5-bfb9-372fc0f04dda
# ╠═743fa7b6-4748-42fa-924a-ada2f7324ded
# ╠═ee6204e5-1110-4dd6-918c-6467becb778a
# ╠═1df44dd1-d83e-497f-907a-cd59f000cf9c
# ╠═8fbc4803-5332-4f01-94e0-20fd7aeb5c7e
# ╠═c0afcea5-ca43-4048-8db0-79266fb98e85
# ╠═cc700ee1-0569-4f71-98ea-170fef8a24c1
# ╠═373ed927-1ff9-46e1-ba31-07d1c3c3ab20
# ╠═3c6ba866-025c-45ec-b9ca-9508facf6e6f
# ╠═723c2b4f-c6e7-4480-b9f6-83b32dcd4c92
# ╠═bb247608-f26b-4d0e-8c40-be4e91a44c46
# ╠═e806f6b3-bc00-4815-9db2-b72691132e1c
# ╠═6b56bd70-c84e-4c7e-86f1-6b69e3c6d65c
# ╠═e205f0d4-e6cb-420e-ba29-a67f5c2575d9
# ╠═654e9565-92fe-4022-80b0-2638a06aef8f
