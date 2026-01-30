
#### PLUTO NOTEBOOK LINK ####
print("Pluto notebook link : https://pluto.land/n/yguc6pup")
#############################

using DataFrames
using Plots

##Define Parameters
α =  0.36 # Capital share in production function
β = 0.99 # Discount factor
δ = 0.025 # Depreciation rate
θ = 1.0                     

##Utility function
function u(c) 
    return log(c)     
end         

##Producution function
function F(k) 
    return θ * k^α    
end       

##Defining the steady state
k_ss = (α/((1/β) - (1-δ)))^(1/(1-α));  # Steady State Capital Stock
y_ss = F(k_ss);                # Steady State Output
c_ss = y_ss - δ*k_ss;           # Steady State Consumption

##Defining the Grid for the Endogenous State Variable: Capital
nk = 501;                            # Number of Grid Points
kmin = 0.2*k_ss; kmax = 1.8*k_ss;    # Bounds for Grid
kg = collect(range(kmin, kmax; length=nk)); # Equally Spaced Grid for Capital

## Build the 2-Dimensional Contemporaneous Utility Grid for the System
U = Array{Float64}(undef, nk, nk)    

for ii = 1:nk       # Loop Over Capital Today
    for jj = 1:nk   # Loop Over Capital Tomorrow
        k  = kg[ii]     # Capital Today
        kp = kg[jj]     # Capital Tomorrow
        # Solve for Consumption at Each Point
        c = F(k) + (1-δ)*k - kp              
        if (kp < 0) || (c <= 0)
            # If Tomorrow's Capital Stock | Today's Consumption is Negative
            U[ii,jj] = -1.0e10               
        else
            U[ii,jj] = u(c)                 
        end
    end
end

##Value Function Iteration

#Initial Guess of the Value Function
V0 = ones(nk)               # nk-vector initial guess

tol = 1e-4;
maxits = 3000; # Define the maximum number of iterations
V_1 = copy(V0);  # The new value function after an iteration
V_0 = copy(V0);  # The value function from which I start each iteration
dif = 1.0;
policy_k_index = Array{Int64}(undef, nk);
conv_time = 0 

function belman_op(U,V_0,V_1, policy_k_index)
    for ii in 1:nk
        vals = U[ii, :] .+ β .* V_0          
        vmax, idx = findmax(vals)            
        V_1[ii] = vmax
        policy_k_index[ii] = idx             
    end
    return policy_k_index, V_1
end

for iter in 1:maxits
    if dif < tol
        println("Converged at iteration: $iter")
        break
    else
        policy_k_index, V_1 = belman_op(U,V_0,V_1,policy_k_index)
        dif = maximum(abs.(V_1 .- V_0))          
        V_0 .= V_1 
        conv_time = conv_time+1                              
    end
end

##Policy function for capital
policy_k = kg[policy_k_index]

##Policy function for consumption
policy_c = Array{Float64}(undef,nk)

for k1 in 1:nk
    policy_c[k1] = F(kg[k1]) + (1-δ)*kg[k1] - policy_k[k1]
end

##Plotting the results
plot(kg,V_1,title="Value Function",label="Value Function")
plot(kg,policy_c,title="Policy function for consumption",label="Optimal Consumption")
plot(kg,policy_k,title="Policy function for capital",label="Optimal Capital")

# BONUS QUESTION. 
V_0_b = [ U[i,i] / (1 - β)  for i in 1:nk ] # we assume that k = k' 
V_1_b = copy(V_0_b);
dif_b = 1.0
policy_k_index_b = Array{Int64}(undef, nk);
conv_time_b = 0

for iter_b in 1:maxits
    if dif_b < tol
        println("Converged at iteration: $iter_b")
        break
    else
        policy_k_index_b, V_1_b = belman_op(U,V_0_b,V_1_b,policy_k_index_b)
        dif_b = maximum(abs.(V_1_b .- V_0_b))          
        V_0_b .= V_1_b  
        conv_time_b = conv_time_b+1                           
    end
end

iteration_time =  DataFrame(Guess = ["intial guess","new guess"], Iteration_time = [conv_time,conv_time_b])
