##############################################################################################################################
## BLOCK 0 ##############################################################################################################################
##############################################################################################################################
using Distributions, LinearAlgebra, Plots

##Define Parameters
struct par_model
    beta::Float64              # discount factor
    mu::Float64            # risk aversion from CRRA parameter
    alpha::Float64            # capital share
    delta::Float64            # depreciation
    A::Float64            # aggregate productivity
    maxits::Int64
    tol::Float64
end

par = par_model(0.95,2.0,0.3,0.03,1.0,3000,1e-6)

##Defining the Grid for the Endogenous State Variable: Capital
Na = 300;
b = -0.2;
amax =  60;
agrid = collect(range(b, length = Na, stop = amax));

##Defining the Grid for the Exogenous State Variable: Technology Shock
rho = 0.9             # persistence of the AR(1) process
sigma = 0.2              # standard deviation of the AR(1) process
Ns = 5



prob = tauchen(0, sigma, rho, Ns)[1];
logs = tauchen(0, sigma, rho, Ns)[2];
sgrid = exp.(logs);


##############################################################################################################################
## BLOCK 1 ##############################################################################################################################
##############################################################################################################################

VFI = function(r, agrid, sgrid, V0, prob, par) 
    # r the interest rate 
    # agrid the asset grid, sgrid the state grid
    # V0 the empty value function matrix
    # par the parameter list
    # prob the transition matrix 

    Ns = length(sgrid)
    Na = length(agrid)
    w = (1-par.alpha)*(par.A*(par.alpha/(r+par.delta))^par.alpha)^(1/(1-par.alpha)); 

    U = zeros(Ns,Na,Na)

    for is in 1:Ns                     # Loop Over skills Today
        for ia in 1:Na                 # Loop Over assets Today
            for ia_p in 1:Na           # Loop Over assets Tomorrow
                a = agrid[ia];     # Technology Today
                a_p = agrid[ia_p];     # Capital Today
                s = sgrid[is];    # Capital Tomorrow
                # Solve for Consumption at Each Point
                c = (1+r)*a + s*w - a_p
                if c .< 0
                    U[is,ia,ia_p] = -10^6;
                else()
                    U[is,ia,ia_p] = c^(1-par.mu)/(1-par.mu);
                end
            end
        end
    end

    Vnew = copy(V0);  # The new value function I obtain after an iteration
    Vguess = copy(V0);  # the  value function from which I start in each new iteration
    policy_a_index = Array{Int64,2}(undef,Ns,Na);
    tv = zeros(Na)

    ### TO FILL ####### 
    for iter in 1:par.maxits
        for is in 1:Ns
            for ia in 1:Na
                tv = U[is,ia,:]'+(par.beta)*prob[is,:]'*Vguess[:,:]
                (Vnew[is,ia], policy_a_index[is,ia]) = findmax(tv[:])
            end
        end
        if maximum(abs,Vguess.-Vnew) < par.tol
            println("Found solution after $iter iterations")
            break
        elseif iter==par.maxits
            println("No solution found after $iter iterations")
            break
        end
        Vguess = copy(Vnew)
    end
    ##################
    return policy_a_index, Vnew
end 




##############################################################################################################################
## BLOCK 2 ##############################################################################################################################
##############################################################################################################################

aiyagari = function(r, par, agrid, sgrid, prob, Vguess)
    Ns = length(sgrid)
    Na = length(agrid)

    # Call the VFI function and get the policy index 
    policy_a_index, Vguess = VFI(r, agrid, sgrid, Vguess, prob, par)

    ##### 1. Building the transition matrix  ####################
    # Build Q as a 4D array
    Q = zeros(Ns,Na,Ns,Na)

    ### TO FILL ####### 
    for is in 1:Ns
        for ia in 1:Na
            ia_p = policy_a_index[is,ia]
            for is_p in 1:Ns
                Q[is,ia,is_p,ia_p] = prob[is,is_p]
            end
        end
    end
    ##################

    # Then reshape it if Q was 4D
    global Q = reshape(Q, Ns*Na, Ns*Na)

    # Check that the rows sum to 1! 
    sum(Q, dims=2)

    ###### 2. Computing the stable distribution  ###################################
    dist = ones(1, Ns*Na) / (Ns*Na);
    dist = get_stable_dist(dist, Q,par)

    # Check that the distribution vector sums to 1! 
    sum(dist)
    
    # Reshape dist as a Ns x Na dimension (more readable) (don't have to)
    global dist = reshape(dist, Ns, Na)

    ###### 3. Computing the aggregate #############################################
    agg_a = 0
    ### TO FILL ####### 
    for is in 1:Ns
        for ia in 1:Na
            next_a_index = policy_a_index[is,ia]
            agg_a += dist[is,ia]*agrid[next_a_index]
        end
    end
    ##################
    global a_idx = policy_a_index       # CHECK
    ##################
    return agg_a, Vguess
end 

function get_stable_dist(invdist, P,par)
    for iter in 1:par.maxits
        invdist2 = invdist * P
        if maximum(abs, invdist2 .- invdist) < 1e-9
            println("Found solution after $iter iterations")
            return invdist2
        elseif iter == par.maxits
            error("No solution found after $iter iterations")
            return invdist
        end
        err = maximum(abs, invdist2 - invdist)
        invdist = invdist2
    end
end

##############################################################################################################################
## BLOCK 3 ##############################################################################################################################
##############################################################################################################################

## Intermediate, do not touch
# we just loop over 3 possible interest rates
r_vec = [0, 0.01, 0.02, 0.03, 0.04, 0.05]
Vguess = zeros(Ns,Na)
agg_a = zeros(length(r_vec))
agg_D = zeros(length(r_vec))

for (ir,r) in enumerate(r_vec)
    agg_a[ir], Vguess = aiyagari(r, par, agrid, sgrid,prob, Vguess)
    agg_D[ir] = (par.alpha/(r+par.delta))^(1/(1-par.alpha))
end 

Plots.plot(agg_a,r_vec[:])
Plots.plot!(agg_D,r_vec[:])

# Bisection to retrieve the values of r, A and K at equilibrium
r_up = 1/par.beta-1
r_low = -par.delta
while r_up-r_low > par.tol
    global r = 0.5*r_up + 0.5*r_low
    global S, Vguess = aiyagari(r, par, agrid, sgrid,prob, Vguess)
    global D = (par.alpha)/(r+par.delta))^(1/(1-par.alpha))
    if S-D > 0
        r_up = r
    else 
        r_low = r
    end
end

plot(agrid, Vguess[:,:]')
plot(agrid, dist[:,:]')
plot(agrid, a_idx[:,:]')
println("The equilibrium interest rate is : $r")
println("The equilibrium aggregate asset supply : $S")
println("The equilibrium aggregate asset demand : $D")