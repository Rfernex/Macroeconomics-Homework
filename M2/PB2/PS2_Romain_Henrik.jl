#### Problem Set 2 - Markov chain process
using Distributions, LinearAlgebra


### Part A - Tauchen method ###############################################################
# Input: tauchen function (do not change)
function tauchen(mean, sd, rho, num_states; q=3)

    uncond_sd = sd/sqrt(1-rho^2)
    y = range(-q*uncond_sd, stop = q*uncond_sd, length = num_states)
    d = y[2]-y[1]

    Pi = zeros(num_states,num_states)

    for row = 1:num_states
      # end points
          Pi[row,1] = cdf(Normal(),(y[1] - rho*y[row] + d/2)/sd)
          Pi[row,num_states] = 1 - cdf(Normal(), (y[num_states] - rho*y[row] - d/2)/sd)

      # middle columns
          for col = 2:num_states-1
              Pi[row, col] = (cdf(Normal(),(y[col] - rho*y[row] + d/2) / sd) -
                             cdf(Normal(),(y[col] - rho*y[row] - d/2) / sd))
          end
    end

  yy = y .+ mean # center process around its mean

  Pi = Pi./sum(Pi, dims = 2) # renormalize

  return Pi, yy
end 

# Parameters value
rho = 0.8
sigma = 0.1225
N = 5
N_iter = 2000;

## Applying the tauchen method (function below) to get the Markov chain process
TransMat,log_y_grid = tauchen( 0, sigma,  rho,  N)
y_grid = exp.(log_y_grid)

## Create a function to get the stable distribution
function get_invdist(pi_old, delta, maxiter, TransMat)
   for iter in 1:N_iter
       pi_next = pi_old*TransMat
       error = maximum(pi_next .- pi_old)
       if error < delta 
           return pi_next
           break
       elseif iter == maxiter
           error("No solution found after $iter iterations")
       end
       pi_old = pi_next
   end
end


## Find the stable distribution for the process generated above 
p = 1/N
init = transpose(repeat([p],5))
delta = 0.01
maxiter = 3000
stable_dist = get_invdist(init, delta, maxiter, TransMat)

## What is the mean income? 
mean_income = sum(transpose(stable_dist).*y_grid)

## What is the share of households with income y_5?
share_top = stable_dist[N]


### Part B - Recessions ###############################################################
P = [0.971 0.029 0.000 
    0.145 0.778 0.077
    0.000 0.508 0.492] 


## Assume we are at normal growth, what is the proba of going in recession in t+1? 
init_growth = transpose([1,0,0])
T1_growth = init_growth*P

## In t+6? 
T6_growth = init_growth*P^6

## Compute the stable distribution
stable_growth = get_invdist(init_growth, delta, maxiter, P)









