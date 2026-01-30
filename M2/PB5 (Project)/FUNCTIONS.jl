############################ FUNCTIONS USED IN THE MAIN ################################

### FUNCTION : Utility of consumption function
function uval(c,par)
    #return c # case 1 : linear utility function (less realistic)
    return (c^(1-par.mu)-1)/(1-par.mu)  # case 2 : CRRA utility function
    #return 
end 

### FUNCTION : Disutility of labor function
function ldis(l,par)
    #return l # case 1 : linear disutility of labor (less realistic)
    return par.psi*((l)^(1+par.gamma))/(1+par.gamma) # case 2 : isoelastic disutility of labor 
end

### FUNCTION : Budget constraint function
function cons(ia,ia_p,il,job,status, next_job, par)
    if job == 1 
        return (1+par.r)*agrid[ia] + lgrid[il]*prod_grid[status,job]-agrid[ia_p] 
    else 
        if next_job == 1 
            return (1+par.r)*agrid[ia] + par.l_bar*prod_grid[status,job]-agrid[ia_p] - par.F 
        else
            return (1+par.r)*agrid[ia] + par.l_bar*prod_grid[status,job]-agrid[ia_p]
        end 
    end 
end 

### FUNCTION : VFI 

VFI = function(V0,P0,prob,par,U)
    Vguess_se = copy(V0)
    Vguess_e = copy(V0)
    Vnew_se = copy(V0)
    Vnew_e = copy(V0)
    policy_idx_se_a = copy(P0)
    policy_idx_e_a = copy(P0)
    policy_idx_se_job = copy(P0)
    policy_idx_e_job = copy(P0)
    policy_idx_se_l = copy(P0)
    for iter = 1:par.maxits
        for ia = 1:Na, status = 1:2
            tv_1_se, idx_1_se = findmax(U[:,ia,:,1,status,1].+par.beta.*(prob[status,:]'*Vguess_se[:,:]))
            tv_1_e, idx_1_e = findmax(U[:,ia,:,1,status,2].+par.beta.*(prob[status,:]'*Vguess_e[:,:]))
            tv_2_se, idx_2_se = findmax(U[idx_lbar,ia,:,2,status,1].+par.beta.*(prob[status,:]'*Vguess_se[:,:]))
            tv_2_e, idx_2_e = findmax(U[idx_lbar,ia,:,2,status,2].+par.beta.*(prob[status,:]'*Vguess_e[:,:]))
            Vnew_se[status,ia], policy_idx_se_job[status,ia] = findmax([tv_1_se,tv_1_e])
            Vnew_e[status,ia], policy_idx_e_job[status,ia] = findmax([tv_2_se,tv_2_e])
            if policy_idx_se_job[status,ia] == 1 
                policy_idx_se_a[status,ia] = idx_1_se[2]
                policy_idx_se_l[status,ia] = idx_1_se[1]
            else 
                policy_idx_se_a[status,ia] = idx_1_e[2]
                policy_idx_se_l[status,ia] = idx_1_e[1]
            end
            if policy_idx_e_job[status,ia] == 1 
                policy_idx_e_a[status,ia] = idx_2_se[2]
            else 
                policy_idx_e_a[status,ia] = idx_2_e[2]
            end
        end
        if maximum(abs,Vguess_se.-Vnew_se) < par.tol && maximum(abs,Vguess_e.-Vnew_e) < par.tol
            println("Found solution after $iter iterations")
            break
        elseif iter==par.maxits
            println("No solution found after $iter iterations")
            break
        end
        Vguess_se = copy(Vnew_se)
        Vguess_e = copy(Vnew_e)
    end
    return Vnew_se, Vnew_e, policy_idx_se_a, policy_idx_e_a,  policy_idx_se_job, policy_idx_e_job, policy_idx_se_l 
end

### FUNCTION : get policy function values 
get_pol_func = function(policy_idx_a,policy_idx_l,policy_idx_job,job) 
    job_policy_a = zeros(2,Na)
    job_policy_l = zeros(2,Na)
    job_policy_job = Array{String,2}(undef, 2, Na)
    for ia = 1:Na, status = 1:2
        job_policy_a[status,ia] = agrid[Int(policy_idx_a[job,status,ia])]
        job_policy_l[status,ia] = lgrid[Int(policy_idx_l[job,status,ia])]
        job_policy_job[status,ia] = empgrid[Int(policy_idx_job[job,status,ia])]
    end 
    return job_policy_l, job_policy_a, job_policy_job
end

### FUNCTION :  Helper function to save result matrix
function save_matrix_to_csv(filename, matrix_data,folder_name)
    df = DataFrame(matrix_data, :auto)
    rename!(df, [Symbol("a$i") for i in 1:size(df, 2)])
    path = joinpath(folder_name, filename)
    CSV.write(path, df)
    println("Saved: $path")
end




