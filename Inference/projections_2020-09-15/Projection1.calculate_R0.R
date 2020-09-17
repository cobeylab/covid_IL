library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)
library(matlib)

jacobian=function(states, elist, pts){
    
    k=0
    jl=list(NULL)
    for(i in 1:length(states)){
        assign(paste("jj", i, sep = "."), lapply(lapply(elist, deriv, states[i]), eval, pts))
        for(j in 1:length(states)){
            k=k+1
            jl[[k]]=attr(eval(as.name(paste("jj", i, sep=".")))[[j]], "gradient")[1,]
        }
    }
    
    J=matrix(as.numeric(as.matrix(jl)[,1]), ncol=length(states))
    return(J)
}

get_beta_from_covar = function(input_params, beta_covariate_file, beta_covariate_column){
    covar_table = read.csv(covid_get_path(beta_covariate_file))
    row.names(covar_table) = covar_table$time
    times = covar_table$time
    foreach(region=1:5, .combine='cbind') %do%{
        input_params[[paste0('beta1_', region)]] * covar_table[, paste0(beta_covariate_column, '_',region)]
    } -> temp
    temp = data.frame(temp)
    row.names(temp) = times
    names(temp) = 1:5
    temp
}
 
get_beta_from_parameters = function(input_params){
    t_sip = 62
    foreach(region=1:11, .combine='cbind') %do%{
        t_phase3_max = input_params[[paste0('t_phase3_max_', region)]]
        t_phase4_max = input_params[[paste0('t_phase4_max_', region)]]
        t_phase5_max = input_params[[paste0('t_phase5_max_', region)]]

        scale_phase3 = input_params[[paste0('scale_phase3_', region)]]
        scale_phase4 = input_params[[paste0('scale_phase4_', region)]]
        scale_phase5 = input_params[[paste0('scale_phase5_', region)]]

        t_phase3 = input_params[[paste0('t_phase3_', region)]]
        t_phase4 = input_params[[paste0('t_phase4_', region)]] 
        t_phase5 = input_params[[paste0('t_phase5_', region)]] 

        beta1 = input_params[[paste0('beta1_', region)]]
        beta2 = input_params[[paste0('beta2_', region)]]

        foreach (t = 1:500, .combine='rbind') %do%{
            if (t < t_sip){
                beta = beta1
            } else if (t < t_phase3){
                t_coord = t - t_sip + 1
                slope = (beta1 - beta2) / 7
                beta = if_else(t_coord < 7, beta1 - slope * t_coord, beta2)
            } else if (t < t_phase4){
                phase3_slope = (scale_phase3 * beta2) / (t_phase3_max - t_phase3)
                transmission_phase3 = beta2 + phase3_slope * (t - t_phase3)
                transmission_phase3_max = beta2 + phase3_slope * (t_phase3_max - t_phase3)
                beta = if_else(t >= t_phase3_max, transmission_phase3_max, transmission_phase3)
            } else if (t<t_phase5){
                phase4_change = (scale_phase4 - scale_phase3) * beta2
                phase4_starting_point = beta2 * (1 + scale_phase3)
                phase4_slope = phase4_change / (t_phase4_max - t_phase4)
                transmission_phase4 = phase4_starting_point + phase4_slope * (t - t_phase4)
                transmission_phase4_max = phase4_starting_point + phase4_slope * (t_phase4_max - t_phase4)
                beta = if_else(t >= t_phase4_max, transmission_phase4_max, transmission_phase4)
            } else{
                phase5_change = (scale_phase5 - scale_phase4)  * beta2
                phase5_starting_point = beta2 * (1 + scale_phase4)
                phase5_slope = phase5_change / (t_phase5_max - t_phase5)
                transmission_phase5 = phase5_starting_point + phase5_slope * (t - t_phase5)
                transmission_phase5_max = phase5_starting_point + phase5_slope * (t_phase5_max - t_phase5)
                beta = if_else(t >= t_phase5_max, transmission_phase5_max, transmission_phase5)
            }
            beta
        } -> new_scalings  
    } ->temp
    temp = data.frame(temp)
    row.names(temp) = 1:500
    names(temp) = 1:11
    return(temp)
}

get_foi_string <- function(scaled_beta,
                t,
                ag, 
                region_num,
                population,
                q_vec,
                contact_covars) {
    foi_string = list()
    for (jj in seq(1, 9)){
        C_val = contact_covars[t, sprintf('C_%s_%s_%s', ag, jj, region_num)]
        infectious = sprintf("(%s+%s+%s+%s)", paste0('P', jj), paste0('A', jj), paste0('IM', jj), paste0('IS', jj))
        foi = sprintf("(%s * %s * %s * %s / %s)", q_vec[ag], scaled_beta, C_val, infectious, population[jj, 'POPULATION'])
        foi_string[jj] = foi
    }
    
    paste(foi_string, collapse="+")
}


calculate_phi <- function(rho, kappa, psi1, psi2, frac_nonhosp_deaths){
    1 - ((kappa * (1 - psi1 - psi2) * frac_nonhosp_deaths) / ((1 - frac_nonhosp_deaths) * (1 - kappa)))
}

get_R0 <- function(region_name, beta_value,
                   pop,
                   paras,
                   t,
                   contacts){
    region_key = c(1, 1, 2, 4, 4, 2, 3, 3, 3, 3, 5)
    region_cons = region_key[region_name]
    paras$phi = calculate_phi(rho = as.numeric(paras$rho),
                                kappa=paras$kappa,
                                psi1=paras$psi1,
                                psi2=paras$psi2,
                                frac_nonhosp_deaths = 0.51)

    # Set up disease-free equilibrium
    deq = list()
    for (age_group in seq(1, 9)){
        deq[paste0('S', age_group)] = pop[age_group, 'POPULATION']
        deq[paste0('A', age_group)] = 0
        deq[paste0('E', age_group)] = 0
        deq[paste0('P', age_group)] = 0
        deq[paste0('IM', age_group)] = 0
        deq[paste0('IS', age_group)] = 0
    }

        # Generate expressions for the F matrix. Each entry describes the flow of NEWLY infected people into each infected class, so for our model, this only includes E.
        F_matrix = c()
        for (ag in seq(1, 9)){
            

             with(paras, {
                foi_string = get_foi_string(scaled_beta=beta_value, t=t, ag=ag, region_num=region_cons, population=pop, q_vec=qq, contact_covars=contacts)
                dE = parse(text=sprintf("(%s) * S%s", foi_string, ag))
                dA = parse(text="0")
                dP = parse(text="0")
                dIM = parse(text="0")
                dIS = parse(text="0")
                l = list(E=dE, A=dA, P=dP, IM=dIM, IS=dIS)
                names(l) = paste0(names(l),ag)
                l
            }) -> elist_temp
           F_matrix = c(F_matrix, elist_temp)
        }
        
        # Generate expressions for the V matrix. Each entry describes the flow between infectected compartments. Includes E. Does NOT include the generation of new infections.
        V_matrix = c()
        for (ag in seq(1, 9)){
             with(paras, {
                foi_string = get_foi_string(scaled_beta=beta_value, t=t, ag=ag, region_num=region_cons, population=pop, q_vec=qq, contact_covars=contacts)
                dE = parse(text=sprintf("-%s * E%s", sigma[ag], ag))
                dA = parse(text=sprintf("%s * %s * E%s - %s * A%s", rho[ag], sigma[ag], ag, eta[ag], ag))
                dP = parse(text=sprintf("(1-%s) * %s * E%s - %s * P%s", rho[ag], sigma[ag], ag, zeta_s[ag], ag))
                dIM = parse(text=sprintf("(1-%s) * %s * P%s - IM%s * (%s * %s + (1-%s) * %s)", kappa[ag], zeta_s[ag], ag, ag, phi[ag], gamma_m[ag], phi[ag], mu_m[ag]))
                dIS = parse(text=sprintf("%s * %s * P%s - %s * IS%s", kappa[ag], zeta_s[ag] , ag, zeta_h[ag], ag))
                
                l = list(E=dE, A=dA, P=dP, IM=dIM, IS=dIS)
                names(l) = paste0(names(l),ag)
                l
            }) -> elist_temp
           V_matrix = c(V_matrix, elist_temp)
        }
        
        # To evaluate at the disease-free equilibrium, i.e. to compute R0.
        JJ_F=jacobian(states=names(F_matrix), elist=F_matrix, pts=deq)
        JJ_V=jacobian(states=names(V_matrix), elist=V_matrix, pts=deq)
        
        R0 = eigen(JJ_F %*% -inv(JJ_V))$values[1]
        R0
}



root <- '../../'
source(file.path(root, '_covid_root.R'))
covid_set_root(root)

source('./input_file_specification_11regions.R')
default_par_file = './default_parameter_values.csv'
output_dir = projection_dir
source(covid_get_path(inference_file))
source('set_up_covariates_and_data_11regions.R')


#pop_totals = as.numeric(lapply(FUN=sum, population_list))
#    chicago_pop = pop_totals[5]
#    pop_totals = pop_totals[1:4]
#    pop_totals = c(pop_totals, chicago_pop, sum(pop_totals))
#    names(pop_totals) = c(region_order[1:4], 'chicago', 'Illinois')

pop_fracs = population_list %>% 
    group_by(covid_region) %>%
    summarize(total=sum(POPULATION)) %>% 
    ungroup() %>%
    mutate(pop_frac = total / sum(total)) %>%
    arrange(as.numeric(covid_region)) %>%
    select(pop_frac) %>%
    unlist(use.names=F)

print(pop_fracs)



par_frame = read.csv(default_par_file)
pars = as.numeric(par_frame$value)
names(pars) = par_frame$param_name
pars = as.list(pars)

parsets = read.csv('./full.final_points.csv')

contacts = pomp_contacts
#fnhd = as.numeric(pars[paste0('region_non_hosp_', seq(1, 5))])

# Get region number
#region_names = c('restore_northcentral', 'restore_central','restore_northeast','restore_southern', 'chicago')



kappa = (as.numeric(pars[paste0('IHR_logit_', seq(1, 9))]) + as.numeric(pars[paste0('IHR_min_', seq(1, 9))])) / (1 + as.numeric(pars[paste0('IHR_min_', seq(1, 9))])) * as.numeric(pars[paste0('IHR_max_', seq(1, 9))]) / (1 - as.numeric(pars[paste0('rho_', seq(1, 9))]))
zeta_h = (pars[['zeta_h_logit']] + pars[['zeta_h_min']]) / (1 + pars[['zeta_h_min']]) * pars[['zeta_h_max']]
print(1/zeta_h)
R0_pars <- list(
        sigma = rep(1/pars$inv_sigma, 9),
        eta = rep(1/pars$inv_eta, 9),
        zeta_s = rep(1/pars$inv_zeta_s, 9),
        zeta_h= rep(zeta_h, 9),
        mu_m = rep(1/pars$inv_mu_m, 9),
        gamma_m = rep(1/pars$inv_gamma_m, 9),
        kappa = kappa,
        rho = as.numeric(pars[paste0('rho_', seq(1, 9))]),
        psi1 = as.numeric(pars[paste0('psi1_', seq(1, 9))]),
        psi2 = as.numeric(pars[paste0('psi2_', seq(1, 9))]),
        qq=as.numeric(pars[paste0('q_', seq(1, 9))]),
        f_nurse=c(0,0,0,0,0,0,0.0117,0.03159,0.11115))



cl <- makeCluster(28)
registerDoParallel(cl)

tmin = simstart
tmax = as.numeric(Sys.Date() - as.Date(t_ref))
tstartloop = Sys.time()
foreach (parset=1:nrow(parsets), .combine='rbind') %do%{
    
    # reassign parameters
    for(n in names(parsets)){
        pars[[n]] = parsets[parset, n]
    }
    
    # Calculate betas
    if (use_changepoint){
       betas = get_beta_from_parameters(pars) 
    } else{
        betas = get_beta_from_covar(pars, beta_covariate_file=beta_covariate_file, beta_covariate_column=beta_covariate_column)
    }

    foreach(rnum=1:11, .combine = 'rbind') %do%{
        print(rnum)
        # Calculate R0 for each timepoint and region
        foreach(time=tmin:tmax, .combine='rbind', 
            .packages=c('dplyr','matlib')) %dopar%{
            scaled_beta = betas[time, rnum]
            print(scaled_beta)
            R0 = Re(get_R0(region_name=rnum,
               pop = population_list %>% filter(covid_region == rnum) %>% arrange(AGE_GROUP_MIN),
               beta_value = scaled_beta,
               paras=R0_pars,
               t=time,
               contacts=contacts))
            c(time, R0, rnum, parset)
        } -> R0_vals
        R0_vals = data.frame(R0_vals)
        names(R0_vals) = c('time','R0','Region','parset')
        last_row = R0_vals[rep(nrow(R0_vals), 300), ]
        rbind(R0_vals, last_row) %>% 
            mutate(time=tmin:(tmin+nrow(.) -1))     
    } -> R0_consolidated
    R0_consolidated = data.frame(R0_consolidated)

    R0_consolidated
} -> final_R0
print(Sys.time()-tstartloop)

IL_R = final_R0 %>% 
    group_by(time, parset) %>% 
    summarize(R0=sum(R0*pop_fracs[Region])) %>% 
    ungroup() %>% 
    mutate(Region=0)

final_R0 = rbind(final_R0, IL_R)
write.csv(final_R0, paste0(projection_dir, 'R0_values.csv'), row.names=F)