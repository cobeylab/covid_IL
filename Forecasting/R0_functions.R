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


get_beta_from_parameters = function(input_params){
    t_sip = 62
    foreach(region=1:5, .combine='cbind') %do%{
        t_phase3_max = input_params[[paste0('t_phase3_max_', region)]]
        t_phase4_max = input_params[[paste0('t_phase4_max_', region)]]

        scale_phase3 = input_params[[paste0('scale_phase3_', region)]]
        scale_phase4 = input_params[[paste0('scale_phase4_', region)]]

        t_phase3 = input_params[[paste0('t_phase3_', region)]]
        t_phase4 = input_params[[paste0('t_phase4_', region)]] 

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
            } else{
                phase4_change = scale_phase4 - scale_phase3
                phase4_starting_point = beta2 * (1 + scale_phase3)
                phase4_slope = phase4_change / (t_phase4_max - t_phase4)
                transmission_phase4 = phase4_starting_point + phase4_slope * (t - t_phase4)
                transmission_phase4_max = phase4_starting_point + phase4_slope * (t_phase4_max - t_phase4)
                beta = if_else(t >= t_phase4_max, transmission_phase4_max, transmission_phase4)
            }
            beta
        } -> new_scalings  
    } ->temp
    temp = data.frame(temp)
    row.names(temp) = 1:500
    names(temp) = 1:5
    return(temp)
}

get_foi_string <- function(scaled_beta,
                t,
                ag, 
                region_num,
                population,
                q_vec,
                f_nurse_vec,
                contact_covars) {
    foi_string = list()
    for (jj in seq(1, 9)){
        C_val = contact_covars[sprintf('C_%s_%s_%s', ag, jj, region_num), t]
        infectious = sprintf("(%s+%s+%s+%s)", paste0('P', jj), paste0('A', jj), paste0('IM', jj), paste0('IS', jj))
        foi = sprintf("(%s * %s * %s * %s / %s)", q_vec[ag], scaled_beta, C_val, infectious, population[jj, 'POPEST2018_CIV'])
        foi_string[jj] = foi
    }
    
    paste(foi_string, collapse="+")
}

get_single_beta_scale = function(date,
                     t_logistic_start,
                     intervention_lift,
                     max_scale){
    logistic = function(x, mscale=max_scale, shift=intervention_lift){
        mean = (mscale-1)/(1+exp(-(as.numeric(x - shift)))) + 1
        mean
    }

    if (date < t_logistic_start){
        scale_beta=1
    } else{
        scale_beta = logistic(date)
    }
    scale_beta
}



calculate_phi <- function(rho, kappa, psi1, psi2, frac_nonhosp_deaths){
    1 - ((kappa * (1 - psi1 - psi2) * frac_nonhosp_deaths) / ((1 - frac_nonhosp_deaths) * (1 - kappa)))
}

get_R0 <- function(region_name, beta_value, b_elderly,b_young, date,
                   paras, fnhd,
                   t,
                   contacts
                   t_intervention_start = as.Date('2020-03-16')){

    # Get region number
    region_names = c('restore_northcentral', 'restore_central','restore_northeast','restore_southern', 'chicago')
    region = match(region_name, region_names)

    # Population files
    population_files = c('../Data/IL_population_north-central.csv',
                         '../Data/IL_population_central.csv',
                         '../Data/IL_population_northeast.csv',
                         '../Data/IL_population_southern.csv',
                         '../Data/IL_population_chicago.csv')
    names(population_files) = region_names
    
    pop = read.csv(population_files[[region_name]])
    paras$phi = calculate_phi(rho = as.numeric(paras$rho),
                                kappa=paras$kappa,
                                psi1=paras$psi1,
                                psi2=paras$psi2,
                                frac_nonhosp_deaths = fnhd[region])

    # Set up disease-free equilibrium
    deq = list()
    for (age_group in seq(1, 9)){
        deq[paste0('S', age_group)] = pop[age_group, 'POPEST2018_CIV']
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
                foi_string = get_foi_string(beta_value, ag, contacts, pop, q_vec=qq)
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
                foi_string = get_foi_string(beta_value, ag, contacts, pop, q_vec=qq)
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