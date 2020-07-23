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

get_foi_string <- function(beta,
                ag, 
                contact_matrix,
                population,
                q_vec,
                f_nurse_vec,
                b_elderly,
                b_scale_young=1) {
    foi_string = list()
    for (jj in seq(1, 9)){
        infectious = sprintf("(%s+%s+%s+%s)", paste0('P', jj), paste0('A', jj), paste0('IM', jj), paste0('IS', jj))
        if (ag %in% seq(7,9) & jj %in% seq(7,9)){
            c_nurse = b_elderly * f_nurse_vec[ag]
        } else{
            c_nurse = 0
        }
        
        if (ag %in% seq(3,6)){
            b_scale = b_scale_young
        } else{
            b_scale = 1
        }

        foi = sprintf("(%s * %s * %s * %s / %s)", q_vec[ag], beta * b_scale, (contact_matrix[ag, jj] + c_nurse), infectious, population[jj, 'POPEST2018_CIV'])
        foi_string[jj] = foi
    }
    
    paste(foi_string, collapse="+")
}

get_contact_matrix <- function(flat_matrices,
                               preSIP=T,
                               region){
    mat_start = 1 + (81 * (region - 1))
    mat_end = 81 + (81 * (region - 1))
    if(preSIP){
        C = flat_matrices$work[mat_start:mat_end] + flat_matrices$school[mat_start:mat_end] + flat_matrices$home[mat_start:mat_end] + flat_matrices$other[mat_start:mat_end]
    }else{
        C = 0.6 * flat_matrices$work[mat_start:mat_end] + flat_matrices$home[mat_start:mat_end] + 0.5 * flat_matrices$other[mat_start:mat_end]
    }
        matrix(C, ncol=9, nrow=9)
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
                   flat_contacts=pomp_contacts,
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

    # Determine if in SIP
    if(date < t_intervention_start){
        preSIP = T
    } else{
        preSIP = F
    }
    
    # Get contact matrix
    contacts = get_contact_matrix(flat_contacts, preSIP=preSIP, region=region)
    
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
                foi_string = get_foi_string(beta_value, ag, contacts, pop, q_vec=qq, f_nurse_vec=f_nurse, b_elderly=b_elderly, b_scale_young=b_young)
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
                foi_string = get_foi_string(beta_value, ag, contacts, pop, q_vec=qq, f_nurse_vec=f_nurse, b_elderly=b_elderly)
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