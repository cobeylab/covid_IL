library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)
library(matlib)

get_par_value = function( p0, pf,t0, tf, t_now ) {
   slope = (pf - p0) / (tf - t0)
   if ( t_now >= tf){
       value = pf
   } else if (t_now <= t0){
       value = p0
   } else{
       value = p0 + slope * (t_now - t0)
   }
  value
}
unlogit = function( logit_value, pmax, pmin) {
  (logit_value * (pmax - pmin)) + pmin
}


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


get_foi_string <- function(scaled_beta,
                t,
                ag, 
                region_num,
                population) {
    foi_string = list()
    for (jj in seq(1, 1)){
        infectious = paste(c(paste0('P', 1:3, jj), paste0('IM', 1:3, jj), paste0('IS', 1:3, jj), paste0('IMdead', 1:3, jj)), collapse ="+")
        foi = sprintf("(%s * (%s) / %s)", scaled_beta, infectious, population[jj, 'POPULATION'])
        foi_string[jj] = foi
    }
    
    paste(foi_string, collapse="+")
}


calculate_phi = function(f, HFR, IHR){
    (f * HFR * IHR) /(1 - IHR - f + IHR * f) 
}


get_R0 <- function(region_cons, 
                   beta_value,
                   pop,
                   paras,
                   t,
                   contacts,
                   IHR,
                   HFR){

    paras$phi = calculate_phi(f = paras$frac_nonhosp_deaths,
                              HFR = HFR,
                              IHR = IHR
                              )

    # Set up disease-free equilibrium
    deq = list()
    for (age_group in seq(1, 1)){
        deq[paste0('S', age_group)] = pop[age_group, 'POPULATION']
        deq[paste0('E', 1:3, age_group)] = 0
        deq[paste0('P', 1:3, age_group)] = 0
        deq[paste0('IMdead', 1:3, age_group)] = 0
        deq[paste0('IM', 1:3, age_group)] = 0
        deq[paste0('IS', 1:3, age_group)] = 0
    }

        # Generate expressions for the F matrix. Each entry describes the flow of NEWLY infected people into each infected class, so for our model, this only includes E.
        F_matrix = c()
        for (ag in seq(1, 1)){
             with(paras, {
                foi_string = get_foi_string(scaled_beta=beta_value, t=t, ag=ag, region_num=region_cons, population=pop)
                dE1 = parse(text=sprintf("(%s) * S%s", foi_string, ag))
                dE2 = parse(text="0")
                dE3 = parse(text="0")
                dP1 = parse(text="0")
                dP2 = parse(text="0")
                dP3 = parse(text="0")
                dIMdead1 = parse(text="0")
                dIMdead2 = parse(text="0")
                dIMdead3 = parse(text="0")
                dIM1 = parse(text="0")
                dIM2 = parse(text="0")
                dIM3 = parse(text="0")
                dIS1 = parse(text="0")
                dIS2 = parse(text="0")
                dIS3 = parse(text="0")
                l = list(E1=dE1, E2=dE2, E3=dE3, 
                         P1=dP1, P2=dP2, P3=dP3,
                         IMdead1=dIMdead1, IMdead2=dIMdead2, IMdead3=dIMdead3,
                         IM1=dIM1, IM2=dIM2, IM3=dIM3, 
                         IS1=dIS1, IS2=dIS2, IS3=dIS3)
                names(l) = paste0(names(l),ag)
                l
            }) -> elist_temp
           F_matrix = c(F_matrix, elist_temp)
        }
        
        # Generate expressions for the V matrix. Each entry describes the flow between infectected compartments. Includes E. Does NOT include the generation of new infections.
        V_matrix = c()
        for (ag in seq(1, 1)){
             with(paras, {
                foi_string = get_foi_string(scaled_beta=beta_value, t=t, ag=ag, region_num=region_cons, population=pop)
                dE1 = parse(text=sprintf("-%s * E1%s", sigma *3 , ag))
                dE2 = parse(text=sprintf("%s * E1%s - %s * E2%s",sigma*3, ag, sigma *3 , ag))
                dE3 = parse(text=sprintf("%s * E2%s -%s * E3%s", sigma*3, ag, sigma *3 , ag))
                
                dP1 = parse(text=sprintf("%s * E3%s - %s * P1%s", sigma * 3, ag, zeta_s * 3, ag))
                dP2 = parse(text=sprintf("%s * P1%s - %s * P2%s", zeta_s * 3, ag, zeta_s * 3, ag))
                dP3 = parse(text=sprintf("%s * P2%s - %s * P3%s", zeta_s * 3, ag, zeta_s * 3, ag))
                
                dIMdead1 = parse(text=sprintf("(1-%s) * %s * %s * P3%s - IMdead1%s * %s", IHR, zeta_s *3, phi[ag], ag, ag, mu_m *3))
                dIMdead2 = parse(text=sprintf("%s * IMdead1%s - IMdead2%s * %s", mu_m *3, ag, ag, mu_m *3))
                dIMdead3 = parse(text=sprintf("%s * IMdead2%s - IMdead3%s * %s", mu_m*3, ag, ag, mu_m *3))
                
                dIM1 = parse(text=sprintf("(1-%s) * %s * (1-%s) * P3%s - IM1%s * %s", IHR, zeta_s * 3, phi[ag], ag, ag, gamma_m *3))
                dIM2 = parse(text=sprintf("%s * IM1%s - IM2%s * %s",  gamma_m*3, ag, ag, gamma_m *3))
                dIM3 = parse(text=sprintf("%s * IM2%s - IM3%s * %s", gamma_m*3, ag, ag, gamma_m *3))
                
                dIS1 = parse(text=sprintf("%s * %s * P3%s - %s * IS1%s", IHR, zeta_s *3, ag, zeta_h *3, ag))
                dIS2 = parse(text=sprintf("%s * IS1%s - %s * IS2%s", zeta_h *3, ag, zeta_h *3, ag))
                dIS3 = parse(text=sprintf("%s * IS2%s - %s * IS3%s", zeta_h *3, ag, zeta_h *3, ag))
                
                l = list(E1=dE1, E2=dE2, E3=dE3, 
                         P1=dP1, P2=dP2, P3=dP3,
                         IMdead1=dIMdead1, IMdead2=dIMdead2, IMdead3=dIMdead3,
                         IM1=dIM1, IM2=dIM2, IM3=dIM3, 
                         IS1=dIS1, IS2=dIS2, IS3=dIS3)
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
