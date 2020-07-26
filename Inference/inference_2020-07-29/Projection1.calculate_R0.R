library(dplyr)
library(tidyr)
library(foreach)

dir.create(file.path(projection_dir), showWarnings = FALSE)

root <- '../../'
source(file.path(root, '_covid_root.R'))
covid_set_root(root)
source(covid_get_path('Forecasting/simulation_statewide.R'))
source(covid_get_path('Forecasting/R0_functions.R'))
source('./input_file_specification.R')
region_order = c('northcentral','central','northeast','southern', 'chicago')

default_par_file = './final_mle_pars.csv'
fit_df = read.csv('./final_points.csv')

output_file = paste0(projetion_dir,'/R0_over_time.csv')

load(contact_filename)
### User-specification of parameters and interventions for model 
par_frame = read.csv(default_par_file)
pars = as.numeric(par_frame$value)
names(pars) = par_frame$param_name
betamax = pars[['beta1_max_1']]

pars = as.list(pars)
population1 = read.csv(population_filename_1)
colnames(population1) <- c("AGE_GROUP_MIN", "POPULATION")
population2 = read.csv(population_filename_2)
colnames(population2) <- c("AGE_GROUP_MIN", "POPULATION")
population3 = read.csv(population_filename_3)
colnames(population3) <- c("AGE_GROUP_MIN", "POPULATION")
population4 = read.csv(population_filename_4)
colnames(population4) <- c("AGE_GROUP_MIN", "POPULATION")
population5 = read.csv(population_filename_5)
colnames(population5) <- c("AGE_GROUP_MIN", "POPULATION")
n_age_groups = nrow(population1)

# ASSUME FOLLOWING ORDER: NORTH-CENTRAL, CENTRAL, NORTHEAST, SOUTHERN 
population_list=list(population1, population2, population3, population4, population5)   

 pop_totals = as.numeric(lapply(FUN=sum, population_list))
    chicago_pop = pop_totals[5]
    pop_totals = pop_totals[1:4]
    pop_totals = c(pop_totals, chicago_pop, sum(pop_totals))
    names(pop_totals) = c(region_order[1:4], 'chicago', 'Illinois')

foreach (parset_id=1:max(fit_df$parset), .combine='rbind') %do%{
    pars_to_use = fit_df %>% filter(parset == parset_id) %>% unlist(use.names=T)
    pars[names(pars_to_use)] = pars_to_use
    
    # Read in parameters
    R0_pars <- list(
                sigma = rep(1/pars$inv_sigma, 9),
                eta = rep(1/pars$inv_eta, 9),
                zeta_s = rep(1/pars$inv_zeta_s, 9),
                zeta_h= rep(1/pars$inv_zeta_h, 9),
                mu_m = rep(1/pars$inv_mu_m, 9),
                gamma_m = rep(1/pars$inv_gamma_m, 9),
                kappa = as.numeric(pars[paste0('kappa_', seq(1, 9))]),
                rho = as.numeric(pars[paste0('rho_', seq(1, 9))]),
                psi1 = as.numeric(pars[paste0('psi1_', seq(1, 9))]),
                psi2 = as.numeric(pars[paste0('psi2_', seq(1, 9))]),
                qq=as.numeric(pars[paste0('q_', seq(1, 9))]),
                f_nurse=c(0,0,0,0,0,0,0.0117,0.03159,0.11115))


    R0_regions = c('restore_northcentral', 'restore_central','restore_northeast','restore_southern', 'chicago')
    fnhd = as.numeric(pars[paste0('region_non_hosp_', seq(1, 5))])
    R0_scales_reg = c(0.1, 0.2, 0.05, 0.2, 0.05)

    foreach (regionnum=1:length(R0_regions), .combine='rbind') %do%{
        region = R0_regions[regionnum]
        b_pre = pars[[paste0('beta1_logit_', regionnum)]] * betamax
        b_post = pars[[paste0('beta2_',regionnum)]]
        b_eld = pars$b_elderly
        pre_int = get_R0(region, b_pre, b_eld, b_young=1,as.Date('2020-03-01'), R0_pars, fnhd)
        
        R0_diff = ((pre_int - post_int) * R0_scales_reg[regionnum] + post_int) / post_int

        c(pre_int_R0=as.numeric(pre_int),
            post_int_R0=as.numeric(post_int),
            phase3_scale = 1 + pars[[paste0('scale_phase3_', regionnum)]],
            phase4_scale=R0_diff)

    } -> R0_values

    R0_values = as_tibble(R0_values) %>% 
        mutate(restore_region = str_replace(R0_regions, "restore_",""),
               pop_frac = as.numeric(pop_totals[restore_region]/pop_totals['Illinois']))

    R0_no_chicago = R0_values %>% filter(restore_region != 'chicago')

    IL_phase4= sum(as.numeric(R0_no_chicago$phase4_scale) * as.numeric(R0_no_chicago$pop_frac))
    IL_phase3 = sum(as.numeric(R0_no_chicago$phase3_scale) * as.numeric(R0_no_chicago$pop_frac))
    IL_pre = sum(R0_no_chicago$pre_int_R0 * R0_no_chicago$pop_frac)
    IL_post = sum(R0_no_chicago$post_int_R0 * R0_no_chicago$pop_frac)

    R0_values = rbind(R0_values, c(pre_int_R0=IL_pre, post_int_R0=IL_post,  phase3_scale=IL_phase3, phase4_scale=IL_phase4, restore_region='Illinois',  pop_frac=1))
    R0_values = as.data.frame(R0_values)
    R0_values$parset = parset_id
    R0_values
} ->final_R0

write.csv(final_R0, output_file, row.names=F)