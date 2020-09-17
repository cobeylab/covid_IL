library(dplyr)
library(data.table)
library(foreach)
library(pomp)

select <- dplyr::select
# Read in arguments
args = commandArgs(trailingOnly=TRUE)
output_dir = args[1]
source('./input_file_specification.R')

# Aggregate points
print('Aggregating points')

if (use_changepoint){
  parnames = c('beta1', 'beta2', 'num_init', 'scale_phase3', 'scale_phase4', 'scale_phase5') 
} else{
  parnames = c('beta1', 'num_init')
}

output_filename <- paste0(output_dir, "consolidated.final.pfilter.csv")
file_path = output_dir #Replace this with a file path to the directory storing the csv files 
file_list <- Sys.glob(paste0(output_dir, "/output*.csv"))
print(length(file_list))
data <- rbindlist(lapply(file_list, fread))
df_output <- as.data.frame(data) %>% arrange(desc(loglik)) 
write.csv(df_output, file = paste0(output_filename))

# Find MLE for each region
mles = df_output %>% 
    group_by(region_to_test) %>%
    filter(loglik == max(loglik))

print('Getting grid')
# Get a grid of points around each mle
foreach(r=1:nrow(mles), .combine='rbind') %do%{
    reg = as.numeric(mles[r, 'region_to_test'])
    pars = mles[r, paste0(parnames,'_',reg)] %>% unlist()
    names(pars) = parnames
    pars$region_to_test = reg
    as.data.frame(pars)
} -> final_design


## Points for grid search
saveRDS(final_design, 'final.mle_grid.rds')

all_pars = unique(do.call(paste0, expand.grid(expand.grid(parnames,'_',1:6))))
# Find MLE for each region
mles = df_output %>% 
    group_by(region_to_test) %>%
    filter(loglik > max(loglik) - 5) %>%
    ungroup() %>%
    select(region_to_test, loglik, loglik_se, all_pars) %>%
    arrange(region_to_test, loglik)

global_mle = mles %>% 
    group_by(region_to_test) %>%
    filter(loglik == max(loglik)) %>%
    ungroup() %>%
    select(region_to_test, loglik, loglik_se, all_pars) %>%
    arrange(region_to_test, loglik)

global_ll = sum(global_mle$loglik)

# produce mle parfile for easy simulation
default_pars = read.csv('./default_parameter_values.csv')
rownames(default_pars) = default_pars$param_name

for(reg in 1:6){
    cols = c(paste0(parnames, '_', reg))
    pars = global_mle %>%
        filter(region_to_test == reg) %>%
        select(cols) %>% 
        unlist()
    default_pars[names(pars), 'value'] = pars
}

default_pars['loglik', 'value'] = global_ll
default_pars['npars_estimated', 'value'] = length(parnames) * 5

write.csv(default_pars, './final_mle_pars.csv', row.names=F) # final parameters of MLE

# Produce sampling df
foreach (covidreg=1:11, .combine='cbind') %do%{
    reg = aggregated_regions[covidreg]
    cols = c(paste0(parnames, '_', reg))

    mles %>% filter(region_to_test==reg) -> df_output
    
    num_sims = 100
    if (sum(exp(df_output$loglik)) == 0){
        probs = df_output$loglik / max(df_output$loglik)
    } else{
        probs = exp(df_output$loglik)
    }
    tot_sims = rmultinom(1, num_sims, probs)
    
    df_output$num_sims = tot_sims
    df_output = df_output %>% select(-loglik, -loglik_se, -region_to_test)

    
     foreach(row=1:nrow(df_output), .combine='rbind') %do%{
            nsim = df_output[row, 'num_sims']
            df_output[rep(row, nsim), c(cols)]
        } -> expanded

     names(expanded) = paste0(parnames, '_', covidreg)
     expanded
} -> concatenated

final = concatenated %>% group_by_all() %>% 
    summarize(num_sims=length(beta1_1)) %>% 
    ungroup() %>% 
    mutate(parset = 1:nrow(.))
write.csv(final, file="final_points.csv") # File that specifies how many simulations to do