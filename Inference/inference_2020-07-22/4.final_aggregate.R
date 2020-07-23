library(dplyr)
library(data.table)
library(foreach)
library(pomp)

select <- dplyr::select
# Read in arguments
args = commandArgs(trailingOnly=TRUE)
#output_dir = args[1]

# Aggregate points
print('Aggregating points')
parnames = c('beta1_logit', 'beta2', 'num_init', 'scale_phase3')
output_filename <- paste0(output_dir, "consolidated.final.pfilter.csv")
file_path = output_dir #Replace this with a file path to the directory storing the csv files 
file_list <- Sys.glob(paste0(output_dir, "/output*.csv"))
print(length(file_list))
data <- rbindlist(lapply(file_list, fread))
df_output <- as.data.frame(data) %>% arrange(desc(loglik)) 
write.csv(df_output, file = paste0(output_filename))


all_pars = unique(do.call(paste0, expand.grid(expand.grid(parnames,'_',1:5))))
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

# produce mle parfile for easy simulation
default_pars = read.csv('./default_parameter_values.csv')
rownames(default_pars) = default_pars$param_name

for(reg in 1:5){
    cols = c(paste0(c('beta1_logit', 'beta2', 'num_init', 'scale_phase3'), '_', reg))
    pars = global_mle %>%
        filter(region_to_test == reg) %>%
        select(cols) %>% 
        unlist()
    default_pars[names(pars), 'value'] = pars
}

write.csv(default_pars, './final_mle_pars.csv', row.names=F)

# Produce sampling df
foreach (reg=1:5, .combine='cbind') %do%{
    cols = c(paste0(c('beta1_logit', 'beta2', 'num_init', 'scale_phase3'), '_', reg))

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
     expanded
} -> concatenated

final = concatenated %>% group_by_all() %>% summarize(num_sims=length(beta2_1)) %>% ungroup() %>% mutate(parset = 1:nrow(.))
write.csv(final, file="final_points.csv")