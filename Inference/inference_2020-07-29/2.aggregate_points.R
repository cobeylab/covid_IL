library(dplyr)
library(data.table)
library(foreach)
library(pomp)

select <- dplyr::select
# Read in arguments
args = commandArgs(trailingOnly=TRUE)
output_dir = args[1]
num_points=2500 # number of points to search in the grid
source('./input_file_specification.R')
# Aggregate points
print('Aggregating points')
if (use_changepoint){
  parnames = c('beta1', 'beta2', 'num_init', 'scale_phase3') 
} else{
  parnames = c('beta1', 'num_init')
}

output_filename <- paste0(output_dir, "consolidated.mif.pfilter.csv")
file_path = output_dir #Replace this with a file path to the directory storing the csv files 
file_list <- Sys.glob(paste0(output_dir, "/output_mif*.csv"))
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

    lower_pars = pars * 0.1
    upper_pars = pars * 1.1

    names(lower_pars) = names(upper_pars) = parnames
    design = sobolDesign(upper=upper_pars,
        lower=lower_pars,
        num_points)
    design$region_to_test = reg
    design
} -> final_design


## Points for grid search
saveRDS(final_design, 'mle_grid.rds')

