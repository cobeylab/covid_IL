library(dplyr)
library(foreach)

source('./input_file_specification_no_age.R')

file_list = list.files(path=output_dir, pattern='.*.csv')
file_list = paste0(output_dir, file_list)
outdf = foreach(i=1:length(file_list), .combine='rbind') %do%{
    read.csv(file_list[i])
}

write.csv(outdf %>% arrange(-loglik), sprintf('agg_mif_%s.csv', model_name), row.names=F)

best = outdf %>% 
    group_by(region_to_test) %>%
    filter(loglik==max(loglik)) %>%
    ungroup() %>%
    select(region_to_test, loglik)


ll = sum(best$loglik)
npar = sum(unlist(outdf[1, paste0('n_changepoints_', 1:11)])) + sum(unlist(outdf[1, paste0('n_HFR_changepoints_', 1:11)])) + 2
print(npar)
print(2 * npar - 2 * ll)

