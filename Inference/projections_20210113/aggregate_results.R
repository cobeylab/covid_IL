library(dplyr)
library(foreach)

source('./input_file_specification_no_age.R')

file_list = list.files(path=output_dir, pattern='.*.csv')
file_list = paste0(output_dir, file_list)
outdf = foreach(i=1:length(file_list), .combine='rbind') %do%{
    read.csv(file_list[i])
}

write.csv(outdf %>% arrange(-loglik), sprintf('agg_mif_%s.csv', model_name), row.names=F)
