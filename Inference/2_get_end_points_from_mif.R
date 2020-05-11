library(dplyr)
library(foreach)
library(tidyverse)

select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise
contains <- dplyr::contains

setwd('./results/')
file_list = list.files(pattern='*.rds')

foreach(i=1:length(file_list), .combine='rbind') %do%{
    mf <- readRDS(file_list[i])
    mf_params <- as.data.frame(mf@traces) %>% mutate(iter=1:nrow(mf@traces)) %>% select(beta1, beta2_1, beta2_2, beta2_3, num_init_1, num_init_2, num_init_3, loglik, iter)  
    mf_params
    } -> plotliks


plotliks %>% drop_na() -> outdf
best_points <- outdf %>% filter(loglik > max(loglik) - 50) %>% arrange(desc(loglik))
write.csv(best_points, '../mif_50iter_best_points.csv')