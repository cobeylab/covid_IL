library(dplyr)
library(foreach)
library(ggplot2)
library(tidyverse)

select <- dplyr::select
rename <- dplyr::rename
summarize <- dplyr::summarise
contains <- dplyr::contains

# columns to gather for plotting
columns = c('beta1', 'beta2_1', 'beta2_2', 'beta2_3', 'num_init_1', 'num_init_2', 'num_init_3', 'loglik', 'nfail')

args = commandArgs(trailingOnly=TRUE)
result_dir = args[1]

setwd(result_dir)
file_list = list.files(pattern='*.rds')

foreach(i=1:length(file_list), .combine='rbind') %do%{
    mf <- readRDS(file_list[i])
    lastiter = nrow(mf@traces)
    mf_params <- as.data.frame(mf@traces) %>% mutate(iter=1:lastiter, chain=i)  
    mf_params
    } -> plotliks


plotliks %>% drop_na() -> outdf
last_iter <- plotliks %>% filter(iter==lastiter)
write.csv(last_iter, '../mif_50iter_end_points.csv')

plotliks %>% group_by(chain) %>% ungroup() %>% drop_na() -> outdf
outdf.long <- outdf %>% gather(param, value, columns)
ggplot(outdf.long, aes(x=iter, y=value, group=chain)) + geom_line(alpha=0.05) + facet_wrap(~param,scales="free_y",ncol=2)
ggsave('../chains.png', width=7, height=12)