library(dplyr)
library(foreach)
library(tidyr)
library(pomp)

source('./input_file_specification_no_age.R')

file_list = Sys.glob(sprintf("%s/*_pfilter_*.rds", output_dir))

outdf = foreach(i=1:length(file_list), .combine='bind_rows') %do%{
    pf = readRDS(file_list[i])
    df = data.frame(t(pf@params))
    df %>% mutate(loglik.pf = pf@loglik, filename=file_list[i])
}

outdf %>% 
  group_by_at(vars(-filename, -loglik.pf)) %>% 
  summarize(n=length(loglik.pf),
            filename = filename[1],
            loglik.pf = logmeanexp(loglik.pf)) %>% 
  ungroup() %>%
  filter(!is.na(loglik.pf)) -> outdf

write.csv(outdf %>% arrange(-loglik.pf), sprintf('agg_mif_%s.csv', model_name), row.names=F)
