## Estimate with epinow2
outpath = '../epinow2_cli_estimates'

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(cowplot)
#library(EpiEstim)
library(optparse)
theme_set(theme_bw())
source('../code/util.R')
source('../code/load_timeseries.R')
dt <- max(load_cli()$date) ## Use the last date in the timeseries to set the output folder name.
tooday <- Sys.Date()


## Read in options from midway
option_list = list(make_option("--var", type = "numeric", default=NULL, help="array_task_number"),
                   make_option("--debug", type = 'logical', default = TRUE, help='if debug=T, run very short chains'),
                   make_option("--midway", type = 'logical', default = TRUE, help='T is running on midway'),
                   make_option("--outpath", type = "character", default = NULL, help = 'optional outpath spec for testing. Default is epinow2_cli_estimates/')); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); # Now you have a list called "opt" with elements opt$var and opt$out
midway = ifelse(length(opt$midway)>0, TRUE, FALSE) ## Set whether running on midway
outpath = ifelse(length(opt$outpath)>0, opt$outpath, outpath) ## Set outpath, letting user input override default if specified
dir_check(outpath)
if(opt$debug == TRUE){
  tooday='debug'
  outpath = paste0(outpath, '/debug')
}

dir_check(outpath)
dir_check(sprintf('%s/%s', outpath, dt))


## Load data --------------------------------------------
dat <- load_cli()
## Plot data -------
dat %>%
  pivot_longer(c(nadmit, smoothed)) %>%
  ggplot()+
  geom_line(aes(x = date, y = value, color = name))+
  facet_wrap(.~region, scales = 'free_y')+
  scale_color_manual("", values = c('black', 'red', 'yellow'))+
  theme(legend.position = 'bottom')+
  ggtitle(sprintf('CLI up to %s', dt))
ggsave(sprintf('%s/%s/cli_input.png', outpath, dt))

## Estimate rt using epinow2
library('EpiNow2')
source('../code/epinow2.R')
get_region <- function(rr){
  dat %>% filter(region == rr) %>%
    group_by(date) %>%
    summarise(nadmit = sum(nadmit),
              smoothed = sum(smoothed))
}


regions = unique(dat$region)
if(midway){ ## If running on midway
  cat('midway is true\n')
  regions = regions[opt$var] ## Only run for the current slurm array task id
  print(regions)
} ## Else, run for all regions.

for(region.in in regions){
  cat(sprintf('Running for region %s', region.in))
  dir_check(sprintf('%s/%s/%s', outpath, dt, region.in))
  run_epinow2(dat_df = get_region(region.in),
              obs_colname = 'nadmit',
              prior_smoothing_window = 1,
              output_folder = sprintf('%s/%s', outpath, region.in), 
              dbug = opt$debug)
}

cat('sampler has run')
runpar_file = sprintf('run_params_cli_%s.rds', as.numeric(Sys.getenv("SLURM_ARRAY_JOB_ID")))
if(!file.exists(runpar_file)){
  write_rds(list(outpath=outpath,
                 dt=dt,
                 today=tooday),
            path = runpar_file)
  cat('runpars written')
}


