library(dplyr)
library(data.table)
#library(RSQLite)
select <- dplyr::select
# Read in arguments
args = commandArgs(trailingOnly=TRUE)
output_dir = args[1]

# Aggregate points
print('Aggregating points')

#output_filename <- paste0(output_dir, "consolidated.all.pfilter.db")
#conn <- dbConnect(SQLite(), output_filename)

#file_path = output_dir #Replace this with a file path to the directory storing the csv files 
#setwd(file_path)
#file_list <- list.files(pattern = "output_pfilter_.*.csv$")

#for (f in file_list){
#    df = read.csv(f)
#    df = df %>% select(-one_of(names(df)[grepl('C_.*', names(df))])) # Exclude contact matrix parameters
#    dbWriteTable(conn, "pfilter_results", df, append=T)
#}

output_filename <- paste0(output_dir, "consolidated.all.pfilter.csv")
file_path = output_dir #Replace this with a file path to the directory storing the csv files 
setwd(file_path)
file_list <- list.files(pattern = "output_pfilter_.*.csv$")
data <- rbindlist(lapply(file_list, fread))
df_output <- as.data.frame(data) %>% arrange(desc(loglik)) 
write.csv(df_output, file = paste0('../', output_filename))
