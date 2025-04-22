suppressPackageStartupMessages( library(tidyverse, quietly = TRUE) )

# load data
data_sv = read_delim("output/gt1kb.cnv.bed", col_names = FALSE) %>% suppressMessages()
colnames(data_sv) = c("chr","start","end","sv_type")

# calculate length
data_sv = data_sv %>% mutate(length=end-start)
cat("======== SV lengths ========\n")
print(data_sv)

# generate summary table
cat("\n======== SV length summary stats ========\n")
options(echo=TRUE)
data_sv %>% summarise(
    avg_len = mean(length, na.rm=TRUE)
) %>% print()

data_sv %>% group_by(sv_type) %>% summarise(
    avg_len = mean(length, na.rm=TRUE)
) %>% print()