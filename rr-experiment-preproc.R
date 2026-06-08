library(tidyverse)

rr0 <- read_csv("data-raw/Resistance.Reversal_ALL.kdr_2025-07-24.csv")

# TODO: we'll need to work out (if possible...) the remaining genotype to haplotypes,
# since you and Silvie only did the starting six...
# -> looks like some imply recombination!!!!

rr0 |> 
    select(-1) |> 
    mutate(
        gen_410=str_replace_all(gen_410, c("S"="V", "R"="L")),
        gen_1016=str_replace_all(gen_1016, c("S"="V", "R"="I")),
        genotype=str_c(gen_410, gen_1016, gen_1534, sep="-")
    ) |> 
    count(genotype, sort=TRUE)
