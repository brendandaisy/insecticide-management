library(tidyverse)

cols <- c("VV/FF", "VV/FC",	"VV/CC", "IV/FF", "IV/FC", "IV/CC",	"II/FF", "II/FC", "II/CC")
dat <- read_csv("Acp_site.csv")

res_geno <- model_summary(res, prop_genotype)

colnames(dat)[5:13] <- c("VFVF", "VFVC", "VCVC", "VFIF", "VFIC", "VCIC", "IFIF", "IFIC", "ICIC")

dat_proc <- dat |> 
    pivot_longer(VFVF:ICIC, names_to="genotype") |> 
    filter(genotype %in% unique(res_geno$genotype))



gg +
    geom_line(aes(Generation, value/N, group=interaction(genotype, Rep)), data=dat_proc, col="gray50")
