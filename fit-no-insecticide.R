library(tidyverse)

source("model-helpers.R")
source("base-model.R")

dat <- read_csv("Acp_site.csv") |> 
    rename(site=Site, t=Generation, rep=Rep)

colnames(dat)[5:13] <- c("VFVF", "VFVC", "VCVC", "IFVF", "VFIC", "VCIC", "IFIF", "IFIC", "ICIC")

dat_summ <- dat |> 
    pivot_longer(VFVF:ICIC, names_to="genotype") |> 
    filter(genotype %in% unique(res_geno$genotype)) |> 
    group_by(genotype, t) |> 
    summarise(prop=sum(value) / sum(N), .groups="drop")

R <- 1
M <- 1
num_adults <- 150
Ainit <- init_adults(dat_summ)
Ainit[1, 1] <- 1
Ainit[3, 3] <- Ainit[3, 3] - 1

## heatmap of SSE over s and h
params <- expand_grid(s=seq(0, 1, 0.1), h=seq(0, 1, 0.1))

param_sse <- params |> 
    mutate(
        sse=map2_dbl(s, h, \(s, h) {
            W <- fitness_weights(s, h, type=1)
            res <- run_model(Ainit, 8, list(R=R, W=W, M=M))
            res_geno <- model_summary(res, prop_genotype)
            comb <- left_join(dat_summ, res_geno, by=c("genotype", "t"))
            sum((comb$prop.x - comb$prop.y)^2)
        })
    )
    
ggplot(param_sse, aes(s, h, fill=sse)) +
    geom_raster()

## compare best results with data
slice_min(param_sse, sse, n=1)
W <- fitness_weights(0.8, 0.5, type=1)
res <- run_model(Ainit, 10, list(R=R, W=W, M=M))
res_geno <- model_summary(res, prop_genotype)

ggplot(res_geno, aes(t, prop, col=genotype)) +
    geom_line(alpha=0.5, linewidth=1.2) +
    # geom_point() +
    geom_line(aes(group=genotype), data=dat_summ, col="white", linewidth=1.1) +
    geom_line(data=dat_summ, linetype="dashed", linewidth=1.1, alpha=0.8) +
    theme_bw()
