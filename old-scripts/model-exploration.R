library(tidyverse)

source("run-model.R")
source("model-helpers.R")

# TODO: Nandita, try to play around with values of s1, s2, and R[3, 3] to find settings
# for all three scenarios where VF, VC, or IC dominate

#  parameter settings-------------------------------------------------------------
h <- 0.05 # h -> 1 means more penalty for heterozygotes, compared to corresponding homozygote
s1 <- 0.4 # s1 -> 1 more penalty for VC haplotype
s2 <- 0.8 # s2 -> 1 more penalty for IC haplotype
W <- fitness_weights(h, s1, s2) # fitness matrix

R <- matrix(1, 3, 3) # resistance matrix
R[3, 3] <- 1 # for now, assume only homozygotes in IC enjoy any immunity from insecticides

M <- 1 # maturation matrix (we've been assuming this constant throughout and probably won't use it)

#  initial condition--------------------------------------------------------------
num_gen <- 20 # number of generations to run experiment
num_adults <- 100
# number of adults of each genotype - contributions: female (rows) and males (cols)
A <- matrix(0, 3, 3)
colnames(A) <- rownames(A) <- haplotypes
diag(A) <- num_adults * c(0.1, 0.3, 0.6) # assume there are more resistant mosquitoes, and only homozygotes initially

# run the deterministic and stochastic versions for `num_gen` generations----------------

res_det <- run_model(A, R, W, M, num_gen, num_adults=num_adults, gfun=run_gen)

res_sto <- map_dfr(1:500, \(rep) { # run for 500 replications
    r <- run_model(A, R, W, M, num_gen, gfun=run_gen_stochastic, num_adults=num_adults)
    tibble_row(res=list(r), rep=rep)
})

# plot the results by haplotype---------------------------------------------------
res_sto_hap <- res_sto |> 
    unnest(res) |> 
    group_by(rep) |> 
    reframe(model_summary(res, mat2haplotype)) |> 
    group_by(haplotype, generation) |> 
    summarize(ymin=quantile(value, 0.025), ymax=quantile(value, 0.975), value=mean(value), .groups="drop")

mean_hap <- res_sto_hap |> 
    select(-ymin, -ymax) |> 
    mutate(type="stochastic avg.") |> 
    bind_rows(
        mutate(model_summary(res_det, mat2haplotype), type="deterministic")
    )

ggplot(mean_hap, aes(generation, value, col=haplotype)) +
    # geom_line(aes(group=interaction(haplotype, rep)), data=filter(res_sto, rep %in% sample(50, 5)), alpha=0.3) +
    geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=haplotype), res_sto_hap, alpha=0.25, col=NA) +
    geom_line(aes(linetype=type), linewidth=1.2) +
    labs(y="count", linetype=NULL) +
    theme_bw()

#  plot the results by genotype---------------------------------------------------
res_sto_gen <- res_sto |> 
    unnest(res) |> 
    group_by(rep) |> 
    reframe(model_summary(res, mat2genotype)) |> 
    group_by(genotype, generation) |> 
    summarize(ymin=quantile(value, 0.025), ymax=quantile(value, 0.975), value=mean(value), .groups="drop")

mean_gen <- res_sto_gen |> 
    select(-ymin, -ymax) |> 
    mutate(type="stochastic avg.") |> 
    bind_rows(
        mutate(model_summary(res_det, mat2genotype), type="deterministic")
    )

ggplot(mean_gen, aes(generation, value)) +
    geom_ribbon(aes(ymin=ymin, ymax=ymax), res_sto_gen, alpha=0.25, fill="gray60", col=NA) +
    geom_line(aes(linetype=type), linewidth=1.2) +
    facet_wrap(~genotype, nrow=3) +
    labs(y="count", linetype=NULL) +
    theme_bw()
