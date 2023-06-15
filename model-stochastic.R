library(tidyverse)

num_adults <- 600
s <- 0.4
h <- 0.05
H <- matrix(c(
    0, 0, h,
    0, 0, h,
    h, h, 1
), nrow=3, byrow=TRUE)
W <- 1 - s*H
R <- matrix(0.55, nh, nh)
R[1, 1] <- 0.3
R[3, 3] <- 0.9
M <- matrix(1, nh, nh)

A <- matrix(0, nh, nh) # number of adults of each genotype - contributions: female (rows) and males (cols)
colnames(A) <- rownames(A) <- haplotypes
diag(A) <- c(500, 50, 50)

res_det <- run_model(A, 50, fun=run_generation)
res_sto <- map_dfr(1:50, \(rep) {
    r <- run_model(A, 50, fun=run_gen_stochastic)
    mutate(model_summary(r), rep=rep)
})

res <- res_sto |> 
    group_by(t, haplotype) |> 
    summarize(prop=mean(prop), .groups="drop") |> 
    relocate(haplotype, prop, t) |> 
    mutate(type="stochastic") |> 
    bind_rows(mutate(model_summary(res_det), type="deterministic"))

ggplot(res, aes(t, prop, col=haplotype)) +
    geom_line(aes(group=interaction(haplotype, rep)), data=filter(res_sto, rep %in% sample(50, 5)), alpha=0.3) +
    stat_smooth(
        aes(ymin=after_stat(ymin), ymax=after_stat(ymax), group=haplotype), data=res_sto, 
        geom="ribbon", col=NA, fill="gray50", alpha=0.3, linewidth=1.2
    ) +
    geom_line(aes(linetype=type), linewidth=1.2) +
    theme_bw()

ggplot(model_summary(res, prop_genotype), aes(t, prop, col=genotype)) +
    geom_line() +
    geom_point() +
    theme_bw()

ggplot(model_summary(res, prop_pairs), aes(t, prop, col=pair)) +
    geom_line() +
    geom_point() +
    theme_bw()