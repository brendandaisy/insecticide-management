library(tidyverse)
library(MASS, exclude=c("select"))
library(tikzDevice)

source("model-helpers.R")
source("base-model.R")

check_bounds <- function(params) {
    all(c(0, 0, 0, 0) <= params & params <= c(0.1, 1, 1, 1))
}

# log_posterior <- function(params, inits, data, num_gen=8) {
#     if (!check_bounds(params))
#         return(-Inf)
#     d <- params[1]
#     h <- params[2]
#     s1 <- params[3]
#     s2 <- params[4]
#     W <- fitness_weights(h, s1, s2)
#     
#     lpdf_reps <- map2_dbl(inits, data, \(irep, drep) {
#         A1 <- irep$A[[1]] + d
#         A1 <- A1 / sum(A1)
#         
#         res_sto <- map_dfr(1:50, ~{
#             res <- run_model(A1, 1, W, 1, num_gen)
#             model_summary(res, mat2genotype)
#         })
#         
#         res_geno <- res_sto |> 
#             group_by(genotype, generation) |> 
#             summarise(prop=mean(prop), .groups="drop") |> 
#             pull(prop)
#         
#         sum(map_dbl(1:(num_gen-1), \(g) {
#             slice <- 6*(g-1)+1:6
#             # TODO: the stochastic model could have more specific likelihood of choose 50 out of 500 (rather than proportions)
#             dmultinom(drep$count[slice], prob=res_geno[slice], log=TRUE)
#         }))
#     })
#     print(params)
#     print(sum(lpdf_reps))
#     return(sum(lpdf_reps)) # assuming "uniform" priors
# }

log_posterior <- function(params, inits, data, num_gen=8) {
    d <- params[1]
    h <- params[2]
    s1 <- params[3]
    s2 <- params[4]
    W <- fitness_weights(h, s1, s2)
    
    lpdf_reps <- map2_dbl(inits, data, \(irep, drep) {
        A1 <- irep$A[[1]] + d
        A1 <- A1 / sum(A1)
        
        res <- run_model(A1, 1, W, 1, num_gen)
        res_geno <- model_summary(res, mat2genotype) |> 
            pull(prop)
        
        sum(map_dbl(1:(num_gen-1), \(g) {
            slice <- 6*(g-1)+1:6
            dmultinom(drep$count[slice], prob=res_geno[slice], log=TRUE)
        }))
    })
    return(sum(lpdf_reps)) # assuming "uniform" priors
}

vm <- read_csv("data-proc/vera-maloof.csv") |> 
    pivot_longer(VFVF:ICIC, names_to="genotype", values_to="count") |> 
    filter(!str_detect(genotype, "IF"))

vm_init <- prop_init_adults(vm) |>
    group_by(site, rep) |>
    group_split()

vm_mod_data <- vm |> 
    group_by(site, rep, generation) |>
    mutate(N=sum(count)) |> # need to recompute since IF removed
    ungroup(generation) |> 
    group_split()

map_opt <- optim(
    c(0.01, 0.1, 0.9, 0.9),
    \(params) log_posterior(params, vm_init[1:6], vm_mod_data[1:6]),
    method="L-BFGS-B",
    lower=c(0.00001, 0, 0, 0), upper=c(1, 1, 1, 1), 
    hessian=TRUE, control=list(trace=3, fnscale=-1)
)

# add confidence intervals
mu <- map_opt$par
cov_mat <- solve(-map_opt$hessian)

post_samples <- mvrnorm(10000, mu, cov_mat)
colnames(post_samples) <- c("init. noise", "dominance coeff.", "VC selection coeff.", "IC selection coeff.")

post_summ <- post_samples |> 
    as_tibble() |> 
    pivot_longer(everything()) |> 
    mutate(name=fct_inorder(name)) |> 
    filter(value >= 0, value <= 1) |> 
    group_by(name) |> 
    summarise(ymin=quantile(value, 0.025), ymax=quantile(value, 0.975)) |> 
    mutate(mean=mu)

ggplot(post_summ, aes(name)) +
    geom_linerange(aes(ymin=ymin, ymax=ymax), col="#00798c", linewidth=1.03) +
    geom_point(aes(y=mean), col="#00798c", size=1.35) +
    scale_x_discrete(guide=guide_axis(angle=30)) +
    labs(x=NULL, y="posterior marginal") +
    theme_gray() +
    theme(text=element_text(color="gray30"), axis.ticks=element_line(color="gray92"))

ggsave("figs/post-two-sites.pdf", width=3.3, height=3.3)

## compare best results with data
colnames(post_samples) <- c("d", "h", "s1", "s2")

# TODO: why did you resample again after the stochastic run?
sample_prop_haplotype <- function(Alist, num_adults=50) {
    Alist_samp <- map(Alist, \(A) {
        As <- matrix(rmultinom(1, num_adults, c(A)), 3, 3)
        As / sum(As)
    })
    model_summary(Alist_samp, mat2haplotype)
}

summ_post_sims <- function(inits) {
    post_sims <- post_samples |> 
        as_tibble() |> 
        slice_sample(n=200) |> 
        rowwise() |> 
        mutate(
            A1=list(inits$A[[1]] + d),
            W=list(fitness_weights(h, s1, s2)),
            res=list(run_model(A1, 1, W, 1, 8, run_gen_stochastic)),
            res_hap=list(sample_prop_haplotype(res, num_adults=50)) # TODO: should be actual num_adults in data
        )
    
    post_sims |> 
        dplyr::select(-A1, -W) |> 
        unnest(res_hap) |> 
        group_by(generation, haplotype) |> 
        summarise(mean=mean(prop), ymin=quantile(prop, 0.025), ymax=quantile(prop, 0.975)) |> 
        ungroup()
}

post_summ <- map_dfr(vm_init[1:6], summ_post_sims)

ps_dat <- read_csv("data-proc/vera-maloof.csv") |> 
    filter(site %in% c("Ac", "Acp")) |>  
    pivot_longer(VF:IC, names_to="haplotype", values_to="count") |> 
    group_by(site, rep, generation) |> 
    reframe(haplotype=haplotype, prop=count/sum(count))

post_summ |> 
    mutate(site=ps_dat$site, rep=ps_dat$rep) |> 
    ggplot(aes(generation, mean, col=haplotype, group=haplotype)) +
    # geom_line(alpha=0.8, linewidth=0.95, linetype="14") +
    geom_ribbon(aes(ymin=ymin, ymax=ymax), fill="gray70", col=NA, alpha=0.4) +
    # geom_line(aes(y=prop), data=ps_dat, alpha=0.7, linetype="13", linewidth=1) +
    geom_point(aes(y=prop), data=ps_dat, size=1.23, alpha=0.9) +
    geom_line(linewidth=1, alpha=0.9) +
    facet_grid(site~rep) +
    scale_color_manual(values=c("#d1495b", "#edae49", "#66a182")) +
    labs(y="proportion") +
    theme_bw() +
    theme(
        text=element_text(color="gray30"), 
        axis.ticks=element_line(color="gray92"), 
        panel.grid.minor=element_blank(),
        strip.text=element_blank(),
        strip.background=element_blank(),
        panel.spacing.y=unit(0.8, "cm")
    )

ggsave("figs/post-sims-2site.pdf", width=6.5, height=3)
