library(tidyverse)
library(ggtern)

source("run-model.R")
source("model-helpers.R")

get_dom_haplo <- function(h, s, r, A, num_gen) {
    W <- fitness_weights(h, s, s)
    R <- matrix(r, 3, 3)
    R[3, 3] <- 1
    res <- run_model(A, R, W, 1, num_gen, num_adults=1, gfun=run_gen)
    if (any(is.nan(last(res))))
        print(res[1:10])
    hap_fin <- mat2haplotype(last(res))
    haplotypes[which.max(hap_fin$value)]
}

A <- matrix(0, 3, 3)
colnames(A) <- rownames(A) <- haplotypes
diag(A) <- c(1/3, 1/3, 1/3) # assume there are more resistant mosquitoes, and only homozygotes initially

run_grid <- expand_grid(
    h=seq(0, 1, length.out=15), 
    s=seq(0.001, 1, length.out=15), 
    r=seq(0.001, 1, length.out=15) # can't really be 0 because then everything dies
) |> 
    rowwise() |> 
    mutate(dom_haplo=get_dom_haplo(h, s, r, A, num_gen=300))

# tenary plot!
ggtern(run_grid, aes(h, s, r)) +
    geom_point(aes(col=dom_haplo), shape=1) +
    geom_Risoprop(value=0.6, alpha=0.5) +
    theme_rgbg()

ggsave("figs/sens-ternary.pdf", width=10, height=9)

# also can try to look at things with 2d plots
ggplot(run_grid, aes(h, s, fill=dom_haplo)) +
    geom_tile(alpha=0.3) + # alpha < 1 allows seeing different outcomes for the missing parameter (i.e. s2 in this example)
    theme_bw()

ggplot(run_grid, aes(s, r, fill=dom_haplo)) +
    geom_tile(alpha=0.3) +
    theme_bw()

ggplot(run_grid, aes(r, h, fill=dom_haplo)) +
    geom_tile(alpha=0.3) +
    theme_bw()
