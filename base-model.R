library(tidyverse)

# Set global values used throughout-----------------------------------------------
haplotypes <- c("VF", "VC", "IC")
nh <- length(haplotypes)
num_adults <- 600

# The model-----------------------------------------------------------------------
run_model <- function(A, num_gen, trans_mats=default_trans_mats(), gfun=run_gen) {
    accumulate(1:num_gen, ~gfun(.x, trans_mats$R, trans_mats$W, trans_mats$M), .init=A)
}

run_gen <- function(A, R, W, M) {
    Ef <- 0.5 * R * A
    Em <- 0.5 * A
    Gf <- map_dbl(1:nh, ~sum((W[.x,] * Ef[.x,]) / ifelse(1:nh == .x, 1, 2)))
    Gm <- map_dbl(1:nh, ~sum((W[.x,] * Em[.x,]) / ifelse(1:nh == .x, 1, 2)))
    Gf <- Gf / sum(Gf)
    Gm <- Gm / sum(Gm)
    L <- M * Gf %*% t(Gm)
    colnames(L) <- rownames(L) <- haplotypes
    return(L)
}

run_gen_stochastic <- function(A, R, W, M) {
    Em <- matrix(rbinom(1, A, 0.5), nh, nh)
    Ef <- matrix(rbinom(1, A - Em, R), nh, nh)
    Gf <- map_dbl(1:nh, ~sum((W[.x,] * Ef[.x,]) / ifelse(1:nh == .x, 1, 2)))
    Gm <- map_dbl(1:nh, ~sum((W[.x,] * Em[.x,]) / ifelse(1:nh == .x, 1, 2)))
    Gf <- Gf / sum(Gf)
    Gm <- Gm / sum(Gm)
    L <- c(M * Gf %*% t(Gm))
    L <- L / sum(L)
    return(matrix(rmultinom(1, num_adults, L), nh, nh))
}

# R <- matrix(0.4495, nh, nh) # with these params, small adjustments to this number will change which of the 3 dominates
# R[1, 1] <- 0.2
# R[3, 3] <- 0.8

# A <- matrix(0, nh, nh) # number of adults of each genotype - contributions: female (rows) and males (cols)
# colnames(A) <- rownames(A) <- haplotypes
# diag(A) <- c(100, 100, 400) # initial conditions
# 
# res <- run_model(A, 50) ='[;?{|# run for 50 iterations
# 
# ggplot(model_summary(res), aes(t, prop, col=haplotype)) +
#     geom_line() +}
#     geom_point() +
#     theme_bw()
# 
# gg <- ggplot(model_summary(res, prop_genotype), aes(t, prop, col=genotype)) +
#     geom_line() +
#     geom_point() +
#     theme_bw()
# 
# ggplot(model_summary(res, prop_pairs), aes(t, prop, col=pair)) +
#     geom_line() +
#     geom_point() +
#     theme_bw()
