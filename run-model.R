library(tidyverse)

# Set global values used throughout-----------------------------------------------
haplotypes <- c("VF", "VC", "IC")
nh <- length(haplotypes)

# The model-----------------------------------------------------------------------
run_model <- function(A1, R, W, M, num_gen, num_adults, gfun=run_gen) {
    C <- matrix(0.5, 3, 3)
    diag(C) <- 1
    accumulate(1:(num_gen-1), ~gfun(.x, R, W, M, C, num_adults), .init=A1)
}

# num_adults unused
run_gen <- function(A, R, W, M, C, num_adults) {
    Ef <- 0.5 * R * A
    Em <- 0.5 * A
    Gf <- colSums(C * W * Ef)
    Gm <- colSums(C * W * Em)
    Gf <- Gf / sum(Gf)
    Gm <- Gm / sum(Gm)
    L <- M * Gf %*% t(Gm)
    return(L / sum(L) * num_adults)
}

run_gen_stochastic <- function(A, R, W, M, C, num_adults) {
    Ef <- 0.5 * R * A # VM maintained exactly 50% for each sex
    Em <- 0.5 * A
    Gf <- colSums(C * W * Ef)
    Gm <- colSums(C * W * Em)
    Gf <- Gf / sum(Gf)
    Gm <- Gm / sum(Gm)
    L <- c(M * Gf %*% t(Gm))
    A <- matrix(rmultinom(1, num_adults, L), 3, 3)
    return(A)
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
