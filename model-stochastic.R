library(tidyverse)

run_gen_stochastic <- function(A) {
    nh <- nrow(A)
    Em <- Ef <- matrix(0, nh, nh) 
    for (i in 1:nh) {
        for (j in 1:nh) {
            Em[i, j] <- rbinom(1, A[i, j], 0.5)
            Ef[i, j] <- rbinom(1, A[i, j] - Em[i, j], R[i, j])
        }
    }
    Gf <- map_dbl(1:nh, ~sum((W[.x,] * Ef[.x,]) / ifelse(1:nh == .x, 1, 2)))
    Gm <- map_dbl(1:nh, ~sum((W[.x,] * Em[.x,]) / ifelse(1:nh == .x, 1, 2)))
    Gf <- Gf / sum(Gf)
    Gm <- Gm / sum(Gm)
    Avec <- rmultinom(1, num_adults, c(M * Gf %*% t(Gm)))
    return(matrix(Avec, nh, nh))
}

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
res_sto <- run_model(A, 50, fun=run_gen_stochastic)

res <- bind_rows(
    mutate(model_summary(res_det), type="deterministic"),
    mutate(model_summary(res_sto), type="stochastic")
)

ggplot(res, aes(t, prop, col=haplotype, linetype=type)) +
    geom_line() +
    theme_bw()

ggplot(model_summary(res, prop_genotype), aes(t, prop, col=genotype)) +
    geom_line() +
    geom_point() +
    theme_bw()

ggplot(model_summary(res, prop_pairs), aes(t, prop, col=pair)) +
    geom_line() +
    geom_point() +
    theme_bw()