library(tidyverse)

run_model <- function(A, num_gen) {
    accumulate(1:num_gen, ~run_generation(.x), .init=A)
}

run_generation <- function(A) {
    nh <- nrow(A)
    Ef <- 0.5 * R * A
    Em <- 0.5 * A
    Gf <- map_dbl(1:nh, ~sum((W[.x,] * Ef[.x,]) / ifelse(1:nh == .x, 1, 2)))
    Gm <- map_dbl(1:nh, ~sum((W[.x,] * Em[.x,]) / ifelse(1:nh == .x, 1, 2)))
    Gf <- Gf / sum(Gf)
    Gm <- Gm / sum(Gm)
    L <- M * Gf %*% t(Gm)
    L <- L / sum(L)
    colnames(L) <- rownames(L) <- haplotypes
    return(num_adults * L)
} 

model_summary <- function(result, fun=prop_haplotype) {
    imap_dfr(result, ~mutate(fun(.x), t=.y))
}

prop_pairs <- function(A) {
    props <- c(t(A)) / num_adults
    names(props) <- pairs
    enframe(props, "pair", "prop")
}

prop_haplotype <- function(A) {
    props <- map_dbl(1:nrow(A), ~sum(A[.x,] + A[,.x]) / (2*num_adults))
    names(props) <- haplotypes
    enframe(props, "haplotype", "prop")
}

prop_genotype <- function(A) {
    expand_grid(i=1:nrow(A), j=1:nrow(A)) |> 
        filter(i >= j) |> 
        mutate(
            genotype=str_c(haplotypes[j], haplotypes[i]),
            prop=map2_dbl(i, j, \(i, j) if (i == j) A[i, j] else A[i, j] + A[j, i]) / num_adults
        )
}

haplotypes <- c("VF", "VC", "IC")
pairs <- expand_grid(h1=haplotypes, h2=haplotypes) |> mutate(g=str_c(h1, "_", h2)) |> pull(g)
nh <- length(haplotypes)

num_adults <- 600
s <- 0.4
h <- 0.02
H <- matrix(c(
    0, 0, h,
    0, h, h,
    h, h, 1
), nrow=3, byrow=TRUE)
W <- 1 - s*H
R <- matrix(0.4495, nh, nh) # with these params, small adustments to this number will change which of the 3 dominates
R[1, 1] <- 0.2
R[3, 3] <- 0.8
M <- matrix(1, nh, nh)

A <- matrix(0, nh, nh) # number of adults of each genotype - contributions: female (rows) and males (cols)
colnames(A) <- rownames(A) <- haplotypes
diag(A) <- c(580, 10, 10)

res <- run_model(A, 50) # run for 50 iterations

ggplot(model_summary(res), aes(t, prop, col=haplotype)) +
    geom_line() +
    geom_point() +
    theme_bw()

ggplot(model_summary(res, prop_genotype), aes(t, prop, col=genotype)) +
    geom_line() +
    geom_point() +
    theme_bw()

ggplot(model_summary(res, prop_pairs), aes(t, prop, col=pair)) +
    geom_line() +
    geom_point() +
    theme_bw()
