library(tidyverse)

geno2index <- list( # upper triangle of A
    "VFVF"=c(1, 1), "VFVC"=c(1, 2), "VCVC"=c(2, 2), 
    "VFIC"=c(1, 3), "VCIC"=c(2, 3), "ICIC"=c(3, 3))
)

index2geno <- expand_grid(i=1:nh, j=1:nh) |> 
    mutate(genotype=ifelse(i >= j, str_c(haplotypes[j], haplotypes[i]), str_c(haplotypes[i], haplotypes[j]))) |> 
    group_by(genotype)

init_adults <- function(dat) {
    A <- matrix(0, nh, nh)
    iwalk(geno2index, \(i, g) {
        genos <- filter(dat, genotype == g, t == 1)$prop[1]
        if (i[1] == i[2]) {
            A[i[1], i[2]] <<- genos
        } else {
            A[i[1], i[2]] <<- genos / 2
            A[i[2], i[1]] <<- genos / 2
        }
    })
    return(num_adults * A)
}

# a simple form of fitness depending only on a "selection coefficient" s and "dominance coef." h
fitness_weights <- function(s, h, type=1:2) {
    x <- if (type == 1) 1 else if (type == 2) h
    H <- matrix(c(
        0, h, h,
        h, x, x,
        h, x, 1
    ), nrow=3, byrow=TRUE)
    return(1 - s*H)
}

# Model summary functions----------------------------------------------------------
model_summary <- function(result, fun=prop_haplotype) {
    imap_dfr(result, ~mutate(fun(.x), t=.y))
}

prop_pairs <- function(A) {
    pairs <- expand_grid(h1=haplotypes, h2=haplotypes) |> mutate(g=str_c(h1, "_", h2)) |> pull(g)
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
    index2geno |> 
        summarize(prop=sum(map2_dbl(i, j, \(i, j) A[i, j])) / num_adults, .groups="drop")
    # expand_grid(i=1:nrow(A), j=1:nrow(A)) |> 
    #     filter(i >= j) |> 
    #     mutate(
    #         genotype=str_c(haplotypes[j], haplotypes[i]),
    #         prop=map2_dbl(i, j, \(i, j) if (i == j) A[i, j] else A[i, j] + A[j, i]) / num_adults
    #     )
}
