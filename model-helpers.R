library(tidyverse)

#  some globals-------------------------------------------------------------------
haplotypes <- fct_inorder(c("VF", "VC", "IC"))
genotypes <- fct_inorder(c("VFVF", "VFVC", "VCVC", "VFIC", "VCIC", "ICIC"))
nh <- length(haplotypes)

geno2index <- list( # upper triangle of A
    "VFVF"=c(1, 1), "VFVC"=c(1, 2), "VCVC"=c(2, 2), 
    "VFIC"=c(1, 3), "VCIC"=c(2, 3), "ICIC"=c(3, 3)
)

index2geno <- expand_grid(i=1:nh, j=1:nh) |> 
    mutate(
        genotype=factor(
            ifelse(i >= j, str_c(haplotypes[j], haplotypes[i]), str_c(haplotypes[i], haplotypes[j])),
            levels=levels(genotypes)
        )
    ) |> 
    arrange(genotype) |> 
    group_by(genotype)

#  functions for processing Vera-Maloof data--------------------------------------
prop_adults_rep <- function(rep_df, key) {
    A <- matrix(0, nh, nh)
    imap(geno2index, \(i, g) {
        genos <- filter(rep_df, genotype == g)$count[1] # get the counts for this genotype
        if (i[1] == i[2]) {
            A[i[1], i[2]] <<- genos
        } else {
            A[i[1], i[2]] <<- genos / 2
            A[i[2], i[1]] <<- genos / 2
        }
    })
    n <- sum(A)
    return(tibble_row(A=list(A / n), N=n))
}

prop_init_adults <- function(vm_long) {
    vm_long |> 
        filter(generation == 1) |> 
        group_by(site, rep, generation) |> 
        group_modify(prop_adults_rep) |> 
        ungroup()
}

# Model helper functions----------------------------------------------------------

# a simple form of fitness depending only on a "selection coefficient" s and "dominance coef." h
fitness_weights <- function(h, s1, s2) {
    1 - matrix(c(
        0, h*s1, h*s2,
        h*s1, s1, h*s2,
        h*s2, h*s2, s2
    ), nrow=3, byrow=TRUE)
}

# Model summary functions----------------------------------------------------------
model_summary <- function(result, fun) {
    imap_dfr(result, ~mutate(fun(.x), generation=.y))
}

# prop_pairs <- function(A) {
#     pairs <- expand_grid(h1=haplotypes, h2=haplotypes) |> mutate(g=str_c(h1, "_", h2)) |> pull(g)
#     props <- c(t(A)) / num_adults
#     names(props) <- pairs
#     enframe(props, "pair", "prop")
# }

mat2haplotype <- function(A) {
    props <- map_dbl(1:nrow(A), ~sum(A[.x,] + A[,.x]) / 2)
    names(props) <- haplotypes
    enframe(props, "haplotype", "value")
}

mat2genotype <- function(A) {
    index2geno |> 
        summarize(value=sum(map2_dbl(i, j, \(i, j) A[i, j])), .groups="drop")
}
