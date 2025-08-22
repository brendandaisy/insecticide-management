library(tidyverse)

vm <- read_csv("data-proc/vera-maloof.csv")

num_reps <- distinct(vm, site, rep) |> nrow()

#  How prevalent was the IF haplotype, at a minimum?------------------------------

# The VI-FF, II-FF, and II-FC genotypes must contain IF, while VI-FC may or may not
# so we get the percentage of genotypes that definitely contained IF, as a minimum:
vm_if <- vm |> 
    group_by(site, generation, rep) |> 
    mutate(prop_IF = IF / (2*N), num_w_IF=`VF/IF`+`IF/IF`+`IF/IC`) |> 
    ungroup()

vm_if |> 
    ggplot(aes(generation, prop_IF, group=rep)) +
    geom_line() +
    geom_point() +
    facet_wrap(~site, nrow=4) +
    theme_bw() +
    labs(y="proportion IF")

ggsave("figs/vm-prop-IF.pdf", width=4.5, height=5.2)

# How many replicates had more than one IF in the first generation?
nrow(filter(vm_if, generation == 1, num_w_IF > 0))

# Percentage of haplotypes that were IF, including/excluding VI-FC genotype:
sum(vm_if$IF) / sum(2*vm_if$N)
sum(vm_if$IF + vm_if$`VF/IC`) / sum(2*vm_if$N)

# vm_proc |> 
#     pivot_longer(VF:IF, names_to="haplotype", values_to="count") |> 
#     ggplot(aes(generation, count, col=haplotype, group=haplotype)) +
#     geom_line() +
#     geom_point() +
#     facet_grid(site ~ rep, labeller=labeller(site=label_value, rep=label_both)) +
#     theme_bw()
# 
# ggsave("figs/vera-maloof-haplo.pdf", width=8, height=8)

vm |> 
    pivot_longer(VF:IF, names_to="haplotype", values_to="count") |> 
    group_by(generation, haplotype) |> 
    summarise(count=sum(count), .groups="drop") |> 
    ggplot(aes(generation, count, col=haplotype, group=haplotype)) +
    geom_line() +
    geom_point() +
    theme_bw()

ggsave("figs/vera-maloof-totals.pdf", width=4, height=3.6)

vmi <- vm |> 
    mutate(I=(IC+IF)/(2*N), .after=site)

lm(I~generation, vmi)

# TODO: Silvie 5/21:
# appear to be two groups of slopes for the I frequency. Maybe these two groups 
# correspond to presence absence of the resistant locus at 410 (i.e. the steeper
# slope has the R allele)

ggplot(vmi, aes(generation, I, group=site, col=site)) +
    geom_point(alpha=0.3) +
    geom_smooth(method="lm", se=FALSE)
