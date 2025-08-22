library(tidyverse)

col_names <- c("site", "generation", "rep", "N", "VFVF", "VFVC", "VCVC", "VFIF", "VFIC", "VCIC", "IFIF", "IFIC", "ICIC")

vm_raw <- readxl::read_xlsx("data-raw/Vera-Maloof tableS2.xlsx", col_names=col_names) |> 
    slice(3:n())

# remove contaminate rows containing group totals, seperator text, etc.
bad_rows <- filter(vm_raw, str_detect(site, "Contin|Site") | `VFVF` == "VVFF" | rep == "Total")

vm_proc <- vm_raw |> 
    anti_join(bad_rows) |> 
    fill(site, generation) |> 
    filter(!is.na(N)) |> 
    mutate(generation=as.double(str_extract(generation, "\\d+")), across(rep:`ICIC`, as.double)) |> 
    # the NUMBER of haplotypes floating around, i.e. 2 * num_adults
    mutate(
        VF=2*`VFVF`+`VFVC`+`VFIF`+`VFIC`,
        VC=`VFVC`+2*`VCVC`+`VCIC`,
        IC=`VFIC`+`VCIC`+`IFIC`+2*`ICIC`,
        IF=`VFIF`+2*`IFIF`+`IFIC`
    )

write_csv(vm_proc, "data-proc/vera-maloof.csv")

vm_proc |> 
    pivot_longer(`VFVF`:`ICIC`, names_to="genotype", values_to="count") |> 
    ggplot(aes(generation, count, col=genotype, group=genotype)) +
    geom_line() +
    geom_point() +
    facet_grid(site ~ rep, labeller=labeller(site=label_value, rep=label_both)) +
    theme_bw()

ggsave("figs/vera-maloof-geno.pdf", width=8, height=8)
