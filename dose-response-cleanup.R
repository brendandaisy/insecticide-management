library(tidyverse)

genotypes <- fct_inorder(c("VFVF", "VFVC", "VCVC", "VFIC", "VCIC", "ICIC"))
dose0 <- read_csv("data-raw/Silvie data - dose-response insecticide mortality data.csv")

dose <- dose0 |> 
    mutate(genotype=factor(case_match(
        strain,
        "VVFF" ~ "VFVF",
        "VVFC" ~ "VFVC",
        "VIFC" ~ "VFIC",
        "VVCC" ~ "VCVC",
        "VICC" ~ "VCIC",
        "IICC" ~ "ICIC"
    ), levels=levels(genotypes))
    ) |> 
    dplyr::select(
        genotype, generation, exposure.date, mortality.date, 
        sex, concentration, dose.per.mosq, mortality
    )

write_csv(dose, "data-proc/dose-response.csv")

ggplot(dose, aes(dose.per.mosq+1, mortality, col=genotype)) +
    geom_point(shape=1, alpha=0.5) +
    geom_smooth(se=FALSE, method="glm", method.args=list(family="binomial"), linewidth=1.1) +
    scale_x_continuous(trans="log") +
    scale_color_manual(values=c("#00798c","#d1495b","#edae49","#66a182","#2e4057","#a274ab")) +
    coord_cartesian(ylim=c(0, 1)) +
    labs(x="dose") +
    theme_gray() +
    theme(
        text=element_text(color="gray30"), 
        axis.ticks=element_line(color="gray92"),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=rel(1.2)),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
    )

ggsave("figs/dose-response.pdf", width=3.9, height=3.1)
