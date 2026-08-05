library(ggplot2)
library(reshape2)
library(patchwork)
library(ggsignif)
library(patchwork)
library(scales)
library("ggbeeswarm")


p1 <- read.table("merge_num.txt",header=T)

p_pop <- read.table("pop.txt",header=T)

data <- merge.data.frame(p1,p_pop,by="ID")


subset_data<-subset(data,lowdepth=="no")

names<-c(
  'MIX_CAP','SC_CAP','SC_WLP','QLI','SC_WSP'
)

subset_data$POP <- factor(subset_data$POP,
                          levels = names,
                          ordered = T)

SC_WLP_subset<-subset(subset_data,POP=='SC_WLP')

my_theme <- function() {
  theme_classic() +
    theme(axis.title.x = element_text(size = 15, colour = "black"),
          axis.title.y = element_text(size = 15, colour = "black"),
          axis.text.x = element_text( size = 15,colour = "black",angle = 45,hjust=1),
          axis.text.y = element_text( size = 15, colour = "black"),
          axis.ticks.length = unit(0.2, "cm"),  
          axis.ticks = element_line(size = 1),  
          axis.line = element_line(size = 1),  
          legend.position="none",
          text = element_text(family = "Times"),
          strip.background = element_rect(fill = "grey90", color = "black", size = 1),  
          strip.text = element_text(color = "black",size =15)
    )
}

top_theme <- function() {
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank())
}

ratio_theme <- function() {
  theme_classic()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5), 
          axis.title.x = element_text(size = 15, colour = "black"),
          axis.title.y = element_text(size = 15, colour = "black"),
          axis.text.x = element_text(size = 15,colour = "black",angle = 45,hjust=1),
          axis.text.y = element_text(size = 15, colour = "black"),
          axis.ticks.length = unit(0.2, "cm"),  
          axis.ticks = element_line(size = 1),  
          axis.line = element_line(size = 1),  
          legend.position="none",
          legend.text = element_text(size=15),
          text = element_text(family = "Times"))
}

# Color settings
fill_cols <- c(
  SC_CAP = "#F2CC8E",
  QLI = "#82B29A",
  SC_WLP = "#DF7A5E",
  MIX_CAP = "#D2BFA5",
  SC_WSP = "#8e7cc3"
)
# Function for unified del_num style
make_num_plot <- function(data, y_col, comparisons, y_position, y_limits, y_breaks) {
  ggplot(data, aes(x = POP)) +
    geom_quasirandom(
      aes(y = .data[[y_col]], color = POP),
      width = 0.4,
      shape = 16,
      size = 2,
      alpha = 0.8
    ) +
    geom_boxplot(
      aes(y = .data[[y_col]]),
      fill = NA,
      color = "black",
      width = 0.6,
      outlier.shape = NA
    ) +
    geom_signif(
      aes(y = .data[[y_col]], group = POP),
      comparisons = comparisons,
      textsize = 3,
      y_position = y_position,
      map_signif_level = FALSE
    ) +
    scale_y_continuous(
      limits = y_limits,
      labels = y_breaks,
      breaks = y_breaks,
      expand = c(0, 0)
    ) +
    scale_color_manual(values = fill_cols)
}

# Heterozygous and Homozygous
LOF_het <- make_num_plot(
  data = subset_data,
  y_col = "snpLOFNum_het",
  comparisons = list(c("MIX_CAP", "SC_CAP")),
  y_position = c(185),
  y_limits = c(60, 200),
  y_breaks = seq(60, 200, 70)
)

LOF_hom <- make_num_plot(
  data = subset_data,
  y_col = "snpLOFNum_hom",
  comparisons = list(c("SC_WLP", "QLI")),
  y_position = c(120),
  y_limits = c(60, 140),
  y_breaks = seq(60, 140, 40)
)


del_het <- make_num_plot(
  data = subset_data,
  y_col = "snpMisDelNum_het",
  comparisons = list(
    c("MIX_CAP", "SC_CAP"),
    c("SC_WLP", "QLI")
  ),
  y_position = c(1875, 1800),
  y_limits = c(700, 2000),
  y_breaks = seq(700, 2000, 650)
)

del_hom <- make_num_plot(
  data = subset_data,
  y_col = "snpMisDelNum_hom",
  comparisons = list(
    c("MIX_CAP", "SC_CAP"),
    c("SC_CAP", "SC_WLP"),
    c("SC_WLP", "QLI")
  ),
  y_position = c(1100, 1150, 1200),
  y_limits = c(700, 1300),
  y_breaks = seq(700, 1400, 300)
)


Tol_het <- make_num_plot(
  data = subset_data,
  y_col = "snpMisTolNum_het",
  comparisons = list(
    c("MIX_CAP", "SC_CAP"),
    c("SC_CAP", "SC_WLP"),
    c("SC_WLP", "QLI")
  ),
  y_position = c(6100, 5800, 5500),
  y_limits = c(3000, 6500),
  y_breaks = seq(3000, 6500, 1750)
)

Tol_hom <- make_num_plot(
  data = subset_data,
  y_col = "snpMisTolNum_hom",
  comparisons = list(
    c("MIX_CAP", "SC_CAP"),
    c("SC_CAP", "SC_WLP"),
    c("SC_WLP", "QLI")
  ),
  y_position = c(4300, 4500, 4700),
  y_limits = c(3000, 5000),
  y_breaks = seq(3000, 5000, 1000)
)

Syn_het <- make_num_plot(
  data = subset_data,
  y_col = "snpSynNum_het",
  comparisons = list(
    c("MIX_CAP", "SC_CAP"),
    c("SC_CAP", "SC_WLP"),
    c("SC_WLP", "QLI")
  ),
  y_position = c(12200, 11500, 10800),
  y_limits = c(6000, 13000),
  y_breaks = seq(6000, 13000, 3500)
)

Syn_hom <- make_num_plot(
  data = subset_data,
  y_col = "snpSynNum_hom",
  comparisons = list(
    c("MIX_CAP", "SC_CAP"),
    c("SC_CAP", "SC_WLP"),
    c("SC_WLP", "QLI")
  ),
  y_position = c(9500, 9750, 10000),
  y_limits = c(7000, 10500),
  y_breaks = seq(7000, 10500, 1750)
)


CNE_het <- make_num_plot(
  data = subset_data,
  y_col = "snpCNENum_het",
  comparisons = list(
    c("MIX_CAP", "SC_CAP"),
    c("SC_CAP", "SC_WLP"),
    c("SC_WLP", "QLI")
  ),
  y_position = c(34000, 32000, 30000),
  y_limits = c(13000, 36000),
  y_breaks = seq(13000, 36000, 11500)
)

CNE_hom <- make_num_plot(
  data = subset_data,
  y_col = "snpCNENum_hom",
  comparisons = list(
    c("MIX_CAP", "SC_CAP"),
    c("SC_CAP", "SC_WLP"),
    c("SC_WLP", "QLI")
  ),
  y_position = c(21000, 22000, 23000),
  y_limits = c(13000, 25000),
  y_breaks = seq(13000, 25000, 6000)
)


## =========================
## Add facet labels
## =========================

## =========================
## Theme helpers for facet strips
## =========================
remove_col_strip <- function() {
  theme(
    strip.text.x = element_blank(),
    strip.background.x = element_blank()
  )
}
remove_row_strip <- function() {
  theme(
    strip.text.y = element_blank(),
    strip.background.y = element_blank()
  )
}
remove_all_strips <- function() {
  theme(
    strip.text.x = element_blank(),
    strip.background.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background.y = element_blank()
  )
}

## =========================
## Homozygous row
## Keep column strips only in this row
## Keep row strip only in the last column
## =========================
Syn_hom <- Syn_hom +
  facet_grid("Homozygous" ~ "Synonymous") +
  labs(x = '', y = '') +
  my_theme() +
  top_theme() +
  remove_row_strip()
Tol_hom <- Tol_hom +
  facet_grid("Homozygous" ~ "Tolerated Mis.") +
  labs(x = '', y = '') +
  my_theme() +
  top_theme() +
  remove_row_strip()
del_hom <- del_hom +
  facet_grid("Homozygous" ~ "Deleterious Mis.") +
  labs(x = '', y = '') +
  my_theme() +
  top_theme() +
  remove_row_strip()
CNE_hom <- CNE_hom +
  facet_grid("Homozygous" ~ "CNE") +
  labs(x = '', y = '') +
  my_theme() +
  top_theme() +
  remove_row_strip()
LOF_hom <- LOF_hom +
  facet_grid("Homozygous" ~ "LOF") +
  labs(x = '', y = '') +
  my_theme() +
  top_theme() +
  remove_row_strip()

## =========================
## Heterozygous row
## Remove column strips in this row
## Keep row strip only in the last column
## =========================
Syn_het <- Syn_het +
  facet_grid("Heterozygous" ~ "Synonymous") +
  labs(x = '', y = '') +
  my_theme() +
  top_theme() +
  remove_all_strips()
Tol_het <- Tol_het +
  facet_grid("Heterozygous" ~ "Tolerated Mis.") +
  labs(x = '', y = '') +
  my_theme() +
  top_theme() +
  remove_all_strips()
del_het <- del_het +
  facet_grid("Heterozygous" ~ "Deleterious Mis.") +
  labs(x = '', y = '') +
  my_theme() +
  top_theme() +
  remove_all_strips()
CNE_het <- CNE_het +
  facet_grid("Heterozygous" ~ "CNE") +
  labs(x = '', y = '') +
  my_theme() +
  top_theme() +
  remove_all_strips()
LOF_het <- LOF_het +
  facet_grid("Heterozygous" ~ "LOF") +
  labs(x = '', y = '') +
  my_theme() +
  top_theme() +
  remove_all_strips()



hom_num <- Syn_hom | Tol_hom | del_hom | CNE_hom | LOF_hom
het_num <- Syn_het | Tol_het | del_het | CNE_het | LOF_het
del_num <- hom_num / het_num
del_num <- del_num & theme(
  plot.margin = margin(2, 2, 2, 2)
)
del_num

ggsave(del_num, file='compare_variant.pdf', width=8, height=20)

# Function for unified plot style using expression
make_expr_plot <- function(data, y_expr, comparisons = NULL, y_position = NULL,
                           y_limits, y_breaks) {
  y_expr <- rlang::enquo(y_expr)
  
  p <- ggplot(data, aes(x = POP)) +
    geom_quasirandom(
      aes(y = !!y_expr, color = POP),
      width = 0.4,
      shape = 16,
      size = 2,
      alpha = 0.8
    ) +
    geom_boxplot(
      aes(y = !!y_expr),
      fill = NA,
      color = "black",
      width = 0.6,
      outlier.shape = NA
    ) +
    scale_y_continuous(
      limits = y_limits,
      labels = y_breaks,
      breaks = y_breaks,
      expand = c(0, 0)
    ) +
    scale_color_manual(values = fill_cols)
  
  if (!is.null(comparisons)) {
    p <- p +
      geom_signif(
        aes(y = !!y_expr, group = POP),
        comparisons = comparisons,
        textsize = 3,
        y_position = y_position,
        map_signif_level = FALSE
      )
  }
  
  return(p)
}

make_ratio_plot <- function(data, y_expr, comparisons = NULL, y_position = NULL,
                            ref_mean_data = syn_ref_mean_data,
                            y_limits = c(0.4, 0.8),
                            y_breaks = seq(0.4, 0.8, 0.1)) {
  y_expr <- rlang::enquo(y_expr)
  
  p <- ggplot(data, aes(x = POP)) +
    geom_quasirandom(
      aes(y = !!y_expr, color = POP),
      width = 0.4,
      shape = 16,
      size = 2,
      alpha = 0.8
    ) +
    geom_boxplot(
      aes(y = !!y_expr),
      fill = NA,
      color = "black",
      width = 0.6,
      outlier.shape = NA
    ) +
    scale_y_continuous(
      limits = y_limits,
      labels = y_breaks,
      breaks = y_breaks,
      expand = c(0, 0)
    ) +
    scale_color_manual(values = fill_cols)
  
  if (!is.null(comparisons)) {
    p <- p +
      geom_signif(
        aes(y = !!y_expr, group = POP),
        comparisons = comparisons,
        textsize = 3,
        y_position = y_position,
        map_signif_level = FALSE
      )
  }
  
  return(p)
}

ratio_comparisons <- list(
  c("MIX_CAP", "SC_CAP"),
  c("SC_CAP", "SC_WLP"),
  c("SC_WLP", "QLI")
)

ratio_y_position <- c(0.675, 0.7, 0.725)

Syn_ratio <- make_ratio_plot(
  data = subset_data,
  y_expr = snpSynNum_hom * 2 / (snpSynNum_hom * 2 + snpSynNum_het),
  comparisons = ratio_comparisons,
  y_position = ratio_y_position
)

LOF_ratio <- make_ratio_plot(
  data = subset_data,
  y_expr = snpLOFNum_hom * 2 / (snpLOFNum_hom * 2 + snpLOFNum_het),
  comparisons = ratio_comparisons,
  y_position = ratio_y_position
)

del_ratio <- make_ratio_plot(
  data = subset_data,
  y_expr = snpMisDelNum_hom * 2 / (snpMisDelNum_hom * 2 + snpMisDelNum_het),
  comparisons = ratio_comparisons,
  y_position = ratio_y_position
)

Tol_ratio <- make_ratio_plot(
  data = subset_data,
  y_expr = snpMisTolNum_hom * 2 / (snpMisTolNum_hom * 2 + snpMisTolNum_het),
  comparisons = ratio_comparisons,
  y_position = ratio_y_position
)

Mis_ratio <- make_ratio_plot(
  data = subset_data,
  y_expr = snpMisNum_hom * 2 / (snpMisNum_hom * 2 + snpMisNum_het),
  comparisons = ratio_comparisons,
  y_position = ratio_y_position
)

CNE_ratio <- make_ratio_plot(
  data = subset_data,
  y_expr = snpCNENum_hom * 2 / (snpCNENum_hom * 2 + snpCNENum_het),
  comparisons = ratio_comparisons,
  y_position = ratio_y_position
)

LOF_all <- make_expr_plot(
  data = subset_data,
  y_expr = snpLOFNum_het + snpLOFNum_hom,
  comparisons = list(c("MIX_CAP", "SC_CAP")),
  y_position = c(275),
  y_limits = c(150, 300),
  y_breaks = seq(150, 300, 75)
)
del_all <- make_expr_plot(
  data = subset_data,
  y_expr = snpMisDelNum_het + snpMisDelNum_hom,
  comparisons = list(
    c("MIX_CAP", "SC_CAP"),
    c("SC_WLP", "QLI")
  ),
  y_position = c(2900, 2800),
  y_limits = c(2000, 3000),
  y_breaks = seq(2000, 3000, 1000)
)
Tol_all <- make_expr_plot(
  data = subset_data,
  y_expr = snpMisTolNum_het + snpMisTolNum_hom,
  comparisons = list(
    c("MIX_CAP", "SC_CAP"),
    c("SC_WLP", "QLI")
  ),
  y_position = c(9600, 9500),
  y_limits = c(8000, 10000),
  y_breaks = seq(8000, 10000, 1000)
)
Mis_all <- make_expr_plot(
  data = subset_data,
  y_expr = snpMisNum_het + snpMisNum_hom,
  comparisons = list(
    c("MIX_CAP", "SC_CAP"),
    c("SC_WLP", "QLI")
  ),
  y_position = c(14000, 13700),
  y_limits = c(11000, 15000),
  y_breaks = seq(11000, 15000, 1000)
)
Syn_all <- make_expr_plot(
  data = subset_data,
  y_expr = snpSynNum_het + snpSynNum_hom,
  comparisons = list(
    c("MIX_CAP", "SC_CAP"),
    c("SC_WLP", "QLI")
  ),
  y_position = c(20000, 19500),
  y_limits = c(16000, 21000),
  y_breaks = seq(16000, 21000, 2500)
)
CNE_all <- make_expr_plot(
  data = subset_data,
  y_expr = snpCNENum_hom + snpCNENum_het,
  comparisons = list(
    c("MIX_CAP", "SC_CAP"),
    c("SC_CAP", "SC_WLP"),
    c("SC_WLP", "QLI")
  ),
  y_position = c(50000, 49000, 48000),
  y_limits = c(40000, 52000),
  y_breaks = seq(40000, 52000, 6000)
)

LOF_all_derived <- make_expr_plot(
  data = subset_data,
  y_expr = snpLOFNum_het + snpLOFNum_hom * 2,
  comparisons = list(
    c("MIX_CAP", "SC_CAP")
  ),
  y_position = c(380),
  y_limits = c(250, 400),
  y_breaks = seq(250, 400, 75)
)

del_all_derived <- make_expr_plot(
  data = subset_data,
  y_expr = snpMisDelNum_het + snpMisDelNum_hom * 2,
  comparisons = NULL,
  y_position = NULL,
  y_limits = c(3000, 3750),
  y_breaks = seq(3000, 3750, 375)
)

Tol_all_derived <- make_expr_plot(
  data = subset_data,
  y_expr = snpMisTolNum_het + snpMisTolNum_hom * 2,
  comparisons = list(
    c("SC_CAP", "SC_WLP")
  ),
  y_position = c(13500),
  y_limits = c(11000, 14000),
  y_breaks = seq(11000, 14000, 1500)
)

Mis_all_derived <- make_expr_plot(
  data = subset_data,
  y_expr = snpMisNum_het + snpMisNum_hom * 2,
  comparisons = NULL,
  y_position = NULL,
  y_limits = c(16000, 20000),
  y_breaks = seq(16000, 20000, 1500)
)

Syn_all_derived <- make_expr_plot(
  data = subset_data,
  y_expr = snpSynNum_het + snpSynNum_hom * 2,
  comparisons = list(
    c("MIX_CAP", "SC_CAP")
  ),
  y_position = c(27500),
  y_limits = c(24000, 28000),
  y_breaks = seq(24000, 28000, 2000)
)

CNE_all_derived <- make_expr_plot(
  data = subset_data,
  y_expr = snpCNENum_hom * 2 + snpCNENum_het,
  comparisons = list(
    c("SC_WLP", "QLI")
  ),
  y_position = c(68000),
  y_limits = c(55000, 70000),
  y_breaks = seq(55000, 70000, 7500)
)


Syn_ratio <- Syn_ratio +
  facet_grid(. ~ "Synonymous") +
  labs(x = "", y = "2*hom/(2*hom+hete)") +
  my_theme() +
  remove_col_strip()
Tol_ratio <- Tol_ratio +
  facet_grid(. ~ "Tolerated Missense") +
  labs(x = "", y = "") +
  my_theme() +
  remove_col_strip()
Mis_ratio <- Mis_ratio +
  facet_grid(. ~ "Missense") +
  labs(x = "", y = "") +
  my_theme() +
  remove_col_strip()
del_ratio <- del_ratio +
  facet_grid(. ~ "Deleterious Missense") +
  labs(x = "", y = "") +
  my_theme() +
  remove_col_strip()
CNE_ratio <- CNE_ratio +
  facet_grid(. ~ "CNE") +
  labs(x = "", y = "") +
  my_theme() +
  remove_col_strip()
LOF_ratio <- LOF_ratio +
  facet_grid(. ~ "LOF") +
  labs(x = "", y = "") +
  my_theme() +
  remove_col_strip()

Syn_all_derived <- Syn_all_derived +
  facet_grid(. ~ "Synonymous") +
  labs(x = "", y = "Derived alleles") +
  my_theme() +
  top_theme()
Tol_all_derived <- Tol_all_derived +
  facet_grid(. ~ "Tolerated Missense") +
  labs(x = "", y = "") +
  my_theme() +
  top_theme()
Mis_all_derived <- Mis_all_derived +
  facet_grid(. ~ "Missense") +
  labs(x = "", y = "") +
  my_theme() +
  top_theme()
del_all_derived <- del_all_derived +
  facet_grid(. ~ "Deleterious Missense") +
  labs(x = "", y = "") +
  my_theme() +
  top_theme()
CNE_all_derived <- CNE_all_derived +
  facet_grid(. ~ "CNE") +
  labs(x = "", y = "") +
  my_theme() +
  top_theme()
LOF_all_derived <- LOF_all_derived +
  facet_grid(. ~ "LOF") +
  labs(x = "", y = "") +
  my_theme() +
  top_theme()


Syn_all <- Syn_all +
  labs(x = "", y = "Derived genotypes") +
  my_theme()
Tol_all <- Tol_all +
  labs(x = "", y = "") +
  my_theme()
Mis_all <- Mis_all +
  labs(x = "", y = "") +
  my_theme()
del_all <- del_all +
  labs(x = "", y = "") +
  my_theme()
CNE_all <- CNE_all +
  labs(x = "", y = "") +
  my_theme()
LOF_all <- LOF_all +
  labs(x = "", y = "") +
  my_theme()


hom_num <- Syn_hom | Tol_hom | del_hom | CNE_hom | LOF_hom
het_num <- Syn_het | Tol_het | del_het | CNE_het | LOF_het
ratio <- Syn_ratio | Tol_ratio | del_ratio | CNE_ratio | LOF_ratio

del_num <- hom_num / het_num / ratio
del_num <- del_num & theme(
  plot.margin = margin(2, 2, 2, 2)
)

# del_num
del_num <- wrap_plots(
  Syn_hom, Tol_hom, del_hom, CNE_hom, LOF_hom,
  Syn_het, Tol_het, del_het, CNE_het, LOF_het,
  Syn_ratio, Tol_ratio, del_ratio, CNE_ratio, LOF_ratio,
  ncol = 5,
  byrow = TRUE
) &
  theme(plot.margin = margin(2, 2, 2, 2))
del_num

ggsave(del_num, file='plot_2C.pdf', width=16, height=9)

Syn <- Syn_all_derived / Syn_all
Tol <- Tol_all_derived / Tol_all
Mis <- Mis_all_derived / Mis_all
Del <- del_all_derived / del_all
CNE <- CNE_all_derived / CNE_all
LOF <- LOF_all_derived / LOF_all

all_S <- Syn | Tol | Del | CNE | LOF

all_S <- all_S & theme(
  plot.margin = margin(2, 2, 2, 2)
)
all_S

ggsave(all_S, file='plot_S_del.pdf', width=16, height=8)


# /syn ratio 
make_ratio_plot <- function(data, num_prefix, ylab, limits, breaks, bar_width = 0.4) {
  plot_data <- data %>%
    transmute(
      POP,
      homozygous   = .data[[paste0(num_prefix, "_hom")]] / snpSynNum_hom,
      heterozygous = .data[[paste0(num_prefix, "_het")]] / snpSynNum_het
    ) %>%
    pivot_longer(
      cols = c(homozygous, heterozygous),
      names_to = "zygosity",
      values_to = "ratio"
    ) %>%
    mutate(zygosity = factor(zygosity, levels = c("homozygous", "heterozygous")))
  
  ggplot(plot_data, aes(x = POP, y = ratio, fill = POP)) +
    geom_point(
      aes(shape = zygosity, group = zygosity),
      size = 3, color = "NA", alpha = 0.8,
      position = position_dodge(width = 0.55)
    ) +
    stat_summary(
      aes(group = zygosity),
      fun = mean, geom = "crossbar",
      width = bar_width,
      position = position_dodge(width = 0.55)
    ) +
    scale_y_continuous(
      limits = limits,
      breaks = breaks,
      labels = breaks,
      expand = c(0, 0)
    ) +
    scale_shape_manual(values = c(homozygous = 21, heterozygous = 23)) +
    scale_fill_manual(values = fill_cols) +
    guides(fill = "none") +
    xlab(NULL) +
    ylab(ylab) +
    ratio_theme()
}

Del_syn <- make_ratio_plot(
  subset_data, "snpMisDelNum",
  "Deleterious missense/Synonymous variant ratio",
  c(0.09, 0.18), seq(0.09, 0.18, 0.03)
)

LOF_syn <- make_ratio_plot(
  subset_data, "snpLOFNum",
  "LOF/Synonymous variant ratio",
  c(0.009, 0.018), seq(0.009, 0.018, 0.003)
)

CNE_syn <- make_ratio_plot(
  subset_data, "snpCNENum",
  "CNE/Synonymous variant ratio",
  c(1.5, 3), seq(1.5, 3, 0.5)
)

Tol_syn <- make_ratio_plot(
  subset_data, "snpMisTolNum",
  "Tolerated missense/Synonymous variant ratio",
  c(0.3, 0.6), seq(0.3, 0.6, 0.1)
)

Mis_syn <- make_ratio_plot(
  subset_data, "snpMisNum",
  "Missense/Synonymous variant ratio",
  c(0, 0.6), seq(0, 0.6, 0.1)
)

Ratio_plot <- (Del_syn/LOF_syn)|(CNE_syn/Tol_syn)
Ratio_plot
ggsave(Ratio_plot, file='plot_2B.pdf', width=10, height=8)
######################### MIX_CAP ##########################
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(ggrepel)
library(patchwork)
library(ggsignif)


# Read data
p1 <- read.table("merge_num.txt", header = TRUE, sep = "\t", check.names = FALSE)
p_pop <- read.table("pop.txt", header = TRUE, sep = "\t", check.names = FALSE)

data <- merge(p1, p_pop, by = "ID")

# Keep SC_CAP and MIX_CAP only
plot_dat <- data %>%
  filter(lowdepth == "no", POP %in% c("SC_CAP", "MIX_CAP")) %>%
  mutate(
    generation = ifelse(POP == "MIX_CAP", `Pedigree generation`, NA),
    comparison_group = ifelse(POP == "SC_CAP", "SC_CAP", generation),
    comparison_group = factor(
      comparison_group,
      levels = c("F1", "F2", "F3", "F4","SC_CAP"),
      ordered = TRUE
    ),
    het_del = snpMisDelNum_het,
    het_del_syn = snpMisDelNum_het / snpSynNum_het,
    label_id = ifelse(POP == "MIX_CAP", ID, NA)
  )

# MIX_CAP-only data for generation trend test
mix_dat <- plot_dat %>%
  filter(POP == "MIX_CAP") %>%
  mutate(
    generation_num = as.integer(factor(generation, levels = c("F1", "F2", "F3", "F4"), ordered = TRUE))
  )

group_cols <- c(
  SC_CAP = "#F2CC8E",
  F1 = "#D2BFA5",
  F2 = "#D2BFA5",
  F3 = "#D2BFA5",
  F4 = "#D2BFA5"
)

make_generation_plot <- function(df, y_col, y_lab, show_signif = TRUE) {
  jitter_pos <- position_jitter(width = 0.08, height = 0, seed = 123)
  
  y_values <- df[[y_col]]
  y_min <- min(y_values, na.rm = TRUE)
  y_max <- max(y_values, na.rm = TRUE)
  y_range <- y_max - y_min
  
  signif_comparisons <- list(
    c("SC_CAP", "F1"),
    c("SC_CAP", "F3"),
    c("F1", "F3")
  )
  
  signif_y <- c(
    y_max + y_range * 0.08,
    y_max + y_range * 0.18,
    y_max + y_range * 0.28
  )
  
  p <- ggplot(df, aes(x = comparison_group, y = .data[[y_col]], color = comparison_group)) +
    geom_vline(
      xintercept = 4.5,
      linetype = "dashed",
      color = "grey60",
      linewidth = 0.5
    ) +
    geom_boxplot(
      fill = NA,
      color = "black",
      width = 0.55,
      outlier.shape = NA,
      linewidth = 0.8
    ) +
    geom_point(
      position = jitter_pos,
      size = 3,
      alpha = 0.8,
      shape = 16
    ) +
    scale_x_discrete(
      labels = c(
        SC_CAP = "SC_CAP",
        F1 = "MIX_CAP F1",
        F2 = "MIX_CAP F2",
        F3 = "MIX_CAP F3",
        F4 = "MIX_CAP F4"
      )
    ) +
    scale_color_manual(values = group_cols) +
    scale_y_continuous(
      limits = c(y_min - y_range * 0.05, y_max + y_range * 0.40),
      breaks = pretty(y_values),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(x = NULL, y = y_lab) +
    my_theme()
  
  if (show_signif) {
    p <- p +
      geom_signif(
        comparisons = signif_comparisons,
        test = "wilcox.test",
        map_signif_level = FALSE,
        textsize = 3,
        y_position = signif_y,
        tip_length = 0.01,
        inherit.aes = FALSE,
        aes(x = comparison_group, y = .data[[y_col]]),
        color = "black"
      )
  }
  p
}

p_raw <- make_generation_plot(
  df = plot_dat,
  y_col = "het_del",
  y_lab = "Heterozygous deleterious missense count",
  show_signif = TRUE
  ) + 
  scale_y_continuous(limits=c(1300,2100),
                     labels = seq(1300,2100,200),
                     breaks = seq(1300,2100,200),
                     expand = c(0,0)) 

p_ratio <- make_generation_plot(
  df = plot_dat,
  y_col = "het_del_syn",
  y_lab = "Heterozygous deleterious missense / synonymous variant ratio",
  show_signif = TRUE
  )+ ylim(0.14,0.18) +
  scale_y_continuous(limits=c(0.14,0.18),
                     labels = seq(0.14,0.18,0.02),
                     breaks = seq(0.14,0.18,0.02),
                     expand = c(0,0)) 

supp_mix_sc <- p_raw | p_ratio +
  plot_annotation(tag_levels = "A") &
  theme(plot.margin = margin(3, 3, 3, 3))

print(supp_mix_sc)

ggsave(
  filename = "Supp_MIXCAP_generation_SC_CAP.pdf",
  plot = supp_mix_sc,
  width = 10,
  height = 5 
)