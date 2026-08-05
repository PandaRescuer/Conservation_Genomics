library(ggplot2)
library(reshape2)
library(dplyr)


SFS_raw <- read.table("merge_SFS.txt", header = TRUE, check.names = FALSE)


# Convert wide format to long format
SFS_melt <- melt(
  SFS_raw,
  id.vars = c("POP", "DAF"),
  measure.vars = c("Synonymous", "Tolerated", "Missense", "Deleterous", "CNE", "LOF"),
  variable.name = "Variant_type",
  value.name = "Count"
)

# Population order
SFS_melt$POP <- factor(
  SFS_melt$POP,
  levels = c("MIX_CAP", "SC_CAP", "SC_WLP", "QLI")
)

# Keep representative variant classes
SFS_summary_data <- SFS_melt %>%
  filter(Variant_type %in% c("Synonymous", "Tolerated", "Deleterous", "CNE", "LOF")) %>%
  mutate(
    Variant_type = factor(
      Variant_type,
      levels = c("Synonymous", "Tolerated", "Deleterous", "CNE", "LOF"),
      labels = c(
        "Synonymous",
        "Tolerated\nmissense",
        "Deleterious\nmissense",
        "CNE",
        "LOF"
      )
    )
  )
# Calculate proportions of low-frequency and high-frequency derived alleles
SFS_summary <- SFS_summary_data %>%
  group_by(POP, Variant_type) %>%
  summarise(
    total_count = sum(Count, na.rm = TRUE),
    low_count = sum(Count[DAF <= 0.2], na.rm = TRUE),
    high_count = sum(Count[DAF >= 0.8], na.rm = TRUE),
    low_prop = low_count / total_count,
    high_prop = high_count / total_count,
    .groups = "drop"
  )

# Convert to long format
SFS_summary_long <- SFS_summary %>%
  select(POP, Variant_type, low_prop, high_prop) %>%
  melt(
    id.vars = c("POP", "Variant_type"),
    variable.name = "Metric",
    value.name = "Proportion"
  )

# Rename metrics
SFS_summary_long$Metric <- factor(
  SFS_summary_long$Metric,
  levels = c("low_prop", "high_prop"),
  labels = c(
    "Low-frequency derived alleles: DAF ≤ 0.2",
    "High-frequency derived alleles: DAF ≥ 0.8"
  )
)

# Population colors
pop_colors <- c(
  MIX_CAP = "#D2BFA5",
  SC_CAP  = "#F2CC8E",
  SC_WLP  = "#DF7A5E",
  QLI     = "#82B29A"
)

pop_shapes <- c(
  SC_CAP  = 21,
  QLI     = 22,
  SC_WLP  = 24,
  MIX_CAP = 23
)

# Line types for frequency classes
metric_linetypes <- c(
  "Low-frequency derived alleles: DAF ≤ 0.2" = "solid",
  "High-frequency derived alleles: DAF ≥ 0.8" = "31"
)

# Point shapes for frequency classes
# metric_shapes <- c(
#   "Low-frequency derived alleles: DAF ≤ 0.2" = 16,
#   "High-frequency derived alleles: DAF ≥ 0.8" = 17
# )

p_sfs_combined <-
  ggplot(
    SFS_summary_long,
    aes(
      x = Variant_type,
      y = Proportion,
      group = interaction(POP, Metric)
    )
  ) +
  geom_line(
    aes(
      color = POP,
      linetype = Metric
    ),
    linewidth = 1.0,
    alpha = 0.9
  ) +
  geom_point(
    aes(
      fill = POP,
      shape = POP
    ),
    color = "black",
    size = 4,
    stroke = 0.5,
    alpha = 0.8
  ) +
  scale_color_manual(values = pop_colors) +
  scale_fill_manual(values = pop_colors) +
  scale_shape_manual(values = pop_shapes) +
  scale_linetype_manual(values = metric_linetypes) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Proportion of variants",
    color = "Population",
    linetype = NULL,
    shape = NULL
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    axis.title.y = element_text(size = 20, colour = "black"),
    axis.text.x = element_text(size = 15, colour = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 15, colour = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1),
    legend.position = "right",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    text = element_text(family = "Times")
  )

p_sfs_combined

ggsave(
  p_sfs_combined,
  file = "SFS_summary_combined_Fig2AS.pdf",
  width = 10,
  height = 8,
  device = cairo_pdf
)