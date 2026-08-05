library(ggplot2)
library(dplyr)
library(tidyr)
library(ggbeeswarm)
library(ggsignif)
library(ggh4x)

setwd("bcftools_roh")

my_theme <- function() {
  theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
      axis.title.x = element_text(face = "bold", size = 15, colour = "black"),
      axis.title.y = element_text(face = "bold", size = 15, colour = "black"),
      axis.text.x = element_text(face = "bold", size = 15, colour = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(face = "bold", size = 15, colour = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.ticks = element_line(size = 1),
      axis.line = element_line(size = 1),
      legend.position = "none",
      legend.text = element_text(size = 15),
      text = element_text(family = "Times")
    )
}


###############################################################################
# Load bcftools ROH length bins
ROH_data <- read.table("bcftools_roh.length.tsv", header = TRUE, sep = "\t", check.names = FALSE)
# Make sure column names are consistent
colnames(ROH_data) <- c("ID", "0.1~1MB", "1~5MB", "5+MB")
# Convert to long format
mydata_long <- pivot_longer(
  ROH_data,
  cols = c("0.1~1MB", "1~5MB", "5+MB"),
  names_to = "Size",
  values_to = "lengthofROH"
)
# Load population and heterozygosity metadata
data_pop <- read.table("pop.txt", header = TRUE)
data_het <- read.table("het.txt", header = TRUE)
# Merge metadata
add_PoP_set <- merge.data.frame(mydata_long, data_pop, by = "ID")
het_ROH_data <- merge.data.frame(add_PoP_set, data_het, by = "NAME")
# Remove low-depth samples
notfsgset <- subset(het_ROH_data, lowdepth != "yes")
# Population order
desired_order <- c("SC_WSP", "QLI", "SC_WLP", "SC_CAP", "MIX_CAP")
notfsgset$POP <- factor(
  notfsgset$POP,
  levels = desired_order,
  ordered = TRUE
)
# Facet strip color
strip_theme <- strip_themed(
  background_y = elem_list_rect(
    fill = c("#8e7cc3", "#82B29A", "#DF7A5E", "#F2CC8E", "#D2BFA5")
  )
)
# Optional: use your manually defined sample order
sample_order_bottom_to_top <- c(
  "GP51", "GP52", "GP53", "GP50", "GP48", "GP49", "GP39", "GP45", "GP41", "GP42",
  "GP47", "GP25", "GP17", "GP13", "GP35", "GP28", "GP30", "GP19", "GP12", "GP44",
  "GP29", "GP26", "GP40", "GP21", "GP37", "GP43", "GP24", "GP16", "GP46", "GP27",
  "GP22", "GP34", "GP04", "GP06", "GP15", "GP18", "GP31", "GP32", "GP03", "GP14",
  "GP07", "GP38", "GP01", "GP33", "GP10", "GP02", "GP09", "GP23", "GP05", "GP08",
  "GP20", "GP36", "GP55", "QLI07c", "QLI06c", "QLI05c", "GP54", "QLI04c", "GP11"
)
notfsgset$ID <- factor(
  notfsgset$ID,
  levels = sample_order_bottom_to_top
)
# Make ROH length plot
plot_ROH_het_bcftools <-
  ggplot(notfsgset, aes(y = ID)) +
  geom_col(
    aes(x = lengthofROH / 1000, fill = Size),
    position = position_stack(reverse = TRUE),
    width = 0.8,
    color = "black",
    linewidth = 0.05
  ) +
  geom_point(
    aes(x = het_whole_genome * 800000, fill = POP4, shape = POP4),
    size = 4,
    alpha = 0.8
  ) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(
    text = element_text(size = 15, family = "Times"),
    strip.text = element_text(size = 12.5),
    legend.position = "none"
  ) +
  facet_grid2(
    POP ~ .,
    scales = "free",
    space = "free",
    strip = strip_theme
  ) +
  scale_shape_manual(
    values = c(
      SC_CAP = 21,
      QLI = 22,
      SC_QLA = 24,
      SC_MSH = 25,
      MIX_CAP = 23,
      SC_XXL = 10
    )
  ) +
  theme(strip.text.y = element_text(angle = 0)) +
  scale_y_discrete(expand = expansion(add = 1.5)) +
  scale_x_continuous(
    sec.axis = sec_axis(
      ~ . / 800,
      breaks = seq(0, 1.5, by = 0.3),
      name = "Heterozygosity per kbp"
    ),
    limits = c(0, 1200)
  ) +
  scale_fill_manual(
    values = c(
      "0.1~1MB" = "#ddbfa3",
      "1~5MB" = "#9c7554",
      "5+MB" = "#5a2d18",
      SC_CAP = "#F2CC8E",
      QLI = "#82B29A",
      SC_QLA = "#df995e",
      SC_MSH = "#DF7A5E",
      MIX_CAP = "#D2BFA5",
      SC_WSP = "#8e7cc3"
    )
  ) +
  labs(x = "Length of ROH (Mb)", y = "")
plot_ROH_het_bcftools
ggsave(
  plot_ROH_het_bcftools,
  file = "ROH_bcftools.pdf",
  device = cairo_pdf,
  width = 8.5,
  height = 13
)
