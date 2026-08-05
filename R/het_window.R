library(dplyr)
library(ggplot2)
library(grid)
windowsFonts(Times = windowsFont("Times New Roman"))


# Global files
het_mean_file <- "het_mean.txt"
roh_file <- "ROH.renamed.hom.indiv"

# Population colors
pop_colors <- c(
  SC_CAP = "#F2CC8E",
  QLI = "#82B29A",
  SC_QLA = "#df995e",
  SC_MSH = "#DF7A5E",
  SC_WLP = "#DF7A5E",
  MIX_CAP = "#D2BFA5",
  SC_XXL = "#8e7cc3",
  SC_WSP = "#8e7cc3"
)

# Read all window files in one folder
read_window_data <- function(path) {
  files <- list.files(path, pattern = "\\.txt$", full.names = TRUE)
  
  do.call(rbind, lapply(files, function(file) {
    sample_name <- sub("\\..*", "", basename(file))
    message(sample_name)
    
    df <- read.table(file, header = FALSE)
    colnames(df) <- c("chr", "bp", "het")
    df$sample <- sample_name
    df
  }))
}

# Prepare merged data
prepare_data <- function(path, pop_col = "POP2", genome_size = 2173986137) {
  all_data <- read_window_data(path)
  
  het_data <- read.table(het_mean_file, header = TRUE)
  roh_data <- read.table(roh_file, header = TRUE)
  roh_data$FROH <- (roh_data$KB * 1000) / genome_size
  
  sample_means <- het_data %>%
    select(sample, all_of(pop_col), het_whole_genome) %>%
    mutate(het_mean = round(het_whole_genome, 8)) %>%
    select(sample, all_of(pop_col), het_mean)
  
  colnames(sample_means)[2] <- "POP_GROUP"
  
  all_data <- merge(all_data, sample_means, by = "sample")
  sample_means <- merge(
    sample_means,
    roh_data[, c("sample", "FROH")],
    by = "sample"
  )
  
  list(all_data = all_data, sample_means = sample_means)
}

# Build facet labels
make_facet_labels <- function(sample_means, froh_digits = 4) {
  facet_labels <- paste0(
    as.character(sample_means$sample),
    "   mean het. = ",
    sprintf("%.4f", sample_means$het_mean * 1000),
    "   FROH = ",
    sprintf(paste0("%.", froh_digits, "f"), sample_means$FROH * 100),
    " %"
  )
  
  names(facet_labels) <- as.character(sample_means$sample)
  facet_labels
}

# Plot function
plot_window_het <- function(path,
                            output_file,
                            pop_col = "POP2",
                            ncol = 4,
                            width = 17,
                            height = 13,
                            sample_order = NULL,
                            froh_digits = 4) {
  
  dat <- prepare_data(path, pop_col = pop_col)
  all_data <- dat$all_data
  sample_means <- dat$sample_means
  
  if (!is.null(sample_order)) {
    all_data$sample <- factor(all_data$sample, levels = sample_order)
    sample_means$sample <- factor(sample_means$sample, levels = sample_order)
    
    all_data <- all_data %>% arrange(sample)
    sample_means <- sample_means %>% arrange(sample)
  }
  
  facet_labels <- make_facet_labels(sample_means, froh_digits = froh_digits)
  
  p <- ggplot(all_data, aes(x = bp / 1000000, y = het * 1000,
                            fill = POP_GROUP, color = POP_GROUP)) +
    geom_col(linewidth = 0.2) +
    facet_wrap(
      ~sample,
      ncol = ncol,
      labeller = labeller(sample = facet_labels)
    ) +
    scale_fill_manual(values = pop_colors) +
    scale_color_manual(values = pop_colors) +
    scale_y_continuous(
      limits = c(0, 6),
      labels = seq(0, 6, 2),
      breaks = seq(0, 6, 2),
      expand = c(0, 0)
    ) +
    labs(x = "Chromosome (Mb)", y = "Heterozygosity per kbp") +
    theme_classic() +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 1),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x.bottom = element_text(size = 12.5, colour = "black"),
      axis.text.y.left = element_text(size = 12.5, colour = "black"),
      axis.ticks.x.bottom = element_line(size = 1),
      axis.ticks.y.left = element_line(size = 1),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position = "none",
      text = element_text(family = "Times"),
      axis.title.y = element_text(size = 12.5, colour = "black"),
      axis.title.x = element_text(size = 12.5, colour = "black"),
      strip.text = element_text(size = 12.5, colour = "black"),
      strip.background = element_blank()
    )
  
  ggsave(output_file, p, width = width, height = height, limitsize = FALSE)
  return("over")
}
plot_window_het(
  path = "window_plot_no_XXL",
  output_file = "combined_facet_plot_no_XXL.pdf",
  pop_col = "POP2",
  ncol = 4,
  width = 17,
  height = 13
)
plot_window_het(
  path = "window_plot_fang",
  output_file = "combined_facet_plot_fang.pdf",
  pop_col = "POP2",
  ncol = 4,
  width = 17,
  height = 6
)
plot_window_het(
  path = "window_plot_only_XXL",
  output_file = "combined_facet_plot_only_XXL2.pdf",
  pop_col = "POP2",
  ncol = 1,
  width = 5,
  height = 5
)
sample_order <- c("GP11", "QLI06c", "GP01", "GP28", "GP50")
plot_window_het(
  path = "figure1",
  output_file = "combined_facet_plot_figure1.pdf",
  pop_col = "POP",
  ncol = 1,
  width = 8,
  height = 9,
  sample_order = sample_order
)
