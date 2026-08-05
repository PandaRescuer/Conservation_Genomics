# ============================================================
# Plot phyloP distributions and perform statistical tests
# ============================================================

library(ggplot2)

# ===================== User input =====================
# Add or remove groups here.
# The name before "=" will be used as the group label in the plot.


files <- c(
  DEL = "mis_snp_DELETERIOUS.phyloP.tsv",
  TOL = "mis_snp_TOLERATED.phyloP.tsv"
)

out_prefix <- "phyloP_DEL_vs_TOL"

# If your input file has a header, set this to TRUE.
# Your current awk output usually has no header, so FALSE is correct.
has_header <- FALSE

# The phyloP score is in which column?
# Your current file: variant_id phyloP, so score_col = 2.
score_col <- 2

# ===================== Read data =====================

read_phyloP <- function(file, group_name, has_header = FALSE, score_col = 2) {
  dat <- read.table(
    file,
    header = has_header,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  
  score <- suppressWarnings(as.numeric(dat[[score_col]]))
  score <- score[!is.na(score)]
  
  if (length(score) == 0) {
    stop(paste("No numeric phyloP scores found in:", file))
  }
  
  data.frame(
    group = group_name,
    phyloP = score,
    stringsAsFactors = FALSE
  )
}

df_list <- lapply(names(files), function(g) {
  read_phyloP(files[[g]], g, has_header = has_header, score_col = score_col)
})

df <- do.call(rbind, df_list)
df$group <- factor(df$group, levels = names(files))

# ===================== Summary statistics =====================

summary_df <- do.call(rbind, lapply(split(df$phyloP, df$group), function(x) {
  data.frame(
    n = length(x),
    mean = mean(x),
    median = median(x),
    q1 = unname(quantile(x, 0.25)),
    q3 = unname(quantile(x, 0.75)),
    min = min(x),
    max = max(x)
  )
}))

summary_df$group <- rownames(summary_df)
rownames(summary_df) <- NULL
summary_df <- summary_df[, c("group", "n", "mean", "median", "q1", "q3", "min", "max")]

write.table(
  summary_df,
  paste0(out_prefix, ".summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(summary_df)

# ===================== Statistical tests =====================

groups <- levels(df$group)

pairwise_results <- list()
k <- 1

for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    g1 <- groups[i]
    g2 <- groups[j]
    
    x1 <- df$phyloP[df$group == g1]
    x2 <- df$phyloP[df$group == g2]
    
    test <- wilcox.test(x1, x2, exact = FALSE)
    
    pairwise_results[[k]] <- data.frame(
      group1 = g1,
      group2 = g2,
      n1 = length(x1),
      n2 = length(x2),
      mean1 = mean(x1),
      mean2 = mean(x2),
      median1 = median(x1),
      median2 = median(x2),
      p_value = test$p.value,
      stringsAsFactors = FALSE
    )
    
    k <- k + 1
  }
}

pairwise_df <- do.call(rbind, pairwise_results)
pairwise_df$p_adj_BH <- p.adjust(pairwise_df$p_value, method = "BH")

write.table(
  pairwise_df,
  paste0(out_prefix, ".pairwise_wilcoxon.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(pairwise_df)

# Kruskal-Wallis test if there are more than two groups
if (length(groups) > 2) {
  kw <- kruskal.test(phyloP ~ group, data = df)
  
  kw_df <- data.frame(
    test = "Kruskal-Wallis",
    statistic = unname(kw$statistic),
    df = unname(kw$parameter),
    p_value = kw$p.value
  )
  
  write.table(
    kw_df,
    paste0(out_prefix, ".kruskal_wallis.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  print(kw_df)
}

# ===================== Plot =====================

p <- ggplot(df, aes(x = group, y = phyloP, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.65, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9, color = "black") +
  geom_signif(
    comparisons = list(c("DEL", "TOL")),
    map_signif_level = FALSE,
    y_position = 10,
    textsize = 5
  ) +
  geom_hline(aes(yintercept=0),linetype=5,linewidth=1,col="#C0C0C0")+
  scale_y_continuous(
    limits = c(-25,20),
    labels = seq(-25,20,5),
    breaks = seq(-25,20,5),
    expand = c(0,0)) +
  scale_fill_manual(
    values = c(
      DEL = "#DF7A5E",
      TOL = "#F2CC8E"
    )
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "phyloP score"
  ) +
  theme(axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text( size = 15,colour = "black"),
        axis.text.y = element_text( size = 15, colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),  
        axis.ticks = element_line(size = 1),  
        axis.line = element_line(size = 1),  
        legend.position="none",
        text = element_text(family = "Times"),
        strip.background = element_rect(fill = "grey90", color = "black", size = 1),  
        strip.text = element_text(color = "black",size =15)
  )

print(p)

ggsave(
  paste0("plot",out_prefix, ".pdf"),
  plot = p,
  width = 4,
  height = 4
)

cat("Done.\n")
cat("Output prefix:", out_prefix, "\n")