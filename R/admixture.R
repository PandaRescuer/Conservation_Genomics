library(ggplot2)
library(dplyr)
library(grid)

setwd("admixture")

# Read CV error results
df <- read.csv("admixture.cv_error.csv")
df$K <- as.numeric(df$K)
df$seed <- as.factor(df$seed)
df$CV_error <- as.numeric(df$CV_error)

summary_df <- df %>%
  group_by(K) %>%
  summarise(
    mean_cv = mean(CV_error, na.rm = TRUE),
    sd_cv = sd(CV_error, na.rm = TRUE),
    .groups = "drop"
  )

# Read ADMIXTURE similarity results
sim_df <- read.delim(
  "admixture_similarity.summary.tsv",
  header = TRUE,
  check.names = FALSE
)

sim_df$K <- as.numeric(sim_df$K)
sim_df$Average_pairwise_similarity <-
  as.numeric(sim_df$Average_pairwise_similarity)

# Define transformation between similarity and CV-error axes
cv_min <- 0.53
cv_max <- 0.55

sim_to_cv <- function(x) {
  cv_min + x * (cv_max - cv_min)
}

cv_to_sim <- function(x) {
  (x - cv_min) / (cv_max - cv_min)
}

sim_df$similarity_scaled <-
  sim_to_cv(sim_df$Average_pairwise_similarity)
p <- ggplot() +
  # Individual CV-error values
  geom_line(
    data = df,
    aes(x = K, y = CV_error, group = seed),
    color = "grey70",
    alpha = 0.8,
    linewidth = 0.8
  ) +
  geom_point(
    data = df,
    aes(x = K, y = CV_error),
    color = "grey70",
    alpha = 0.9,
    size = 2
  ) +
  
  # Mean CV error
  geom_errorbar(
    data = summary_df,
    aes(
      x = K,
      ymin = mean_cv - sd_cv,
      ymax = mean_cv + sd_cv
    ),
    width = 0.1,
    color = "black",
    linewidth = 0.6
  ) +
  geom_line(
    data = summary_df,
    aes(
      x = K,
      y = mean_cv,
      color = "CV error",
      linetype = "CV error"
    ),
    linewidth = 1.2
  ) +
  geom_point(
    data = summary_df,
    aes(x = K, y = mean_cv, color = "CV error"),
    size = 3
  ) +
  
  # Average pairwise similarity
  geom_line(
    data = sim_df,
    aes(
      x = K,
      y = similarity_scaled,
      color = "Average pairwise similarity",
      linetype = "Average pairwise similarity",
      group = 1
    ),
    linewidth = 1.6
  ) +
  geom_point(
    data = sim_df,
    aes(
      x = K,
      y = similarity_scaled,
      color = "Average pairwise similarity"
    ),
    size = 3
  ) +
  
  scale_x_continuous(
    breaks = sort(unique(df$K))
  ) +
  
  scale_y_continuous(
    name = "CV error",
    limits = c(cv_min, cv_max),
    breaks = seq(cv_min, cv_max, 0.01),
    sec.axis = sec_axis(
      trans = ~ cv_to_sim(.),
      name = "Average pairwise similarity",
      breaks = seq(0, 1, 0.2)
    )
  ) +
  
  scale_color_manual(
    values = c(
      "CV error" = "black",
      "Average pairwise similarity" = "#D95F02"
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      "CV error" = "solid",
      "Average pairwise similarity" = "11"
    )
  ) +
  
  labs(
    x = "K",
    color = NULL,
    linetype = NULL
  ) +
  
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      keywidth = unit(1, "cm"),
      keyheight = unit(0.7, "cm")
    ),
    linetype = guide_legend(
      nrow = 1,
      byrow = TRUE,
      keywidth = unit(1, "cm"),
      keyheight = unit(0.7, "cm")
    )
  ) +
  
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    axis.title   = element_text(face = "bold", size = 15, colour = "black"),
    axis.text    = element_text(face = "bold", size = 15, colour = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks   = element_line(linewidth = 1),
    axis.line    = element_line(linewidth = 1),
    legend.position = "top",
    legend.text  = element_text(size = 13),
    legend.key.width  = unit(2.8, "cm"),
    legend.key.height = unit(0.7, "cm"),
    legend.spacing.x  = unit(0.5, "cm"),
    legend.margin     = margin(b = 5),
    text = element_text(family = "Times")
  )


p

ggsave("CV_error.pdf",device = cairo_pdf, p, width = 5, height = 5)


library(ggplot2)
setwd("admixture")
mydat <- read.table("admixture_out.txt",header = T)
names<-c(
  'QLI01', 'QLI02', 'QLI03', 'QLI05c', 'QLI06c', 'QLI07c', 'QLI04c', 'GP54', 'GP55', 'GP52', 'GP51', 'GP50', 'GP48','GP53', 'GP49', 'GP42', 'GP39', 'GP41', 'GP17', 'QLA01c', 'QLA02c', 'GP14', 'GP15', 'GP20', 'GP23', 'GP33', 'GP18', 'GP28', 'GP24', 'GP32', 'GP47', 'GP40', 'GP01', 'GP02', 'GP11', 'GP38', 'QLA05', 'GP08', 'GP27', 'GP05', 'GP03', 'QLA03', 'GP06', 'QLA07', 'QLA06', 'QLA08', 'GP31', 'QLA04', 'GP35', 'GP04', 'GP07', 'GP34', 'GP37', 'GP19', 'GP13', 'GP26', 'GP21', 'GP30', 'GP29', 'GP43', 'GP12', 'GP44', 'GP22', 'GP45', 'GP25', 'GP46', 'GP09', 'MSH09', 'MSH05', 'MSH03', 'MSH01', 'MSH04', 'MSH07', 'MSH02', 'MSH06', 'MSH08', 'GP16', 'GP36', 'GP10'
)
mydat$sample_name <- factor(mydat$sample_name,
                            levels = names,
                            ordered = T)
ggplot(mydat,aes(x=sample_name, weight = percent, fill = group)) +
  geom_bar(position="stack") +theme_classic()+
  theme(axis.text.x = element_text(angle =45,hjust=1.2,vjust=1.6,size=5))+
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())+
  labs(x="",y="")+
  scale_fill_manual(values=c(QLI = "#82B29A",QLA = "#F2CC8E", MSH = "#DF7A5E"))


library(ggplot2)
setwd("admixture")
names<-c(
  'QLI01', 'QLI02', 'QLI03', 'QLI05c', 'QLI06c', 'QLI07c', 'QLI04c', 'GP54', 'GP55','GP52', 'GP51', 'GP50', 'GP48','GP53', 'GP49', 'GP45','GP42', 'GP41','GP47','GP39', 'GP11', 'GP28','GP17','GP14', 'GP15', 'GP20', 'GP23', 'GP33', 'GP18',  'GP24', 'GP32', 'GP40', 'GP38', 'GP27' , 'GP31',  'GP35',  'GP34', 'GP37', 'GP19', 'GP13', 'GP26', 'GP21', 'GP30', 'GP29', 'GP43', 'GP12', 'GP44', 'GP22', 'GP25','GP46','GP16', 'GP36','GP01', 'GP02','GP05', 'GP03','GP06','GP04', 'GP07','GP08','QLA01c', 'QLA02c','QLA05', 'QLA04','QLA03','QLA07', 'QLA06', 'QLA08',  'GP09', 'GP10','MSH09', 'MSH05', 'MSH03', 'MSH01', 'MSH04', 'MSH07', 'MSH02', 'MSH06', 'MSH08'
  
)
#k=2
mydat2 <- read.table("admixture2_out.txt",header = T)
mydat2$sample_name <- factor(mydat2$sample_name,
                             levels = names,
                             ordered = T)
#k=3
mydat3 <- read.table("admixture_out.txt",header = T)

mydat3$sample_name <- factor(mydat3$sample_name,
                             levels = names,
                             ordered = T)

#k=4
mydat4 <- read.table("admixture4_out.txt",header = T)
mydat4$sample_name <- factor(mydat4$sample_name,
                             levels = names,
                             ordered = T)


admix_plot_2<-ggplot(mydat2,aes(x=sample_name, weight = percent, fill = factor(group, levels = c("SC", "QLI")))) +
  geom_bar(position="stack",width=1,linewidth=0.5) +theme_classic()+
  #geom_bar(position="stack",width=1.1) +theme_classic()+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())+
  guides(fill="none")+
  labs(x="",y="")+
  scale_fill_manual(values=c(QLI = "#82B29A",SC="#F2CC8E"))

admix_plot_3 <- ggplot(mydat3,aes(x=sample_name, weight = percent, fill = factor(group, levels = c("MSH", "QLA", "QLI")))) +
  geom_bar(position="stack",width=1,,linewidth=0.5) +theme_classic()+
  #geom_bar(position="stack",width=1.1) +theme_classic()+
  #theme(axis.text.x = element_text(angle =45,face = "bold",hjust=1.2,vjust=1.6,size=7))+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())+
  guides(fill="none")+
  labs(x="",y="")+
  scale_fill_manual(values=c(QLI = "#82B29A",QLA = "#df995e", MSH = "#DF7A5E"))

admix_plot_4<-ggplot(mydat4,aes(x=sample_name, weight = percent, fill = factor(group, levels = c("MSH","Other", "QLA", "QLI")))) +
  geom_bar(position="stack",width=1,,linewidth=0.5) +theme_classic()+
  # geom_bar(position="stack",width=1,color="white",linewidth=0.5) +theme_classic()+
  #geom_bar(position="stack",width=1.1) +theme_classic()+
  theme(axis.text.x = element_text(angle =45,hjust=1.2,vjust=1.6,size=5))+
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())+
  guides(fill="none")+
  labs(x="",y="")+
  scale_fill_manual(values=c(QLI = "#82B29A",QLA = "#df995e", MSH = "#DF7A5E",Other="#F2CC8E"))

merge_plot<-admix_plot_2/admix_plot_3/admix_plot_4
merge_plot

ggsave(merge_plot, file='admixture.pdf',device = cairo_pdf, width=16, height=12)

