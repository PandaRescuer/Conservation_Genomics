library(ggplot2)
library(ggsignif)
library(ggbeeswarm)
library(dplyr)
library(patchwork)
library(rlang) 

setwd("D:\\重测序项目\\分析部分\\slim\\final\\stat")

my_theme <- function() {
  theme_classic()+
    theme(axis.title.x = element_text(face = "bold", size = 15, colour = "black"),
          axis.title.y = element_text(face = "bold", size = 15, colour = "black"),
          axis.text.x = element_text(face = "bold", size = 15,colour = "black",angle = 45,hjust=1),
          axis.text.y = element_text(face = "bold", size = 15, colour = "black"),
          axis.ticks.length = unit(0.2, "cm"), 
          axis.ticks = element_line(size = 1),  
          axis.line = element_line(size = 1),  
          legend.position="none",
          legend.text = element_text(size=15),
          text = element_text(family = "Times")
    )
}

my_theme2 <- function() {
  theme_classic()+
    theme(axis.title.x = element_text(face = "bold", size = 15, colour = "black"),
          axis.title.y = element_text(face = "bold", size = 15, colour = "black"),
          axis.text.x = element_text(face = "bold", size = 15,colour = "black"),
          axis.text.y = element_text(face = "bold", size = 15, colour = "black"),
          axis.ticks.length = unit(0.2, "cm"),  
          axis.ticks = element_line(size = 1),  
          axis.line = element_line(size = 1),  
          legend.position="none",
          legend.text = element_text(size=15),
          text = element_text(family = "Times")
    )
}

top_theme <- function() {
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank())
}

add_significance <- function(
    df,
    group_col,
    y_col,
    comparisons = NULL,
    group_levels = NULL,
    p_cutoff = 0.05,
    y_position = NULL,
    ylim = NULL,
    bracket_start = 0.72,
    bracket_step = 0.06,
    textsize = 3,
    min_n = 2,
    p_digits = 3,
    p_eps = 0.001
) {
  df_test <- df %>%
    dplyr::filter(
      !is.na(.data[[group_col]]),
      !is.na(.data[[y_col]])
    )
  
  # Use specified levels, factor levels, or observed data order.
  if (is.null(group_levels)) {
    group_values <- df_test[[group_col]]
    
    if (is.factor(group_values)) {
      group_levels <- levels(droplevels(group_values))
    } else {
      group_levels <- unique(as.character(group_values))
    }
  }
  
  df_test <- df_test %>%
    dplyr::mutate(
      .signif_group = factor(
        as.character(.data[[group_col]]),
        levels = group_levels
      )
    )
  
  # Compare every pair of groups when comparisons are not supplied.
  if (is.null(comparisons)) {
    if (length(group_levels) < 2) {
      return(NULL)
    }
    
    comparisons <- combn(group_levels, 2, simplify = FALSE)
  }
  
  # Calculate pairwise Wilcoxon rank-sum test P values.
  p_values <- vapply(
    comparisons,
    function(comp) {
      x <- df_test %>%
        dplyr::filter(.signif_group == comp[1]) %>%
        dplyr::pull(.data[[y_col]])
      
      y <- df_test %>%
        dplyr::filter(.signif_group == comp[2]) %>%
        dplyr::pull(.data[[y_col]])
      
      if (length(x) < min_n || length(y) < min_n) {
        return(NA_real_)
      }
      
      suppressWarnings(
        stats::wilcox.test(x, y, exact = FALSE)$p.value
      )
    },
    numeric(1)
  )
  
  significant <- !is.na(p_values) & p_values < p_cutoff
  
  if (!any(significant)) {
    return(NULL)
  }
  
  significant_comparisons <- comparisons[significant]
  significant_pvalues <- p_values[significant]
  n_significant <- length(significant_comparisons)
  
  if (is.null(y_position)) {
    
    # Use plot limits when supplied; otherwise use the observed data range.
    if (is.null(ylim)) {
      y_range <- range(df_test[[y_col]], na.rm = TRUE)
    } else {
      if (
        !is.numeric(ylim) ||
        length(ylim) != 2 ||
        any(!is.finite(ylim)) ||
        ylim[1] >= ylim[2]
      ) {
        stop("`ylim` must be a numeric vector: c(ymin, ymax).")
      }
      
      y_range <- ylim
    }
    
    y_span <- diff(y_range)
    
    if (!is.finite(y_span) || y_span <= 0) {
      stop("The y-axis range must have a positive finite span.")
    }
    
    if (
      !is.numeric(bracket_start) ||
      length(bracket_start) != 1 ||
      !is.finite(bracket_start) ||
      bracket_start < 0 ||
      bracket_start > 1
    ) {
      stop("`bracket_start` must be a single number between 0 and 1.")
    }
    
    if (
      !is.numeric(bracket_step) ||
      length(bracket_step) != 1 ||
      !is.finite(bracket_step) ||
      bracket_step <= 0
    ) {
      stop("`bracket_step` must be a single positive number.")
    }
    
    # Use relative y-axis positions to keep visual spacing consistent across plots.
    bracket_fraction <- bracket_start +
      (seq_len(n_significant) - 1) * bracket_step
    
    if (any(bracket_fraction >= 1)) {
      warning(
        "Some significance brackets are at or above the y-axis upper limit. ",
        "Decrease `bracket_start` or `bracket_step`, or increase the plot y-axis limit."
      )
    }
    
    y_position <- y_range[1] + bracket_fraction * y_span
    
  } else if (length(y_position) != n_significant) {
    stop(
      "`y_position` must have the same length as the number of ",
      "significant comparisons: ",
      n_significant
    )
  }
  
  # Report exact P values unless they are below the reporting threshold.
  p_eps_label <- formatC(
    p_eps,
    format = "fg",
    digits = p_digits
  )
  
  p_labels <- ifelse(
    significant_pvalues < p_eps,
    paste0(" < ", p_eps_label),
    formatC(
        significant_pvalues,
        format = "fg",
        digits = p_digits
      
    )
  )
  
  ggsignif::geom_signif(
    comparisons = significant_comparisons,
    annotations = p_labels,
    y_position = y_position,
    textsize = textsize,
    map_signif_level = FALSE,
    step_increase = 0
  )
}

add_significance_adjoin <- function(
    df,
    group_col,
    y_col,
    comparisons = NULL,
    group_levels = NULL,
    p_cutoff = 0.05,
    y_position = NULL,
    ylim = NULL,
    bracket_start = 0.72,
    bracket_step = 0.06,
    textsize = 3,
    min_n = 2,
    p_digits = 3,
    p_eps = 0.001
) {
  df_test <- df %>%
    dplyr::filter(
      !is.na(.data[[group_col]]),
      !is.na(.data[[y_col]])
    )
  
  # Use specified levels, factor levels, or observed data order.
  if (is.null(group_levels)) {
    group_values <- df_test[[group_col]]
    
    if (is.factor(group_values)) {
      group_levels <- levels(droplevels(group_values))
    } else {
      group_levels <- unique(as.character(group_values))
    }
  }
  
  df_test <- df_test %>%
    dplyr::mutate(
      .signif_group = factor(
        as.character(.data[[group_col]]),
        levels = group_levels
      )
    )
  
  # Compare every pair of groups when comparisons are not supplied.
  if (is.null(comparisons)) {
    if (length(group_levels) < 2) {
      return(NULL)
    }
    
    # Compare adjacent groups only.
    comparisons <- lapply(
      seq_len(length(group_levels) - 1),
      function(i) group_levels[c(i, i + 1)]
    )
  }
  
  # Calculate pairwise Wilcoxon rank-sum test P values.
  p_values <- vapply(
    comparisons,
    function(comp) {
      x <- df_test %>%
        dplyr::filter(.signif_group == comp[1]) %>%
        dplyr::pull(.data[[y_col]])
      
      y <- df_test %>%
        dplyr::filter(.signif_group == comp[2]) %>%
        dplyr::pull(.data[[y_col]])
      
      if (length(x) < min_n || length(y) < min_n) {
        return(NA_real_)
      }
      
      suppressWarnings(
        stats::wilcox.test(x, y, exact = FALSE)$p.value
      )
    },
    numeric(1)
  )
  
  significant <- !is.na(p_values) & p_values < p_cutoff
  
  if (!any(significant)) {
    return(NULL)
  }
  
  significant_comparisons <- comparisons[significant]
  significant_pvalues <- p_values[significant]
  n_significant <- length(significant_comparisons)
  
  if (is.null(y_position)) {
    
    # Use plot limits when supplied; otherwise use the observed data range.
    if (is.null(ylim)) {
      y_range <- range(df_test[[y_col]], na.rm = TRUE)
    } else {
      if (
        !is.numeric(ylim) ||
        length(ylim) != 2 ||
        any(!is.finite(ylim)) ||
        ylim[1] >= ylim[2]
      ) {
        stop("`ylim` must be a numeric vector: c(ymin, ymax).")
      }
      
      y_range <- ylim
    }
    
    y_span <- diff(y_range)
    
    if (!is.finite(y_span) || y_span <= 0) {
      stop("The y-axis range must have a positive finite span.")
    }
    
    if (
      !is.numeric(bracket_start) ||
      length(bracket_start) != 1 ||
      !is.finite(bracket_start) ||
      bracket_start < 0 ||
      bracket_start > 1
    ) {
      stop("`bracket_start` must be a single number between 0 and 1.")
    }
    
    if (
      !is.numeric(bracket_step) ||
      length(bracket_step) != 1 ||
      !is.finite(bracket_step) ||
      bracket_step <= 0
    ) {
      stop("`bracket_step` must be a single positive number.")
    }
    
    # Use relative y-axis positions to keep visual spacing consistent across plots.
    bracket_fraction <- bracket_start +
      (seq_len(n_significant) - 1) * bracket_step
    
    if (any(bracket_fraction >= 1)) {
      warning(
        "Some significance brackets are at or above the y-axis upper limit. ",
        "Decrease `bracket_start` or `bracket_step`, or increase the plot y-axis limit."
      )
    }
    
    y_position <- y_range[1] + bracket_fraction * y_span
    
  } else if (length(y_position) != n_significant) {
    stop(
      "`y_position` must have the same length as the number of ",
      "significant comparisons: ",
      n_significant
    )
  }
  
  # Report exact P values unless they are below the reporting threshold.
  p_eps_label <- formatC(
    p_eps,
    format = "fg",
    digits = p_digits
  )
  
  p_labels <- ifelse(
    significant_pvalues < p_eps,
    paste0(" < ", p_eps_label),
    formatC(
      significant_pvalues,
      format = "fg",
      digits = p_digits
      
    )
  )
  
  ggsignif::geom_signif(
    comparisons = significant_comparisons,
    annotations = p_labels,
    y_position = y_position,
    textsize = textsize,
    map_signif_level = FALSE,
    step_increase = 0
  )
}

# time compare ------------------------------------------------------------


draw_after <- function(picture,scenarios_item){
  picture <- picture +
    geom_line(
      data=subset(slim_result,scenarios==scenarios_item),
      mapping=aes(x=num, y=after, color = scenarios),
      size=1) +
    geom_ribbon(
      data = subset(slim_result, scenarios == scenarios_item),
      mapping = aes(x = num, ymin = after_ci_min, ymax = after_ci_max, fill =scenarios),
      alpha = 0.1
    ) +
    geom_point(
      data = subset(slim_result, scenarios == scenarios_item),
      mapping = aes(x = num, y = after, color = scenarios,shape=scenarios),
      size = 3
    )
  return(picture)
}






# N
slim_result <- read.table('time2ex.txt',header=F,sep='\t',fill = T)
names(slim_result) <- c("scenarios","num","after","after_ci_min","after_ci_max","before","before_ci_min","before_ci_max")
slim_result[slim_result < 0 & !is.na(slim_result)] <- 0 

pic <- ggplot()
pic <- draw_after(pic,'N10')
pic <- draw_after(pic,'N20')
pic <- draw_after(pic,'N30')
pic <- draw_after(pic,'N40')
pic <- draw_after(pic,'N50')
pic <- draw_after(pic,'N100')
N_after_pic <- 
  pic +
  geom_hline(aes(yintercept=20),linetype=5,linewidth=1,col="#C0C0C0")+
  scale_shape_manual(values = c(0,1,2,3,4,5,8))+
  scale_fill_manual(values=c(N10 = "#DF7A5E", N20 = "#F2CC8E", N30 = "#82B29A",N40="#D2BFA5",N50="#8e7cc3",N100="#8FB3DC"))+
  scale_color_manual(values=c(N10 = "#DF7A5E", N20 = "#F2CC8E", N30 = "#82B29A",N40="#D2BFA5",N50="#8e7cc3",N100="#8FB3DC"))+
  labs(x = "Population size", y = "Survival time to extinction")+
  scale_y_continuous(limits=c(0,100),
                     labels = seq(0,100,20),
                     breaks = seq(0,100,20),
  ) +
  xlim(2,8)+
  my_theme2()

N_after_pic

pic <- ggplot()
pic <- draw_after(pic,'K10')
pic <- draw_after(pic,'K20')
pic <- draw_after(pic,'K30')
pic <- draw_after(pic,'K40')
pic <- draw_after(pic,'K50')
pic <- draw_before(pic,'N100')
K_after_pic <- pic +
  geom_hline(aes(yintercept=20),linetype=5,linewidth=1,col="#C0C0C0")+
  scale_shape_manual(values = c(0,1,2,3,4,5,8))+
  scale_fill_manual(values=c(K10 = "#DF7A5E", K20 = "#F2CC8E", K30 = "#82B29A",K40="#D2BFA5",K50="#8e7cc3",N100="#8FB3DC"))+
  scale_color_manual(values=c(K10 = "#DF7A5E", K20 = "#F2CC8E", K30 = "#82B29A",K40="#D2BFA5",K50="#8e7cc3",N100="#8FB3DC"))+
  labs(x = "Population size", y = "Survival time to extinction")+
  scale_y_continuous(limits=c(0,100),
                     labels = seq(0,100,20),
                     breaks = seq(0,100,20),
  ) +
  xlim(2,8)+
  my_theme2()



# no rescue ---------------------------------------------------------------

# Nto4
slim_result <- read.table('time2_4_no_rescue.txt',header=T,sep='\t',fill = T)
colnames(slim_result)<-c("scenarios","num")
names<-c(
  'N10','N20','N30','N40','N50','N100'
)
slim_sub_result <- subset(slim_result, as.character(scenarios) %in% names)
slim_sub_result$scenarios <- factor(slim_sub_result$scenarios,
                                    levels = names,
                                    ordered = T)
scenario_labels <- c(
  "N10"  = "N=10",
  "N20"  = "N=20",
  "N30"  = "N=30",
  "N40"  = "N=40",
  "N50"  = "N=50",
  "N100"  = "N=100"
)
N_time2_4_pic <- 
  ggplot(slim_sub_result,aes(scenarios,num,fill=scenarios)) + 
  geom_hline(aes(yintercept=100),linetype=5,linewidth=1,col="#C0C0C0")+
  geom_hline(aes(yintercept=200),linetype=5,linewidth=1,col="#C0C0C0")+
  geom_quasirandom(
    aes(
      color = scenarios
    ),
    shape = 16,
    width = 0.4,
    size = 3,
    alpha = 0.8
  ) +
  geom_boxplot(
    fill = NA,
    color = "black",
    width = 0.6,
    outlier.shape = NA
  ) +
  geom_point(aes(x = scenarios, y = 500, color = scenarios),shape = 21, size = 10, col="black", alpha = 0.5) +  # 添加大圆点
  scale_fill_manual(values=c(N10 = "#DF7A5E" ,N20 = "#F2CC8E", N30 = "#82B29A", N40 = "#D2BFA5",N50="#8e7cc3",N100="#8FB3DC"))+
  scale_color_manual(values=c(N10 = "#DF7A5E" ,N20 = "#F2CC8E", N30 = "#82B29A", N40 = "#D2BFA5",N50="#8e7cc3",N100="#8FB3DC"))+
  add_significance_adjoin(
    df = slim_sub_result,
    group_col = "scenarios",
    y_col = 'num',                    
    group_levels = c('N10','N20','N30','N40','N50','N100'),
    ylim = c(0,700),
    bracket_start = 0.75,
    bracket_step = 0.05
  ) +
  scale_x_discrete(labels = scenario_labels) +
  scale_y_continuous(limits=c(0,700),
                     labels = seq(0,700,100),
                     breaks = seq(0,700,100),
                     expand = c(0,0))+
  labs(x = "Scenarios", y = "Time to minimum population size")+
  guides(fill="none")+
  my_theme2()
  

# Kto4
slim_result <- read.table('time2_4_no_rescue.txt',header=T,sep='\t',fill = T)
colnames(slim_result)<-c("scenarios","num")
names<-c(
  'K10','K20','K30','K40','K50','N100'
)
slim_sub_result <- subset(slim_result, as.character(scenarios) %in% names)
slim_sub_result$scenarios <- factor(slim_sub_result$scenarios,
                                    levels = names,
                                    ordered = T)
K_time2_4_pic <- 
  ggplot(slim_sub_result,aes(scenarios,num,fill=scenarios)) + 
  geom_quasirandom(aes(color=scenarios),width = 0.2,shape=21,size=4,col="black", alpha = 0.5)+
  geom_hline(aes(yintercept=100),linetype=5,linewidth=1,col="#C0C0C0")+
  geom_hline(aes(yintercept=200),linetype=5,linewidth=1,col="#C0C0C0")+
  geom_boxplot(aes(y=num),width=0.4,outlier.shape = NA)+
  geom_point(aes(x = scenarios, y = 500, color = scenarios),shape = 21, size = 10, col="black", alpha = 0.5) +  # 添加大圆点
  scale_fill_manual(values=c(K10 = "#DF7A5E" ,K20 = "#F2CC8E", K30 = "#82B29A", K40 = "#D2BFA5",K50="#8e7cc3",N100="#8FB3DC"))+
  scale_color_manual(values=c(K10 = "#DF7A5E" ,K20 = "#F2CC8E", K30 = "#82B29A", K40 = "#D2BFA5",K50="#8e7cc3",N100="#8FB3DC"))+
  geom_signif(
    aes(y = num, group = scenarios),
    comparisons = list(c("K10", "K20"),
                       c("K20", "K30"),
                       c("K30", "K40"),
                       c("K50", "N100")
    ),
    textsize = 3,
    y_position=c(550,600,550,600),
    map_signif_level=FALSE)+
  scale_y_continuous(limits=c(0,700),
                     labels = seq(0,700,100),
                     breaks = seq(0,700,100),
                     expand = c(0,0))+
  labs(x = "Scenarios", y = "Time to minimum population size")+
  guides(fill="none")+
  my_theme2()


# rescue ------------------------------------------------------------------

slim_result <- read.table(
  "time2_4_rescue2.txt",
  header = FALSE,
  sep = "\t",
  fill = TRUE
)
colnames(slim_result) <- c("scenarios", "num")
scenario_order <- c(
  "r5", "N10_r5",
  "r10", "N10_r10",
  "r15", "N10_r15",
  "r20", "N10_r20",
  "r30", "r40", "r50"
)
slim_sub_result <- slim_result %>%
  filter(scenarios %in% scenario_order) %>%
  mutate(
    release_num = as.numeric(sub(".*r", "", scenarios)),
    rescue_timing = ifelse(
      grepl("^N10", scenarios),
      "T=10",
      "T=4"
    ),
    rescue_timing = factor(
      rescue_timing,
      levels = c("T=4", "T=10")
    ),
    release_num = factor(
      release_num,
      levels = c(5, 10, 15, 20, 30, 40, 50)
    )
  )
scenario_labels <- c(
  "r5"      = "T=4\nF=5/1/1",
  "N10_r5"  = "T=10\nF=5/1/1",
  "r10"     = "T=4\nF=10/1/1",
  "N10_r10" = "T=10\nF=10/1/1",
  "r15"     = "T=4\nF=15/1/1",
  "N10_r15" = "T=10\nF=15/1/1",
  "r20"     = "T=4\nF=20/1/1",
  "N10_r20" = "T=10\nF=20/1/1",
  "r30"     = "T=4\nF=30/1/1",
  "r40"     = "T=4\nF=40/1/1",
  "r50"     = "T=4\nF=50/1/1"
)
rescue_num_pic <- 
  ggplot(slim_sub_result,aes(scenarios,num,fill=rescue_timing)) + 
  geom_hline(aes(yintercept=100),linetype=5,linewidth=1,col="#C0C0C0")+
  geom_hline(aes(yintercept=200),linetype=5,linewidth=1,col="#C0C0C0")+
  geom_quasirandom(
    aes(
      color = rescue_timing
    ),
    shape = 16,
    width = 0.4,
    size = 3,
    alpha = 0.8
  ) +
  geom_boxplot(
    fill = NA,
    color = "black",
    width = 0.6,
    outlier.shape = NA
  ) +
  geom_point(aes(x = scenarios, y = 500, color = rescue_timing),shape = 21, size = 10, col="black", alpha = 0.5) +  # 添加大圆点
  add_significance_adjoin(
    df = slim_sub_result,
    group_col = "scenarios",
    y_col = 'num',                   
    group_levels = c('r5','N10_r5','r10','N10_r10','r15','N10_r15','r20','N10_r20','r30','r40','r50'),
    ylim = c(0,700),
    bracket_start = 0.8,
    bracket_step = 0.05
  ) +
  scale_y_continuous(limits=c(0,700),
                     labels = seq(0,700,100),
                     breaks = seq(0,700,100),
                     expand = c(0,0))+
  scale_x_discrete(
    limits = scenario_order,
    labels = scenario_labels[scenario_order]
  )+
  scale_shape_manual(values = c(0,1,2,3,4))+
  scale_colour_manual(
    values = c(
      "T=4" = "#8FB3DC",
      "T=10" = "#F2CC8E"
    ),
    name = "Rescue initiated at"
  ) +
  scale_fill_manual(
    values = c(
      "T=4" = "#8FB3DC",
      "T=10" = "#F2CC8E"
    ),
    guide = "none"
  ) +  labs(x = "Scenarios", y = "Time to minimum population size\nafter rescue")+
  guides(fill="none")+
  my_theme2()
rescue_num_pic

# freq
slim_result <- read.table('time2_4_rescue2.txt',header=F,sep='\t',fill = T)
colnames(slim_result)<-c("scenarios","num")
names<-c(
  'r30','r30_3_10_15','r30_6_5_6','r30_6_5_4','r30_6_5_2'
)
slim_sub_result <- subset(slim_result, as.character(scenarios) %in% names)

slim_sub_result$scenarios <- factor(slim_sub_result$scenarios,
                                    levels = names,
                                    ordered = T)
freq_labels <- c(
  "r30"          = "F=30/1/1",
  "r30_3_10_15"  = "F=3/10/15",
  "r30_6_5_6"    = "F=6/5/6",
  "r30_6_5_4"    = "F=6/5/4",
  "r30_6_5_2"    = "F=6/5/2"
)
rescue_freq_pic <- 
  ggplot(slim_sub_result,aes(scenarios,num,fill=scenarios)) + 
  geom_hline(aes(yintercept=100),linetype=5,linewidth=1,col="#C0C0C0")+
  geom_hline(aes(yintercept=200),linetype=5,linewidth=1,col="#C0C0C0")+
  geom_quasirandom(
    aes(
      color = scenarios
    ),
    shape = 16,
    width = 0.4,
    size = 3,
    alpha = 0.8
  ) +
  geom_boxplot(
    fill = NA,
    color = "black",
    width = 0.6,
    outlier.shape = NA
  ) +
  geom_point(aes(x = scenarios, y = 500, fill = scenarios),shape = 21, size = 10, col="black", alpha = 0.5) +  # 添加大圆点
  scale_y_continuous(limits=c(0,700),
                     labels = seq(0,700,100),
                     breaks = seq(0,700,100),
                     expand = c(0,0))+
  scale_x_discrete(
    limits = names,
    labels = freq_labels[names]
  )+
  scale_fill_manual(values=c(r30 = "#82B29A" ,r30_3_10_15 = "#DF7A5E", r30_6_5_6 = "#F2CC8E", r30_6_5_4 = "#D2BFA5",r30_6_5_2="#8e7cc3"))+
  scale_color_manual(values=c(r30 = "#82B29A" ,r30_3_10_15 = "#DF7A5E", r30_6_5_6 = "#F2CC8E", r30_6_5_4 = "#D2BFA5",r30_6_5_2="#8e7cc3"))+
  labs(x = "Scenarios", y = "Time to minimum population size\nafter rescue")+
  guides(fill="none")+
  my_theme2()
rescue_freq_pic

# method 示意图 ---------------------------------------------------------------------

slim_het_del <- read.table("D://重测序项目//分析部分//slim//final//stat//het_del.txt",header=T)
names(slim_het_del) <- c("het","del")
het_del_plot <- 
  ggplot() + 
  geom_point(data = subset(slim_het_del,del >= 2.8974e+07*het-14.8329 & het >= mean(het)),aes(x=het,y=del),fill="#D2BFA5",size=3,shape=24,col="black",alpha = 0.8,stroke=0.5) +
  geom_point(data = subset(slim_het_del,del >= 2.8974e+07*het-14.8329 & het < mean(het)),aes(x=het,y=del),fill="#DF7A5E",size=3,shape=23,col="black",alpha = 0.8,stroke=0.5) +
  geom_point(data = subset(slim_het_del,del < 2.8974e+07*het-14.8329 & het < mean(het)),aes(x=het,y=del),fill="#8e7cc3",size=3,shape=22,col="black",alpha = 0.8,stroke=0.5) + 
  geom_point(data = subset(slim_het_del,del < 2.8974e+07*het-14.8329 & het >= mean(het)),aes(x=het,y=del),fill="#F2CC8E",size=3,shape=25,col="black",alpha = 0.8,stroke=0.5) +
  geom_vline(xintercept = mean(slim_het_del$het), linetype = "dashed", color = "grey",size=1.5) +
  labs(x = "Heterozygosity", y = "Het. deleterious variant\ngenotypes count")+
  geom_abline(slope =  2.8974e+07, intercept = -14.8329, size=1.5, color = "black") +
  my_theme2()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
  )+
  theme(legend.position=c(0.05, 0.95))
#het_del_plot

# method time -------------------------------------------------------------
# selcect
slim_result <- read.table('time2_4_rescue2.txt',header=F,sep='\t',fill = T)
colnames(slim_result)<-c("scenarios","num")
names<-c(
  'r30','r30_high_low','r30_high_high','r30_low_low','r30_low_high'
)
slim_sub_result <- subset(slim_result, as.character(scenarios) %in% names)

slim_sub_result$scenarios <- factor(slim_sub_result$scenarios,
                                    levels = names,
                                    ordered = T)
scenario_labels <- c(
  "r30" = "Random",
  "r30_high_high" = "HH",
  "r30_high_low" = "HL",
  "r30_low_high" = "LH",
  "r30_low_low" = "LL"
)
rescue_select_pic <- 
  ggplot(slim_sub_result,aes(scenarios,num,fill=scenarios)) + 
  geom_hline(aes(yintercept=100),linetype=5,linewidth=1,col="#C0C0C0")+
  geom_hline(aes(yintercept=200),linetype=5,linewidth=1,col="#C0C0C0")+
  #geom_quasirandom(aes(color=scenarios),width = 0.4,shape=21,size=4,col="black", alpha = 0.5)+
  geom_quasirandom(
    aes(
      color = scenarios,
      shape = scenarios
    ),
    width = 0.4,
    size = 2,
    col = "black",
    stroke=0.2,
    alpha = 0.8
  ) +
  geom_boxplot(
    fill = NA,
    color = "black",
    width = 0.6,
    outlier.shape = NA
  ) +
  geom_point(aes(x = scenarios, y = 500, color = scenarios),shape = 21, size = 10, col="black", alpha = 0.5) +  
  scale_shape_manual(values = c(r30 = 21 ,r30_high_low = 25, r30_low_high = 23,r30_high_high=24,r30_low_low=22))+
  add_significance(
    df = slim_sub_result,
    group_col = "scenarios",
    y_col = 'num',                     
    group_levels = c("r30", "r30_high_high", "r30_high_low", "r30_low_high", "r30_low_low"),
    ylim = c(0,700),
    bracket_start = 0.65,
    bracket_step = 0.05
  ) +
  scale_y_continuous(
    limits = c(0,700),
    labels = seq(0,700,100),
    breaks = seq(0,700,100),
    expand = c(0,0)
  )+
  scale_x_discrete(labels = scenario_labels) +
  scale_fill_manual(values=c(r30 = "#82B29A" ,r30_high_low = "#F2CC8E", r30_low_high = "#DF7A5E",r30_high_high="#D2BFA5",r30_low_low="#8e7cc3"))+
  scale_color_manual(values=c(r30 = "#82B29A" ,r30_high_low = "#F2CC8E", r30_low_high = "#DF7A5E",r30_high_high="#D2BFA5",r30_low_low="#8e7cc3"))+
  labs(x = "Scenarios", y = "Time to minimum population size\nafter rescue")+
  guides(fill="none")+
  my_theme2()
# rescue_select_pic


#  method het del fitness-----------------------------------------------------------


### selcet info
selcet_fun <- function(slim_out,arg_item) {
  pic <-
    ggplot(slim_out,aes(scenarios,{{ arg_item }},fill=scenarios)) + 
    geom_quasirandom(
      aes(
        color = scenarios,
        shape = scenarios
      ),
      width = 0.4,
      size = 2,
      col = "black",
      stroke=0.2,
      alpha = 0.8
    ) +
    geom_boxplot(
      aes(y={{ arg_item }}),
      fill = NA,
      color = "black",
      width = 0.6,
      outlier.shape = NA
    ) +
    scale_shape_manual(values = c(r30 = 21 ,r30_high_low = 25, r30_low_high = 23,r30_high_high=24,r30_low_low=22))+
    scale_fill_manual(values=c(r30 = "#82B29A" ,r30_high_low = "#F2CC8E", r30_low_high = "#DF7A5E",r30_high_high="#D2BFA5",r30_low_low="#8e7cc3"))+
    scale_color_manual(values=c(r30 = "#82B29A" ,r30_high_low = "#F2CC8E", r30_low_high = "#DF7A5E",r30_high_high="#D2BFA5",r30_low_low="#8e7cc3"))+
    guides(fill="none")+
    my_theme()
  return (pic)
}

scenario_labels <- c(
  "r30" = "Random",
  "r30_high_high" = "HH",
  "r30_high_low" = "HL",
  "r30_low_high" = "LH",
  "r30_low_low" = "LL"
)

plot_metric <- function(data, y_col, y_limits, y_label) {
  selcet_fun(data, !!sym(y_col)) +
    add_significance(
      df = data,
      group_col = "scenarios",
      y_col = y_col,
      group_levels = c(
        "r30",
        "r30_high_high",
        "r30_high_low",
        "r30_low_high",
        "r30_low_low"
      ),
      ylim = y_limits,
      bracket_start = 0.65,
      bracket_step = 0.05
    ) +
    scale_x_discrete(labels = scenario_labels) +
    scale_y_continuous(
      name = y_label,
      limits = y_limits,
      breaks = seq(y_limits[1], y_limits[2], length.out = 5),
      expand = c(0, 0)
    )
}
process_file <- function(file) {
  slim_result <- read.table(
    file,
    header = FALSE,
    sep = "\t",
    fill = TRUE
  )
  colnames(slim_result) <- c(
    "scenarios", "popSize", "meanFitness", "meanHet",
    "B_gen", "B_year", "FROH_1Mb", "avgStrDel",
    "avgModDel", "avgWkDel", "avgHomDel", "avgHetDel"
  )
  scenarios_keep <- c(
    "r30",
    "r30_high_low",
    "r30_high_high",
    "r30_low_low",
    "r30_low_high"
  )
  slim_sub <- subset(
    slim_result,
    as.character(scenarios) %in% scenarios_keep
  )
  slim_sub$scenarios <- factor(
    slim_sub$scenarios,
    levels = scenarios_keep,
    ordered = TRUE
  )
  # Convert heterozygosity from per bp to per kbp
  slim_sub$meanHet_kbp <- slim_sub$meanHet * 1000
  ranges <- list(
    meanHet_kbp = c(3.8e-2, 5.0e-2),
    avgHetDel = c(500, 1000),
    meanFitness = c(0.95, 1)
  )
  list(
    het = plot_metric(
      slim_sub,
      "meanHet_kbp",
      ranges$meanHet_kbp,
      "Heterozygosity per kbp"
    ),
    hetdel = plot_metric(
      slim_sub,
      "avgHetDel",
      ranges$avgHetDel,
      "Het. deleterious variant\ngenotypes count"
    ),
    fitness = plot_metric(
      slim_sub,
      "meanFitness",
      ranges$meanFitness,
      "Fitness"
    )
  )
}

plots_1g <- process_file("het_del_after_rescue_1ge.txt")
plots_4g <- process_file("het_del_after_rescue_50.txt")

top_row <- (
  plots_1g$het |
    plots_1g$hetdel |
    plots_1g$fitness
) & top_theme()
bottom_row <- (
  plots_4g$het |
    plots_4g$hetdel |
    plots_4g$fitness
)
# top_row / bottom_row

plot_compare_method <- ((het_del_plot / rescue_select_pic) | (top_row / bottom_row)) + plot_layout(widths = c(1, 3))

# sex
slim_result <- read.table('time2_4_rescue2.txt',header=F,sep='\t',fill = T)
colnames(slim_result)<-c("scenarios","num")
names<-c(
 'r30','r30_randomsex'
)
slim_sub_result <- subset(slim_result, as.character(scenarios) %in% names)

slim_sub_result$scenarios <- factor(slim_sub_result$scenarios,
                                    levels = names,
                                    ordered = T)
# all ----------------------------------------------------------------------
null_plot <- ggplot()
slim_1 <- null_plot|N_time2_4_pic|(N_after_pic/K_after_pic)
slim_2 <- (rescue_num_pic|rescue_freq_pic)+ plot_layout(widths = c(2, 1))

slim_main <- (slim_1/slim_2/plot_compare_method)+plot_layout(heights = c(1, 1, 2))
ggsave(slim_main, file='Plot_slim.pdf',device = cairo_pdf, width=16, height=20)

