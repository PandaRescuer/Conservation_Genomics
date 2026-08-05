library(ggplot2)
library(ggbeeswarm)
library(ggpmisc)
library(ggsignif)
library(patchwork)
library(dplyr)
library(ggrepel)
windowsFonts(Times = windowsFont("Times New Roman"))

my_theme <- function() {
  theme_classic()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5))+ 
    theme(axis.title.x = element_text(size = 15, colour = "black"),
          axis.title.y = element_text(size = 15, colour = "black"),
          axis.text.x = element_text(size = 15,colour = "black",angle = 45,hjust=1),
          axis.text.y = element_text(size = 15, colour = "black"),
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
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5))+ 
    theme(axis.title.x = element_text(size = 15, colour = "black"),
          axis.title.y = element_text(size = 15, colour = "black"),
          axis.text.x = element_text(size = 15,colour = "black"),
          axis.text.y = element_text(size = 15, colour = "black"),
          axis.ticks.length = unit(0.2, "cm"),  
          axis.ticks = element_line(size = 1),  
          axis.line = element_line(size = 1),  
          legend.position="none",
          legend.text = element_text(size=15),
          text = element_text(family = "Times")
    )
}

add_significance <- function(
    df,
    group_col,
    y_col,
    comparisons = NULL,
    group_levels = NULL,
    p_cutoff = 0.05,
    y_position = NULL,
    step_increase = 0.2,
    textsize = 3,
    min_n = 2,
    p_digits = 3
) {
  df_test <- df %>%
    dplyr::filter(
      !is.na(.data[[group_col]]),
      !is.na(.data[[y_col]])
    )
  
  # Use specified levels, factor levels, or data order
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
  
  # Compare all groups by default
  if (is.null(comparisons)) {
    if (length(group_levels) < 2) {
      return(NULL)
    }
    
    comparisons <- combn(group_levels, 2, simplify = FALSE)
  }
  
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
  
  if (is.null(y_position)) {
    y_range <- range(df_test[[y_col]], na.rm = TRUE)
    y_span <- diff(y_range)
    
    if (!is.finite(y_span) || y_span == 0) {
      y_span <- max(abs(y_range), na.rm = TRUE)
    }
    
    if (!is.finite(y_span) || y_span == 0) {
      y_span <- 1
    }
    
    y_position <- max(y_range) +
      seq_along(significant_comparisons) * step_increase * y_span
  } else if (length(y_position) != length(significant_comparisons)) {
    stop(
      "`y_position` must have the same length as the number of ",
      "significant comparisons: ",
      length(significant_comparisons)
    )
  }
  
  p_labels <- format(
    significant_pvalues,
    scientific = FALSE,
    digits = p_digits
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


# DEL 1 absent -------------------------------------------------------------------


new_del_LOF_DEL <- read.table("new_del_LOF_DEL.txt",header=T)

new_del_LOF_DEL2 <- new_del_LOF_DEL %>%
  mutate(
    # Set group order on x-axis
    rescue_pop = factor(
      rescue_pop,
      levels = c("MIX_CAP", "SC_CAP", "SC_QLA", "SC_MSH")
    ),
    border_col = ifelse(
      rescue_pop == "MIX_CAP" & Pedigree_generation == "F1",
      "#8FB3DC",
      as.character(rescue_pop)
    ),
    border_size = ifelse(
      rescue_pop == "MIX_CAP" & Pedigree_generation == "F1",
      1.5,
      0.2
    )
  )

het_DEL_absent <- 
  ggplot(new_del_LOF_DEL2, aes(rescue_pop, het_DEL, fill = rescue_pop)) +
  geom_quasirandom(
    aes(
      color = border_col,
      stroke = border_size
    ),
    width = 0.4,
    shape = 21,
    size = 2.5,
    alpha = 0.8
  ) +
  geom_boxplot(
    aes(y = het_DEL),
    fill = NA,
    color = "black",
    width = 0.6,
    outlier.shape = NA
  ) +
  ylab("Het. DEL genotypes absent from SC_WSP") +
  xlab(NULL) +
  guides(fill = "none") +
  scale_y_continuous(
    limits = c(800, 1300),
    labels = seq(800, 1300, 100),
    breaks = seq(800, 1300, 100),
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    limits = c("MIX_CAP", "SC_CAP", "SC_QLA", "SC_MSH")
  ) +
  scale_color_manual(
    values = c(
      MIX_CAP = "#D2BFA5",
      SC_CAP = "#F2CC8E",
      SC_QLA = "#82B29A",
      SC_MSH = "#DF7A5E",
      "#8FB3DC" = "#8FB3DC"
    )
  ) +
  scale_fill_manual(
    values = c(
      MIX_CAP = "#D2BFA5",
      SC_CAP = "#F2CC8E",
      SC_QLA = "#82B29A",
      SC_MSH = "#DF7A5E"
    )
  ) +
  add_significance(
    df = new_del_LOF_DEL2,
    group_col = "rescue_pop",
    y_col = "het_DEL",
    group_levels = c("MIX_CAP", "SC_CAP", "SC_QLA", "SC_MSH"),
    p_cutoff = 0.05,
    step_increase = 0.15
  )

hom_DEL_absent <-
  ggplot(new_del_LOF_DEL2, aes(rescue_pop, hom_DEL, fill = rescue_pop)) + 
  geom_quasirandom(
    aes(
      color = border_col,
      stroke = border_size
    ),
    width = 0.4,
    shape = 21,
    size = 2.5,
    alpha = 0.8
  ) +
  geom_boxplot(
    aes(y = hom_DEL),
    fill = NA,
    color = "black",
    width = 0.6,
    outlier.shape = NA
  ) +
  ylab("Hom. DEL genotypes absent from SC_WSP") +
  xlab(NULL) +
  guides(
    fill = "none",
    color = "none",
    stroke = "none"
  ) +
  scale_x_discrete(
    limits = c("MIX_CAP", "SC_CAP", "SC_QLA", "SC_MSH")
  ) +
  scale_y_continuous(
    limits = c(100, 400),
    labels = seq(100, 400, 100),
    breaks = seq(100, 400, 100),
    expand = c(0, 0)
  ) +
  scale_color_manual(
    values = c(
      MIX_CAP = "#D2BFA5",
      SC_CAP = "#F2CC8E",
      SC_QLA = "#82B29A",
      SC_MSH = "#DF7A5E",
      "#8FB3DC" = "#8FB3DC"
    )
  ) +
  scale_fill_manual(
    values = c(
      MIX_CAP = "#D2BFA5",
      SC_CAP = "#F2CC8E",
      SC_QLA = "#82B29A",
      SC_MSH = "#DF7A5E"
    )
  )+ 
  add_significance(
    df = new_del_LOF_DEL2,
    group_col = "rescue_pop",
    y_col = "hom_DEL",
    group_levels = c("MIX_CAP", "SC_CAP", "SC_QLA", "SC_MSH"),
    p_cutoff = 0.05,
    step_increase = 0.15
  )

het_DEL_absent <- het_DEL_absent + my_theme()
hom_DEL_absent <- hom_DEL_absent + my_theme()

# DEL 2 method 1/2-------------------------------------------------------------------


new_del_plot_ref_strict <- function(x_item,mut_type,xmin,xmax,ymin,ymax){
  pic <- 
    ggplot() + 
    geom_rect(
      aes(
        xmin = {{ xmin }},    
        xmax = {{ xmax }},    
        ymin = {{ ymin }}, 
        ymax = {{ ymax }}
      ),
      fill = "#82B29A",  
      alpha = 0.1          
    ) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_MSH'),aes(x={{ x_item }},y={{ mut_type }}),fill="#DF7A5E",size=4,shape=24,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),aes(x={{ x_item }},y={{ mut_type }}),fill="#82B29A",size=4,shape=25,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP' & Pedigree_generation == 'F1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#D2BFA5',size=4,shape=23,col="#8FB3DC",alpha = 0.8,stroke= 1.5) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP'& Pedigree_generation != 'F1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#D2BFA5',size=4,shape=23,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_1!='1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_1=='1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="#8e7cc3",alpha = 0.8, stroke=1.5) +   
    #geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good=='2'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="#9acfea",alpha = 0.8,stroke= 1.5) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                color = "black",
                alpha = 0.3,
                level = 0.95,
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                color = "#82B29A",
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                se=FALSE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                color = "black",
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                color = "#F2CC8E",
                se=FALSE) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
      # aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..eq.label.., ..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      hjust = -0.2
    ) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
      # aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..eq.label.., ..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      vjust = 2,  
      hjust = -0.2
    ) +
    guides(fill="none")
  return (pic)
}

new_del_plot_ref_loose <- function(x_item,mut_type){
  pic <- 
    ggplot() + 
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_2=='0'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_2!='0'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="#00796B",alpha = 0.8, stroke=1.5) +   
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                color = "black",
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                color = "#F2CC8E",
                se=FALSE) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
      # aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..eq.label.., ..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      vjust = 2,  
      hjust = -0.2
    ) +
    guides(fill="none")
  return (pic)
}

# use func

het_new_het_Del_ref_loose <- 
  new_del_plot_ref_loose(het_whole_genome,het_DEL)+
  ylim(750,1200)+
  geom_vline(xintercept = 1.21, linetype = "dashed", color = "grey",size=1)+  
  labs(x = "Heterozygosity per kbp", y = "Het. DEL genotypes absent from SC_WSP")+
  my_theme2()+  
  scale_x_continuous(
    limits = c(1.0, 1.4),
    labels = seq(1.0, 1.4, 0.1),
    breaks = seq(1.0, 1.4, 0.1),
    expand = c(0, 0)
  ) 

het_new_het_Del_ref_strict <- new_del_plot_ref_strict(het_whole_genome,het_DEL,1.15,Inf,-Inf,928.25)+
  ylim(750,1200)+
  geom_hline(yintercept = 928.25, linetype = "dashed", color = "grey",size=1)+
  geom_vline(xintercept = 1.15, linetype = "dashed", color = "grey",size=1)+ 
  labs(x = "Heterozygosity per kbp", y = "Het. DEL genotypes absent from SC_WSP")+
  my_theme2()+  
  scale_x_continuous(
    limits = c(1.0, 1.4),
    labels = seq(1.0, 1.4, 0.1),
    breaks = seq(1.0, 1.4, 0.1),
    expand = c(0, 0)
  ) 

# DEL 3 method 3/4-------------------------------------------------------------------

CAP_new_del_plot_DEL <- 
  ggplot() +
  geom_rect(
    aes(
      xmin = 1.15,    
      xmax = Inf,   
      ymin = -Inf, 
      ymax = 1569.375
    ),
    fill = "#82B29A",  
    alpha = 0.2          
  ) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_1=='0'),aes(x=het_whole_genome,y=snpMisDelNum_het,fill=het_DEL),size=4,shape=21,col="black",alpha = 0.8) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_1=='1'),aes(x=het_whole_genome,y=snpMisDelNum_het,fill=het_DEL),size=4,shape=21,col="#8e7cc3",alpha = 0.8,stroke= 1.5) +   
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_2=='1' & good_DEL_method_1!='1'),aes(x=het_whole_genome,y=snpMisDelNum_het,fill=het_DEL),size=4,shape=21,col="#00796B",alpha = 0.8,stroke= 1.5) +   
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP' & Pedigree_generation == 'F1'),aes(x=het_whole_genome,y=snpMisDelNum_het,fill=het_DEL),size=4,shape=23,col="#8FB3DC",alpha = 0.8,stroke= 1.5) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP'& Pedigree_generation != 'F1'),aes(x=het_whole_genome,y=snpMisDelNum_het,fill=het_DEL),size=4,shape=23,col="black",alpha = 0.8) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),aes(x=het_whole_genome,y=snpMisDelNum_het,fill=het_DEL),size=4,shape=25,col="black",alpha = 0.8) +   
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_MSH'),aes(x=het_whole_genome,y=snpMisDelNum_het,fill=het_DEL),size=4,shape=24,col="black",alpha = 0.8) +  
  geom_vline(xintercept = 1.15, linetype = "dashed", color = "grey",size=1) +
  geom_hline(yintercept = 1569.375, linetype = "dashed", color = "grey",size=1)+
  ylim(1300,2000)+
  labs(x = "Heterozygosity per kbp", y = "Het. DEL genotypes count")+
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
              aes(x=het_whole_genome,y=snpMisDelNum_het),
              method = "lm",
              alpha = 0.3,
              level = 0.95,
              color = "black",
              se=TRUE) +
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
              aes(x=het_whole_genome,y=snpMisDelNum_het),
              method = "lm",
              color = "#82B29A",
              alpha = 0.3,
              level = 0.95,
              linetype= 'dashed',
              se=FALSE) +
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
              aes(x=het_whole_genome,y=snpMisDelNum_het),
              method = "lm", 
              alpha = 0.3,
              level = 0.95,
              color = "black",
              se=TRUE) +
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
              aes(x=het_whole_genome,y=snpMisDelNum_het),
              method = "lm", 
              color = "#F2CC8E",
              alpha = 0.3,
              level = 0.95,
              linetype= 'dashed',
              se=FALSE) +
  stat_poly_eq(
    data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
    aes(x = het_whole_genome, y = snpMisDelNum_het,  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    family = "Times",
    size = 4,
    vjust = 2,  
    hjust = -0.2
  )+ 
  stat_poly_eq(
    data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
    aes(x = het_whole_genome, y = snpMisDelNum_het,  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    family = "Times",
    size = 4,
    hjust = -0.2
  )+
  scale_fill_gradient2(
    low      = "#82B29A",           
    mid      = "white",         
    high     = "firebrick",      
    midpoint = 928.25  
  ) +
  my_theme2()+
  theme(legend.position=c(0.8, 0.2))+  
  scale_x_continuous(
    limits = c(1.0, 1.4),
    labels = seq(1.0, 1.4, 0.1),
    breaks = seq(1.0, 1.4, 0.1),
    expand = c(0, 0)
  ) 



CAP_new_del_plot_2_DEL <- 
  ggplot() + 
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_1=='0'),aes(x=het_whole_genome,y=snpMisDelNum_het,fill=het_DEL),size=4,shape=21,col="black",alpha = 0.8) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_1=='1'),aes(x=het_whole_genome,y=snpMisDelNum_het,fill=het_DEL),size=4,shape=21,col="#8e7cc3",alpha = 0.8,stroke= 1.5) +   
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_2=='1' & good_DEL_method_1!='1'),aes(x=het_whole_genome,y=snpMisDelNum_het,fill=het_DEL),size=4,shape=21,col="#00796B",alpha = 0.8,stroke= 1.5) +   
  geom_vline(xintercept = 1.21, linetype = "dashed", color = "grey",size=1) +
  ylim(1300,2000)+
  labs(x = "Heterozygosity per kbp", y = "Het. DEL genotypes count")+
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
              aes(x=het_whole_genome,y=snpMisDelNum_het),
              method = "lm", 
              alpha = 0.3,
              level = 0.95,
              color = "black",
              se=TRUE) +
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
              aes(x=het_whole_genome,y=snpMisDelNum_het),
              method = "lm", 
              color = "#F2CC8E",
              alpha = 0.3,
              level = 0.95,
              linetype= 'dashed',
              se=FALSE) +
  stat_poly_eq(
    data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
    aes(x = het_whole_genome, y = snpMisDelNum_het,  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    family = "Times",
    size = 4,
    vjust = 2,  
    hjust = -0.2
  )+ 
  scale_fill_gradient2(
    low      = "#82B29A",          
    mid      = "white",        
    high     = "firebrick",      
    midpoint = 960.37  
  ) +
  my_theme2()+
  theme(legend.position=c(0.85, 0.2))+  
  scale_x_continuous(
    limits = c(1.0, 1.4),
    labels = seq(1.0, 1.4, 0.1),
    breaks = seq(1.0, 1.4, 0.1),
    expand = c(0, 0)
  ) 


# DEL 4 compare method-------------------------------------------------------------------
new_del_LOF_DEL <- new_del_LOF_DEL %>%
  mutate(
    method1_ref_group = case_when(
      rescue_pop != "SC_CAP" ~ rescue_pop,
      het_whole_genome >= 1.15 & het_DEL >= 928.25 ~ "HH",
      het_whole_genome >= 1.15 & het_DEL < 928.25 ~ "HL",
      het_whole_genome < 1.15 & het_DEL >= 928.25 ~ "LH",
      TRUE ~ "LL"
    ),
    method1_group = case_when(
      rescue_pop != "SC_CAP" ~ rescue_pop,
      het_whole_genome >= 1.21 & (het_whole_genome * 761 + 37.9 - het_DEL) <= 0 ~ "HH",
      het_whole_genome >= 1.21 & (het_whole_genome * 761 + 37.9 - het_DEL) > 0 ~ "HL",
      het_whole_genome < 1.21 & (het_whole_genome * 761 + 37.9 - het_DEL) <= 0 ~ "LH",
      TRUE ~ "LL"
    ),
    method2_group = case_when(
      rescue_pop != "SC_CAP" ~ rescue_pop,
      het_whole_genome >= 1.15 & snpMisDelNum_het >= 1569.375 ~ "HH",
      het_whole_genome >= 1.15 & snpMisDelNum_het < 1569.375 ~ "HL",
      het_whole_genome < 1.15 & snpMisDelNum_het >= 1569.375 ~ "LH",
      TRUE ~ "LL"
    ),
    group_misdel = 1210 * het_whole_genome + 146 - snpMisDelNum_het,
    misdel_flag = if_else(group_misdel > 0, "lower", "higher"),
    method3_group = case_when(
      rescue_pop != "SC_CAP" ~ rescue_pop,
      misdel_flag == "higher" & het_whole_genome >= 1.21 ~ "HH",
      misdel_flag == "higher" & het_whole_genome < 1.21 ~ "LH",
      misdel_flag == "lower" & het_whole_genome >= 1.21 ~ "HL",
      TRUE ~ "LL"
    )
  )

library(ggbeeswarm)

plot_group <- function(df, group_col, y_col, y_lab, y_lim) {
  p_cutoff = 0.05
  df_sub <- df %>%
    filter(rescue_pop %in% c("SC_CAP", "SC_QLA")) %>%
    mutate(
      group = factor(.data[[group_col]], levels = c(
        "LL",
        "LH",
        "HL",
        "HH",
        "SC_QLA"
      ))
    )
  group_levels <- levels(df_sub$group)
  comparisons <- combn(group_levels, 2, simplify = FALSE)
  pvals <- sapply(comparisons, function(comp) {
    x <- df_sub %>% filter(group == comp[1]) %>% pull(.data[[y_col]])
    y <- df_sub %>% filter(group == comp[2]) %>% pull(.data[[y_col]])
    if (length(x) < 2 || length(y) < 2) return(NA_real_)
    suppressWarnings(wilcox.test(x, y)$p.value)
  })
  sig_idx <- which(!is.na(pvals) & pvals < p_cutoff)
  p <- ggplot(
    data = df_sub,
    aes(x = .data[[group_col]], y = .data[[y_col]])
  ) +
    geom_quasirandom(
      aes( fill = .data[[group_col]],colour = .data[[group_col]]),
      width = 0.3,
      shape = 16,
      size = 2.5,
      alpha = 0.8
    ) +
    geom_boxplot(
      fill = NA,
      color = "black",
      width = 0.6,
      outlier.shape = NA
    ) +
    scale_y_continuous(
      limits = c(y_lim[1], y_lim[2]),
      breaks = seq(y_lim[1], y_lim[2], y_lim[3]),
      labels = seq(y_lim[1], y_lim[2], y_lim[3]),
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = y_lab) +
    scale_fill_manual(values=c(HL = "#8e7cc3" ,LH = "#F2CC8E", HH = "#DF7A5E",LL="#D2BFA5",SC_QLA = "#82B29A"))+
    scale_color_manual(values=c(HL = "#8e7cc3" ,LH = "#F2CC8E", HH = "#DF7A5E",LL="#D2BFA5",SC_QLA = "#82B29A"))+
    my_theme()+
    theme(panel.border = element_blank())
  
  if (length(sig_idx) > 0) {
    y_max <- max(df_sub[[y_col]], na.rm = TRUE)
    y_rng <- diff(range(df_sub[[y_col]], na.rm = TRUE))
    if (y_rng == 0) y_rng <- abs(y_max)
    if (y_rng == 0) y_rng <- 1
    p <- p + geom_signif(
      comparisons = comparisons[sig_idx],
      y_position = y_max + seq_along(sig_idx) * 0.2 * y_rng,
      map_signif_level = FALSE,
      textsize = 3,
      step_increase = 0
    )
  }
  p
}

top_theme <- function() {
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank())
}

# 3) 生成图
method1_ref_het_DEL_absent_count <- plot_group(
  new_del_LOF_DEL, "method1_ref_group", "het_DEL",
  "Recipient-absent Het. DEL genotype count", c(800, 1400,200)
)
method1_het_DEL_absent_count <- plot_group(
  new_del_LOF_DEL, "method1_group", "het_DEL",
  "Recipient-absent Het. DEL genotype count", c(800, 1400,200)
)
method2_het_DEL_absent_count <- plot_group(
  new_del_LOF_DEL, "method2_group", "het_DEL",
  "Recipient-absent Het. DEL genotype count", c(800, 1400,200)
)+top_theme()
method3_het_DEL_absent_count <- plot_group(
  new_del_LOF_DEL, "method3_group", "het_DEL",
  "Recipient-absent Het. DEL genotype count", c(800, 1400,200)
)+top_theme()
method1_ref_hetkbp <- plot_group(
  new_del_LOF_DEL, "method1_ref_group", "het_whole_genome",
  "Heterozygosity per kbp", c(1.0, 1.8, 0.2)
)
method1_hetkbp <- plot_group(
  new_del_LOF_DEL, "method1_group", "het_whole_genome",
  "Heterozygosity per kbp", c(1.0, 1.8, 0.2)
)
method2_hetkbp <- plot_group(
  new_del_LOF_DEL, "method2_group", "het_whole_genome",
  "Heterozygosity per kbp", c(1.0, 1.8, 0.2)
)
method3_hetkbp <- plot_group(
  new_del_LOF_DEL, "method3_group", "het_whole_genome",
  "Heterozygosity per kbp", c(1.0, 1.8, 0.2)
)

# LOF 1 absent  -----------------------------------------------------------

het_LOF_absent <- 
  ggplot(new_del_LOF_DEL2, aes(rescue_pop, het_LOF, fill = rescue_pop)) +
  geom_quasirandom(
    aes(
      color = border_col,
      stroke = border_size
    ),
    width = 0.4,
    shape = 21,
    size = 2.5,
    alpha = 0.8
  ) +
  geom_boxplot(
    aes(y = het_LOF),
    fill = NA,
    color = "black",
    width = 0.6,
    outlier.shape = NA
  ) +
  ylab("Het. LOF genotypes absent from SC_WSP") +
  xlab(NULL) +
  guides(fill = "none") +
  scale_y_continuous(limits=c(50,130),
                     labels = seq(50,130,20),
                     breaks = seq(50,130,20),
                     expand = c(0,0)) +
  scale_x_discrete(
    limits = c("MIX_CAP", "SC_CAP", "SC_QLA", "SC_MSH")
  ) +
  scale_color_manual(
    values = c(
      MIX_CAP = "#D2BFA5",
      SC_CAP = "#F2CC8E",
      SC_QLA = "#82B29A",
      SC_MSH = "#DF7A5E",
      "#8FB3DC" = "#8FB3DC"
    )
  ) +
  scale_fill_manual(
    values = c(
      MIX_CAP = "#D2BFA5",
      SC_CAP = "#F2CC8E",
      SC_QLA = "#82B29A",
      SC_MSH = "#DF7A5E"
    )
  ) + 
  add_significance(
    df = new_del_LOF_DEL2,
    group_col = "rescue_pop",
    y_col = "het_LOF",
    group_levels = c("MIX_CAP", "SC_CAP", "SC_QLA", "SC_MSH"),
    p_cutoff = 0.05,
    step_increase = 0.15
  )

hom_LOF_absent <-
  ggplot(new_del_LOF_DEL2, aes(rescue_pop, hom_LOF, fill = rescue_pop)) + 
  geom_quasirandom(
    aes(
      color = border_col,
      stroke = border_size
    ),
    width = 0.4,
    shape = 21,
    size = 2.5,
    alpha = 0.8
  ) +
  geom_boxplot(
    aes(y = hom_LOF),
    fill = NA,
    color = "black",
    width = 0.6,
    outlier.shape = NA
  ) +
  scale_y_continuous(limits=c(0,50),
                     labels = seq(0,50,10),
                     breaks = seq(0,50,10),
                     expand = c(0,0)) +
  ylab("Hom. LOF genotypes absent from SC_WSP") +
  xlab(NULL) +
  guides(
    fill = "none",
    color = "none",
    stroke = "none"
  ) +
  scale_x_discrete(
    limits = c("MIX_CAP", "SC_CAP", "SC_QLA", "SC_MSH")
  ) +
  scale_color_manual(
    values = c(
      MIX_CAP = "#D2BFA5",
      SC_CAP = "#F2CC8E",
      SC_QLA = "#82B29A",
      SC_MSH = "#DF7A5E",
      "#8FB3DC" = "#8FB3DC"
    )
  ) +
  scale_fill_manual(
    values = c(
      MIX_CAP = "#D2BFA5",
      SC_CAP = "#F2CC8E",
      SC_QLA = "#82B29A",
      SC_MSH = "#DF7A5E"
    )
  ) + 
  add_significance(
    df = new_del_LOF_DEL2,
    group_col = "rescue_pop",
    y_col = "hom_LOF",
    group_levels = c("MIX_CAP", "SC_CAP", "SC_QLA", "SC_MSH"),
    p_cutoff = 0.05,
    step_increase = 0.15
  )

het_LOF_absent <- het_LOF_absent + my_theme()
hom_LOF_absent <- hom_LOF_absent + my_theme()


# LOF 2 method 1/2 --------------------------------------------------------
new_del_plot_LOF_compare_DEL <- function(x_item,mut_type,xmin,xmax,ymin,ymax){
  pic <- 
    ggplot() + 
    geom_rect(
      aes(
        xmin = {{ xmin }},    
        xmax = {{ xmax }},    
        ymin = {{ ymin }}, 
        ymax = {{ ymax }}
      ),
      fill = "#82B29A",  
      alpha = 0.1          
    ) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_MSH'),aes(x={{ x_item }},y={{ mut_type }}),fill="#DF7A5E",size=4,shape=24,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),aes(x={{ x_item }},y={{ mut_type }}),fill="#82B29A",size=4,shape=25,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP' & Pedigree_generation == 'F1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#D2BFA5',size=4,shape=23,col="#8FB3DC",alpha = 0.8,stroke= 1.5) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP'& Pedigree_generation != 'F1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#D2BFA5',size=4,shape=23,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_1=='0'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_1=='1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="#8e7cc3",alpha = 0.8, stroke=1.5) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                color = "black",
                alpha = 0.3,
                level = 0.95,
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                color = "#82B29A",
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                se=FALSE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                color = "black",
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                color = "#F2CC8E",
                se=FALSE) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
      # aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..eq.label.., ..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      hjust = -0.2
    ) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
      # aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..eq.label.., ..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      vjust = 2,  
      hjust = -0.2
    ) +
    guides(fill="none")
  return (pic)
}

new_del_plot_LOF_ref_loose_compare_DEL <- function(x_item,mut_type){
  pic <- 
    ggplot() + 
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_2=='0'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_2!='0'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="#00796B",alpha = 0.8, stroke=1.5) +   
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                color = "black",
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                color = "#F2CC8E",
                se=FALSE) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
      # aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..eq.label.., ..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      vjust = 2,  
      hjust = -0.2
    ) +
    guides(fill="none")
  return (pic)
}

DEL_LOF_het_new_het_LOF <- 
  new_del_plot_LOF_compare_DEL(het_whole_genome,het_LOF,1.15,Inf,-Inf,72.5)+
  ylim(40,120)+
  geom_hline(yintercept = 72.5, linetype = "dashed", color = "grey",size=1)+
  geom_vline(xintercept = 1.15, linetype = "dashed", color = "grey",size=1)+
  labs(x = "Heterozygosity per kbp", y = "Het. LOF genotypes absent from SC_WSP")+
  my_theme2()+  
  scale_x_continuous(
    limits = c(1.0, 1.4),
    labels = seq(1.0, 1.4, 0.1),
    breaks = seq(1.0, 1.4, 0.1),
    expand = c(0, 0)
  ) 

DEL_LOF_het_new_het_LOF_ref_loose <- 
  new_del_plot_LOF_ref_loose_compare_DEL(het_whole_genome,het_LOF)+
  ylim(40,120)+
  geom_vline(xintercept = 1.21, linetype = "dashed", color = "grey",size=1)+
  labs(x = "Heterozygosity per kbp", y = "Het. LOF genotypes absent from SC_WSP")+
  my_theme2()+  
  scale_x_continuous(
    limits = c(1.0, 1.4),
    labels = seq(1.0, 1.4, 0.1),
    breaks = seq(1.0, 1.4, 0.1),
    expand = c(0, 0)
  ) 


# LOF 3 method 3/4 --------------------------------------------------------


DEL_LOF_CAP_new_del_plot_LOF <- 
  ggplot() + 
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_1=='0'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=21,col="black",alpha = 0.8) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_1=='1'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=21,col="#8e7cc3",alpha = 0.8,stroke= 1.5) +   
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_2=='1' & good_DEL_method_1!='1'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=21,col="#00796B",alpha = 0.8,stroke= 1.5) +   
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP' & Pedigree_generation == 'F1'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=23,col="#8FB3DC",alpha = 0.8,stroke= 1.5) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP'& Pedigree_generation != 'F1'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=23,col="black",alpha = 0.8) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=25,col="black",alpha = 0.8) +   
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_MSH'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=24,col="black",alpha = 0.8) +  
  geom_vline(xintercept = 1.15, linetype = "dashed", color = "grey",size=1) +
  geom_hline(yintercept = 135.625, linetype = "dashed", color = "grey",size=1)+
  labs(x = "Heterozygosity", y = "Het. LOF mutations of each individual")+
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
              aes(x=het_whole_genome,y=snpLOFNum_het),
              method = "lm",
              alpha = 0.3,
              level = 0.95,
              color = "black",
              se=TRUE) +
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
              aes(x=het_whole_genome,y=snpLOFNum_het),
              method = "lm",
              color = "#82B29A",
              alpha = 0.3,
              level = 0.95,
              linetype= 'dashed',
              se=FALSE) +
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
              aes(x=het_whole_genome,y=snpLOFNum_het),
              method = "lm", 
              alpha = 0.3,
              level = 0.95,
              color = "black",
              se=TRUE) +
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
              aes(x=het_whole_genome,y=snpLOFNum_het),
              method = "lm", 
              color = "#F2CC8E",
              alpha = 0.3,
              level = 0.95,
              linetype= 'dashed',
              se=FALSE) +
  stat_poly_eq(
    data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
    aes(x = het_whole_genome, y = snpLOFNum_het,  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    family = "Times",
    size = 4,
    vjust = 2 
  )+ 
  stat_poly_eq(
    data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
    aes(x = het_whole_genome, y = snpLOFNum_het,  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    family = "Times",
    size = 4,
  )+
  scale_fill_gradient2(
    low      = "#82B29A",           
    mid      = "white",          
    high     = "firebrick",      
    midpoint = 72.5  
  ) +    
  my_theme()+
  theme(legend.position=c(0.8, 0.2))+  
  scale_x_continuous(
    limits = c(1.0, 1.4),
    labels = seq(1.0, 1.4, 0.1),
    breaks = seq(1.0, 1.4, 0.1),
    expand = c(0, 0)
  ) 


DEL_LOF_CAP_new_del_plot_2_LOF <- 
  ggplot() + 
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_1=='0'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=21,col="black",alpha = 0.8) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_1=='1'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=21,col="#8e7cc3",alpha = 0.8,stroke= 1.5) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_2=='1' & good_DEL_method_1 != '1'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=21,col="#00796B",alpha = 0.8,stroke= 1.5) +
  geom_vline(xintercept = 1.21, linetype = "dashed", color = "grey",size=1) +
  labs(x = "Heterozygosity", y = "Het. LOF mutations of each individual")+
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
              aes(x=het_whole_genome,y=snpLOFNum_het),
              method = "lm", 
              alpha = 0.3,
              level = 0.95,
              color = "black",
              se=TRUE) +
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
              aes(x=het_whole_genome,y=snpLOFNum_het),
              method = "lm", 
              color = "#F2CC8E",
              alpha = 0.3,
              level = 0.95,
              linetype= 'dashed',
              se=FALSE) +
  stat_poly_eq(
    data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
    aes(x = het_whole_genome, y = snpLOFNum_het,  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    family = "Times",
    size = 4,
    vjust = 2  
  )+ 
  scale_fill_gradient2(
    low      = "#82B29A",           
    mid      = "white",          
    high     = "firebrick",      
    midpoint = 76.033  
  ) +

  my_theme()+
  theme(legend.position=c(0.8, 0.2))+  
  scale_x_continuous(
    limits = c(1.0, 1.4),
    labels = seq(1.0, 1.4, 0.1),
    breaks = seq(1.0, 1.4, 0.1),
    expand = c(0, 0)
  ) 

(DEL_LOF_het_new_het_LOF|DEL_LOF_het_new_het_LOF_ref_loose)/(DEL_LOF_CAP_new_del_plot_LOF|DEL_LOF_CAP_new_del_plot_2_LOF)

# LOF self ----------------------------------------------------------------

new_del_LOF_DEL <- new_del_LOF_DEL %>%
  mutate(
    method1_LOF = case_when(
      rescue_pop != "SC_CAP" ~ rescue_pop,
      het_whole_genome >= 1.15 & het_LOF >= 72.5 ~ "HH",
      het_whole_genome >= 1.15 & het_LOF < 72.5 ~ "HL",
      het_whole_genome < 1.15 & het_LOF >= 72.5 ~ "LH",
      TRUE ~ "LL"
    ),
    method2_LOF = case_when(
      rescue_pop != "SC_CAP" ~ rescue_pop,
      het_whole_genome >= 1.21 & (het_whole_genome * 80.2 - 21.1 - het_LOF) <= 0 ~ "HH",
      het_whole_genome >= 1.21 & (het_whole_genome * 80.2 - 21.1 - het_LOF) > 0 ~ "HL",
      het_whole_genome < 1.21 & (het_whole_genome * 80.2 - 21.1 - het_LOF) <= 0 ~ "LH",
      TRUE ~ "LL"
    ),
    method3_LOF = case_when(
      rescue_pop != "SC_CAP" ~ rescue_pop,
      het_whole_genome >= 1.15 & snpLOFNum_het >= 135.625 ~ "HH",
      het_whole_genome >= 1.15 & snpLOFNum_het < 135.625 ~ "HL",
      het_whole_genome < 1.15 & snpLOFNum_het >= 135.625 ~ "LH",
      TRUE ~ "LL"
    ),
    group_LOF = 119 * het_whole_genome -2.68 - snpLOFNum_het,
    LOF_flag = if_else(group_LOF > 0, "lower", "higher"),
    method4_LOF = case_when(
      rescue_pop != "SC_CAP" ~ rescue_pop,
      LOF_flag == "higher" & het_whole_genome >= 1.21 ~ "HH",
      LOF_flag == "higher" & het_whole_genome < 1.21 ~ "LH",
      LOF_flag == "lower" & het_whole_genome >= 1.21 ~ "HL",
      TRUE ~ "LL"
    )
  )
write.table(new_del_LOF_DEL, 'merge_method_DEL.txt', row.names = TRUE, sep = '\t', quote = FALSE)
library(ggbeeswarm)

plot_group <- function(df, group_col, y_col, y_lab, y_lim = NULL) {
  p_cutoff <- 0.05
  
  df_sub <- df %>%
    filter(rescue_pop %in% c("SC_CAP", "SC_QLA")) %>%
    mutate(
      group = factor(
        .data[[group_col]],
        levels = c("LL", "LH", "HL", "HH", "SC_QLA")
      )
    )
  
  group_levels <- levels(df_sub$group)
  comparisons <- combn(group_levels, 2, simplify = FALSE)
  
  pvals <- sapply(comparisons, function(comp) {
    x <- df_sub %>%
      filter(group == comp[1]) %>%
      pull(.data[[y_col]])
    
    y <- df_sub %>%
      filter(group == comp[2]) %>%
      pull(.data[[y_col]])
    
    if (length(x) < 2 || length(y) < 2) {
      return(NA_real_)
    }
    
    suppressWarnings(wilcox.test(x, y)$p.value)
  })
  
  sig_idx <- which(!is.na(pvals) & pvals < p_cutoff)
  
  p <- ggplot(
    data = df_sub,
    aes(x = .data[[group_col]], y = .data[[y_col]])
  ) +
    geom_quasirandom(
      aes(
        fill = .data[[group_col]],
        colour = .data[[group_col]]
      ),
      width = 0.3,
      shape = 16,
      size = 2.5,
      alpha = 0.8
    ) +
    geom_boxplot(
      fill = NA,
      color = "black",
      width = 0.6,
      outlier.shape = NA
    ) +
    labs(x = NULL, y = y_lab) +
    scale_fill_manual(
      values = c(
        HL = "#8e7cc3",
        LH = "#F2CC8E",
        HH = "#DF7A5E",
        LL = "#D2BFA5",
        SC_QLA = "#82B29A"
      )
    ) +
    scale_color_manual(
      values = c(
        HL = "#8e7cc3",
        LH = "#F2CC8E",
        HH = "#DF7A5E",
        LL = "#D2BFA5",
        SC_QLA = "#82B29A"
      )
    ) +
    my_theme() +
    theme(panel.border = element_blank())
  
  # Add y-axis settings only when y_lim is provided
  if (!is.null(y_lim)) {
    if (length(y_lim) != 3) {
      stop("y_lim must be NULL or a numeric vector: c(min, max, step)")
    }
    
    y_breaks <- seq(y_lim[1], y_lim[2], by = y_lim[3])
    
    p <- p +
      scale_y_continuous(
        limits = y_lim[1:2],
        breaks = y_breaks,
        labels = y_breaks,
        expand = c(0, 0)
      )
  }
  
  # Add significant comparisons
  if (length(sig_idx) > 0) {
    y_max <- max(df_sub[[y_col]], na.rm = TRUE)
    y_rng <- diff(range(df_sub[[y_col]], na.rm = TRUE))
    
    if (y_rng == 0) {
      y_rng <- abs(y_max)
    }
    
    if (y_rng == 0) {
      y_rng <- 1
    }
    
    p <- p +
      geom_signif(
        comparisons = comparisons[sig_idx],
        y_position = y_max + seq_along(sig_idx) * 0.2 * y_rng,
        map_signif_level = FALSE,
        textsize = 3,
        step_increase = 0
      )
  }
  
  return(p)
}

top_theme <- function() {
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank())
}


LOF_method1_ref_het_LOF_absent_count <- plot_group(
  new_del_LOF_DEL, "method1_LOF", "het_LOF",
  "Recipient-absent Het. LOF genotype count",c(50,130,40)
)
LOF_method1_het_LOF_absent_count <- plot_group(
  new_del_LOF_DEL, "method2_LOF", "het_LOF",
  "Recipient-absent Het. LOF genotype count", c(50,130,40)
)
LOF_method2_het_LOF_absent_count <- plot_group(
  new_del_LOF_DEL, "method3_LOF", "het_LOF",
  "Recipient-absent Het. LOF genotype count", c(50,130,40)
)+top_theme()
LOF_method3_het_LOF_absent_count <- plot_group(
  new_del_LOF_DEL, "method4_LOF", "het_LOF",
  "Recipient-absent Het. LOF genotype count", c(50,130,40)
)+top_theme()
LOF_method1_ref_hetkbp <- plot_group(
  new_del_LOF_DEL, "method1_LOF", "het_whole_genome",
  "Heterozygosity per kbp", c(1.0, 1.8, 0.2)
)
LOF_method1_hetkbp <- plot_group(
  new_del_LOF_DEL, "method2_LOF", "het_whole_genome",
  "Heterozygosity per kbp", c(1.0, 1.8, 0.2)
)
LOF_method2_hetkbp <- plot_group(
  new_del_LOF_DEL, "method3_LOF", "het_whole_genome",
  "Heterozygosity per kbp", c(1.0, 1.8, 0.2)
)
LOF_method3_hetkbp <- plot_group(
  new_del_LOF_DEL, "method4_LOF", "het_whole_genome",
  "Heterozygosity per kbp", c(1.0, 1.8, 0.2)
)



# LOF self method 1/2 --------------------------------------------------------
self_new_del_plot_LOF_compare_DEL <- function(x_item,mut_type,xmin,xmax,ymin,ymax){
  pic <- 
    ggplot() + 
    geom_rect(
      aes(
        xmin = {{ xmin }},    
        xmax = {{ xmax }},    
        ymin = {{ ymin }}, 
        ymax = {{ ymax }}
      ),
      fill = "#82B29A",  
      alpha = 0.1          
    ) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_MSH'),aes(x={{ x_item }},y={{ mut_type }}),fill="#DF7A5E",size=4,shape=24,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),aes(x={{ x_item }},y={{ mut_type }}),fill="#82B29A",size=4,shape=25,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP' & Pedigree_generation == 'F1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#D2BFA5',size=4,shape=23,col="#8FB3DC",alpha = 0.8,stroke= 1.5) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP'& Pedigree_generation != 'F1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#D2BFA5',size=4,shape=23,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method1_LOF!='HL'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method1_LOF=='HL'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="#8e7cc3",alpha = 0.8, stroke=1.5) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                color = "black",
                alpha = 0.3,
                level = 0.95,
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                color = "#82B29A",
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                se=FALSE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                color = "black",
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                color = "#F2CC8E",
                se=FALSE) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      hjust = -0.2
    ) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      vjust = 2,  
      hjust = -0.2
    ) +
    guides(fill="none")
  return (pic)
}

self_new_del_plot_LOF_ref_loose_compare_DEL <- function(x_item,mut_type){
  pic <- 
    ggplot() +     
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method2_LOF %in% c("HH", "LH")),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method2_LOF %in% c("HL", "LL")),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="#00796B",alpha = 0.8, stroke=1.5) +      
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                color = "black",
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                color = "#F2CC8E",
                se=FALSE) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
      # aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..eq.label.., ..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      vjust = 2,  
      hjust = -0.2
    ) +
    guides(fill="none")
  return (pic)
}

het_new_het_LOF <- 
  self_new_del_plot_LOF_compare_DEL(het_whole_genome,het_LOF,1.15,Inf,-Inf,72.5)+
  ylim(40,120)+
  geom_hline(yintercept = 72.5, linetype = "dashed", color = "grey",size=1)+
  geom_vline(xintercept = 1.15, linetype = "dashed", color = "grey",size=1)+
  labs(x = "Heterozygosity per kbp", y = "Het. LOF genotypes absent from SC_WSP")+
  my_theme2()+  
  scale_x_continuous(
    limits = c(1.0, 1.4),
    labels = seq(1.0, 1.4, 0.1),
    breaks = seq(1.0, 1.4, 0.1),
    expand = c(0, 0)
  ) 

het_new_het_LOF_ref_loose <- 
  self_new_del_plot_LOF_ref_loose_compare_DEL(het_whole_genome,het_LOF)+
  ylim(40,120)+
  geom_vline(xintercept = 1.21, linetype = "dashed", color = "grey",size=1)+
  labs(x = "Heterozygosity per kbp", y = "Het. LOF genotypes absent from SC_WSP")+
  my_theme2()+  
  scale_x_continuous(
    limits = c(1.0, 1.4),
    labels = seq(1.0, 1.4, 0.1),
    breaks = seq(1.0, 1.4, 0.1),
    expand = c(0, 0)
  ) 



# LOF self method 3/4 --------------------------------------------------------


CAP_new_del_plot_LOF <- 
  ggplot() + 
  geom_rect(
    aes(
      xmin = 1.15,    
      xmax = Inf,    
      ymin = -Inf, 
      ymax = 135.625
    ),
    fill = "#82B29A",  
    alpha = 0.2          
  ) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method1_LOF !='HL' & method2_LOF %in% c("HH", "LH")),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=21,col="black",alpha = 0.8) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method1_LOF=='HL'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=21,col="#8e7cc3",alpha = 0.8,stroke= 1.5) +   
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method2_LOF %in% c("HL", "LL") & method1_LOF!='HL'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=21,col="#00796B",alpha = 0.8,stroke= 1.5) +   
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP' & Pedigree_generation == 'F1'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=23,col="#8FB3DC",alpha = 0.8,stroke= 1.5) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP'& Pedigree_generation != 'F1'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=23,col="black",alpha = 0.8) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=25,col="black",alpha = 0.8) +   
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_MSH'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=24,col="black",alpha = 0.8) +  
  geom_vline(xintercept = 1.15, linetype = "dashed", color = "grey",size=1) +
  geom_hline(yintercept = 135.625, linetype = "dashed", color = "grey",size=1)+
  ylim(100,200)+
  labs(x = "Heterozygosity", y = "Het. LOF mutations of each individual")+
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
              aes(x=het_whole_genome,y=snpLOFNum_het),
              method = "lm",
              alpha = 0.3,
              level = 0.95,
              color = "black",
              se=TRUE) +
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
              aes(x=het_whole_genome,y=snpLOFNum_het),
              method = "lm",
              color = "#82B29A",
              alpha = 0.3,
              level = 0.95,
              linetype= 'dashed',
              se=FALSE) +
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
              aes(x=het_whole_genome,y=snpLOFNum_het),
              method = "lm", 
              alpha = 0.3,
              level = 0.95,
              color = "black",
              se=TRUE) +
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
              aes(x=het_whole_genome,y=snpLOFNum_het),
              method = "lm", 
              color = "#F2CC8E",
              alpha = 0.3,
              level = 0.95,
              linetype= 'dashed',
              se=FALSE) +
  stat_poly_eq(
    data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
    aes(x = het_whole_genome, y = snpLOFNum_het,  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    family = "Times",
    size = 4,
    vjust = 2,  
    hjust = -0.2
  )+ 
  stat_poly_eq(
    data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
    aes(x = het_whole_genome, y = snpLOFNum_het,  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    family = "Times",
    size = 4,
    hjust = -0.2
  )+
  scale_fill_gradient2(
    low      = "#82B29A",           
    mid      = "white",         
    high     = "firebrick",      
    midpoint = 72.5  
  ) +    
  my_theme2()+
  theme(legend.position=c(0.8, 0.2))+  
  scale_x_continuous(
    limits = c(1.0, 1.4),
    labels = seq(1.0, 1.4, 0.1),
    breaks = seq(1.0, 1.4, 0.1),
    expand = c(0, 0)
  ) 


CAP_new_del_plot_2_LOF <- 
  ggplot() + 
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method1_LOF !='HL' & method2_LOF %in% c("HH", "LH")),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=21,col="black",alpha = 0.8) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method1_LOF=='HL'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=21,col="#8e7cc3",alpha = 0.8,stroke= 1.5) +   
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method2_LOF %in% c("HL", "LL") & method1_LOF!='HL'),aes(x=het_whole_genome,y=snpLOFNum_het,fill=het_LOF),size=4,shape=21,col="#00796B",alpha = 0.8,stroke= 1.5) +   
  geom_vline(xintercept = 1.21, linetype = "dashed", color = "grey",size=1) +
  ylim(100,200)+
  labs(x = "Heterozygosity", y = "Het. LOF mutations of each individual")+
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
              aes(x=het_whole_genome,y=snpLOFNum_het),
              method = "lm", 
              alpha = 0.3,
              level = 0.95,
              color = "black",
              se=TRUE) +
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
              aes(x=het_whole_genome,y=snpLOFNum_het),
              method = "lm", 
              color = "#F2CC8E",
              alpha = 0.3,
              level = 0.95,
              linetype= 'dashed',
              se=FALSE) +
  stat_poly_eq(
    data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
    aes(x = het_whole_genome, y = snpLOFNum_het,  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    family = "Times",
    size = 4,
    vjust = 2,  
    hjust = -0.2
  )+ 
  scale_fill_gradient2(
    low      = "#82B29A",           
    mid      = "white",          
    high     = "firebrick",      
    midpoint = 76.033  
  ) +
  my_theme2()+
  theme(legend.position=c(0.8, 0.2))+  
  scale_x_continuous(
    limits = c(1.0, 1.4),
    labels = seq(1.0, 1.4, 0.1),
    breaks = seq(1.0, 1.4, 0.1),
    expand = c(0, 0)
  ) 

(het_new_het_LOF|het_new_het_LOF_ref_loose)/(CAP_new_del_plot_LOF|CAP_new_del_plot_2_LOF)

# LOF self compare method -------------------------------------------------

LOF_merge_ref_strict <- (
  het_new_het_LOF |
    stack_with_gap(
      LOF_method1_ref_het_LOF_absent_count,
      LOF_method1_ref_hetkbp
    )
) + plot_layout(widths = c(2, 1))
LOF_merge_ref_loose <- (
  het_new_het_LOF_ref_loose |
    stack_with_gap(
      LOF_method1_het_LOF_absent_count,
      LOF_method1_hetkbp
    )
) + plot_layout(widths = c(2, 1))
LOF_merge_wild_benchmark <- (
  CAP_new_del_plot_LOF |
    stack_with_gap(
      LOF_method2_het_LOF_absent_count,
      LOF_method2_hetkbp
    )
) + plot_layout(widths = c(2, 1))
LOF_merge_non_ref <- (
  CAP_new_del_plot_2_LOF |
    stack_with_gap(
      LOF_method3_het_LOF_absent_count,
      LOF_method3_hetkbp
    )
) + plot_layout(widths = c(2, 1))

LOF_merge_ref <- wrap_plots(
  LOF_merge_ref_strict,
  LOF_merge_ref_loose,
  LOF_merge_wild_benchmark,
  LOF_merge_non_ref,
  ncol = 2
)
LOF_merge_ref



# Recall DEL from LOF -----------------------------------------------------

new_del_plot_DEL_compare_LOF <- function(x_item,mut_type,xmin,xmax,ymin,ymax){
  pic <- 
    ggplot() + 
    geom_rect(
      aes(
        xmin = {{ xmin }},   
        xmax = {{ xmax }},    
        ymin = {{ ymin }}, 
        ymax = {{ ymax }}
      ),
      fill = "#82B29A",  
      alpha = 0.1         
    ) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_MSH'),aes(x={{ x_item }},y={{ mut_type }}),fill="#DF7A5E",size=4,shape=24,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),aes(x={{ x_item }},y={{ mut_type }}),fill="#82B29A",size=4,shape=25,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP' & Pedigree_generation == 'F1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#D2BFA5',size=4,shape=23,col="#8FB3DC",alpha = 0.8,stroke= 1.5) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP'& Pedigree_generation != 'F1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#D2BFA5',size=4,shape=23,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method1_LOF!='HL'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method1_LOF=='HL'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="#8e7cc3",alpha = 0.8, stroke=1.5) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                color = "black",
                alpha = 0.3,
                level = 0.95,
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                color = "#82B29A",
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                se=FALSE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                color = "black",
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                color = "#F2CC8E",
                se=FALSE) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..eq.label.., ..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      hjust = -0.2
    ) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..eq.label.., ..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      vjust = 2,  
      hjust = -0.2
    ) +
    guides(fill="none")
  return (pic)
}

new_del_plot_DEL_ref_loose_compare_LOF <- function(x_item,mut_type){
  pic <- 
    ggplot() +   
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method2_LOF %in% c("HH", "LH")),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method2_LOF %in% c("HL", "LL")),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="#00796B",alpha = 0.8, stroke=1.5) +   
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                color = "black",
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                color = "#F2CC8E",
                se=FALSE) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..eq.label.., ..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      vjust = 2,  
      hjust = -0.2
    ) +
    guides(fill="none")
  return (pic)
}

LOF_DEL_method1 <- 
  new_del_plot_DEL_compare_LOF(het_whole_genome,het_DEL,1.15,Inf,-Inf,928.25)+
  ylim(750,1200)+
  geom_hline(yintercept = 928.25, linetype = "dashed", color = "grey",size=1)+
  geom_vline(xintercept = 1.15, linetype = "dashed", color = "grey",size=1)+
  labs(x = "Heterozygosity per kbp", y = "Het. LOF genotypes absent from SC_WSP")+
  my_theme2()

LOF_DEL_method2 <- 
  new_del_plot_DEL_ref_loose_compare_LOF(het_whole_genome,het_DEL)+
  ylim(750,1200)+
  geom_vline(xintercept = 1.21, linetype = "dashed", color = "grey",size=1)+
  labs(x = "Heterozygosity per kbp", y = "Het. LOF genotypes absent from SC_WSP")+
  my_theme2()


# DEL+LOF method1 -----------------------------------------------------------------
plot_method1 <- function(x_item,mut_type,xmin,xmax,ymin,ymax){
  pic <- 
    ggplot() + 
    geom_rect(
      aes(
        xmin = {{ xmin }},    
        xmax = {{ xmax }},    
        ymin = {{ ymin }},
        ymax = {{ ymax }}
      ),
      fill = "#82B29A",  
      alpha = 0.1          
    ) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_MSH'),aes(x={{ x_item }},y={{ mut_type }}),fill="#DF7A5E",size=4,shape=24,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),aes(x={{ x_item }},y={{ mut_type }}),fill="#82B29A",size=4,shape=25,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP' & Pedigree_generation == 'F1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#D2BFA5',size=4,shape=23,col="#8FB3DC",alpha = 0.8,stroke= 1.5) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP'& Pedigree_generation != 'F1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#D2BFA5',size=4,shape=23,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method1_DEL_LOF != 'HL'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method1_DEL_LOF == 'HL'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="#8e7cc3",alpha = 0.8, stroke=1.5) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                color = "black",
                alpha = 0.3,
                level = 0.95,
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                color = "#82B29A",
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                se=FALSE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                color = "black",
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                color = "#F2CC8E",
                se=FALSE) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..eq.label.., ..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      hjust = -0.2
    ) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..eq.label.., ..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      vjust = 2,  
      hjust = -0.2
    ) +
    guides(fill="none")
  return (pic)
}

plot_method2 <- function(x_item,mut_type){
  pic <- 
    ggplot() + 
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method2_DEL_LOF %in% c("HH", "LH")),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method2_DEL_LOF %in% c("HL", "LL")),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="#00796B",alpha = 0.8, stroke=1.5) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                color = "black",
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                color = "#F2CC8E",
                se=FALSE) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..eq.label.., ..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      vjust = 2,  
      hjust = -0.2
    ) +
    guides(fill="none")
  return (pic)
}

# use func


method1_DEL_LOF <- 
  plot_method1(het_whole_genome,het_DEL+het_LOF,1.15,Inf,-Inf,928.25+72.5)+
  geom_hline(yintercept = 928.25+72.5, linetype = "dashed", color = "grey",size=1)+
  geom_vline(xintercept = 1.15, linetype = "dashed", color = "grey",size=1)+
  labs(x = "Heterozygosity per kbp", y = "Het. DEL+LOF genotypes absent from SC_WSP")+
  my_theme2()


method2_DEL_LOF <- 
  plot_method2(het_whole_genome,het_DEL+het_LOF)+
  geom_vline(xintercept = 1.21, linetype = "dashed", color = "grey",size=1)+  
  labs(x = "Heterozygosity per kbp", y = "Het. DEL+LOF genotypes absent from SC_WSP")+
  my_theme2()


# DEL+LOF method 2/3 --------------------------------------------------------------


method3_DEL_LOF <- 
  ggplot() +
  geom_rect(
    aes(
      xmin = 1.15,    
      xmax = Inf,    
      ymin = -Inf, 
      ymax = 1705
    ),
    fill = "#82B29A",  
    alpha = 0.2          
  ) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'  & method1_DEL_LOF != "HL"),aes(x=het_whole_genome,y=snpMisDelNum_het+snpLOFNum_het,fill=het_DEL+het_LOF),size=4,shape=21,col="black",alpha = 0.8) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method1_DEL_LOF == "HL"),aes(x=het_whole_genome,y=snpMisDelNum_het+snpLOFNum_het,fill=het_DEL+het_LOF),size=4,shape=21,col="#8e7cc3",alpha = 0.8,stroke= 1.5) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method2_DEL_LOF %in% c("HL", "LL") & method1_DEL_LOF!= "HL"),aes(x=het_whole_genome,y=snpMisDelNum_het+snpLOFNum_het,fill=het_DEL+het_LOF),size=4,shape=21,col="#00796B",alpha = 0.8,stroke= 1.5) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP' & Pedigree_generation == 'F1'),aes(x=het_whole_genome,y=snpMisDelNum_het+snpLOFNum_het,fill=het_DEL+het_LOF),size=4,shape=23,col="#8FB3DC",alpha = 0.8,stroke= 1.5) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP'& Pedigree_generation != 'F1'),aes(x=het_whole_genome,y=snpMisDelNum_het+snpLOFNum_het,fill=het_DEL+het_LOF),size=4,shape=23,col="black",alpha = 0.8) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),aes(x=het_whole_genome,y=snpMisDelNum_het+snpLOFNum_het,fill=het_DEL+het_LOF),size=4,shape=25,col="black",alpha = 0.8) +   
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_MSH'),aes(x=het_whole_genome,y=snpMisDelNum_het+snpLOFNum_het,fill=het_DEL+het_LOF),size=4,shape=24,col="black",alpha = 0.8) +  
  geom_vline(xintercept = 1.15, linetype = "dashed", color = "grey",size=1) +
  geom_hline(yintercept = 1705, linetype = "dashed", color = "grey",size=1)+
  labs(x = "Heterozygosity per kbp", y = "Het. DEL+LOF genotypes count")+
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
              aes(x=het_whole_genome,y=snpMisDelNum_het+snpLOFNum_het),
              method = "lm",
              alpha = 0.3,
              level = 0.95,
              color = "black",
              se=TRUE) +
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
              aes(x=het_whole_genome,y=snpMisDelNum_het+snpLOFNum_het),
              method = "lm",
              color = "#82B29A",
              alpha = 0.3,
              level = 0.95,
              linetype= 'dashed',
              se=FALSE) +
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
              aes(x=het_whole_genome,y=snpMisDelNum_het+snpLOFNum_het),
              method = "lm", 
              alpha = 0.3,
              level = 0.95,
              color = "black",
              se=TRUE) +
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
              aes(x=het_whole_genome,y=snpMisDelNum_het+snpLOFNum_het),
              method = "lm", 
              color = "#F2CC8E",
              alpha = 0.3,
              level = 0.95,
              linetype= 'dashed',
              se=FALSE) +
  stat_poly_eq(
    data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
    aes(x = het_whole_genome, y = snpMisDelNum_het+snpLOFNum_het,  label = paste(..eq.label.., ..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    family = "Times",
    size = 4,
    vjust = 2,  
    hjust = -0.2
  )+ 
  stat_poly_eq(
    data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
    aes(x = het_whole_genome, y = snpMisDelNum_het+snpLOFNum_het,  label = paste(..eq.label.., ..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    family = "Times",
    size = 4,
    hjust = -0.2
  )+
  scale_fill_gradient2(
    low      = "#82B29A",          
    mid      = "white",          
    high     = "firebrick",      
    midpoint = 1000.75  
  ) +
  my_theme2()+
  theme(legend.position=c(0.85, 0.2))


method4_DEL_LOF <- 
  ggplot() + 
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'  & method1_DEL_LOF != "HL"),aes(x=het_whole_genome,y=snpMisDelNum_het+snpLOFNum_het,fill=het_DEL+het_LOF),size=4,shape=21,col="black",alpha = 0.8) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method1_DEL_LOF == "HL"),aes(x=het_whole_genome,y=snpMisDelNum_het+snpLOFNum_het,fill=het_DEL+het_LOF),size=4,shape=21,col="#8e7cc3",alpha = 0.8,stroke= 1.5) +
  geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method2_DEL_LOF %in% c("HL", "LL") & method1_DEL_LOF!= "HL"),aes(x=het_whole_genome,y=snpMisDelNum_het+snpLOFNum_het,fill=het_DEL+het_LOF),size=4,shape=21,col="#00796B",alpha = 0.8,stroke= 1.5) +
  geom_vline(xintercept = 1.21, linetype = "dashed", color = "grey",size=1) +
  labs(x = "Heterozygosity per kbp", y = "Het. DEL genotypes count")+
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
              aes(x=het_whole_genome,y=snpMisDelNum_het+snpLOFNum_het),
              method = "lm", 
              alpha = 0.3,
              level = 0.95,
              color = "black",
              se=TRUE) +
  geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
              aes(x=het_whole_genome,y=snpMisDelNum_het+snpLOFNum_het),
              method = "lm", 
              color = "#F2CC8E",
              alpha = 0.3,
              level = 0.95,
              linetype= 'dashed',
              se=FALSE) +
  stat_poly_eq(
    data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
    aes(x = het_whole_genome, y = snpMisDelNum_het+snpLOFNum_het,  label = paste(..eq.label.., ..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    family = "Times",
    size = 4,
    vjust = 2,  
    hjust = -0.2
  )+ 
  scale_fill_gradient2(
    low      = "#82B29A",           
    mid      = "white",         
    high     = "firebrick",      
    midpoint = 960.37+76.033  
  ) +
  my_theme2()+
  theme(legend.position=c(0.85, 0.2))

(method1_DEL_LOF|method2_DEL_LOF)/(method3_DEL_LOF|method4_DEL_LOF)

# DEL+LOF compare method --------------------------------------------------


new_del_LOF_DEL <- new_del_LOF_DEL %>%
  mutate(
    method1_DEL_LOF = case_when(
      rescue_pop != "SC_CAP" ~ rescue_pop,
      het_whole_genome >= 1.15 & het_LOF+het_DEL >= 928.25+72.5 ~ "HH",
      het_whole_genome >= 1.15 & het_LOF+het_DEL < 928.25+72.5 ~ "HL",
      het_whole_genome < 1.15 & het_LOF+het_DEL >= 928.25+72.5 ~ "LH",
      TRUE ~ "LL"
    ),
    method2_DEL_LOF = case_when(
      rescue_pop != "SC_CAP" ~ rescue_pop,
      het_whole_genome >= 1.21 & (het_whole_genome * 841 + 16.8 - het_LOF -het_DEL) <= 0 ~ "HH",
      het_whole_genome >= 1.21 & (het_whole_genome * 841 + 16.8 - het_LOF-het_DEL) > 0 ~ "HL",
      het_whole_genome < 1.21 & (het_whole_genome * 841 + 16.8 - het_LOF-het_DEL) <= 0 ~ "LH",
      TRUE ~ "LL"
    ),
    method3_DEL_LOF = case_when(
      rescue_pop != "SC_CAP" ~ rescue_pop,
      het_whole_genome >= 1.15 & (snpLOFNum_het+snpMisDelNum_het) >= 1705 ~ "HH",
      het_whole_genome >= 1.15 & (snpLOFNum_het+snpMisDelNum_het) < 1705 ~ "HL",
      het_whole_genome < 1.15 & (snpLOFNum_het+snpMisDelNum_het) >= 1705 ~ "LH",
      TRUE ~ "LL"
    ),
    group_DEL_LOF = 1330 * het_whole_genome +144 - snpLOFNum_het -snpMisDelNum_het,
    DEL_LOF_flag = if_else(group_DEL_LOF > 0, "lower", "higher"),
    method4_DEL_LOF = case_when(
      rescue_pop != "SC_CAP" ~ rescue_pop,
      DEL_LOF_flag == "higher" & het_whole_genome >= 1.21 ~ "HH",
      DEL_LOF_flag == "higher" & het_whole_genome < 1.21 ~ "LH",
      DEL_LOF_flag == "lower" & het_whole_genome >= 1.21 ~ "HL",
      TRUE ~ "LL"
    )
  )

DEL_LOF_method1_ref_het_LOF_absent_count <- plot_group(
  new_del_LOF_DEL, "method1_DEL_LOF", "het_LOF+het_DEL",
  "Het DEL+LOF genotypes\nabsent from SC_WSP"
)+top_theme()
DEL_LOF_method1_het_LOF_absent_count <- plot_group(
  new_del_LOF_DEL, "method2_DEL_LOF", "het_LOF+het_DEL",
  "Het DEL+LOF genotypes\nabsent from SC_WSP"
)+top_theme()
DEL_LOF_method2_het_LOF_absent_count <- plot_group(
  new_del_LOF_DEL, "method3_DEL_LOF", "het_LOF+het_DEL",
  "Het DEL+LOF genotypes\nabsent from SC_WSP"
)+top_theme()
DEL_LOF_method3_het_LOF_absent_count <- plot_group(
  new_del_LOF_DEL, "method4_DEL_LOF", "het_LOF+het_DEL",
  "Het DEL+LOF genotypes\nabsent from SC_WSP"
)+top_theme()


# merge plot --------------------------------------------------------------

combined_plot <- (
  (het_DEL_absent | hom_DEL_absent | het_LOF_absent | hom_LOF_absent) /
  (het_new_het_Del_ref_strict|het_new_het_Del_ref_loose|CAP_new_del_plot_DEL|CAP_new_del_plot_2_DEL)  /
  (het_new_het_LOF|het_new_het_LOF_ref_loose|CAP_new_del_plot_LOF|CAP_new_del_plot_2_LOF)
) + plot_layout(heights = c(0.5, 1, 1))
ggsave(combined_plot, file='D://重测序项目//结果部分//绘图//投稿版本//commb//combined_plot_select.pdf',device = cairo_pdf, width=16, height=16)

# Supplementary -----------------------------------------------------------
new_del_plot <- function(x_item,mut_type,xmin,xmax,ymin,ymax){
  pic <- 
    ggplot() + 
    geom_rect(
      aes(
        xmin = {{ xmin }},    
        xmax = {{ xmax }},    
        ymin = {{ ymin }}, 
        ymax = {{ ymax }}
      ),
      fill = "#82B29A", 
      alpha = 0.1          
    ) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_MSH'),aes(x={{ x_item }},y={{ mut_type }}),fill="#DF7A5E",size=4,shape=24,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),aes(x={{ x_item }},y={{ mut_type }}),fill="#82B29A",size=4,shape=25,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP' & Pedigree_generation == 'F1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#D2BFA5',size=4,shape=23,col="#8FB3DC",alpha = 0.8,stroke= 1.5) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP'& Pedigree_generation != 'F1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#D2BFA5',size=4,shape=23,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_1=='0'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good_DEL_method_1=='1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="#8e7cc3",alpha = 0.8, stroke=1.5) +   
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' &  good_DEL_method_2=='1' & good_DEL_method_1!='1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="#00796B",alpha = 0.8, stroke=1.5) +   
    #geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & good=='2'),aes(x={{ x_item }},y={{ mut_type }}),fill='#F2CC8E',size=4,shape=21,col="#9acfea",alpha = 0.8,stroke= 1.5) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                color = "black",
                alpha = 0.3,
                level = 0.95,
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                color = "#82B29A",
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                se=FALSE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                color = "black",
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                color = "#F2CC8E",
                se=FALSE) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
      # aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..eq.label.., ..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      hjust = -0.2
    ) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      vjust = 2,  
      hjust = -0.2
    ) +
    guides(fill="none")
  return (pic)
}

new_del_plot_LOF <- function(x_item,mut_type,xmin,xmax,ymin,ymax){
  pic <- 
    ggplot() + 
    geom_rect(
      aes(
        xmin = {{ xmin }},    
        xmax = {{ xmax }},    
        ymin = {{ ymin }}, 
        ymax = {{ ymax }}
      ),
      fill = "#82B29A",  
      alpha = 0.1          
    ) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_MSH'),aes(x={{ x_item }},y={{ mut_type }}),fill="#DF7A5E",size=4,shape=24,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),aes(x={{ x_item }},y={{ mut_type }}),fill="#82B29A",size=4,shape=25,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP' & Pedigree_generation == 'F1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#D2BFA5',size=4,shape=23,col="#8FB3DC",alpha = 0.8,stroke= 1.5) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='MIX_CAP'& Pedigree_generation != 'F1'),aes(x={{ x_item }},y={{ mut_type }}),fill='#D2BFA5',size=4,shape=23,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method1_LOF !='HL' & method2_LOF %in% c("HH", "LH")),aes(x={{ x_item }},y={{ mut_type }}),fill="#F2CC8E",size=4,shape=21,col="black",alpha = 0.8) +
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method1_LOF=='HL'),aes(x={{ x_item }},y={{ mut_type }}),fill="#F2CC8E",size=4,shape=21,col="#8e7cc3",alpha = 0.8,stroke= 1.5) +   
    geom_point(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP' & method2_LOF %in% c("HL", "LL") & method1_LOF!='HL'),aes(x={{ x_item }},y={{ mut_type }}),fill="#F2CC8E",size=4,shape=21,col="#00796B",alpha = 0.8,stroke= 1.5) +   
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                color = "black",
                alpha = 0.3,
                level = 0.95,
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                color = "#82B29A",
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                se=FALSE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                color = "black",
                se=TRUE) +
    geom_smooth(data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
                aes(x={{ x_item }},y={{ mut_type }}),
                method = "lm", 
                alpha = 0.3,
                level = 0.95,
                linetype= 'dashed',
                color = "#F2CC8E",
                se=FALSE) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_QLA'),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      hjust = -0.2
    ) +
    stat_poly_eq(
      data = subset(new_del_LOF_DEL, rescue_pop =='SC_CAP'),
      aes(x = {{ x_item }}, y = {{ mut_type }},  label = paste(..rr.label.., "italic(p) == ", signif(..p.value.., digits = 3), sep = "~~~")),
      formula = y ~ x,
      parse = TRUE,
      family = "Times",
      size = 4,
      vjust = 2,  
      hjust = -0.2
    ) +
    guides(fill="none")
  return (pic)
}

het_new_hom_Del <- 
  new_del_plot(het_whole_genome,hom_DEL,1.15,Inf,-Inf,218.5)+
  ylim(120,300)+
  geom_hline(yintercept = 218.5, linetype = "dashed", color = "grey",size=1)+
  geom_vline(xintercept = 1.15, linetype = "dashed", color = "grey",size=1)+
  labs(x = "Heterozygosity per kbp", y = "Recipient-absent Hom. DEL genotype count")+
  my_theme2()

het_new_hom_LOF <-
  new_del_plot_LOF(het_whole_genome,hom_LOF,1.15,Inf,-Inf,21.375)+
  geom_hline(yintercept = 21.375, linetype = "dashed", color = "grey",size=1)+
  geom_vline(xintercept = 1.15, linetype = "dashed", color = "grey",size=1)+
  labs(x = "Heterozygosity per kbp", y = "Recipient-absent Hom. LOF genotype count")+
  my_theme2()


stack_with_gap <- function(top_plot, bottom_plot, gap = 0.15) {
  wrap_plots(
    top_plot,
    plot_spacer(),
    bottom_plot,
    ncol = 1,
    heights = c(1, gap, 1)
  )
}

S_absent <- stack_with_gap(
    (method2_het_DEL_absent_count|method3_het_DEL_absent_count|LOF_method2_het_LOF_absent_count|LOF_method3_het_LOF_absent_count), 
    (method2_hetkbp|method3_hetkbp|LOF_method2_hetkbp|LOF_method3_hetkbp)
  )

combined_plot_S <- (
   (het_new_hom_Del| het_new_hom_LOF) /
   S_absent
) + plot_layout(heights = c(1, 2))
combined_plot_S

ggsave(combined_plot_S, file='combined_plot_select_S.pdf',device = cairo_pdf, width=16, height=12)
