library(ggplot2)
slim_result <- read.table('panda_slim_50000_keep.txt',header=T,sep=',',fill = T)

my_theme <- function() {
  theme_classic()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5))+ 
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

fitness_pic <- ggplot()+
  geom_line(
    data=slim_result,
    mapping=aes(x=gen, y=meanFitness),
    color = '#DF7A5E',
    size=1) +
  labs(x = "Year", y = "Fitness")+
  guides(fill="none")+
  scale_y_continuous(limits=c(0.96,1),
                     labels = seq(0.96,1,0.01),
                     breaks = seq(0.96,1,0.01),
  ) +
  my_theme()
fitness_pic

hete_pic <- ggplot()+
  geom_line(
    data=slim_result,
    mapping=aes(x=gen, y=meanHet),
    color = '#DF7A5E',
    size=1) +
  labs(x = "Year", y = "Heterozygosity")+
  guides(fill="none")+
  scale_y_continuous(limits=c(0,1e-04),
                     labels = seq(0,1e-04,0.25e-04),
                     breaks = seq(0,1e-04,0.25e-04),
  ) +
  my_theme()
hete_pic

FROH_pic <- ggplot()+
  geom_line(
    data=slim_result,
    mapping=aes(x=gen, y=FROH_1Mb),
    color = '#DF7A5E',
    size=1) +
  labs(x = "Year", y = "FROH")+
  guides(fill="none")+
  my_theme()
FROH_pic

load_pic <- ggplot()+
  geom_line(
    data=slim_result,
    mapping=aes(x=gen, y=B_gen),
    color = '#DF7A5E',
    size=1) +
  labs(x = "Year", y = "Inbreeding load")+
  guides(fill="none")+
  scale_y_continuous(limits=c(0,4),
                     labels = seq(0,4,1),
                     breaks = seq(0,4,1),
  ) +
  my_theme()
load_pic

het_pic <- ggplot()+
  geom_line(
    data=slim_result,
    mapping=aes(x=gen, y=avgHetDel),
    color = '#DF7A5E',
    size=1) +
  labs(x = "Year", y = "Average heterozygous Deleteriou mutations")+
  guides(fill="none")+
  my_theme()
het_pic

hom_pic <- ggplot()+
  geom_line(
    data=slim_result,
    mapping=aes(x=gen, y=avgHomDel ),
    color = '#DF7A5E',
    size=1) +
  labs(x = "Year", y = "Average homozygous Deleteriou mutations")+
  guides(fill="none")+
  my_theme()
hom_pic

burn_in_plot <- (fitness_pic|FROH_pic)/(hete_pic|load_pic)
ggsave(burn_in_plot, file='S_burn_in.pdf',device = cairo_pdf, width=16, height=8)
