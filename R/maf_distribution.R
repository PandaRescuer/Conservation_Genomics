library(ggplot2)
library(reshape2)

p1 <- read.table("LOF_snp_RMA_MAF_distribution.txt",header=F)
p2 <- read.table("synonymous_snp_RMA_MAF_distribution.txt",header=F)
p3 <- read.table("CNE_snp_RMA_MAF_distribution.txt",header=F)
p4 <- read.table("missense_snp_RMA_MAF_distribution.txt",header=F)
p5 <- read.table("mis_snp_DELETERIOUS_MAF_distribution.txt",header=F)
p6 <- read.table("mis_snp_TOLERATED_MAF_distribution.txt",header=F)


sumV2 <- sum(p1$V2)

p1$V2 <- p1$V2 / sumV2

sumV2 <- sum(p2$V2)
p2$V2 <- p2$V2 / sumV2

sumV2 <- sum(p3$V2)
p3$V2 <- p3$V2 / sumV2

sumV2 <- sum(p4$V2)
p4$V2 <- p4$V2 / sumV2

sumV2 <- sum(p5$V2)
p5$V2 <- p5$V2 / sumV2

sumV2 <- sum(p6$V2)
p6$V2 <- p6$V2 / sumV2

names(p1) <- c("maf", "LOF")
names(p2) <- c("maf", "Syn")
names(p3) <- c("maf", "CNE")
names(p4) <- c("maf", "Mis")
names(p5) <- c("maf", "Del")
names(p6) <- c("maf", "Tol")

data_merge <- merge.data.frame(p1,p5,by="maf")
data_merge <- merge.data.frame(data_merge,p3,by="maf")
data_merge <- merge.data.frame(data_merge,p4,by="maf")
data_merge <- merge.data.frame(data_merge,p6,by="maf")
data_merge <- merge.data.frame(data_merge,p2,by="maf")
data_merge <- melt(data_merge, id="maf")

MAF <- ggplot(data=data_merge,aes(x=maf+0.025, y=value, fill=variable))+
  geom_bar(stat="identity",position = 'dodge', colour = "black",linewidth=0.25)+
  theme_bw()+
  theme_classic()+
  scale_x_continuous(
    limits=c(0,1),
    breaks=seq(0,1,0.05))+
  scale_y_continuous(
    limits=c(0,0.6),
    labels = seq(0,0.6,0.1),
    breaks = seq(0,0.6,0.1),
    expand = c(0,0)) +
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
  ) +
  labs(x = "Derived allele frequency", y = "Propotion",fill = "Type") +
  scale_fill_manual(values =c("#DF7A5E","#F2CC8E","#D2BFA5","#82B29A","#8e7cc3","#7fdae9"))


ggsave(
  paste0("S1B.pdf"),
  plot = MAF,
  width = 8,
  height = 4.5
)

data_RMA_prop <- read.table("MAF_RMA.txt",header=T)
names<-c('LOF','Mis_DEL','CNE','Missense','Mis_TOL','Synonymous')
data_RMA_prop$variant <- factor(data_RMA_prop$variant,
                              levels = names,
                              ordered = T)
data_RMA_prop$type <- factor(data_RMA_prop$type, 
                              levels = c("RMA", "notRMA"),
                              ordered = T)

S2B <- ggplot(data_RMA_prop,aes(x=variant,y=prop,fill=type))+
  geom_bar(stat = 'identity',position = 'dodge',width = 0.5,colour='black')+
  theme_bw()+
  theme_classic()+
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
  ) +
  scale_x_discrete(labels = c(
    "Mis_DEL" = "DEL",
    "Missense" = "MIS",
    "Mis_TOL" = "TOL",
    "Synonymous" = "SYN"
  ))+
  scale_y_continuous(
    limits = c(0,1),
    labels = seq(0,1,0.25),
    breaks = seq(0,1,0.25),
    expand = c(0,0)) +
  labs(x = "", y = "Propotion") +
  scale_fill_manual(values =c("#DF7A5E","#F2CC8E"))
S2B
ggsave(
  paste0("S2B.pdf"),
  plot = S2B,
  width = 4,
  height = 4.5
)
