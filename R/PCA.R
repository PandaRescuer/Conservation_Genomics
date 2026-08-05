library("ggplot2")
library("ggstar")
eigenvec <- read.table("all_pca.eigenvec", quote="\"", comment.char="")
colnames(eigenvec)<-c("FID","sample",paste0("PC",1:5))
write.table(eigenvec[2:ncol(eigenvec)],file = "pca.eigenvector.xls",sep = "\t",row.names = F,col.names = T,quote = F)

eigenval <- read.table("all_pca.eigenval", quote="\"", comment.char="")
pcs<-paste0("pc",1:nrow(eigenval))
eigenval[nrow(eigenval),1]<-0

percentage<-eigenval$V1/sum(eigenval$V1)*100
eigenval_df<-as.data.frame(cbind(pcs,eigenval[,1],percentage),stringsAsFactors = F)
names(eigenval_df)<-c("pcs","variance","proportation")
eigenval_df$variance<-as.numeric(eigenval_df$variance)
eigenval_df$proportation<-as.numeric(eigenval_df$proportation)
write.table(eigenval_df,file = "pca.eigenvalue.xls",sep = "\t",quote = F,row.names = F,col.names = T)
head(eigenval_df)

eigvec <-read.table("pca.eigenvector.txt",header = T)
eigval <- read.table("pca.eigenvalue.xls",header = T)

p_pop <- read.table("pop.txt",header=T)
p1 <- merge.data.frame(eigvec,p_pop,by="ID")

wildset<-subset(p1,G2=="0")

PCA <- 
  ggplot(data = p1,aes(x = PC1, y = PC2,fill=POP4,shape=POP4)) +
  geom_point(size=4, alpha = 0.8)+
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey",size=0.5)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey",size=0.5)+
  theme(legend.position = "right")+
  theme(panel.grid = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5))+ 
  #theme(legend.text = element_text(size = 16, face = "bold.italic"),legend.title = element_blank())+ 
  xlab(paste0("PC1 (", round(eigval[eigval$pcs == "pc1", 3], 2), "%)"))+ ylab(paste0("PC2 (", round(eigval[eigval$pcs == "pc2", 3], 2), "%)")) +
  theme(axis.title.x = element_text(face = "bold", size = 15, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 15, colour = "black"),
        axis.text.x = element_text(size = 15,face = "bold", colour = "black"),
        axis.text.y = element_text(face = "bold", size = 15, colour = "black"),
        strip.text = element_text(color = "black", face = "bold",size = 15),
        legend.text = element_text(color = "black", face = "bold",size = 15),
        legend.title = element_text(color = "black", face = "bold",size = 15),)+
  scale_shape_manual(values = c(SC_CAP = 21,
                                QLI = 22,
                                SC_QLA = 25,
                                SC_MSH = 24,
                                MIX_CAP = 23,
                                SC_XXL= 10 
  ))+
  scale_fill_manual(values=c(SC_CAP = "#F2CC8E",
                             QLI = "#82B29A",
                             SC_QLA = "#df995e", 
                             SC_MSH = "#DF7A5E" , 
                             MIX_CAP = "#D2BFA5",
                             SC_XXL="#8e7cc3"))
PCA
ggsave(PCA, file='PCA.pdf', width=9.5, height=5.5)

