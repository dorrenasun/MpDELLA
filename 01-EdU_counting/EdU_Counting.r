#Modified by Rui Sun, 2020-02-21
rm(list = ls())
if(!require(spatstat)) install.packages("spatstat")
if(!require(viridis)){install.packages("viridis")}

# library(raster)
library(viridis)
library(spatstat)

#read data
Path<-"./Coordinates"
Cors<-list.files(path = Path,pattern = ".csv")
Cors<-Cors[Cors!="Summary.csv"]

#Create output folder
outPATH<-"./DensityMaps"
if (dir.exists(outPATH)){
  unlink(outPATH,recursive = T)
}
dir.create(outPATH)
outPATH<-paste(outPATH,"/",sep = "")
maskPATH<-"./DensityMasks"
if (dir.exists(maskPATH)){
  unlink(maskPATH,recursive = T)
}
dir.create(maskPATH)
maskPATH<-paste(maskPATH,"/",sep = "")

table<-data.frame(name=Cors,width=NA,height=NA,count=NA,MaxDen=NA,MaskArea=NA,MaskPoints=NA)

#Color maps
Cm<-colourmap(plasma(256),range=c(-1e-10,0.015))
par(fg="white")
par(mar=rep(0,4))
dpi<-300
scale<-2.35 # pixel/um

#Analyze each pic and generate density maps
for (i in 1:length(Cors)){
# for (i in 1:5){
  # i<-19
  data<-read.csv(paste0(Path,"/",Cors[i]))
  
  table$name[i]<-Cors[i]
  table$width[i]<-max(data$Width)
  table$height[i]<-max(data$Height)
  table$count[i]<-nrow(data)-1

  data<-data[-nrow(data),]
  data$X<-as.numeric(data$X)
  data$Y<-as.numeric(data$Y)
  
  Bw<-10
  
  data.ppp<-ppp(data$X,data$Y,c(0,table$width[i]),c(0, table$height[i]))
  #change the unit to positive nuclei/(um^2) Scale: 235 px = 100 um
  data.ppp<-rescale(data.ppp,scale,"um")
  # if(i==19) data.ppp<-rescale(data.ppp,1.32,"um")
  den<-density(data.ppp,sigma=20,edge=F,dimyx=c(ceiling(table$height[i]/Bw),ceiling(table$width[i]/Bw)))
  table$MaxDen[i]<-max(den$v)
  den.mask<-cut(den,c(0.001,0.1))
  den.polygon<-as.polygonal(as.owin(den.mask))
  data.sub<-data.ppp[den.polygon]
  table$MaskArea[i]<-summary(data.sub)$window$area
  table$MaskPoints[i]<-summary(data.sub)$n
  jpeg(filename = paste0(outPATH,Cors[i],".jpg"),width = table$width[i]/0.8325,height = table$height[i]/0.8325,res = dpi)
  par(mar=rep(0,4),
      xaxs="i",
      yaxs="i",
      pty="m",
      fg="white",
      mgp=c(0, 0, 0)#,
      )
  plot(den,col=Cm,main=NULL,asp=table$height[i]/table$width[i],legend=FALSE,scale=FALSE)
  plot(den.polygon,lwd=1,add=T)
  points(data.ppp,pch=".")

  dev.off()
  
  jpeg(filename = paste0(maskPATH,Cors[i],".mask.jpg"),width = table$width[i]/0.8325,height = table$height[i]/0.8325,res = dpi)
  par(mar=rep(0,4),
      xaxs="i",
      yaxs="i",
      pty="m",
      fg="white",
      mgp=c(0, 0, 0)#,
  )
  plot(den,col=Cm,main=NULL,asp=table$height[i]/table$width[i],legend=FALSE,scale=FALSE)
  plot(den.polygon,lwd=1,col="white",add=T)
  dev.off()  
}
pdf(file = "colormap.pdf",width = 5,height = 5)
plot(den,col=Cm,main=NULL,asp=table$height[i]/table$width[i])
dev.off()

while (length(dev.list()>0)) dev.off()

table$MaskRatio<-table$MaskPoints/table$count
ListName<-"list.csv"
if (!file.exists(ListName)){
  emptylist<-subset(table, select = "name")
  emptylist$genotype<-""
  emptylist$rank<""
  write.csv(emptylist,file = ListName,row.names = F)
  stop(paste0("Please fill genotypes in ",ListName))
} else {
  list<-read.csv(ListName)
  table$genotype<-list$genotype[match(table$name,list$name)]
  if (sum(is.na(table$genotype))>0) stop(paste0("Check your ",ListName))
}
table$genotype<-factor(table$genotype,levels=unique(list$genotype[order(list$rank)]))
table$MaskArea<-table$MaskArea/1e6

write.table(table,file="table.tsv",sep = "\t",row.names = F)

#Calculate Tukey groups
if(!require(multcomp)){install.packages("multcomp")}
if(!require(emmeans)){install.packages("emmeans")}

model.ct = lm(count ~ genotype,
               data = table)
marginal.ct = emmeans(model.ct,~ genotype)
CLD.ct = cld(marginal.ct,
             alpha=0.05,
             sort=F,
             Letters=letters,
             adjust="tukey")
CLD.ct$.group<-gsub("[[:space:]]", "", CLD.ct$.group)
model.area = lm(MaskArea ~ genotype,
              data = table)
marginal.area = emmeans(model.area,~ genotype)
CLD.area = cld(marginal.area,
             alpha=0.05,
             sort=F,
             Letters=letters,
             adjust="tukey")
CLD.area$.group<-gsub("[[:space:]]", "", CLD.area$.group)


#plot boxes
if(!require(ggplot2)){install.packages("ggplot2")}
library(ggplot2)

#theme for ggplot
Lwd=0.625
theme_typical<-theme_bw()+theme(
  aspect.ratio = 1,
  plot.margin = unit(c(2,1,1,1), "cm"),
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour="black",size=1.5),
  # axis.line.x = element_line(colour="black",size=1),
  # axis.line.y=element_line(colour="black",size=1),
  # strip.background = element_rect(fill="white",colour="black",size=1.5),
  strip.background = element_blank(),
  strip.text=element_text(face="bold",size = 14, colour = "black"),
  # axis.title.x = element_text(face="bold",size = 14, colour = "black"),
  axis.title.x = element_blank(),
  axis.title.y = element_text(face="bold",size = 14, colour = "black"),
  axis.text.x = element_text(face="bold",size = 12, angle = 315, hjust=0, colour = "black"),
  axis.text.y = element_text(face="bold",size = 12, colour = "black"),
  legend.title = element_text(face = "bold",size = 12),
  legend.text = element_text(face = "bold",size = 12),
  legend.text.align = 0
  #legend.position = "top"
)

p.ct<-ggplot(table, aes(x=genotype, y=count))+
  geom_boxplot(alpha=1)+
  geom_jitter(position=position_jitter(0.2),
              color="darkgreen",
              alpha=0.6,
              stroke=0,
              # width = 0.1,
              size=2) +
  geom_text(data = CLD.ct, aes(x =genotype, y = upper.CL, label = .group),
              position = position_dodge(width = 0.75),
              vjust=-1,hjust=0.5,
              size=2)
p.area<-ggplot(table, aes(x=genotype, y=MaskArea))+
  geom_boxplot(alpha=1)+
  geom_jitter(position=position_jitter(0.2),
              color="darkgreen",
              alpha=0.6,
              stroke=0,
              # width = 0.1,
              size=2) +
  geom_text(data = CLD.area, aes(x =genotype, y = upper.CL, label = .group),
            position = position_dodge(width = 0.75),
            vjust=-1,hjust=0.5,
            size=2)
# 
#Titles
xtitle="genotype"
ytitle.ct="Number of EdU-labelled nuclei"
ytitle.area="Actively dividing area (mm2)"

p.ct<-p.ct+labs(x=xtitle,y=ytitle.ct)+theme_typical
p.area<-p.area+labs(x=xtitle,y=ytitle.area)+theme_typical

Max=round(max(table$count)*1.1,digits = -2)
p.ct<-p.ct+scale_y_continuous(expand= c(0,0),limits = c(0,Max))

library(gridExtra)
p<-grid.arrange(p.ct, p.area, nrow = 2)

print(p)

Fig.wd<-7
Fig.ht<-15

outFile.svg<-"Boxplot.svg"
ggsave(file=outFile.svg, plot=p, width=Fig.wd, height=Fig.ht)

