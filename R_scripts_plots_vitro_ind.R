# Figures for individually in vitro assayed plerixafor sample experiments
# Barplots show means +/-SEM of fold-change combared to paired baseline samples. 
# Paul Miller, paulhmiller@gmail.com


library(reshape2)
library(ggplot2)
library(grid)
library(plyr)
source("C:/Users/paulm/Documents/R/source/functions.R")
source("C:/Users/paulm/Documents/R/source/plotting_themes.R")

setwd('C:/Users/paulm/CRC Paul/PROJECTS/Plerixafor/figures/in vitro individual postthaw')

dat <- read.csv("plerixafor_vitro_individual3.csv",header=T,check.names=F,colClasses=c("Experiment"="factor")) 

# Anonymize donor names 
levels(dat$Donor) <- 1:10
dat <- dat[,1:18] 
assays <- names(dat[7:18])


## Determine fold changes

# Global BL log10 geometric mean:
BLm <- 10^(colMeans(log10(dat[dat$Timepoint=="BL", 7:18])))
# Global BL arithmetic mean:
#BLm <- colMeans(dat[dat$Timepoint=="BL", 7:18])

# Fold values:
# Donor for-loop (need to split into 2 separate loops because of unequal sample
# numbers for each donor group)
foldA <- NULL
for(i in c(1,3,4,6,7)){
	D0 <- dat[dat$Donor==i & dat$Timepoint=="4hrs", ] 
	BL <- dat[dat$Donor==i & dat$Timepoint=="BL", ] 
	D0[,7:18] <- D0[, 7:18] / BL[, 7:18] 
	BL[,7:18] <- BL[, 7:18] / BLm
	foldA <- rbind(foldA,BL,D0)
}
foldB <- NULL
for(i in c(2,5,8,9,10)){
  D0 <- dat[dat$Donor==i & dat$Timepoint=="4hrs", ] 
  Dm1 <- dat[dat$Donor==i & dat$Timepoint=="-24hrs", ] 
  BL <- dat[dat$Donor==i & dat$Timepoint=="BL", ] 
  D0[,7:18] <- D0[, 7:18] / BL[, 7:18] 
  Dm1[,7:18] <- Dm1[, 7:18] / BL[, 7:18] 
  BL[,7:18] <- BL[, 7:18] / BLm
  foldB <- rbind(foldB,BL,D0,Dm1)
}
fold <- rbind(foldA, foldB)

	
mean(fold[fold$Timepoint=="BL",10])

dat <- fold

# Log10 transform
dat[7:ncol(dat)] <- log10(dat[7:ncol(dat)])

# Rename group factors (had to separate because couldn't rename factor based 
# on two columns. Consequence is that BL values now removed, so don't need to 
# remove BL separately)
P <- dat[dat$Group=="A" & dat$Timepoint=="4hrs",]
GP <- dat[dat$Group=="B" & dat$Timepoint=="4hrs",]
G <- dat[dat$Group=="B" & dat$Timepoint=="-24hrs",]
levels(G$Group)[levels(G$Group)=="B"] <- "G"
dat <- rbind(P,GP,G)
levels(dat$Group)[levels(dat$Group)=="A"] <- "P"   
levels(dat$Group)[levels(dat$Group)=="B"] <- "G+P"


# Check for normality for each group
out <- NULL
for (i in c('P','G+P','G')){
  tmp <- dat[dat$Group==i ,7:ncol(dat)]
  tmp <- lapply(tmp, shapiro.test)
  tmp <- sapply(tmp, '[', c("statistic", "p.value"))
  tmp <- cbind(rep(i,2), tmp)
  out <- rbind(out, tmp)
}
print(out)
# The above indicates that the log data is not not-normal. 
# (i.e. can't reject that it isn't normal)

# Convert to tall format
dat<-melt(dat,id.vars=c("Experiment","ID","Donor","Sample","Timepoint","Group"),variable.name="assay",value.name="value",na.rm=TRUE)  

# Re-order factor levels:
dat$Timepoint<-factor(dat$Timepoint,levels=c("BL","4hrs","-24hrs")) 
dat$Group<-factor(dat$Group,levels=c("G","P","G+P")) 

PB<-dat 


#*** add 2 log to all values so barplots look OK. MUST change y-axis breaks. 
PB[, 8] <- PB[, 8] +2

# Exclude BL
PB <- PB[PB$Timepoint!="BL",]

PB<-summarySE(PB, measurevar="value", groupvars=c("Timepoint","Group","assay"),
na.rm=TRUE) 
# SummarySE provides std, SEM, and (default 95%) CI. #measurevar is the x-axis

CFC <- PB[PB$assay=="total.CFC",]
w3LTCTNC <- PB[PB$assay=="TNC.wk3",]
w6LTCTNC <- PB[PB$assay=="TNC.wk6",]
w3LTCIC <- PB[PB$assay=="LTC-IC.wk3",]
w6LTCIC <- PB[PB$assay=="LTC-IC.wk6",]
w3TNC <- PB[PB$assay=="TNC.wk3",]
w6TNC <- PB[PB$assay=="TNC.wk6",]
w3LTC15 <- PB[PB$assay=="CD15.wk3",]
w6LTC15 <- PB[PB$assay=="CD15.wk6",]
w3LTC34 <- PB[PB$assay=="CD34.wk3",]
w6LTC34 <- PB[PB$assay=="CD34.wk6",]


# To look at individual dots. Needs to hash out summarySE above.
#library(beeswarm) 
#beeswarm(CFC$value ~ CFC$Group) 


## Plot by time-point

# Combine into time points
foldw2 <- CFC
foldw3 <- rbind(w3TNC, w3LTC15, w3LTC34, w3LTCIC)
foldw6 <- rbind(w6TNC, w6LTC15, w6LTC34, w6LTCIC)


plot <- ggplot(data=foldw2,aes(x=assay, y=value, fill=Group))+ 
  annotation_logticks(sides = "l", size=0.3) +
  geom_bar(stat="identity", position=position_dodge())+ 
  scale_fill_manual(values=p1) + 
  scale_color_manual(values=p1) +		
  geom_errorbar(aes(ymin=value-se, ymax=value+se, color=Group), lty=1, width=0.5, size=0.75, position=position_dodge(width=0.9)) +
  ggtitle("") + 
  guides(fill=F,color=F)+
  scale_y_continuous(breaks=c(1, 2, 3, 4, 5, 6),labels = c(expression(10^-1), expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4))) + 
  coord_cartesian(ylim=c(2.0, 4.0), expand = FALSE) +
  ylab("fold-increase") + 
  scale_x_discrete("", labels=("CFC")) +
  themePM1()
plot + 	theme(aspect.ratio=NULL) 
ggsave(filename="mPB CFC ind.pdf",width=2.5,height=5.15, units="cm")

# Make x-labels for the following plots
labels <- c("TNC",expression("33/15"^"+"),expression("34"^"+"),"CFC") 

plot %+% foldw3 +
  ggtitle("Week 3 LTC") + 
  scale_x_discrete("", labels=labels) +
  coord_cartesian(ylim=c(2.18,5.9)) +
  theme() 
ggsave(filename="mPB w3LTC ind.pdf",width=5.25,height=5.3, units="cm")


plot %+% foldw6 +
  ggtitle("Week 6 LTC") + 
  scale_x_discrete("", labels=labels) +
  coord_cartesian(ylim=c(2.18,5.9)) +
  theme()  
ggsave(filename="mPB w6LTC ind.pdf",width=5.25,height=5.3, units="cm")



 


# t.test (assumes equal variance) - FDR corrected for multiple testing

# Comparison between arms:
tmp1<-dat[dat$Timepoint=="4hrs",]
tmp2<-dat[dat$Timepoint=="-24hrs",]
PB <- rbind(tmp1,tmp2)

sink('ttest FDR.txt')
stats <- NULL
for(i in 4:12){ 
  tmp <- PB[PB$assay==assays[[i]], c(6,8)]
  print(assays[[i]])
  out <- pairwise.t.test(tmp$value, tmp$Group, p.adjust.method="fdr")
  print(out)
}
sink()

