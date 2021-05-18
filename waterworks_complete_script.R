#Modificeret Rscript til PCS/NMDS/multivariat statistik
#This script contains the bulk of the code used for the papers on waterworks microbiology
#contained in Albers et al. (2015) Environmental Science & Technology, Vol. 49, No. 2, 20.01.2015, p. 839-846.
#- and Harder et al. (2019) FEMS Microbiology Ecology, Vol. 95, No. 11, fiz148, 17.10.2019.

#multiplot function for ggplot2. Note the layout=NULL and the option to apply n number of columns (cols=n) 
multiplot <- function(..., plotlist=NULL, file, cols=2, layout=NULL) {
  require(grid)
  require(gridExtra)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot 
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


setwd("~/Desktop/Final_euk_NN") #definer arbejdsfolder

library('mgcv') #indl?ser de n?dvendige pakker til analyse
library('picante')
library('ade4')
library('ggplot2')
library('RColorBrewer')
library('vegan')
library('ape')
library('MBI')
library('iNEXT')
library("devtools")
library("asbio")
library("grid")
library("MASS")
library("xlsx")
library("ncf")

#Environmental parameters
env<-read.table("vandv_env.txt", header=TRUE, row.names=1)

#Geographical matrices
#all 21 filter separately
geo.dist21<-read.table("geodist_21.txt",header=TRUE,row.names=1)
#11 plants/waterworks with pre- and postfilters combined
geo.dist_mantel<-read.table("geodist_mantel.txt",header=TRUE,row.names=1)

#Reading OTU matrices
#Bacterial - full 21 sample dataset at 3,5 and 10% clustering thresholds
BAKT03<- read.table("Bakt03_sidst.txt",header=TRUE,row.names=1)
BAKT05<- read.table("Bakt05_sidst.txt",header=TRUE,row.names=1)
BAKT10<- read.table("Bakt10_sidst.txt",header=TRUE,row.names=1)

#11 - with pre- and postfilter tests combined
BAKT03_mantel<- read.table("Bakt03_sidst_mantel.txt",header=TRUE,row.names=1)
BAKT05_mantel<- read.table("Bakt05_sidst_mantel.txt",header=TRUE,row.names=1)
BAKT10_mantel<- read.table("Bakt10_sidst_mantel.txt",header=TRUE,row.names=1)

#check rowsums
rowSums(BAKT03)
sum(rowSums(BAKT03))

#Eukaryotic
#full 21 filter dataset at 3,5 and 10% clustering thresholds
EUK03<- read.table("FINAL_0.03_euk_astT_shared.txt",header=TRUE,row.names=1)
EUK05<- read.table("FINAL_0.05_euk_astT_shared.txt",header=TRUE,row.names=1)
EUK10<- read.table("FINAL_0.10_euk_astT_shared",header=TRUE,row.names=1)

#11 - with pre- and postfilter tests combined
EUK03_mantel<- read.table("FINAL_0.03_euk_astT_shared_mantel.txt",header=TRUE,row.names=1)
EUK05_mantel<- read.table("FINAL_0.05_euk_astT_shared_mantel.txt",header=TRUE,row.names=1)
EUK10_mantel<- read.table("FINAL_0.10_euk_astT_shared_mantel.txt",header=TRUE,row.names=1)

#check rowsums
rowSums(EUK03)
sum(rowSums(EUK03))

#BC-distances
BAKT03_mantel.dist<-vegdist(BAKT03_mantel,method="bray")
BAKT05_mantel.dist<-vegdist(BAKT05_mantel,method="bray")
BAKT10_mantel.dist<-vegdist(BAKT10_mantel,method="bray")

EUK03_mantel.dist<-vegdist(EUK03_mantel,method="bray")
EUK05_mantel.dist<-vegdist(EUK05_mantel,method="bray")
EUK10_mantel.dist<-vegdist(EUK10_mantel,method="bray")

#Mantel tests on 11 plants/combined filters
mantel(BAKT03_mantel.dist, geo.dist_mantel)
mantel(BAKT05_mantel.dist, geo.dist_mantel)
mantel(BAKT10_mantel.dist, geo.dist_mantel)

mantel(EUK03_mantel.dist, geo.dist_mantel)
mantel(EUK05_mantel.dist, geo.dist_mantel)
mantel(EUK10_mantel.dist, geo.dist_mantel)

#check coverage
#transpose table to vegan/iNEXT-compatible format
b03<-t(BAKT03)
b05<-t(BAKT05)
b10<-t(BAKT10)

e03<-t(EUK03)
e05<-t(EUK05)
e10<-t(EUK10)

#iNEXT sample coverage calculations (Chao & Jost 2012 - coverage-based rarefaction)
rich_b3_cov_all<-estimateD(b03, datatype="abundance", base ="coverage")


rich_e03_cov_all<-estimateD(e03, datatype="abundance", base ="coverage")
rich_e3<-estimateD(e03, datatype="abundance", base ="coverage", level=0.973)
rich_e5<-estimateD(e05, datatype="abundance", base ="coverage", level=0.973)
rich_e10<-estimateD(e10, datatype="abundance", base ="coverage", level=0.973)

#It turns out that OdsD + LinD in bacteria and SveD+T in eukaryotes pulls down coverage (to <90) at all thresholds, left out

EUK03_u_sve<-EUK03[-c(12,13),]
BAKT03_u_OdsDLinD<-BAKT03[-c(7,10),] 

#Finding coverage levels in reduced datasets
#bacteria
b03_uOdsDLinD<-t(BAKT03_u_OdsDLinD)
bakt03_cov_uOdsDLinD<-estimateD(b03_uOdsDLinD, datatype="abundance", base ="coverage")
#eukaryotes
e03_uSve<-t(EUK03_u_sve)
euk03_cov_uSve<-estimateD(e03_uSve, datatype="abundance", base ="coverage")

#apply the coverage threshold for the lowest sample after removal of below pars (96.9% - acceptable)
rich_b3<-estimateD(b03, datatype="abundance", base ="coverage", level=0.969)
rich_b5<-estimateD(b05, datatype="abundance", base ="coverage", level=0.969)
rich_b10<-estimateD(b10, datatype="abundance", base ="coverage", level=0.969)

#reducing the different geographic matrices to match 97 coverage
#for bacteria 
geo.dist19_u_OdsDLinD<-geo.dist21[-c(7, 10), -c(7, 10)]
#for eukaryotes
geo.dist19_u_sve<-geo.dist21[-c(12,13), -c(12,13)]

#reading in environnemental variables
envB<-env[-c(7,10),]
envE<-env[-c(12,13),]

#Calculating BC distances in reduced datasets
#bacteria
BAKT03_u_OdsDLinD.dist<-vegdist(BAKT03_u_OdsDLinD,method="bray")
#eukaryotes
EUK03_u_sve.dist<-vegdist(EUK03_u_sve,method="bray")

#Mantel tests on reduced data for all 21 (19) filters
#bacteria
mantel(BAKT03_u_OdsDLinD.dist, geo.dist19_u_OdsDLinD)
#eukaryotes
mantel(EUK03_u_sve.dist, geo.dist19_u_sve)

#mantel mildly significant for bacteria (p=0,039), trying to leave out Astrup

BAKT03_u_OdsDLinDAst<-BAKT03[-c(3,4,7,10),] 
bakt03_u_OdsDLinDAst.dist<-vegdist(BAKT03_u_OdsDLinDAst,method="bray")

geo.dist19_u_OdsDLinDAst<-geo.dist21[-c(1,2,12,13), -c(1,2,12,13)]

mantel(bakt03_u_OdsDLinDAst.dist, geo.dist19_u_OdsDLinDAst)
#Mantel no longer significant (p=0.08), geographic effects a result of the deviant Astrup works. Proceeding with this reduced dataset of 19

#testing presence-absence, replacing abundance with 1:0
BAKT03_u_OdsDLinD[BAKT03_u_OdsDLinD > 0] <- 1
BAKT03_u_OdsDLinD

EUK03_preab<-EUK03_u_sve[EUK03_u_sve > 0] <- 1
EUK03_preab<-EUK03_u_sve

#Bakt03_preab.dist<-vegdist(Bakt03_preab,method="bray")
#EUK03_preab.dist<-vegdist(EUK03_preab,method="bray")

#mantel(Bakt03_preab.dist, geo.dist19_u_OdsDLinD)
#mantel(EUK03_preab.dist, geo.dist19_u_sve)

#rereading original files
EUK03_u_sve<-EUK03[-c(12,13),]
BAKT03_u_OdsDLinD<-BAKT03[-c(7,10),] 

#calculating alpha diversity

alpha_div_EUK03<-alpha.div(EUK03_u_sve,index="shan")
alpha_div_BAKT03<-alpha.div(BAKT03_u_OdsDLinD,index="shan")

alpha_div_EUK03<-as.data.frame(alpha_div_EUK03)
alpha_div_BAKT03<-as.data.frame(alpha_div_BAKT03)

#make a combined diversity file

rich_b3<-rich_b3[-c(7,10),]
rich_b5<-rich_b5[-c(7,10),]
rich_b10<-rich_b10[-c(7,10),]

rich_e3<-rich_e3[-c(12,13),]
rich_e5<-rich_e5[-c(12,13),]
rich_e10<-rich_e10[-c(12,13),]

BAKT03_rich_div<-cbind(row.names(alpha_div_BAKT03), envB$age,alpha_div_BAKT03$alpha_div_BAKT03,bakt03_cov_uOdsDLinD$"q = 0",rich_b3$"q = 0",rich_b5$"q = 0",rich_b10$"q = 0")
EUK03_rich_div<-cbind(row.names(alpha_div_EUK03),envE$age,alpha_div_EUK03$alpha_div_EUK03,euk03_cov_uSve$"q = 0",rich_e3$"q = 0",rich_e5$"q = 0",rich_e10$"q = 0")

colnames(BAKT03_rich_div)<-c("Group","age","BAKT03_alphadiv",	"Bakt_rich_all","Bakt03_rich","Bakt05_rich","BAKT03_rich")
colnames(EUK03_rich_div)<-c("Group","age","EUK03_alphadiv",	"Euk_rich_all","Euk03_rich","Euk05_rich","EUK03_rich")

#NMDS analysis
BAKT03_NMDS<-metaMDS(BAKT03_u_OdsDLinD, distance = "bray", k = 2, wascores = TRUE) #laver en Bray-Curtis-matrix og NMDS-analysen p? ?n gang
EUK03_NMDS<-metaMDS(EUK03_u_sve, distance = "bray", k = 2, wascores = TRUE) #laver en Bray-Curtis-matrix og NMDS-analysen p? ?n gang

#define a Regions vector of length 19
Region=c(rep("West Jutland",4),rep("Fyn/East Jutland",9),rep("Sjaelland",6))

# Use function chull() to derive hull vertices
hulls1.df = do.call(rbind, lapply(unique(Region), function(i) {
  data = BAKT03_NMDS$points[which(Region == i), ]
  data.frame(
    Region = i,
    data[chull(data[, 1], data[, 2]),],
    row.names = NULL
  )
} ) )

hulls2.df = do.call(rbind, lapply(unique(Region), function(i) {
  data = EUK03_NMDS$points[which(Region == i), ]
  data.frame(
    Region = i,
    data[chull(data[, 1], data[, 2]),],
    row.names = NULL
  )
} ) )

#reading in age matrix 
agedist<-read.table("agedist_21.txt",header=TRUE,row.names=1)
agedist_mantel<-read.table("agedist_mantel.txt",header=TRUE,row.names=1)
agedist_mantel<-as.dist(agedist_mantel)

agedistB<-as.dist(agedist[-c(7,10), -c(7,10),])
agedistE<-as.dist(agedist[-c(12,13), -c(12,13)])

mantel(EUK03_u_sve.dist, agedistE)
mantel(BAKT03_u_OdsDLinD.dist, agedistB)

#Mantel tests on compound AGE data
mantel(EUK03_mantel.dist, agedist_mantel)
mantel(BAKT03_mantel.dist, agedist_mantel)


#standard envfit of environnemental variables to the NMDS 
BAKT03_NMDS.envfit <- envfit(BAKT03_NMDS, env = envB, perm = 999) 
BAKT03_NMDS.envfit
EUK03_NMDS.envfit <- envfit(EUK03_NMDS, env = envE, perm = 999) 
EUK03_NMDS.envfit

Richness<-estimateR(EUK03)

#Plotting only significant env.var.
envB_sf<-
BAKT03sf_NMDS.envfit <- envfit(BAKT03_NMDS, env = env_sf, perm = 999) 
BAKT03sf_NMDS.envfit

envE_sf<-
EUK03sf_NMDS.envfit <- envfit(EUK03_NMDS, env = envE_sf, perm = 999) 
EUK03sf_NMDS.envfit

# data for the envfit arrows
env.scores.BAKT03 <- as.data.frame(scores(BAKT03_NMDS.envfit, display = "vectors")) #extracts relevant scores from envifit
env.scores.BAKT03 <- cbind(env.scores.BAKT03, env.variables = rownames(env.scores.BAKT03)) #and then gives them their names

env.scores.EUK03 <- as.data.frame(scores(EUK03_NMDS.envfit, display = "vectors")) #extracts relevant scores from envifit
env.scores.EUK03 <- cbind(env.scores.EUK03, env.variables = rownames(env.scores.EUK03)) #and then gives them their names

mult <- 2 #multiplier for the arrows and text for envfit below. You can change this and then rerun the plot command.

p1<-ggplot(
  # Create data object including treatment as a column
  data = cbind(comm = rownames(BAKT03_NMDS$points), treat = Region, as.data.frame(BAKT03_NMDS$points)),
  aes(x = MDS1, y = MDS2, xlab=""))+
  # Create polygon from hull dataframe
  geom_polygon(data = hulls1.df, aes(x = MDS1, y = MDS2, group = Region, fill = Region), alpha = 0.5, lwd = 0.8) +
  scale_fill_manual(values=c("darkorchid4","goldenrod","skyblue1")) +
  # Use aesthetics to color points by treatment
  theme_bw() +
  geom_point(aes(guide=FALSE), size = 2) +
  geom_segment(data = env.scores.BAKT03,
               aes(x = 0, xend = mult*NMDS1, y = 0, yend = mult*NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") + #arrows for envfit.  doubled the length for similarity to the plot() function. NB check ?envfit regarding arrow length if not familiar with lengths
  geom_text(data = env.scores.BAKT03, #labels the environmental variable arrows * "mult" as for the arrows
            aes(x = mult*NMDS1, y = mult*NMDS2, label=env.variables),
            size = 4,fontface="italic",
            hjust = 0.5) +
  geom_text(aes(label = comm), vjust = 1.5, size=5) +
  coord_cartesian(xlim=c(-2,2), ylim=c(-1.3,1.3)) +
  labs(title="Prokaryotes (0.10)", x="") +
  theme(plot.margin=unit(c(0,25,0,0), "mm"))

p2<-ggplot(
  # Create data object including treatment as a column
  data = cbind(comm = rownames(EUK03_NMDS$points), treat = Region, as.data.frame(EUK03_NMDS$points)),
  aes(x = MDS1, y = MDS2, xlab=""))+
  # Create polygon from hull dataframe
  geom_polygon(data = hulls2.df, aes(x = MDS1, y = MDS2, group = Region, fill = Region), alpha = 0.5, lwd = 0.8) +
  scale_fill_manual(values=c("darkorchid4","goldenrod","skyblue1")) + 
  # Use aesthetics to color points by treatment
  theme_bw() +
  geom_point(aes(guide=FALSE), size = 2) +
  geom_segment(data = env.scores.EUK03,
               aes(x = 0, xend = mult*NMDS1, y = 0, yend = mult*NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") + #arrows for envfit.  doubled the length for similarity to the plot() function. NB check ?envfit regarding arrow length if not familiar with lengths
  geom_text(data = env.scores.EUK03, #labels the environmental variable arrows * "mult" as for the arrows
            aes(x = mult*NMDS1, y = mult*NMDS2, label=env.variables),
            size = 4, fontface="italic",
            vjust = 0.5) +
  geom_text(aes(label = comm), vjust = 1.5, size=5) +
  coord_cartesian(xlim=c(-2,2), ylim=c(-1.3,1.3)) +
  labs(title="Eukaryotes (0.10)", x="MDS1", y="") +
  theme(plot.margin=unit(c(0,-5,0,-20), "mm"))

multiplot(p1+guides(fill=FALSE,col=FALSE),p2+guides(col=FALSE))

#Without vectors
p3<-ggplot(
  # Create data object including treatment as a column
  data = cbind(comm = rownames(BAKT03_NMDS$points), treat = Region, as.data.frame(BAKT03_NMDS$points)),
  aes(x = MDS1, y = MDS2, xlab=""))+
  # Create polygon from hull dataframe
  geom_polygon(data = hulls1.df, aes(x = MDS1, y = MDS2, group = Region, fill = Region), alpha = 0.5, lwd = 0.8) +
  scale_fill_manual(values=c("darkorchid4","goldenrod","skyblue1")) +
  # Use aesthetics to color points by treatment
  theme_bw() +
  geom_point(aes(guide=FALSE), size = 2) +
  geom_text(aes(label = comm), vjust = 1.5, size=5) +
  coord_cartesian(xlim=c(-2,2), ylim=c(-1.3,1.3)) +
  labs(title="Prokaryotes (0.03)", x="") +
  theme(plot.margin=unit(c(0,25,0,0), "mm"))

#p4<-
ggplot(
  # Create data object including treatment as a column
  data = cbind(comm = rownames(EUK03_NMDS$points), treat = Region, as.data.frame(EUK03_NMDS$points)),
  aes(x = MDS1, y = MDS2, xlab=""))+
  # Create polygon from hull dataframe
  geom_polygon(data = hulls2.df, aes(x = MDS1, y = MDS2, group = Region, fill = Region), alpha = 0.5, lwd = 0.8) +
  scale_fill_manual(values=c("darkorchid4","goldenrod","skyblue1")) + 
  # Use aesthetics to color points by treatment
  theme_bw() +
  geom_point(aes(guide=FALSE), size = 2) +
  geom_text(aes(label = comm), vjust = 1.5, size=5) +
  coord_cartesian(xlim=c(-2.25,2.45), ylim=c(-1.5,1.5)) +
  labs(title="NMDS(Eukaryotes) ", x="MDS1", y="MDS2") #+
  #theme(plot.margin=unit(c(0,-5,0,-20), "mm"))

#multiplot(p3+guides(fill=FALSE,col=FALSE),p4+guides(col=FALSE))

multiplot(p4+guides(col=FALSE))

#Check without Western Jutland
BAKT03_uVj<-BAKT03_u_OdsDLinD[-c(1,2,3,4),]
  
#BAKT03_uVj_preab<-read.table("lea_0.10_fn.shared_uVj_preab.txt", header=T, row.names=1)
geo_uVj_uOdsDLinD<-geo.dist19_u_OdsDLinD[-c(1,2,3,4), -c(1,2,3,4)]
  
BAKT03_uVj.dist<-vegdist(BAKT03_uVj, method="bray")
#BAKT_uVj_preab.dist<-vegdist(BAKT03_uVj_preab, method="bray")
mantel(BAKT03_uVj.dist, geo_uVj_uOdsDLinD)
#mantel(BAKT_uVj_preab.dist, geo_uVj_uOdsDLinD)
envB_uVj<-envB[-c(1,2,3,4),]
BAKT03_uVj_NMDS<-metaMDS(BAKT03_uVj, distance = "bray", k = 2, wascores = TRUE) #laver en Bray-Curtis-matrix og NMDS-analysen p? ?n gang

EUK03_uVj<-EUK03_u_sve[-c(1,2,3,4), ]
geo_uVj_u_Sve<-geo.dist19_u_sve[-c(1,2,3,4), -c(1,2,3,4)]
EUK03_uVj.dist<-vegdist(EUK03_uVj, method="bray")
mantel(EUK03_uVj.dist,geo_uVj_u_Sve)
envE_uVj<-envE[-c(1,2,3,4),]
EUK03_uVj_NMDS<-metaMDS(EUK03_uVj, distance = "bray", k = 2, wascores = TRUE) #laver en Bray-Curtis-matrix og NMDS-analysen p? ?n gang

#standard envfit of environnemental variables to the NMDS 
BAKT03_uVj_NMDS.envfit <- envfit(BAKT03_uVj_NMDS, env = envB_uVj, perm = 999) 
BAKT03_uVj_NMDS.envfit
EUK03_uVj_NMDS.envfit <- envfit(EUK03_uVj_NMDS, env = envE_uVj, perm = 999) 
EUK03_uVj_NMDS.envfit
EUK03_uVj_NMDS.envfit<-as.matrix(EUK03_uVj_NMDS.envfit)

#Check without Sj??lland
BAKT03_uSj<-BAKT03_u_OdsDLinD[-c(14,15,16,17,18,19),]
#BAKT03_uSj[BAKT03_uSj > 0] <- 1
#BAKT03_uSj_preab<-BAKT03_uSj
#BAKT03_uSj<-BAKT03_u_OdsDLinD[-c(14,15,16,17,18,19),]

geo_uSj_uOdsDLinD<-geo.dist19_u_OdsDLinD[-c(14,15,16,17,18,19), -c(14,15,16,17,18,19)]
BAKT_uSj.dist<-vegdist(BAKT03_uSj, method="bray")
#BAKT_uSj_preab.dist<-vegdist(BAKT03_uSj_preab, method="bray")
mantel(BAKT_uSj.dist, geo_uSj_uOdsDLinD)
#mantel(BAKT_uSj_preab.dist, geo_uSj_uOdsDLinD)
envB_uSj<-envB[-c(14,15,16,17,18,19),]
BAKT03_uSj_NMDS<-metaMDS(BAKT03_uSj, distance = "bray", k = 2, wascores = TRUE) #laver en Bray-Curtis-matrix og NMDS-analysen p? ?n gang

EUK03_uSj<-EUK03_u_sve[-c(14,15,16,17,18,19),]
geo_uSj_u_Sve<-geo.dist19_u_sve[-c(14,15,16,17,18,19), -c(14,15,16,17,18,19)]
EUK03_uSj.dist<-vegdist(EUK03_uSj, method="bray")
mantel(EUK03_uSj.dist,geo_uSj_u_Sve)
envE_uSj<-envE[-c(14,15,16,17,18,19),]
EUK03_uSj_NMDS<-metaMDS(EUK03_uSj, distance = "bray", k = 2, wascores = TRUE) #laver en Bray-Curtis-matrix og NMDS-analysen p? ?n gang

#standard envfit of environnemental variables to the NMDS 
BAKT03_uSj_NMDS.envfit <- envfit(BAKT03_uSj_NMDS, env = envB_uSj, perm = 999) 
BAKT03_uSj_NMDS.envfit
EUK03_uSj_NMDS.envfit <- envfit(EUK03_uSj_NMDS, env = envE_uSj, perm = 999) 
EUK03_uSj_NMDS.envfit

#Check without Fyn/EJylland
BAKT03_uF<-BAKT03_u_OdsDLinD[-c(5,6,7,8,9,10,11,12,13),]
geo_uF<-geo.dist19_u_sve[-c(5,6,7,8,9,10,11,12,13),-c(5,6,7,8,9,10,11,12,13)]
BAKT03_uF.dist<-vegdist(BAKT03_uF, method="bray")
mantel(BAKT03_uF.dist, geo_uF)
env_uF<-env[-c(5,6,7,8,9,10,11,12,13,14,15),]
BAKT03_uF_NMDS<-metaMDS(BAKT03_uF, distance = "bray", k = 2, wascores = TRUE) #laver en Bray-Curtis-matrix og NMDS-analysen p? ?n gang

EUK03_uF<-EUK03_u_sve[-c(5,6,7,8,9,10,11,12,13),]
EUK03_uF.dist<-vegdist(EUK03_uF, method="bray")
mantel(EUK03_uF.dist,geo_uF)
EUK03_uF_NMDS<-metaMDS(EUK03_uF, distance = "bray", k = 2, wascores = TRUE) #laver en Bray-Curtis-matrix og NMDS-analysen p? ?n gang

#standard envfit of environnemental variables to the NMDS 
BAKT03_uF_NMDS.envfit <- envfit(BAKT03_uF_NMDS, env = env_uF, perm = 999) 
BAKT03_uF_NMDS.envfit
EUK03_uF_NMDS.envfit <- envfit(EUK03_uF_NMDS, env = env_uF, perm = 999) 
EUK03_uF_NMDS.envfit

#merge envfit-lists and compare factors correlating significantly
#merge list function

BAKT03_NMDS.list<-mapply(BAKT03_NMDS.envfit,BAKT03_uVj_NMDS.envfit,BAKT03_uF_NMDS.envfit,BAKT03_uSj_NMDS.envfit,FUN=list,SIMPLIFY=FALSE)
EUK03_NMDS.list<-mapply(EUK03_NMDS.envfit,EUK03_uVj_NMDS.envfit,EUK03_uF_NMDS.envfit,EUK03_uSj_NMDS.envfit,FUN=list,SIMPLIFY=FALSE)

#check whether the significant NMDS vectors in the eukaryotes without Western Jutland correlate with distance (longitude)
model_calcium<-lm(envE_uVj$Caw ~ envE_uVj$pH)
summary(model_calcium)

model_age_Long<-lm(envE_uSj$age ~ envE_uSj$Longitude)
summary(model_age_Long)

model_BW_Long<-lm(envE_uVj$BW ~ envE_uVj$Longitude)
summary(model_BW_Long)

model_BET_Long<-lm(envE_uVj$BET ~ envE_uVj$Longitude)
summary(model_BET_Long)

model_Cond_Long<-lm(envE_uVj$Cond ~ envE_uVj$Longitude)
summary(model_Cond_Long)

model_Mgw_Long<-lm(envE_uVj$Mgw ~ envE_uVj$Longitude)
summary(model_Mgw_Long)

model_Naw_Long<-lm(envE_uVj$Naw ~ envE_uVj$Longitude)
summary(model_Naw_Long)

model_NO3w_Long<-lm(envE_uVj$NO3w ~ envE_uVj$Longitude)
summary(model_NO3w_Long)

model_SO4w_Long<-lm(envE_uVj$SO4w ~ envE_uVj$Longitude)
summary(model_SO4w_Long)

model_Few_Long<-lm(envE_uVj$Few ~ envE_uVj$Longitude)
summary(model_Few_Long)

model_NH4_Long<-lm(envE_uVj$NH4 ~ envE_uVj$Longitude)
summary(model_NH4_Long)

model_NVOC_Long<-lm(envE_uVj$NVOC ~ envE_uVj$Longitude)
summary(model_NVOC_Long)

#check whether the significant NMDS vectors in the eukaryotes without Sj??lland correlate w. longitude
model_calcium<-lm(envE_uSj$Caw ~ envE_uSj$pH)
summary(model_calcium)

model_BW_Long<-lm(envE_uSj$BW ~ envE_uSj$Longitude)
summary(model_BW_Long)

model_BET_Long<-lm(envE_uSj$BET ~ envE_uSj$Longitude)
summary(model_BET_Long)

model_Cond_Long<-lm(envE_uSj$Cond ~ envE_uSj$Longitude)
summary(model_Cond_Long)

model_Mgw_Long<-lm(envE_uSj$Mgw ~ envE_uSj$Longitude)
summary(model_Mgw_Long)

model_Naw_Long<-lm(envE_uSj$Naw ~ envE_uSj$Longitude)
summary(model_Naw_Long)

model_NO3w_Long<-lm(envE_uSj$NO3w ~ envE_uSj$Longitude)
summary(model_NO3w_Long)

model_NVOC_Long<-lm(envE_uSj$NVOC ~ envE_uSj$Longitude)
summary(model_NVOC_Long)

model_SO4w_Long<-lm(envE_uSj$SO4w ~ envE_uSj$Longitude)
summary(model_SO4w_Long)

model_Few_Long<-lm(envE_uSj$Few ~ envE_uSj$Longitude)
summary(model_Few_Long)

model_NH4_Long<-lm(envE_uSj$NH4 ~ envE_uSj$Longitude)
summary(model_NH4_Long)


#Plot_age (from models 9 +  13)
x1<-as.numeric(BAKT03_rich_div[,2])
y1<-as.numeric(BAKT03_rich_div[,5])
x2<-as.numeric(EUK03_rich_div[,2])
y2<-as.numeric(EUK03_rich_div[,5])
is.numeric(x1)
is.numeric(y1)
xy1 <- data.frame(x1, y1)
xy2 <- data.frame(x2, y2)

m1<-lm(y1~x1, data=xy1)
m2<-lm(y2~x2, data=xy2)

summary(m1)
summary(m2)

#y3<-as.numeric(BAKT03_rich_div[,3])
#y4<-as.numeric(EUK03_rich_div[,3])

#xy3 <- data.frame(x1, y3)
#xy4 <- data.frame(x2, y4)

#m3<-lm(y3~x1, data=xy3)
#m4<-lm(y4~x2, data=xy4)

#summary(m3)
#summary(m4)


# Plot age/richness correlations

p1 <- ggplot(xy1, aes(x=x1, y=y1)) + theme_bw() +
  annotate("text",x=33, y=50,label="R2=0.06") + 
  annotate("text", x=33, y=37,label="p=0.17") +
  geom_point() + geom_smooth(aes(x=x1, y=y1), method="lm", colour="red") +
  coord_cartesian(ylim=c(0,310), xlim=c(0,41)) + 
  theme(plot.margin=unit(c(1,0,1,1), "mm")) 

p2 <- ggplot(xy1, aes(x=x2, y=y2)) + theme_bw() +
  annotate("text",x=33, y=15,label="R2=0.37") + 
  annotate("text", x=33, y=10,label="p=0.003") +
  geom_point() + geom_smooth(aes(x=x2, y=y2), method="lm", colour="red") +
  coord_cartesian(ylim=c(0,110), xlim=c(0,41)) + 
  theme(plot.margin=unit(c(1,2,1,1), "mm")) 

multiplot(p1 + labs(title = "Prokaryotes", x = "", y = "Richness of OTUs(3%)"), p2 + labs(title = "Eukaryotes",x = "Age of filter material", y = ""))

#p3 <- ggplot(xy3, aes(x=x1, y=y3)) + theme_bw() +
#  annotate("text",x=5, y=5,label="R2=0.05") + 
#  annotate("text", x=5, y=4.65,label="p=0.35") +
#  geom_point() + geom_smooth(aes(x=x1, y=y3), method="lm") +
#  coord_cartesian(ylim=c(0,5), xlim=c(0,41)) + 
#  theme(plot.margin=unit(c(2.5,1,0,4), "mm")) 

#p4 <- ggplot(xy4, aes(x=x2, y=y4)) + theme_bw() +
#  annotate("text",x=5, y=4,label="R2=0.44") + 
#  annotate("text", x=5, y=3.65,label="p<0.01") +
#  geom_point() + geom_smooth(aes(x=x2, y=y4), method="lm") +
#  coord_cartesian(ylim=c(0,4), xlim=c(0,41)) + 
#  theme(plot.margin=unit(c(2.5,1,0,2), "mm")) 

#multiplot(p1 + labs(title = "Prokaryotes", x = "", y = "Richness of OTUs(10%)"), p3 + labs(x = "Age of filter material", y = "alpha diversity(10%)"),p2 +labs(title = "Eukaryotes", x = "", y = ""),p4 +labs(x = "", y = ""))

#Plot distance-decay relationship: Get the lower triangle of each distance matrix and combine
BAKT03.distance.decay = data.frame(
  geo.dist19_u_OdsDLinD = as.matrix(geo.dist19_u_OdsDLinD)[lower.tri(geo.dist19_u_OdsDLinD)],
  BAKT03_u_OdsDLinD.dist = as.matrix(BAKT03_u_OdsDLinD.dist)[lower.tri(BAKT03_u_OdsDLinD.dist)]
)

EUK03.distance.decay = data.frame(
  geo.dist19_u_sve = as.matrix(geo.dist19_u_sve)[lower.tri(geo.dist19_u_sve)],
  EUK03_u_sve.dist = as.matrix(EUK03_u_sve.dist)[lower.tri(EUK03_u_sve.dist)]
)

#define variables for ggplots
x5<-log((BAKT03.distance.decay$geo.dist19_u_OdsDLinD)+1)
y5<-BAKT03.distance.decay$BAKT03_u_OdsDLinD.dist
x6<-log((EUK03.distance.decay$geo.dist19_u_sve)+1)
y6<-EUK03.distance.decay$EUK03_u_sve.dist

#create data frames
xy5 <- data.frame(x5, y5)
xy6 <- data.frame(x6, y6)

m5<-lm(y5~x5, data=xy5)
m6<-lm(y6~x6, data=xy6)

summary(m5)
summary(m6)

#plots
p5 <- ggplot(xy5, aes(x=x5, y=y5)) + theme_bw() +
  annotate("text",x=5, y=0.15,label="R2=0.03") + 
  annotate("text", x=5, y=0.10,label="p<0.05") +
  geom_point() + geom_smooth(aes(x=x5, y=y5), method="lm",colour="darkorchid") +
  coord_cartesian(ylim=c(0,1,1.01), xlim=c(0,6)) 

p6 <- ggplot(xy6, aes(x=x6, y=y6)) + theme_bw() +
  annotate("text",x=5, y=0.15,label="R2=0.33") + 
  annotate("text", x=5, y=0.10,label="p<0.00001") +
  geom_point() + geom_smooth(aes(x=x6, y=y6), method="lm",colour="darkorchid") +
  coord_cartesian(ylim=c(0,1.01), xlim=c(0,6))

multiplot(p5 + labs(title = "Prokaryotes", x = "log(Distance(km)+1)", y = "Bray-Curtis dissimilarities"),p6 + labs(title = "Eukaryotes", x = "", y = ""))

#Delete functions with 0 - correct for autocorrelation in pre/postfilters from same plants
BAKT03.distance_0.decay<-BAKT03.distance.decay[!(apply(BAKT03.distance.decay, 1, function(y) any(y == 0))),]
EUK03.distance_0.decay<-EUK03.distance.decay[!(apply(EUK03.distance.decay, 1, function(y) any(y == 0))),]

#new corrected defined variables for ggplots
x7<-((BAKT03.distance_0.decay$geo.dist19_u_OdsDLinD))
y7<-BAKT03.distance_0.decay$BAKT03_u_OdsDLinD.dist
x8<-((EUK03.distance_0.decay$geo.dist19_u_sve))
y8<-EUK03.distance_0.decay$EUK03_u_sve.dist

xy7 <- data.frame(x7, y7)
xy8 <- data.frame(x8, y8)

m7<-lm(y7~x7, data=xy7)
m8<-lm(y8~x8, data=xy8)

summary(m7)
summary(m8)

#ggplots of corrected variables
p7 <- ggplot(xy7, aes(x=x7, y=y7)) + theme_bw() +
  annotate("text",x=225, y=0.40,label="R2=0.01") + 
  annotate("text", x=215, y=0.37,label="p=0.12") +
  geom_point() + geom_smooth(aes(x=x7, y=y7), method="lm",colour="darkorchid") +
  coord_cartesian(ylim=c(0.35,1), xlim=c(1.9,252)) 

p8 <- ggplot(xy8, aes(x=x8, y=y8)) + theme_bw() +
  annotate("text",x=225, y=0.40,label="R2=0.23") + 
  annotate("text", x=215, y=0.37,label="p<0.00001") +
  geom_point() + geom_smooth(aes(x=x8, y=y8), method="lm",colour="darkorchid") +
  coord_cartesian(ylim=c(0.35,1), xlim=c(1.9,252))

multiplot(p7 + labs(title = "Prokaryotes", x = "Distance(km)", y = "Bray-Curtis dissimilarities"),p8 + labs(title = "Eukaryotes", x = "", y = ""))

#experiments w/ multiple regressions

model1 <-lm(env$BAKT03_rich ~ env$age) 
summary(model1)

model2 <-lm(env$Bakt05_rich ~ env$age)
summary(model2)

model3 <-lm(env$Bakt03_rich ~ env$age)
summary(model3)

model4 <-lm(env$Bakt_rich_all ~ env$age)
summary(model4)

model5 <-lm(env$EUK03_rich ~ env$age)
summary(model5)

model6 <-lm(env$Euk05_rich ~ env$age)
summary(model6)

model7 <-lm(env$Euk03_rich ~ env$age)
summary(model7)

model8 <-lm(env$Euk_rich_all ~ env$age)
summary(model8)

model9 <-lm(y1 ~ envB$age) 
summary(model9)

model10 <-lm(envB$Bakt05_rich ~ envB$age)
summary(model10)

model11 <-lm(envB$Bakt03_rich ~ envB$age)
summary(model11)

model12 <-lm(envB$Bakt_rich_all ~ envB$age)
summary(model12)

model13 <-lm(y2 ~ envE$age)
summary(model13)

model14 <-lm(envE$Euk05_rich ~ envE$age)
summary(model14)

model15 <-lm(envE$Euk03_rich ~ envE$age)
summary(model15)

model16 <-lm(envE$Euk_rich_all ~ envE$age)
summary(model16)

model_Euk_alfa <-lm(envE$EUK03_alphadiv ~ envE$age)
summary(model_Euk_alfa)

model_Bakt_alfa <-lm(envB$BAKT03_alphadiv ~ envE$age)
summary(model_Bakt_alfa)

#Eukaryote models

model_euk_age_long <-lm(envB$Longitude ~ envB$age)
summary(model_euk_age_long)

model_euk_age_bw <-lm(y2 ~ envE$age +envE$BW)
summary(model_euk_age_bw)

model_euk_age_bw <-lm(y2 ~ envE$age +envE$BW)
summary(model_euk_age_bw)

model_euk_age_bw_ph <-lm(y2 ~ envE$age +envE$BW + envE$pH)
summary(model_euk_age_bw_ph)

model_euk_age_bw_ph_long <-lm(y2 ~ envE$age +envE$BW + envE$pH + envE$Longitude)
summary(model_euk_age_bw_ph_long)

model_euk_age_bw_ph_long <-lm(y2 ~ envE$age +envE$BW + envE$pH + envE$Longitude)
summary(model_euk_age_bw_ph_long)

model_euk_all <-lm(y2 ~ envE$age +envE$BW + envE$pH + envE$Longitude + envE$Fe + envE$Mn + envE$BET + envE$dH + envE$Cond + envE$Caw + envE$NVOC + envE$Clw + envE$color + envE$O2w + envE$Kw + envE$Mgw + envE$Naw + envE$NO3w + envE$Pw + envE$SO4w + envE$Few + envE$Mnw + envE$NH4 + envE$CH4 + envE$Latitude)
summary(model_euk_all)

#final stepwise regression
step <- stepAIC(model_euk_all, direction="both")
step$anova # display results 

#Bacteria
model_Bakt_bw <-lm(y1 ~ envB$BW)
summary(model_Bakt_bw)

model_Bakt_color_bw <-lm(y1 ~ envB$BW+envB$color)
summary(model_Bakt_color_bw)

model_bakt_bedst <-lm(y1 ~ envB$age +envB$BW + envB$pH + envB$Mn + envB$NVOC)
summary(model_Bakt_bedst)

model_Bakt_age_bw_ph_long <-lm(y1 ~ envB$age +envB$BW + envB$pH + envB$Longitude)
summary(model_Bakt_age_bw_ph_long)

model_Bakt_age_ph_long <-lm(y1 ~ envB$age + envB$pH + envB$Longitude)
summary(model_Bakt_age_ph_long)

model_bakt_all <-lm(y1 ~ envB$age +envB$BW + envB$pH + envB$Longitude + envB$Fe + envB$Mn + envB$BET + envB$dH + envB$Cond + envB$Caw + envB$NVOC + envB$Clw + envB$color + envB$O2w + envB$Kw + envB$Mgw + envB$Naw + envB$NO3w + envB$Pw + envB$SO4w + envB$Few + envB$Mnw + envB$NH4 + envB$CH4 + envB$Latitude)
summary(model_bakt_all)

#final stepwise regression
step <- stepAIC(model_bakt_all, direction="both")
step$anova # display results 

#Check other possible variables explaining the OTU BC distances
#Plot age-decay relationship: Get the lower triangle of each distance matrix and combine
agedistB1<-as.matrix(dist(scale(envB[, c(1)])))
agedistE1<-as.matrix(dist(scale(envE[, c(1)])))

BAKT03.dissimilarity.age = data.frame(
  agedistB1 = as.matrix(agedistB1)[lower.tri(agedistB1)],
  BAKT03_u_OdsDLinD.dist = as.matrix(BAKT03_u_OdsDLinD.dist)[lower.tri(BAKT03_u_OdsDLinD.dist)]
)

EUK03.dissimilarity.age = data.frame(
  agedistE1 = as.matrix(agedistE1)[lower.tri(agedistE1)],
  EUK03_u_sve.dist = as.matrix(EUK03_u_sve.dist)[lower.tri(EUK03_u_sve.dist)]
)

#Delete functions with 0
BAKT03.age_0.decay<-BAKT03.dissimilarity.age[!(apply(BAKT03.dissimilarity.age, 1, function(y) any(y == 0))),]
EUK03.age_0.decay<-EUK03.dissimilarity.age[!(apply(EUK03.dissimilarity.age, 1, function(y) any(y == 0))),]


xbage<-((BAKT03.age_0.decay$agedistB1))
ybage<-BAKT03.age_0.decay$BAKT03_u_OdsDLinD.dist
xeage<-((EUK03.age_0.decay$agedistE1))
yeage<-EUK03.age_0.decay$EUK03_u_sve.dist

xybage <- data.frame(xbage, ybage)
xyeage <- data.frame(xeage, yeage)

mbage<-lm(ybage~log(xbage), data=xybage)
meage<-lm(yeage~xeage, data=xyeage)

summary(mbage)
summary(meage)

#Plot backwasching(bw)-decay relationship: Get the lower triangle of each distance matrix and combine
bwdistB1<-as.matrix(dist(scale(envB[, c(2)])))
bwdistE1<-as.matrix(dist(scale(envE[, c(2)])))

BAKT03.dissimilarity.bw = data.frame(
  bwdistB1 = as.matrix(bwdistB1)[lower.tri(bwdistB1)],
  BAKT03_u_OdsDLinD.dist = as.matrix(BAKT03_u_OdsDLinD.dist)[lower.tri(BAKT03_u_OdsDLinD.dist)]
)

EUK03.dissimilarity.bw = data.frame(
  bwdistE1 = as.matrix(bwdistE1)[lower.tri(bwdistE1)],
  EUK03_u_sve.dist = as.matrix(EUK03_u_sve.dist)[lower.tri(EUK03_u_sve.dist)]
)

#Delete functions with 0
BAKT03.bw_0.decay<-BAKT03.dissimilarity.bw[!(apply(BAKT03.dissimilarity.bw, 1, function(y) any(y == 0))),]
EUK03.bw_0.decay<-EUK03.dissimilarity.bw[!(apply(EUK03.dissimilarity.bw, 1, function(y) any(y == 0))),]

xbbw<-((BAKT03.bw_0.decay$bwdistB1))
ybbw<-BAKT03.bw_0.decay$BAKT03_u_OdsDLinD.dist
xebw<-((EUK03.bw_0.decay$bwdistE1))
yebw<-EUK03.bw_0.decay$EUK03_u_sve.dist

xybbw <- data.frame(xbbw, ybbw)
xyebw <- data.frame(xebw, yebw)

mbbw<-lm(ybbw~log(xbbw), data=xybbw)
mebw<-lm(yebw~xebw, data=xyebw)

summary(mbbw)
summary(mebw)#*

#Plot pH-decay relationship: Get the lower triangle of each distance matrix and combine
phdistB1<-as.matrix(dist(scale(envB[, c(6)])))
phdistE1<-as.matrix(dist(scale(envE[, c(6)])))

BAKT03.dissimilarity.ph = data.frame(
  phdistB1 = as.matrix(phdistB1)[lower.tri(phdistB1)],
  BAKT03_u_OdsDLinD.dist = as.matrix(BAKT03_u_OdsDLinD.dist)[lower.tri(BAKT03_u_OdsDLinD.dist)]
)

EUK03.dissimilarity.ph = data.frame(
  phdistE1 = as.matrix(phdistE1)[lower.tri(phdistE1)],
  EUK03_u_sve.dist = as.matrix(EUK03_u_sve.dist)[lower.tri(EUK03_u_sve.dist)]
)

#Delete functions with 0
BAKT03.ph_0.decay<-BAKT03.dissimilarity.ph[!(apply(BAKT03.dissimilarity.ph, 1, function(y) any(y == 0))),]
EUK03.ph_0.decay<-EUK03.dissimilarity.ph[!(apply(EUK03.dissimilarity.ph, 1, function(y) any(y == 0))),]

xbph<-((BAKT03.ph_0.decay$phdistB1))
ybph<-BAKT03.ph_0.decay$BAKT03_u_OdsDLinD.dist
xeph<-((EUK03.ph_0.decay$phdistE1))
yeph<-EUK03.ph_0.decay$EUK03_u_sve.dist

xybph <- data.frame(xbph, ybph)
xyeph <- data.frame(xeph, yeph)

mbph<-lm(ybph~log(xbph), data=xybph)
meph<-lm(yeph~xeph, data=xyeph)

summary(mbph)
summary(meph)

pbph <- ggplot(xybph, aes(x=xbph, y=ybph)) + theme_bw() +
  geom_point() + geom_smooth(aes(x=xbph, y=ybph), method="lm") + theme(plot.margin=unit(c(-2,0,1,0), "mm"))

peph <- ggplot(xyeph, aes(x=xeph, y=yeph)) + theme_bw() +
  geom_point() + geom_smooth(aes(x=xeph, y=yeph), method="lm") + theme(plot.margin=unit(c(1,0,-2,-1), "mm"))

#Plot water hardness(dh)-decay relationship: Get the lower triangle of each distance matrix and combine
dhdistB1<-as.matrix(dist(scale(envB[, c(7)])))
dhdistE1<-as.matrix(dist(scale(envE[, c(7)])))

BAKT03.dissimilarity.dh = data.frame(
  dhdistB1 = as.matrix(dhdistB1)[lower.tri(dhdistB1)],
  BAKT03_u_OdsDLinD.dist = as.matrix(BAKT03_u_OdsDLinD.dist)[lower.tri(BAKT03_u_OdsDLinD.dist)]
)

EUK03.dissimilarity.dh = data.frame(
  dhdistE1 = as.matrix(dhdistE1)[lower.tri(dhdistE1)],
  EUK03_u_sve.dist = as.matrix(EUK03_u_sve.dist)[lower.tri(EUK03_u_sve.dist)]
)

#Delete functions with 0
BAKT03.dh_0.decay<-BAKT03.dissimilarity.dh[!(apply(BAKT03.dissimilarity.dh, 1, function(y) any(y == 0))),]
EUK03.dh_0.decay<-EUK03.dissimilarity.dh[!(apply(EUK03.dissimilarity.dh, 1, function(y) any(y == 0))),]

xbdh<-((BAKT03.dh_0.decay$dhdistB1))
ybdh<-BAKT03.dh_0.decay$BAKT03_u_OdsDLinD.dist
xedh<-((EUK03.dh_0.decay$dhdistE1))
yedh<-EUK03.dh_0.decay$EUK03_u_sve.dist

xybdh <- data.frame(xbdh, ybdh)
xyedh <- data.frame(xedh, yedh)

mbdh<-lm(ybdh~log(xbdh), data=xybdh)
medh<-lm(yedh~xedh, data=xyedh)

summary(mbdh)
summary(medh)#*

pbdh <- ggplot(xybdh, aes(x=xbdh, y=ybdh)) + theme_bw() +
  geom_point() + geom_smooth(aes(x=xbdh, y=ybdh), method="lm") + theme(plot.margin=unit(c(-2,0,1,0), "mm"))

pedh <- ggplot(xyedh, aes(x=xedh, y=yedh)) + theme_bw() +
  geom_point() + geom_smooth(aes(x=xedh, y=yedh), method="lm") + theme(plot.margin=unit(c(1,0,-2,-1), "mm"))

#Plot O2-decay relationship: Get the lower triangle of each distance matrix and combine
o2distB1<-as.matrix(dist(scale(envB[, c(13)])))
o2distE1<-as.matrix(dist(scale(envE[, c(13)])))

BAKT03.dissimilarity.o2 = data.frame(
  o2distB1 = as.matrix(o2distB1)[lower.tri(o2distB1)],
  BAKT03_u_OdsDLinD.dist = as.matrix(BAKT03_u_OdsDLinD.dist)[lower.tri(BAKT03_u_OdsDLinD.dist)]
)

EUK03.dissimilarity.o2 = data.frame(
  o2distE1 = as.matrix(o2distE1)[lower.tri(o2distE1)],
  EUK03_u_sve.dist = as.matrix(EUK03_u_sve.dist)[lower.tri(EUK03_u_sve.dist)]
)

#Delete functions with 0
BAKT03.o2_0.decay<-BAKT03.dissimilarity.o2[!(apply(BAKT03.dissimilarity.o2, 1, function(y) any(y == 0))),]
EUK03.o2_0.decay<-EUK03.dissimilarity.o2[!(apply(EUK03.dissimilarity.o2, 1, function(y) any(y == 0))),]

xbo2<-((BAKT03.o2_0.decay$o2distB1))
ybo2<-BAKT03.o2_0.decay$BAKT03_u_OdsDLinD.dist
xeo2<-((EUK03.o2_0.decay$o2distE1))
yeo2<-EUK03.o2_0.decay$EUK03_u_sve.dist

xybo2 <- data.frame(xbo2, ybo2)
xyeo2 <- data.frame(xeo2, yeo2)

mbo2<-lm(ybo2~log(xbo2), data=xybo2)
meo2<-lm(yeo2~xeo2, data=xyeo2)

summary(mbo2)
summary(meo2)

pbo2 <- ggplot(xybo2, aes(x=xbo2, y=ybo2)) + theme_bw() +
  geom_point() + geom_smooth(aes(x=xbo2, y=ybo2), method="lm") + theme(plot.margin=unit(c(-2,0,1,0), "mm"))

peo2 <- ggplot(xyeo2, aes(x=xeo2, y=yeo2)) + theme_bw() +
  geom_point() + geom_smooth(aes(x=xeo2, y=yeo2), method="lm") + theme(plot.margin=unit(c(1,0,-2,-1), "mm"))

#Plot no3-decay relationship: Get the lower triangle of each distance matrix and combine
no3distB1<-as.matrix(dist(scale(envB[, c(17)])))
no3distE1<-as.matrix(dist(scale(envE[, c(17)])))

BAKT03.dissimilarity.no3 = data.frame(
  no3distB1 = as.matrix(no3distB1)[lower.tri(no3distB1)],
  BAKT03_u_OdsDLinD.dist = as.matrix(BAKT03_u_OdsDLinD.dist)[lower.tri(BAKT03_u_OdsDLinD.dist)]
)

EUK03.dissimilarity.no3 = data.frame(
  no3distE1 = as.matrix(no3distE1)[lower.tri(no3distE1)],
  EUK03_u_sve.dist = as.matrix(EUK03_u_sve.dist)[lower.tri(EUK03_u_sve.dist)]
)

#Delete functions with 0
BAKT03.no3_0.decay<-BAKT03.dissimilarity.no3[!(apply(BAKT03.dissimilarity.no3, 1, function(y) any(y == 0))),]
EUK03.no3_0.decay<-EUK03.dissimilarity.no3[!(apply(EUK03.dissimilarity.no3, 1, function(y) any(y == 0))),]

xbno3<-((BAKT03.no3_0.decay$no3distB1))
ybno3<-BAKT03.no3_0.decay$BAKT03_u_OdsDLinD.dist
xeno3<-((EUK03.no3_0.decay$no3distE1))
yeno3<-EUK03.no3_0.decay$EUK03_u_sve.dist

xybno3 <- data.frame(xbno3, ybno3)
xyeno3 <- data.frame(xeno3, yeno3)

mbno3<-lm(ybno3~log(xbno3), data=xybno3)
meno3<-lm(yeno3~xeno3, data=xyeno3)

summary(mbno3)
summary(meno3)#*

pbno3 <- ggplot(xybno3, aes(x=xbno3, y=ybno3)) + theme_bw() +
  geom_point() + geom_smooth(aes(x=xbno3, y=ybno3), method="lm") + theme(plot.margin=unit(c(-2,0,1,0), "mm"))

peno3 <- ggplot(xyeno3, aes(x=xeno3, y=yeno3)) + theme_bw() +
  geom_point() + geom_smooth(aes(x=xeno3, y=yeno3), method="lm") + theme(plot.margin=unit(c(1,0,-2,-1), "mm"))


#Plot so4-decay relationship: Get the lower triangle of each distance matrix and combine
so4distB1<-as.matrix(dist(scale(envB[, c(19)])))
so4distE1<-as.matrix(dist(scale(envE[, c(19)])))

BAKT03.dissimilarity.so4 = data.frame(
  so4distB1 = as.matrix(so4distB1)[lower.tri(so4distB1)],
  BAKT03_u_OdsDLinD.dist = as.matrix(BAKT03_u_OdsDLinD.dist)[lower.tri(BAKT03_u_OdsDLinD.dist)]
)

EUK03.dissimilarity.so4 = data.frame(
  so4distE1 = as.matrix(so4distE1)[lower.tri(so4distE1)],
  EUK03_u_sve.dist = as.matrix(EUK03_u_sve.dist)[lower.tri(EUK03_u_sve.dist)]
)

#Delete functions with 0
BAKT03.so4_0.decay<-BAKT03.dissimilarity.so4[!(apply(BAKT03.dissimilarity.so4, 1, function(y) any(y == 0))),]
EUK03.so4_0.decay<-EUK03.dissimilarity.so4[!(apply(EUK03.dissimilarity.so4, 1, function(y) any(y == 0))),]

xbso4<-((BAKT03.so4_0.decay$so4distB1))
ybso4<-BAKT03.so4_0.decay$BAKT03_u_OdsDLinD.dist
xeso4<-((EUK03.so4_0.decay$so4distE1))
yeso4<-EUK03.so4_0.decay$EUK03_u_sve.dist

xybso4 <- data.frame(xbso4, ybso4)
xyeso4 <- data.frame(xeso4, yeso4)

mbso4<-lm(ybso4~log(xbso4), data=xybso4)
meso4<-lm(yeso4~xeso4, data=xyeso4)

summary(mbso4)
summary(meso4)

pbso4 <- ggplot(xybso4, aes(x=xbso4, y=ybso4)) + theme_bw() +
  geom_point() + geom_smooth(aes(x=xbso4, y=ybso4), method="lm") + theme(plot.margin=unit(c(-2,0,1,0), "mm"))

peso4 <- ggplot(xyeso4, aes(x=xeso4, y=yeso4)) + theme_bw() +
  geom_point() + geom_smooth(aes(x=xeso4, y=yeso4), method="lm") + theme(plot.margin=unit(c(1,0,-2,-1), "mm"))

#plot few
fewdistB1<-as.matrix(dist(scale(envB[, c(20)])))
fewdistE1<-as.matrix(dist(scale(envE[, c(20)])))

BAKT03.dissimilarity.few = data.frame(
  fewdistB1 = as.matrix(fewdistB1)[lower.tri(fewdistB1)],
  BAKT03_u_OdsDLinD.dist = as.matrix(BAKT03_u_OdsDLinD.dist)[lower.tri(BAKT03_u_OdsDLinD.dist)]
)

EUK03.dissimilarity.few = data.frame(
  fewdistE1 = as.matrix(fewdistE1)[lower.tri(fewdistE1)],
  EUK03_u_sve.dist = as.matrix(EUK03_u_sve.dist)[lower.tri(EUK03_u_sve.dist)]
)

#Delete functions with 0
BAKT03.few_0.decay<-BAKT03.dissimilarity.few[!(apply(BAKT03.dissimilarity.few, 1, function(y) any(y == 0))),]
EUK03.few_0.decay<-EUK03.dissimilarity.few[!(apply(EUK03.dissimilarity.few, 1, function(y) any(y == 0))),]

xbfew<-((BAKT03.few_0.decay$fewdistB1))
ybfew<-BAKT03.few_0.decay$BAKT03_u_OdsDLinD.dist
xefew<-((EUK03.few_0.decay$fewdistE1))
yefew<-EUK03.few_0.decay$EUK03_u_sve.dist

xybfew <- data.frame(xbfew, ybfew)
xyefew <- data.frame(xefew, yefew)

mbfew<-lm(ybfew~log(xbfew), data=xybfew)
mefew<-lm(yefew~xefew, data=xyefew)

summary(mbfew)
summary(mefew)#*

#plot mnw  
mnwdistB1<-as.matrix(dist(scale(envB[, c(21)])))
mnwdistE1<-as.matrix(dist(scale(envE[, c(21)])))

BAKT03.dissimilarity.mnw = data.frame(
  mnwdistB1 = as.matrix(mnwdistB1)[lower.tri(mnwdistB1)],
  BAKT03_u_OdsDLinD.dist = as.matrix(BAKT03_u_OdsDLinD.dist)[lower.tri(BAKT03_u_OdsDLinD.dist)]
)

EUK03.dissimilarity.mnw = data.frame(
  mnwdistE1 = as.matrix(mnwdistE1)[lower.tri(mnwdistE1)],
  EUK03_u_sve.dist = as.matrix(EUK03_u_sve.dist)[lower.tri(EUK03_u_sve.dist)]
)

#Delete functions with 0
BAKT03.mnw_0.decay<-BAKT03.dissimilarity.mnw[!(apply(BAKT03.dissimilarity.mnw, 1, function(y) any(y == 0))),]
EUK03.mnw_0.decay<-EUK03.dissimilarity.mnw[!(apply(EUK03.dissimilarity.mnw, 1, function(y) any(y == 0))),]

xbmnw<-((BAKT03.mnw_0.decay$mnwdistB1))
ybmnw<-BAKT03.mnw_0.decay$BAKT03_u_OdsDLinD.dist
xemnw<-((EUK03.mnw_0.decay$mnwdistE1))
yemnw<-EUK03.mnw_0.decay$EUK03_u_sve.dist

xybmnw <- data.frame(xbmnw, ybmnw)
xyemnw <- data.frame(xemnw, yemnw)

mbmnw<-lm(ybmnw~log(xbmnw), data=xybmnw)
memnw<-lm(yemnw~xemnw, data=xyemnw)

summary(mbmnw)
summary(memnw)#*

#plot nh4 
nh4distB1<-as.matrix(dist(scale(envB[, c(22)])))
nh4distE1<-as.matrix(dist(scale(envE[, c(22)])))

BAKT03.dissimilarity.nh4 = data.frame(
  nh4distB1 = as.matrix(nh4distB1)[lower.tri(nh4distB1)],
  BAKT03_u_OdsDLinD.dist = as.matrix(BAKT03_u_OdsDLinD.dist)[lower.tri(BAKT03_u_OdsDLinD.dist)]
)

EUK03.dissimilarity.nh4 = data.frame(
  nh4distE1 = as.matrix(nh4distE1)[lower.tri(nh4distE1)],
  EUK03_u_sve.dist = as.matrix(EUK03_u_sve.dist)[lower.tri(EUK03_u_sve.dist)]
)

#Delete functions with 0
BAKT03.nh4_0.decay<-BAKT03.dissimilarity.nh4[!(apply(BAKT03.dissimilarity.nh4, 1, function(y) any(y == 0))),]
EUK03.nh4_0.decay<-EUK03.dissimilarity.nh4[!(apply(EUK03.dissimilarity.nh4, 1, function(y) any(y == 0))),]

xbnh4<-((BAKT03.nh4_0.decay$nh4distB1))
ybnh4<-BAKT03.nh4_0.decay$BAKT03_u_OdsDLinD.dist
xenh4<-((EUK03.nh4_0.decay$nh4distE1))
yenh4<-EUK03.nh4_0.decay$EUK03_u_sve.dist

xybnh4 <- data.frame(xbnh4, ybnh4)
xyenh4 <- data.frame(xenh4, yenh4)

mbnh4<-lm(ybnh4~log(xbnh4), data=xybnh4)
menh4<-lm(yenh4~xenh4, data=xyenh4)

summary(mbnh4)
summary(menh4)#*

#plot latitude
latdistB1<-as.matrix(dist(scale(envB[, c(24)])))
latdistE1<-as.matrix(dist(scale(envE[, c(24)])))

BAKT03.dissimilarity.lat = data.frame(
  latdistB1 = as.matrix(latdistB1)[lower.tri(latdistB1)],
  BAKT03_u_OdsDLinD.dist = as.matrix(BAKT03_u_OdsDLinD.dist)[lower.tri(BAKT03_u_OdsDLinD.dist)]
)

EUK03.dissimilarity.lat = data.frame(
  latdistE1 = as.matrix(latdistE1)[lower.tri(latdistE1)],
  EUK03_u_sve.dist = as.matrix(EUK03_u_sve.dist)[lower.tri(EUK03_u_sve.dist)]
)

#Delete functions with 0
BAKT03.lat_0.decay<-BAKT03.dissimilarity.lat[!(apply(BAKT03.dissimilarity.lat, 1, function(y) any(y == 0))),]
EUK03.lat_0.decay<-EUK03.dissimilarity.lat[!(apply(EUK03.dissimilarity.lat, 1, function(y) any(y == 0))),]

xblat<-((BAKT03.lat_0.decay$latdistB1))
yblat<-BAKT03.lat_0.decay$BAKT03_u_OdsDLinD.dist
xelat<-((EUK03.lat_0.decay$latdistE1))
yelat<-EUK03.lat_0.decay$EUK03_u_sve.dist

xyblat <- data.frame(xblat, yblat)
xyelat <- data.frame(xelat, yelat)

mblat<-lm(yblat~log(xblat), data=xyblat)
melat<-lm(yelat~xelat, data=xyelat)

summary(mblat)
summary(melat)

#plot longitude
londistB1<-as.matrix(dist(scale(envB[, c(25)])))
londistE1<-as.matrix(dist(scale(envE[, c(25)])))

BAKT03.dissimilarity.lon = data.frame(
  londistB1 = as.matrix(londistB1)[lower.tri(londistB1)],
  BAKT03_u_OdsDLinD.dist = as.matrix(BAKT03_u_OdsDLinD.dist)[lower.tri(BAKT03_u_OdsDLinD.dist)]
)

EUK03.dissimilarity.lon = data.frame(
  londistE1 = as.matrix(londistE1)[lower.tri(londistE1)],
  EUK03_u_sve.dist = as.matrix(EUK03_u_sve.dist)[lower.tri(EUK03_u_sve.dist)]
)

#Delete functions with 0
BAKT03.lon_0.decay<-BAKT03.dissimilarity.lon[!(apply(BAKT03.dissimilarity.lon, 1, function(y) any(y == 0))),]
EUK03.lon_0.decay<-EUK03.dissimilarity.lon[!(apply(EUK03.dissimilarity.lon, 1, function(y) any(y == 0))),]

xblon<-((BAKT03.lon_0.decay$londistB1))
yblon<-BAKT03.lon_0.decay$BAKT03_u_OdsDLinD.dist
xelon<-((EUK03.lon_0.decay$londistE1))
yelon<-EUK03.lon_0.decay$EUK03_u_sve.dist

xyblon <- data.frame(xblon, yblon)
xyelon <- data.frame(xelon, yelon)

mblon<-lm(yblon~log(xblon), data=xyblon)
melon<-lm(yelon~xelon, data=xyelon)

summary(mblon)#*
summary(melon)#*

boxplot(BAKT03_u_OdsDLinD.dist, EUK03_u_sve.dist)

#Check bw relative to geographywithout Islev+S??nders?? (outliers BW)
BAKT03_uIs_Soen<-BAKT03_u_OdsDLinD[-c(16,17,18,19),]
geoB_uIs_Soen<-geo.distB[-c(16,17,18,19), -c(16,17,18,19)]
BAKT_uIs_Soen.dist<-vegdist(BAKT03_uIs_Soen, method="bray")
mantel(BAKT_uIslSoen.dist, geo_uSj_uIslSoen)
envB_uIs_Soen<-envB[-c(16,17,18,19),]
BAKT03__uIslSoen_NMDS<-metaMDS(BAKT03__uIslSoen, distance = "bray", k = 2, wascores = TRUE) #laver en Bray-Curtis-matrix og NMDS-analysen p? ?n gang

#eukaryotes
EUK03_uIs_Soen<-EUK03_u_sve[-c(16,17,18,19),]
geoE_uIs_Soen<-geo.distE[-c(16,17,18,19), -c(16,17,18,19)]
EUK03_uIs_Soen.dist<-vegdist(EUK03_uIs_Soen, method="bray")
mantel(EUK03_uIs_Soen.dist,geo_uIs_Soen)
envE_uIs_Soen<-envE[-c(16,17,18,19),]
EUK03_uIs_Soen_NMDS<-metaMDS(EUK03_uIs_Soen, distance = "bray", k = 2, wascores = TRUE) #laver en Bray-Curtis-matrix og NMDS-analysen p? ?n gang

envB_uIs_Soen.dist <- dist(scale(envB_uIs_Soen[, c(2)])) # here i select backwashing
mantel.partial(BAKT_uIs_Soen.dist, geoB_uIs_Soen, envB_uIs_Soen.dist, permutations=999) # influence of env on OTU matrix after partialing out space

envE_uIs_Soen.dist <- dist(scale(envE_uIs_Soen[, c(2)])) 
mantel.partial(EUK03_uIs_Soen.dist, geo_uIs_Soen,envE_uIs_Soen.dist,  permutations=999) # in

#Now check backwashing w/o deviant S??n/Isl
#first reading in backwashing matrix with +/- numbers and doing mantel tests
bwdist<-read.table("BW_dist_u_forf.txt",header=TRUE,row.names=1)

bwdistB_uIs_Soen<-as.dist(bwdist[-c(7,10,18,19,20,21), -c(7,10,18,19,20,21)])
bwdistE_uIs_Soen<-as.dist(bwdist[-c(12,13,18,19,20,21), -c(12,13,18,19,20,21)])

mantel(EUK03_uIs_Soen.dist, bwdistE_uIs_Soen)
mantel(BAKT_uIs_Soen.dist, bwdistB_uIs_Soen)

#Plot reduced isl/s?? dissimilarity-bw relationship: Get the lower triangle of each distance matrix and combine
BAKT03_uIS_Soen.dissimilariy.bw = data.frame(
  bwdistB_uIs_Soen = as.matrix(bwdistB_uIs_Soen)[lower.tri(bwdistB_uIs_Soen)],
  BAKT_uIs_Soen.dist = as.matrix(BAKT_uIs_Soen.dist)[lower.tri(BAKT_uIs_Soen.dist)]
)

EUK03_uIS_Soen.dissimilariy.bw = data.frame(
  bwdistE_uIs_Soen = as.matrix(bwdistE_uIs_Soen)[lower.tri(bwdistE_uIs_Soen)],
  EUK03_uIs_Soen.dist = as.matrix(EUK03_uIs_Soen.dist)[lower.tri(EUK03_uIs_Soen.dist)]
)

plot(BAKT03_uIS_Soen.dissimilariy.bw)
plot(EUK03_uIS_Soen.dissimilariy.bw) 

#again delete functions with 0
BAKT03_uIS_Soen.bw_0.decay<-BAKT03_uIS_Soen.dissimilariy.bw[!(apply(BAKT03_uIS_Soen.dissimilariy.bw, 1, function(y) any(y == 0))),]
EUK03_uIS_Soen.bw_0.decay<-EUK03_uIS_Soen.dissimilariy.bw[!(apply(EUK03_uIS_Soen.dissimilariy.bw, 1, function(y) any(y == 0))),]

plot(BAKT03_uIS_Soen.bw_0.decay)
plot(EUK03_uIS_Soen.bw_0.decay)

xbbw0<-((BAKT03_uIS_Soen.bw_0.decay$bwdistB_uIs_Soen))
ybbw0<-BAKT03_uIS_Soen.bw_0.decay$BAKT_uIs_Soen.dist
xebw0<-((EUK03_uIS_Soen.bw_0.decay$bwdistE_uIs_Soen))
yebw0<-EUK03_uIS_Soen.bw_0.decay$EUK03_uIs_Soen.dist

xybbw0 <- data.frame(xbbw0, ybbw0)
xyebw0 <- data.frame(xebw0, yebw0)

mbbw0<-lm(ybbw0~xbbw0, data=xybbw0)
mebw0<-lm(yebw0~xebw0, data=xyebw0)

summary(mbbw0)
summary(mebw0)

p_bbw0 <- ggplot(xybbw0, aes(x=xbbw0, y=ybbw0)) + theme_bw() +
  #annotate("text",x=225, y=0.45,label="R2=0.01") + 
  #annotate("text", x=215, y=0.37,label="p=0.12") +
  geom_point() + geom_smooth(aes(x=xbbw0, y=ybbw0), method="lm",colour="darkorchid") +
  theme(plot.margin=unit(c(0,0,0,0), "mm"))

p_ebw0 <- ggplot(xyebw0, aes(x=xebw0, y=yebw0)) + theme_bw() +
  #annotate("text",x=225, y=0.45,label="R2=0.23") + 
  #annotate("text", x=215, y=0.37,label="p<0.00001") +
  geom_point() + geom_smooth(aes(x=xebw0, y=yebw0), method="lm",colour="darkorchid") +
  theme(plot.margin=unit(c(0,0,0,0), "mm"))

multiplot(p_bbw0,p_ebw0)

#partial mantel tests correcting geography for the six "significant" environmental variables (bw,dh,no3,few,mnw,nh4) 
geo.distE<-as.matrix(geo.distE)
#geo.distE[upper.tri(geo.distE, diag=TRUE)] <- NA
#bwdistE1[upper.tri(bwdistE1, diag=TRUE)] <- NA

EUK03_u_sve.dist

mantel.partial(EUK03_u_sve.dist, bwdistE1, geo.distE, permutations=999) # in
mantel.partial(EUK03_u_sve.dist, dhdistE1, geo.distE, permutations=999) # in
mantel.partial(EUK03_u_sve.dist, no3distE1, geo.distE, permutations=999) # in
mantel.partial(EUK03_u_sve.dist, fewdistE1, geo.distE, permutations=999) # in
mantel.partial(EUK03_u_sve.dist, mnwdistE1, geo.distE, permutations=999) # in
mantel.partial(EUK03_u_sve.dist, nh4distE1, geo.distE, permutations=999) # in

EUK03_u_sve.dist<-as.matrix(vegdist(EUK03_u_sve,method="bray"))

#testing morans I
long.dists <- as.matrix(dist(EUK03.lon_0.decay$londistE1))

long.dists.inv <- 1/long.dists
long.dists.inv

long.dists.inv[which(long.dists.inv==Inf)] <- 0 
long.dists.inv[1:5, 1:5]

#latitudes
lat.dists <- as.matrix(dist(EUK03.lat_0.decay$latdistE1))

lat.dists.inv <- 1/lat.dists
lat.dists.inv

lat.dists.inv[which(lat.dists.inv==Inf)] <- 0 
lat.dists.inv[1:5, 1:5]

#in the matrix, each off-diagonal entry (i,j) is equal to 1/(distance between i and j). 
#Calculate Moran´s I:

Moran.I(EUK03.lon_0.decay$EUK03_u_sve.dist, long.dists.inv)

lonlatE<-correlog.nc(EUK03_u_sve.dist, geo.distE, latdistE1, increment=50)

plot(lonlatE)

#plot mantel correlogram:

plot(mantel.correlog(EUK03_u_sve.dist, D.geo=geo.distE, cutoff=FALSE,
                     break.pts=quantile(as.numeric(geo.distE), probs=seq(0, 1, 0.1))))

#create more beautiful ggplot correlogram
mcorrE<-mantel.correlog(EUK03_u_sve.dist, D.geo=geo.distE, cutoff=FALSE, mult="holm", n.class=3)

mcorrB<-mantel.correlog(BAKT03_u_OdsDLinD.dist, D.geo=geo.distB,cutoff=FALSE,mult="holm", n.class=3)

plot(mcorrE)
plot(mcorrB)

#create dataframes for the correlograms
e.df <- as.data.frame(mcorrE$mantel.res)
b.df <- as.data.frame(mcorrB$mantel.res)

e.df$Method <- "Pearson"
e.df$Microbiome <- "Eukaryotes"

b.df$Method <- "Pearson"
b.df$Microbiome <- "Prokaryotes"
#make a single data frame
comb.df <- rbind(e.df, b.df)

#make a significance factor
comb.df$Significance <- ifelse(comb.df$"Pr(corrected)" < 0.05, "Significant", "Insignificant")
#make a distance class factor
comb.df$Distance.Class <- rownames(comb.df)

#plot the correlograms i ggplot
p<-ggplot(comb.df, aes(class.index, Mantel.cor, group=Microbiome, shape=Microbiome)) + geom_line(size=1) +
  geom_point(size=4) + scale_shape_manual(values=c(16,1)) + theme(text = element_text(size=12), axis.text.x = element_text(angle=90, vjust=1))  +
  theme_bw()+ ylab("Mantel Correlation") + xlab("Distance Class")  + geom_hline(yintercept=0.00, color="black") +  theme(text = element_text(size=20)) + xlim(0, 250)

p+guides(colour=FALSE, size=FALSE, shape=FALSE)

# where OTU.dist is community matrix (after vegdist(OTU)), and sp.dist geodistance matrix.

#Testing reduced datasets

EUK03_80<- read.table("shared_tax80_reduced.txt",header=TRUE,row.names=1)
sum(rowSums(EUK03_80))

geo.dist21<-read.table("geodist_21.txt",header=TRUE,row.names=1)

EUK03_80.dist<-vegdist(EUK03_80,method="bray")
mantel(EUK03_80.dist, geo.dist21)

e03_80<-t(EUK03_80)
e03_80t<-estimateD(e03_80, datatype="abundance", base ="coverage", level=0.95)

EUK03_proto<- read.table("FINAL_0.03_euk_astT_shared_tax80_PROTO_U_tax.txt",header=TRUE,row.names=1)
sum(rowSums(EUK03_proto))

EUK03_80.dist<-vegdist(EUK03_proto,method="bray")
mantel(EUK03_80.dist, geo.dist21)


geo.dist21[geo.dist21 == 0] <- NA
EUK03_80.dist[is.na(as.dist(geo.dist21))] <- NA

mantel(EUK03_80.dist, geo.dist21, na.rm=TRUE)

EUK03_proto1<-EUK03_proto[-c(12,13),]
geo.dist.proto1<-geo.dist21[-c(12,13), -c(12,13)]
EUK03_proto1.dist<-vegdist(EUK03_proto1,method="bray")
mantel(EUK03_proto1.dist, geo.dist.proto1)

sum(rowSums(EUK03_proto1))

e03_80_proto<-t(EUK03_proto)
e03_80t<-estimateD(e03_80, datatype="abundance", base ="coverage")

EUK03_proto_cov<-EUK03_proto[-c(3,4,12,13,21),]
geo.dist.proto<-geo.dist21[-c(3,4,12,13,21), -c(3,4,12,13,21)]

EUK03_proto_cov.dist<-vegdist(EUK03_proto_cov,method="bray")
mantel(EUK03_proto_cov.dist, geo.dist.proto)

#age test >15 years

EUK03_age15<-EUK03[-c(3,4,5,6,9,12,13,14,15,21),]
geo.dist.age15<-geo.dist21[-c(3,4,5,6,9,12,13,14,15,21), -c(4,5,6,9,12,13,14,15,21)]

EUK03_age15.dist<-vegdist(EUK03_age20,method="bray")
mantel(EUK03_age15.dist, geo.dist.age20)

#bacteria
env$age
BAKT03_age15<-BAKT03[-c(4,5,6,9,12,13,14,15),]

BAKT03_age15.dist<-vegdist(BAKT03_age20,method="bray")
mantel(BAKT03_age15.dist, geo.dist.age20)

env_age20<-env[-c(3,4,5,6,9,12,13,14,15,21),]

EUK03_age20t<-t(EUK03_age20)
EUK03_age20t_rich<-estimateD(EUK03_age20t, datatype="abundance", base ="coverage")

#9000 test

EUK03_plus9000<-EUK03[-c(3,4,12,13,14,15,16,18,19,20,21),]
geo.dist.plus9000<-geo.dist21[-c(3,4,12,13,14,15,16,18,19,20,21), -c(3,4,12,13,14,15,16,18,19,20,21)]

EUK03_plus9000.dist<-vegdist(EUK03_plus9000,method="bray")
mantel(EUK03_plus9000.dist, geo.dist.plus9000)

geo.dist.plus9000[geo.dist.plus9000 == 0] <- NA
EUK03_plus9000.dist[is.na(as.dist(geo.dist.plus9000))] <- NA

mantel(EUK03_plus9000.dist, geo.dist.plus9000, na.rm=TRUE)

env_plus9000<-env[-c(3,4,12,13,14,15,16,18,19,20,21),]

EUK03_plus9000t<-t(EUK03_plus9000)
EUK03_plus9000t_rich<-estimateD(EUK03_plus9000t, datatype="abundance", base ="coverage")

EUK03_plus9000_rich_div<-cbind(row.names(EUK03_plus9000),env_plus9000$age,EUK03_plus9000_rich$"q = 0")

EUK03_plus9000.distance.decay = data.frame(
  geo.dist.plus9000 = as.matrix(geo.dist.plus9000)[lower.tri(geo.dist.plus9000)],
  EUK03_plus9000.dist = as.matrix(EUK03_plus9000.dist)[lower.tri(EUK03_plus9000.dist)]
)

rowSums(EUK03_plus9000)

H_Euk03_plus9000 <- diversity(EUK03_plus9000)

J_EUK03_plus9000 <- H_Euk03_plus9000/(log(specnumber(EUK03_plus9000)))

write.table(J_EUK03_plus9000 , "~/Desktop/Final_euk_NN/J_EUK03_plus9000.txt",sep="\t") 

jeuk9000<-read.table("~/Desktop/Final_euk_NN/J_EUK03_plus9000.txt",header=TRUE)

mPilEuk<-lm(env_plus9000$age ~ jeuk9000$x)
summary(mPilEuk)

EUK03_plus9000.distance.decay<-EUK03_plus9000.distance.decay[!(apply(EUK03_plus9000.distance.decay, 1, function(y) any(y == 0))),]

x8<-((EUK03_plus9000.distance.decay$geo.dist.plus9000))
y8<-EUK03_plus9000.distance.decay$EUK03_plus9000.dist

xy8 <- data.frame(x8, y8)

m8<-lm(y8~x8, data=xy8)

summary(m8)

#test antal sekvenser ~artsrigdom

numspec9000<-specnumber(EUK03_plus9000)
numseq9000<-rowSums(EUK03_plus9000)
check<-lm(numseq9000~numspec9000)
plot(numseq9000~numspec9000)

summary(check)


ggplot(xy8, aes(x=x8, y=y8)) + theme_bw() +
  annotate("text",x=225, y=0.40,label="R2=0.31") + 
  annotate("text", x=215, y=0.37,label="p<0.00001") +
  geom_point() + geom_smooth(aes(x=x8, y=y8), method="lm",colour="darkorchid") +
  coord_cartesian(ylim=c(0.35,1), xlim=c(1.9,252)) +
  theme(plot.margin=unit(c(1,2,1,-1), "mm"))

#9000 test med resampling
rowSums(EUK03_plus9000)
EUK03_plus9000_ss <- rrarefy(EUK03_plus9000, 9125) 

EUK03_plus9000_ss.dist<-vegdist(EUK03_plus9000_ss,method="bray")
mantel(EUK03_plus9000_ss.dist, geo.dist.plus9000)

#evenness
H_Bakt03_ss <- diversity(BAKT03_plus8600_ss)
H_Euk03_ss <- diversity(EUK03_plus9000_ss)

#Pielou's evenness J=H0=log(S):

J_Bakt03_ss <- H_Bakt03_ss/(log(specnumber(BAKT03_plus8600_ss)))
J_Euk03_ss <- H_Euk03_ss/(log(specnumber(EUK03_plus9000_ss)))


#8600 bakteria
rowSums(BAKT03)
BAKT03_plus8600<-BAKT03[-c(1,2,5,6,7,8,9,10,15,18,20,21),]
geo.dist.plus8600<-geo.dist21[-c(1,2,5,6,7,8,9,10,15,18,20,21), -c(1,2,5,6,7,8,9,10,15,18,20,21)]

BAKT03_plus8600.dist<-vegdist(BAKT03_plus8600,method="bray")
mantel(BAKT03_plus8600.dist, geo.dist.plus8600)

geo.dist.plus8600[geo.dist.plus8600 == 0] <- NA
BAKT03_plus8600.dist[is.na(as.dist(geo.dist.plus8600))] <- NA

mantel(BAKT03_plus8600.dist, geo.dist.plus8600, na.rm=TRUE)



#m. resampling
BAKT03_plus8600_ss<-rrarefy(BAKT03_plus8600, 8624)
BAKT03_plus8600_ss.dist<-vegdist(BAKT03_plus8600_ss,method="bray")
mantel(BAKT03_plus8600_ss.dist, geo.dist.plus8600)

#tjek dist-decay for plus9000 uden resampling

x1<-as.numeric(EUK03_plus9000_rich_div[,2])
y1<-as.numeric(EUK03_plus9000_rich_div[,3])

is.numeric(x1)
is.numeric(y1)
xy1 <- data.frame(x1, y1)

m1<-lm(y1~x1, data=xy1)
summary(m1)

ggplot(xy1, aes(x=x1, y=y1)) + theme_bw() +
  annotate("text",x=33, y=50,label="R2=0.1") + 
  annotate("text", x=33, y=37,label="p=0.19") +
  geom_point() + geom_smooth(aes(x=x1, y=y1), method="lm", colour="red") +
  coord_cartesian(ylim=c(0,225), xlim=c(0,41)) + 
  theme(plot.margin=unit(c(1,0,1,1), "mm"))

#Plot distance-decay relationship: Get the lower triangle of each distance matrix and combine
EUK03_plus9000.distance.decay = data.frame(
  geo.dist.plus9000 = as.matrix(geo.dist.plus9000)[lower.tri(geo.dist.plus9000)],
  EUK03_plus9000.dist = as.matrix(EUK03_plus9000.dist)[lower.tri(EUK03_plus9000.dist)]
)

#Delete functions with 0
EUK03_plus9000.distance_0.decay<-EUK03_plus9000.distance.decay[!(apply(EUK03_plus9000.distance.decay, 1, function(y) any(y == 0))),]

x2<-EUK03_plus9000.distance_0.decay$geo.dist.plus9000
y2<-EUK03_plus9000.distance_0.decay$EUK03_plus9000.dist

xy2 <- data.frame(x2, y2)
m2<-lm(y2~x2, data=xy2)
summary(m2)


#protozo_coverage
e03_80_proto_covt<-t(EUK03_proto_cov)
proto_cov1<-estimateD(e03_80_proto_covt, datatype="abundance", base ="coverage")

EUK03_proto_cov_uVj<-EUK03_proto[-c(1,2,3,4,12,13,21),]
geo.dist.proto_uVj<-geo.dist21[-c(1,2,3,4,12,13,21), -c(1,2,3,4,12,13,21)]

EUK03_proto_cov_uVj.dist<-vegdist(EUK03_proto_cov_uVj,method="bray")
mantel(EUK03_proto_cov_uVj.dist, geo.dist.proto_uVj)

#Check if sequence count correlates with richness/diversity
#For the whole dataset

seq_countall<-rowSums(EUK03)
species_all<-specnumber(EUK03)

plot(seq_countall,species_all)

m_seqspecall<-lm(species_all~seq_countall)
summary(m_seqspecall)
#p=0,00043, st??rkt signifikant

seq_count9000<-rowSums(EUK03_plus9000)
species_9000<-specnumber(EUK03_plus9000)

plot(seq_count9000,species_9000)

m_seqspec9000<-lm(species_9000~seq_count9000)
summary(m_seqspec9000)
#p=0,17, ikke signifikant

#Bacteria

seq_countallB<-rowSums(BAKT03)
species_allB<-specnumber(BAKT03)

plot(seq_countallB,species_allB)

m_seqspecallB<-lm(species_allB~seq_countallB)
summary(m_seqspecall)
#p=0,00043, st??rkt signifikant

seq_count8600<-rowSums(BAKT03_plus8600)
species_8600<-specnumber(BAKT03_plus8600)

plot(seq_count8600,species_8600)

m_seqspec8600<-lm(species_8600~seq_count8600)
summary(m_seqspec8600)
#p=0,073, ikke signifikant