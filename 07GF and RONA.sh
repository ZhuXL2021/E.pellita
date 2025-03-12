#GF offset
#in R
library(gradientForest)
Ent.data <- read.table("CPA.Ent_facter.txt", header = TRUE,row.names = "pop") 
SNP.MAF <- read.table("BIOSNP.maf", header = TRUE,row.names = "pop") 
present <- data.frame(t(Ent.data))
SNP.data <- data.frame(t(SNP.MAF))
bioclimatic <- paste("BIO",1:19,sep = "")  # 
maxLevel <- log2(0.368*nrow(SNP.data)/2   )# 
gf_SNP.data <- gradientForest(cbind(present[,bioclimatic], SNP.data),  
                              predictor.vars=colnames(present[,bioclimatic]), 
                              response.vars=colnames(SNP.data), ntree=500,    
                              compact=T,
                              maxLevel=maxLevel, trace=T, corr.threshold=0.70)  
pdf(file="gf_SNP.data01.pdf")
plot(gf_SNP.data, plot.type = "Overall.Importance", col=c(rep("grey",15) ),
     las=2,cex.names=0.8) 
dev.off()  
###Import current climate
greengrid <- read.csv("CPA_pc.19BIO.csv", header = TRUE)
head(greengrid) 
dim(greengrid)
all_tgrid=cbind(greengrid[,c("Lon","Lat")], predict(gf_SNP.data, 
               greengrid[,3:21]))

n<-sum(is.na(all_tgrid))
Trns_grid <- na.omit(all_tgrid) 
n<-sum(is.na(Trns_grid)) 
n<-sum(is.na(all_tgrid))
Trns_grid <- na.omit(all_tgrid)
n<-sum(is.na(Trns_grid))
all_PCs <- prcomp(Trns_grid[,c("BIO2","BIO4","BIO10","BIO12","BIO15")],
                  center=TRUE, scale.=FALSE)
summary(all_PCs)
#set up a colour palette for the mapping
a1 <- all_PCs$x[,1]
a2 <- all_PCs$x[,2]
a3 <- all_PCs$x[,3]
r <- a1+a2
g <- -a2
b <- a3+a2-a1
r <- (r-min(r)) / (max(r)-min(r)) * 255
g <- (g-min(g)) / (max(g)-min(g)) * 255
b <- (b-min(b)) / (max(b)-min(b)) * 255
grid <- greengrid[,c("Lon","Lat")]
grid$R=r
grid$G=g
grid$B=b
nvs <- dim(all_PCs$rotation)[1]
vec <- c("BIO2","BIO4","BIO10","BIO12","BIO15")
lv <- length(vec)
vind <- rownames(all_PCs$rotation) %in% vec
scal <- 60
xrng <- range(all_PCs$x[,1], all_PCs$rotation[,1]/scal)*1.1
yrng <- range(all_PCs$x[,2], all_PCs$rotation[,2]/scal)*1.1
pdf(file="all_PCplot02.pdf")
plot((all_PCs$x[,1:2]), xlim=xrng, ylim=yrng, pch=".", cex=7, col=rgb(r,g,b, max = 255), asp=1)
arrows(rep(0,lv), rep(0,lv), all_PCs$rotation[,1]/scal, all_PCs$rotation[,2]/scal, length = 0.1)
jit <- 0.0015
text(all_PCs$rotation[,1]/scal+jit*sign(all_PCs$rotation[,1]), all_PCs$rotation[,2]/scal+jit*sign(all_PCs$rotation[,2]), labels = vec)
dev.off()

pdf("Map2.pdf")
green.pred <- predict(gf_SNP.data, greengrid[,c("BIO2","BIO4","BIO10","BIO12","BIO15")])
plot(Trns_grid[, c("Lon", "Lat")], pch=15,cex=1.0,asp=1,col=rgb(r,g,b, max=255))
dev.off()
#export map for use in ArcGIS
greencols=rgb(r,g,b,max=255)
greencols2=col2rgb(greencols)
greencols3=t(greencols2)
gradients=cbind(Trns_grid[,1:2],greencols3)
gradients$color=greencols
write.csv(gradients,file="all_gradients4arcgis03.csv",row.names=F,quote=F)
all_tgrid=cbind(greengrid[,c("Lon","Lat")], predict(gf_SNP.data,
                greengrid[,3:21])) 
##### offset
fut <- read.table("CPA_ssp585.fut.txt", header = TRUE) #Import future climate
head(fut)
future_all=cbind(fut[,c("Lon","Lat")], predict(gf_SNP.data,
                fut[,3:21]))
genOffsetAll <-sqrt((future_all[,4]-all_tgrid[,4])^2+(future_all[,6]-all_tgrid[,6])^2+(future_all[,12]-all_tgrid[,12])^2+(future_all[,14]-all_tgrid[,14])^2+(future_all[,17]-all_tgrid[,17])^2)
Offset=cbind(future_all[,c("Lon","Lat")],genOffsetAll) 
colnames(Offset)[3]<-"offset"  
write.csv(Offset, "CPA_ssp585.fut.5BIO.offset.csv", quote=F, row.names=T) 



#RONA
pyRona lfmm -covars 5 -pc CPA_RONA.pc.5BIO.txt -fc CPA_RONA.SSP585_2090.5BIO.txt -out out/CPA_RONA.SSP585_2090.5BIO.pdf -P 0.05 -assoc CPA.RONA.csv -geno LFMM/CPA_LFMM.lfmm >CPA_RONA.SSP585_2090.5BIO.log
