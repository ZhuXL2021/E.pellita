#LMFF
plink --allow-extra-chr --vcf 1286304SNP.vcf.gz  --recode --out CPA_LFMM 
#ped2lfmm
ped2lfmm CPA_LFMM.ped CPA_LFMM.lfmm

#in R
library(lfmm)
X <- read.table("CPA.Env.txt", h = T, stringsAsFactors = F) #environmental data
X <- as.matrix(X[, c(4:8)]) 
Y <- data.table::fread(paste("CPA_LFMM", ".lfmm",
                             sep = ""), header = F)
mod.lfmm <- lfmm::lfmm_ridge(Y = Y,
                             X = X,
                             K = 3) 
pv <- lfmm::lfmm_test(Y = Y, X = X, lfmm = mod.lfmm,
                      calibrate = "gif")
adjpvalues <- pv$calibrated.pvalue 
hist(pv$pvalue[,1], main="Unadjusted p-values")
hist(pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")
write.csv(adjpvalues,"CPA.adjpvalues.csv") 
filtered_data <- adjpvalues[rowSums(adjpvalues < 0.05) > 0, ]
write.csv(filtered_data,"CPA.filtered_data.csv") 
save(list = ls(), file = "CPA.LFMM.RData")


#RDA
#in R
#vcf2ped
#Run in R
library(vegan)
library(psych)
gen <- read.csv("CPA_RDA.lfmm",sep = "", header = F)
env <- read.table("CPA.Env.txt", h = T, stringsAsFactors = F)
env <- env[,c(4:8)]
colnames(env) <- c('BIO2','BIO4','BIO10','BIO12','BIO15') 
wolf.rda <- rda(gen ~ ., data=env, scale=T)
wolf.rda
RsquareAdj(wolf.rda) 
summary(eigenvals(wolf.rda, model = "constrained"))
screeplot(wolf.rda)
signif.full <- anova.cca(wolf.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full 
vif.cca(wolf.rda)
plot(wolf.rda, scaling=3) 
plot(wolf.rda, choices = c(1, 3), scaling=3) 
#SNPS were identified from the three RDAs
load.rda <- scores(wolf.rda, choices=c(1:3), display="species")
pdf("CPA.RDA1.pdf") 
opar<-par(no.readonly=TRUE)
par(mfrow=c(3,2)) 
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3")
par(opar)
dev.off()
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}
cand1 <- outliers(load.rda[,1],3) 
cand2 <- outliers(load.rda[,2],3)
cand3 <- outliers(load.rda[,3],3)
ncand <- length(cand1) + length(cand2) + length(cand3)
ncand
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3)    <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)
foo <- matrix(nrow=(ncand), ncol=5) 
colnames(foo) <- c("BIO2","BIO4","BIO10","BIO12","BIO15")
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen[,nam]
  foo[i,] <- apply(env,2,function(x) cor(x,snp.gen))
}
cand <- cbind.data.frame(cand,foo)  
head(cand)
length(cand$snp[duplicated(cand$snp)])  
foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # n duplicates on axis 1
table(foo[foo[,1]==2,2]) # n duplicates on axis 2
table(foo[foo[,1]==3,2]) # n duplicates on axis 3
cand <- cand[!duplicated(cand$snp),] # remove duplicate detections
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,9] <- names(which.max(abs(bar[4:8]))) 
  cand[i,10] <- max(abs(bar[4:8]))             
}
colnames(cand)[9] <- "predictor"
colnames(cand)[10] <- "correlation"
table(cand$predictor) 
write.csv(cand,"CPA.3RDA_result.csv")


