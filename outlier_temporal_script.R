require(adegenet)
require(diveRsity)
require(hierfstat)
require(rapport)

# set your own working directory here:
setwd("../Desktop/dorte_geneticCline/")
setwd("C:/Users/bmen/Desktop/dorte_geneticCline/")

#GP <- readGenepop(infile = "287SNPS_GP.gen", gp = 2, bootstrap = T)

#### 287 SNPs file ####
GP <- read.genepop("287SNPS_GP_19112020.gen")

# if error appears Error in dimnames(x) <- dn : 
#length of 'dimnames' [1] not equal to array extent
#In addition: Warning message:
#  In matrix(unlist(strsplit(vec.genot, "[[:space:]]+")), ncol = nloc,  :
#              data length [28913] is not a sub-multiple or multiple of the number of rows [129]
# Checks: NA=>00; there is a heading on the first line ""; there is at least one POP after all locus list;
# there is space in between genotypes and not tabs

locus <- read.table("287SNPs.txt",stringsAsFactors = F)
head(locus)

# include the pop map
GP_pop <- read.table("287SNPS_GP_19112020_popmap.txt",stringsAsFactors = F)
GP_pop <- as.data.frame(GP_pop)
GP_pop_fac <- as.factor(GP_pop$V1)
table(GP_pop_fac)

# convert to hierfstat
GP_hf <- genind2hierfstat(GP,pop=GP_pop_fac)
head(GP_hf[1:3,1:3])

GP_stats <- basic.stats(GP_hf,diploid=TRUE,digits=4)
str(GP_stats$pop.freq)

# Explore the dataset
GP_stats$pop.freq
write.table(GP_stats$pop.freq, "GP_stats")
GP_stats$pop.freq[1:2]
GP_stats$pop.freq$Gdist_S284585_2743[1:2, 1:5]

# Save dataset as R object for Henrik's tests
saveRDS(GP_stats, "GP_stats")

#little script by Henrik Baktoft to extract allele freqs per locus
#and generate dataframe according to signasel's output
dat <- readRDS('GP_stats')
dat_out_277 <- c()
for(i in 1:length(dat$pop.freq)){
  locus <- names(dat$pop.freq)[[i]]
  pop <- colnames(dat$pop.freq[[i]])
  if(nrow(dat$pop.freq[[i]]) == 2){
    x1 <- dat$pop.freq[[i]][1,]
    x2 <- dat$pop.freq[[i]][2,]
    
    dat_i <- data.table::data.table(locus, pop, x1,x2, i)
    dat_out_277 <- rbind(dat_out_277, dat_i)
  } else {
    print(paste0("NOTE: Something is missing in ", locus, "! Not included in final output"))
  }
}

# there is a locus that raises an error
# We eliminate it so it does not raise the error
dat$pop.freq[[160]]

# explorations/dimensions of the df
dim(dat_out)
tail(dat_out)

# Include pop sizes
table(GP_pop_fac)
pop_sizes <- c(68,6,6,5,40,31,3,6,21,3,31,5,39,49)
length(pop_sizes)
dat_out_277$N <- rep(c(68,6,6,5,40,31,3,6,21,3,31,5,39,49), 
                 times=length(dat_out_277$pop)/length(pop_sizes))

# Include Ne in the df
#dat_out$Ne <- c(100)
dat_out_277$Ne <- rep(c(0,0,0,0,200,200,0,0,200,0,200,0,200,200),
                  times=length(dat_out_277$pop)/length(pop_sizes))
# Include sample size in the df
dat_out_277$gen <- rep(c(0,0,0,0,0,5,0,0,0,0,2,0,3,4),
                   times=length(dat_out_277$pop)/length(pop_sizes))
dat_out_277$no_copies <- round((dat_out_277$N*2)*(dat_out_277$x1))
tail(dat_out_277,14)

# Signasel:
# signasel_code.R should be downloaded from the developers repository
# https://github.com/hubert-pop/signasel/blob/master/signasel3s-gh1.R
source("../salmon/signasel/signasel_code.R")
example.mat <- matrix(data=c(0, 15, 30, 50, 5, 20, 30, 50, 10, 25, 30, 50), 
                      ncol=4, byrow=T)
example.mat
signaseltest(example.mat)

## For each sample (row) the columns are: g, i, S, N:
## g0 i0 S0 N0
## g1 i1 S1 N1
## g2 i2 S2 N2
## etc.
## with g: generation, i: number of allele copies, S: sample size,
## N: (effective) population size

# initializing storage df
# SK 70 and 280
signa_287output <- data.frame(matrix(ncol = 6, nrow = 286))
colnames(signa_287output) <- c("L0","Lmax","smax","LRT","-log10pvalue", "warn")

# VA:130,0,130,0,280,280
signa_287output_VA <- data.frame(matrix(ncol = 6, nrow = 286))
colnames(signa_287output_VA) <- c("L0","Lmax","smax","LRT","-log10pvalue", "warn")

# SK2011=280; SK2019=280
signa_287_SK_output_2 <- data.frame(matrix(ncol = 6, nrow = 286))
colnames(signa_287_SK_output_2) <- c("L0","Lmax","smax","LRT","-log10pvalue", "warn")

# SK2011=100; SK2019=100
signa_287_SK_output_3 <- data.frame(matrix(ncol = 6, nrow = 286))
colnames(signa_287_SK_output_3) <- c("L0","Lmax","smax","LRT","-log10pvalue", "warn")

# SK2011=200; SK2019=200
signa_287_SK_output_4 <- data.frame(matrix(ncol = 6, nrow = 286))
colnames(signa_287_SK_output_4) <- c("L0","Lmax","smax","LRT","-log10pvalue", "warn")

# VA2011=100; VA2012=100; VA2016=100)
signa_287_VA_output_2 <- data.frame(matrix(ncol = 6, nrow = 286))
colnames(signa_287_VA_output_2) <- c("L0","Lmax","smax","LRT","-log10pvalue", "warn")

# VA2011=200; VA2012=200; VA2016=200)
signa_287_VA_output_3 <- data.frame(matrix(ncol = 6, nrow = 286))
colnames(signa_287_VA_output_3) <- c("L0","Lmax","smax","LRT","-log10pvalue", "warn")



for (loc in c(161:287)){
  print(loc)
  # re-order the columns to obtain the input for signasel
  signasel.mat <- dat_out[dat_out$i==loc,c(8,9,6,7)]
  #in case we only want specific generations
  #signasel.mat <- signasel.mat[c(5,6),] # SK2001 vs SK2011
  signasel.mat <- signasel.mat[c(9,11,13,14),] # VA2004-2010-2012-2016
  signasel.mat <- as.matrix(signasel.mat)
  signasel.mat[,4] <- as.integer(signasel.mat[,4])
  #signa_287output[loc,] <- signaseltest(signasel.mat)
  #signa_287output_VA[loc,] <- signaseltest(signasel.mat)
  #signa_287_SK_output_2[loc,] <- signaseltest(signasel.mat)
  #signa_287_VA_output_2[loc,] <- signaseltest(signasel.mat)
  #signa_287_SK_output_3[loc,] <- signaseltest(signasel.mat)
  #signa_287_SK_output_4[loc,] <- signaseltest(signasel.mat)
  signa_287_VA_output_3[loc,] <- signaseltest(signasel.mat)
}

# SKJERN
head(signa_287output)
tail(signa_287output)

hist(signa_287output$smax, breaks=100)
hist(signa_287output$LRT,breaks=100)

output287_pval1 <- which(signa_287output$`-log10pvalue`> 1.3,) # pval 0.05
output287_pval2 <- which(signa_287output$`-log10pvalue`> 2,) # pval 0.01
output287_pval3 <- which(signa_287output$`-log10pvalue`> 3,) # pval 0.001

signa_287output[c(output287_pval2),]
dat_out[dat_out$i==147,]
#16  25  51  88  95  97 100 107 146 147 157 161 167 210 243 257

# VARDE
output287_VA_pval1 <- which(signa_287output_VA$`-log10pvalue`> 1.3,) # pval 0.05
output287_VA_pval2 <- which(signa_287output_VA$`-log10pvalue`> 2,) # pval 0.01
output287_VA_pval3 <- which(signa_287output_VA$`-log10pvalue`> 3,) # pval 0.001

signa_287output[c(output287_pval2),]
dat_out[dat_out$i==147,]
#16  25  51  88  95  97 100 107 146 147 157 161 167 210 243 257

# save the datasets
save(signa_287output, file="signa_287_SK_output.RData")
save(signa_287_SK_output_2, file="signa_287_SK_output_2.RData")
save(signa_287_SK_output_3, file="signa_287_SK_output3.RData")
save(signa_287_SK_output_4, file="signa_287_SK_output4.RData")

save(signa_287output_VA, file="signa_287output_VA.RData")
save(signa_287_VA_output_2, file="signa_287_VA_output2.RData")
save(signa_287_VA_output_3, file="signa_287_VA_output3.RData")

load("signa_287_VA_output2.RData")
load("signa_287_VA_output3.RData")


#### 3656 SNPs file ####
GP_3656 <- read.genepop("3656SNPS_GP_20112020.gen")

# include the pop map
GP_3656_pop <- read.table("3656SNPS_GP_20112020_popmap.txt",stringsAsFactors = F)
GP_3656_pop <- as.data.frame(GP_3656_pop)
GP_3656_pop_fac <- as.factor(GP_3656_pop$V1)
table(GP_3656_pop_fac)

# convert to hierfstat
GP_3656_hf <- genind2hierfstat(GP_3656,pop=GP_3656_pop_fac)
head(GP_3656_hf[1:3,1:3])

GP_3656_stats <- basic.stats(GP_3656_hf,diploid=TRUE,digits=4)
str(GP_3656_stats$pop.freq)

# Explore the dataset
GP_3656_stats$pop.freq
#write.table(GP_3656_stats$pop.freq, "GP_3656_stats")
GP_3656_stats$pop.freq[1:2]
GP_3656_stats$pop.freq$Gdist_S284585_2743[1:2, 1:5]

# Save dataset as R object for Henrik's tests
saveRDS(GP_3656_stats, "GP_3656_stats")

#little script by Henrik Baktoft to extract allele freqs per locus
#and generate dataframe according to signasel's output
dat_3656 <- readRDS('GP_3656_stats')
dat_3656_out <- c()
for(i in 1:length(dat_3656$pop.freq)){
  locus <- names(dat_3656$pop.freq)[[i]]
  pop <- colnames(dat_3656$pop.freq[[i]])
  if(nrow(dat_3656$pop.freq[[i]]) == 2){
    x1 <- dat_3656$pop.freq[[i]][1,]
    x2 <- dat_3656$pop.freq[[i]][2,]
    
    dat_i <- data.table::data.table(locus, pop, x1,x2, i)
    dat_3656_out <- rbind(dat_3656_out, dat_i)
  } else {
    print(paste0("NOTE: Something is missing in ", locus, "! Not included in final output"))
  }
}

# explorations/dimensions of the df
dim(dat_3656_out)
tail(dat_3656_out)

# Include pop sizes
table(GP_3656_pop_fac)
pop_3656_sizes <- c(68,6,6,5,9,31,18,16,3,6,21,3,31,5,39)
length(pop_3656_sizes)
dat_3656_out$N <- rep(c(68,6,6,5,9,31,18,16,3,6,21,3,31,5,39), 
                 times=length(dat_3656_out$pop)/length(pop_3656_sizes))

# Include Ne in the df
#dat_out$Ne <- c(100)
dat_3656_out$Ne <- rep(c(0,0,0,0,200,200,0,0,0,0,200,0,200,0,200),
                  times=length(dat_3656_out$pop)/length(pop_3656_sizes))
# Include sample size in the df
dat_3656_out$gen <- rep(c(0,0,0,0,0,5,0,0,0,0,0,0,2,0,3),
                   times=length(dat_3656_out$pop)/length(pop_3656_sizes))
dat_3656_out$no_copies <- round((dat_3656_out$N*2)*(dat_3656_out$x1))
tail(dat_3656_out,14)

# Signasel:
# signasel_code.R should be downloaded from the developers repository
# https://github.com/hubert-pop/signasel/blob/master/signasel3s-gh1.R
source("../salmon/signasel/signasel_code.R")

## For each sample (row) the columns are: g, i, S, N:
## g0 i0 S0 N0
## g1 i1 S1 N1
## g2 i2 S2 N2
## etc.
## with g: generation, i: number of allele copies, S: sample size,
## N: (effective) population size

# initializing storage df
# SK2011=70; SK2019=280
signa_3656_SK_output <- data.frame(matrix(ncol = 6, nrow = 3656))
colnames(signa_3656_SK_output) <- c("L0","Lmax","smax","LRT","-log10pvalue", "warn")

# VA2011=130; VA2012=130; VA2016=280)
signa_3656_VA_output <- data.frame(matrix(ncol = 6, nrow = 3656))
colnames(signa_3656_VA_output) <- c("L0","Lmax","smax","LRT","-log10pvalue", "warn")

# initializing storage df
# SK2011=280; SK2019=280
signa_3656_SK_output_2 <- data.frame(matrix(ncol = 6, nrow = 3656))
colnames(signa_3656_SK_output_2) <- c("L0","Lmax","smax","LRT","-log10pvalue", "warn")

# SK2011=100; SK2019=100
signa_3656_SK_output_3 <- data.frame(matrix(ncol = 6, nrow = 3656))
colnames(signa_3656_SK_output_3) <- c("L0","Lmax","smax","LRT","-log10pvalue", "warn")

# SK2011=200; SK2019=200
signa_3656_SK_output_4 <- data.frame(matrix(ncol = 6, nrow = 3656))
colnames(signa_3656_SK_output_4) <- c("L0","Lmax","smax","LRT","-log10pvalue", "warn")

# VA2011=100; VA2012=100; VA2016=100)
signa_3656_VA_output_2 <- data.frame(matrix(ncol = 6, nrow = 3656))
colnames(signa_3656_VA_output_2) <- c("L0","Lmax","smax","LRT","-log10pvalue", "warn")

# VA2011=200; VA2012=200; VA2016=200)
signa_3656_VA_output_3 <- data.frame(matrix(ncol = 6, nrow = 3656))
colnames(signa_3656_VA_output_3) <- c("L0","Lmax","smax","LRT","-log10pvalue", "warn")


for (loc in c(1:3656)){
  print(loc)
  # re-order the columns to obtain the input for signasel
  signasel.mat <- dat_3656_out[dat_3656_out$i==loc,c(8,9,6,7)]
  #in case we only want specific generations
  #signasel.mat <- signasel.mat[c(5,6),] # SK2001 vs SK2011
  signasel.mat <- signasel.mat[c(11,13,15),] # VA2004-2010-2012
  signasel.mat <- as.matrix(signasel.mat)
  signasel.mat[,4] <- as.integer(signasel.mat[,4])
  #signa_3656_SK_output[loc,] <- signaseltest(signasel.mat)# with Ne=70,280
  #signa_3656_VA_output[loc,] <- signaseltest(signasel.mat)# with Ne=130,130,280
  #signa_3656_SK_output_2[loc,] <- signaseltest(signasel.mat)# with Ne=280,280
  #signa_3656_VA_output_2[loc,] <- signaseltest(signasel.mat)# with Ne=100,100,100
  #signa_3656_SK_output_3[loc,] <- signaseltest(signasel.mat)# with Ne=100,100
  #signa_3656_SK_output_4[loc,] <- signaseltest(signasel.mat)# with Ne=200,200
  signa_3656_VA_output_3[loc,] <- signaseltest(signasel.mat)# with Ne=200,200
}

# SKJERN
head(signa_287output)
tail(signa_287output)

hist(signa_287output$smax, breaks=100)
hist(signa_287output$LRT,breaks=100)

output287_pval1 <- which(signa_287output$`-log10pvalue`> 1.3,) # pval 0.05
output287_pval2 <- which(signa_287output$`-log10pvalue`> 2,) # pval 0.01
output287_pval3 <- which(signa_287output$`-log10pvalue`> 3,) # pval 0.001

signa_287output[c(output287_pval2),]
dat_out[dat_out$i==147,]
#16  25  51  88  95  97 100 107 146 147 157 161 167 210 243 257

# VARDE
output287_VA_pval1 <- which(signa_287output_VA$`-log10pvalue`> 1.3,) # pval 0.05
output287_VA_pval2 <- which(signa_287output_VA$`-log10pvalue`> 2,) # pval 0.01
output287_VA_pval3 <- which(signa_287output_VA$`-log10pvalue`> 3,) # pval 0.001

signa_287output[c(output287_pval2),]
dat_out[dat_out$i==147,]
#16  25  51  88  95  97 100 107 146 147 157 161 167 210 243 257

# save the datasets
save(signa_3656_SK_output, file="signa_3656_SK_output.RData")# with Ne=70,280
save(signa_3656_SK_output_2, file="signa_3656_SK_output_2.RData")# with Ne=280,280
save(signa_3656_SK_output_3, file="signa_3656_SK_output_3.RData")# with Ne=100,100,100
save(signa_3656_SK_output_4, file="signa_3656_SK_output_4.RData")# with Ne=200,200,200

load("signa_3656_SK_output_3.RData")# with Ne=100,100,100
load("signa_3656_SK_output_4.RData")# with Ne=200,200,200

save(signa_3656_VA_output, file="signa_3656_VA_output.RData")# with Ne=130,130,280
save(signa_3656_VA_output_2, file="signa_3656_VA_output_2.RData") # with Ne=100,100
save(signa_3656_VA_output_3, file="signa_3656_VA_output_3.RData") # with Ne=200,200

load("signa_3656_VA_output_2.RData")
load("signa_3656_VA_output_3.RData")

#### 288 SNPs file for Skjern river ####
GP <- read.genepop("Skjern_288SNPS_Genepop.gen") 

locus <- read.table("locus_288skjern_list.txt",stringsAsFactors = F)
head(locus)

# include the pop map
GP_pop <- read.table("Skjern_288SNPS_Genepop_popmap.txt",stringsAsFactors = F, header=T)
head(GP_pop)
GP_pop <- as.data.frame(GP_pop)
GP_pop_fac <- as.factor(GP_pop[,4])
table(GP_pop_fac)

# convert to hierfstat
GP_hf <- genind2hierfstat(GP,pop=GP_pop_fac)
head(GP_hf[1:3,1:3])

GP_stats <- basic.stats(GP_hf,diploid=TRUE,digits=4)
str(GP_stats$pop.freq)

# Explore the dataset
GP_stats$pop.freq
write.table(GP_stats$pop.freq, "GP_288_skjern_stats")
GP_stats$pop.freq[1:2]
GP_stats$pop.freq$Gdist_S284585_2743[1:2, 1:5]

# Save dataset as R object for Henrik's tests
saveRDS(GP_stats, "GP_288_skjern_stats")

#little script by Henrik Baktoft to extract allele freqs per locus
#and generate dataframe according to signasel's output
dat <- readRDS('GP_288_skjern_stats')
dat_out_288 <- c()
for(i in 1:length(dat$pop.freq)){
  locus <- names(dat$pop.freq)[[i]]
  pop <- colnames(dat$pop.freq[[i]])
  if(nrow(dat$pop.freq[[i]]) == 2){
    x1 <- dat$pop.freq[[i]][1,]
    x2 <- dat$pop.freq[[i]][2,]
    
    dat_i <- data.table::data.table(locus, pop, x1,x2, i)
    dat_out_288 <- rbind(dat_out_288, dat_i)
  } else {
    print(paste0("NOTE: Something is missing in ", locus, "! Not included in final output"))
  }
}

# there is a locus that raises an error
# We eliminate it so it does not raise the error
dat$pop.freq[[288]]

# explorations/dimensions of the df
dim(dat_out)
tail(dat_out)

# Include pop sizes
table(GP_pop_fac)
pop_sizes <- c(40,4,8,8,2,4,24,2,30,31)
length(pop_sizes)
dat_out_288$N <- rep(c(40,4,8,8,2,4,24,2,30,31), 
                 times=length(dat_out_288$pop)/length(pop_sizes))

# Include Ne in the df
#dat_out$Ne <- c(100)
dat_out_288$Ne <- rep(c(0,0,0,200,0,0,200,0,200,200),
                  times=length(dat_out_288$pop)/length(pop_sizes))
# Include sample size in the df
# 2001, 2007, 2010, 2019
# 2001 -> 2007: 1.7 gen = 2 gen
# 2007 -> 2010: 0.85 = 1 gen
# 2010 -> 2019: 2.57 = 3 gen
dat_out_288$gen <- rep(c(0,0,0,0,0,0,2,0,3,6),
                   times=length(dat_out_288$pop)/length(pop_sizes))
dat_out_288$no_copies <- round((dat_out_288$N*2)*(dat_out_288$x1))
tail(dat_out_288,14)

# Signasel:
# signasel_code.R should be downloaded from the developers repository
# https://github.com/hubert-pop/signasel/blob/master/signasel3s-gh1.R
source("../salmon/signasel/signasel_code.R")
example.mat <- matrix(data=c(0, 15, 30, 50, 5, 20, 30, 50, 10, 25, 30, 50), 
                      ncol=4, byrow=T)
example.mat
signaseltest(example.mat)

## For each sample (row) the columns are: g, i, S, N:
## g0 i0 S0 N0
## g1 i1 S1 N1
## g2 i2 S2 N2
## etc.
## with g: generation, i: number of allele copies, S: sample size,
## N: (effective) population size

# initializing storage df
# SK -> Constant Ne of 100
signa_288_Skjern4pop_Ne100 <- data.frame(matrix(ncol = 6, nrow = 288))
colnames(signa_288_Skjern4pop_Ne100) <- c("L0","Lmax","smax","LRT","-log10pvalue", "warn")

# SK -> Constant Ne of 200
signa_288_Skjern4pop_Ne200 <- data.frame(matrix(ncol = 6, nrow = 288))
colnames(signa_288_Skjern4pop_Ne200) <- c("L0","Lmax","smax","LRT","-log10pvalue", "warn")

for (loc in c(1:288)){
  print(loc)
  # re-order the columns to obtain the input for signasel
  signasel.mat <- dat_out[dat_out$i==loc,c(8,9,6,7)]
  head(signasel.mat)
  #in case we only want specific generations
  signasel.mat <- signasel.mat[c(4,7,9,10),] # SK2001, 2007, 2011 and 2019
  #signasel.mat <- signasel.mat[c(9,11,13,14),] # VA2004-2010-2012-2016
  signasel.mat <- as.matrix(signasel.mat)
  signasel.mat[,4] <- as.integer(signasel.mat[,4])
  #signa_288_Skjern4pop_Ne100[loc,] <- signaseltest(signasel.mat) # With constant Ne of 100
  signa_288_Skjern4pop_Ne200[loc,] <- signaseltest(signasel.mat) # With constant Ne of 200
}

# Save and load the results if needed
save(signa_288_Skjern4pop_Ne100, file="signa_288_Skjern4pop_Ne100.RData") # with Ne=100
save(signa_288_Skjern4pop_Ne200, file="signa_288_Skjern4pop_Ne200.RData") # with Ne=200

load("signa_288_Skjern4pop_Ne100.RData")
load("signa_288_Skjern4pop_Ne200.RData")


# SKJERN
head(signa_288_Skjern4pop_Ne200)
tail(signa_288_Skjern4pop_Ne200)

##### Results ####
# SK 287 SNPs
signa_287output
signa_287_SK_output_2
signa_287_SK_output_3
# SK
signa_287output_pval1 <- which(signa_287output$`-log10pvalue`> 1.3,) # pval 0.05
signa_287output_pval2 <- which(signa_287output$`-log10pvalue`> 2,) # pval 0.01
signa_287output_pval3 <- which(signa_287output$`-log10pvalue`> 3,) # pval 0.001
signa_287output[c(signa_287output_pval1),]
dat_out[dat_out$i==257,]
# SK 280;280
signa_287_SK_output_2_pval1 <- which(signa_287_SK_output_2$`-log10pvalue`> 1.3,) # pval 0.05
signa_287_SK_output_2_pval2 <- which(signa_287_SK_output_2$`-log10pvalue`> 2,) # pval 0.01
signa_287_SK_output_2_pval3 <- which(signa_287_SK_output_2$`-log10pvalue`> 3,) # pval 0.001
signa_287_SK_output_2[c(signa_287_SK_output_2_pval3),]
dat_out[dat_out$i==100,]
# SK
signa_287_SK_output_3_pval1 <- which(signa_287_SK_output_3$`-log10pvalue`> 1.3,) # pval 0.05
signa_287_SK_output_3_pval2 <- which(signa_287_SK_output_3$`-log10pvalue`> 2,) # pval 0.01
signa_287_SK_output_3_pval3 <- which(signa_287_SK_output_3$`-log10pvalue`> 3,) # pval 0.001
signa_287_SK_output_3[c(signa_287_SK_output_3_pval3),]
dat_out[dat_out$i==16,]
# SK 200, 200
signa_287_SK_output_4_pval1 <- which(signa_287_SK_output_4$`-log10pvalue`> 1.3,) # pval 0.05
signa_287_SK_output_4_pval2 <- which(signa_287_SK_output_4$`-log10pvalue`> 2,) # pval 0.01
signa_287_SK_output_4_pval3 <- which(signa_287_SK_output_4$`-log10pvalue`> 3,) # pval 0.001
signa_287_SK_output_4[c(signa_287_SK_output_4_pval2),]
dat_out[dat_out$i==16,]

# SK 3656 SNPs
hist(signa_3656_SK_output_2$`-log10pvalue`)
signa_3656_SK_output_pval1 <- which(signa_3656_SK_output_2$`-log10pvalue`> 1.3,) # pval 0.05
signa_3656_SK_output_pval2 <- which(signa_3656_SK_output_2$`-log10pvalue`> 2,) # pval 0.01
signa_3656_SK_output_pval3 <- which(signa_3656_SK_output_2$`-log10pvalue`> 7,) # pval 0.000001
signa_3656_SK_output_2[c(signa_3656_SK_output_pval3),]
dat_3656_out[dat_3656_out$i==3479,]

hist(signa_3656_SK_output_3$`-log10pvalue`)
signa_3656_SK_output3_pval1 <- which(signa_3656_SK_output_3$`-log10pvalue`> 1.3,) # pval 0.05
signa_3656_SK_output3_pval2 <- which(signa_3656_SK_output_3$`-log10pvalue`> 2,) # pval 0.01
signa_3656_SK_output3_pval3 <- which(signa_3656_SK_output_3$`-log10pvalue`> 4,) # pval 0.001
signa_3656_SK_output_3[c(signa_3656_SK_output3_pval3),]
dat_3656_out[dat_3656_out$i==659,]
# with FDR correction
signa_3656_SK_output_3$pval <- 10^(-signa_3656_SK_output_3$`-log10pvalue`)
signa_3656_SK_output_3$pval_FDR <- round(p.adjust(signa_3656_SK_output_3$pval, "fdr"), 3)
signa_3656_SK_output_3_pvalFDR <- which(signa_3656_SK_output_3$pval_FDR < 0.05)

hist(signa_3656_SK_output_4$`-log10pvalue`)
signa_3656_SK_output4_pval1 <- which(signa_3656_SK_output_4$`-log10pvalue`> 1.3,) # pval 0.05
signa_3656_SK_output4_pval2 <- which(signa_3656_SK_output_4$`-log10pvalue`> 2,) # pval 0.01
signa_3656_SK_output4_pval3 <- which(signa_3656_SK_output_4$`-log10pvalue`> 3,) # pval 0.001
signa_3656_SK_output_4[c(signa_3656_SK_output4_pval3),]
dat_3656_out[dat_3656_out$i==659,]
# with FDR correction
signa_3656_SK_output_4$pval <- 10^(-signa_3656_SK_output_4$`-log10pvalue`)
signa_3656_SK_output_4$pval_FDR <- round(p.adjust(signa_3656_SK_output_4$pval, "fdr"), 3)
signa_3656_SK_output_4_pvalFDR <- which(signa_3656_SK_output_4$pval_FDR < 0.05)
which(dat_3656_out[dat_3656_out$i==signa_3656_SK_output_4_pvalFDR,])

x <- NULL
for(vec in c(unique(dat_3656_out$i))){
  print(vec)
  df_to_look <- dat_3656_out[dat_3656_out$i==vec,]
  x[vec] <- df_to_look[1,1]
  signa_3656_SK_output_3$name_snp[vec] <- df_to_look$locus[1]
  if (vec %in% signa_3656_SK_output_3_pvalFDR){
    signa_3656_SK_output_3$FDRoutlier[vec] <- "Yes"
  }
  else{
    signa_3656_SK_output_3$FDRoutlier[vec] <- "No"
  }
}
write.table(signa_3656_SK_output_4, "results_signasel_3656_SK_200Ne")
write.table(signa_3656_SK_output_3, "results_signasel_3656_SK_100Ne")

# VA 287 SNPs 
hist(signa_287output_VA$`-log10pvalue`)
signa_287_VA_output_pval1 <- which(signa_287output_VA$`-log10pvalue`> 1.3,) # pval 0.05
signa_287_VA_output_pval2 <- which(signa_287output_VA$`-log10pvalue`> 2,) # pval 0.01
signa_287_VA_output_pval3 <- which(signa_287output_VA$`-log10pvalue`> 3,) # pval 0.001
signa_287output_VA[c(signa_287_VA_output_pval3),]
dat_out[dat_out$i==277,]

signa_287_VA_output_2_pval1 <- which(signa_287_VA_output_2$`-log10pvalue`> 1.3,) # pval 0.05
signa_287_VA_output_2_pval2 <- which(signa_287_VA_output_2$`-log10pvalue`> 2,) # pval 0.01
signa_287_VA_output_2_pval3 <- which(signa_287_VA_output_2$`-log10pvalue`> 3,) # pval 0.001
signa_287_VA_output_2[c(signa_287_VA_output_2_pval3),]
# with FDR correction
signa_287_VA_output_2$pval <- 10^(-signa_287_VA_output_2$`-log10pvalue`)
signa_287_VA_output_2$pval_FDR <- round(p.adjust(signa_287_VA_output_2$pval, "fdr"), 3)
signa_287_VA_output_2_pvalFDR <- which(signa_287_VA_output_2$pval_FDR < 0.05)

signa_287_VA_output_3_pval1 <- which(signa_287_VA_output_3$`-log10pvalue`> 1.3,) # pval 0.05
signa_287_VA_output_3_pval2 <- which(signa_287_VA_output_3$`-log10pvalue`> 2,) # pval 0.01
signa_287_VA_output_3_pval3 <- which(signa_287_VA_output_3$`-log10pvalue`> 3,) # pval 0.001
signa_287_VA_output_3[c(signa_287_VA_output_3_pval2),]
dat_out[dat_out$i==200,]
# with FDR correction
signa_287_VA_output_3$pval <- 10^(-signa_287_VA_output_3$`-log10pvalue`)
signa_287_VA_output_3$pval_FDR <- round(p.adjust(signa_287_VA_output_3$pval, "fdr"), 3)
signa_287_VA_output_3_pvalFDR <- which(signa_287_VA_output_3$pval_FDR < 0.05)

x <- NULL
for(vec in c(unique(dat_out_277$i))){
  print(vec)
  df_to_look <- dat_out_277[dat_out_277$i==vec,]
  x[vec] <- df_to_look[1,1]
  signa_287_VA_output_2$name_snp[vec] <- df_to_look$locus[1]
  if (vec %in% signa_287_VA_output_2_pvalFDR){
    signa_287_VA_output_2$FDRoutlier[vec] <- "Yes"
  }
  else{
    signa_287_VA_output_2$FDRoutlier[vec] <- "No"
  }
}
write.table(signa_287_VA_output_2, "results_signasel_287_VA_100Ne")
write.table(signa_287_VA_output_3, "results_signasel_287_VA_200Ne")



# VA 3656 SNPs
signa_3656_VA_output# with Ne=130,130,280

hist(signa_3656_VA_output$`-log10pvalue`)
signa_3656_VA_output_pval1 <- which(signa_3656_VA_output$`-log10pvalue`> 1.3,) # pval 0.05
signa_3656_VA_output_pval2 <- which(signa_3656_VA_output$`-log10pvalue`> 2,) # pval 0.01
signa_3656_VA_output_pval3 <- which(signa_3656_VA_output$`-log10pvalue`> 3,) # pval 0.001
signa_3656_VA_output[c(signa_3656_VA_output_pval3),]
dat_3656_out[dat_3656_out$i==3589,]

signa_3656_VA_output_2 # with Ne=100,100
signa_3656_VA_output_2_pval1 <- which(signa_3656_VA_output_2$`-log10pvalue`> 1.3,) # pval 0.05
signa_3656_VA_output_2_pval2 <- which(signa_3656_VA_output_2$`-log10pvalue`> 2,) # pval 0.01
signa_3656_VA_output_2_pval3 <- which(signa_3656_VA_output_2$`-log10pvalue`> 3,) # pval 0.001
signa_3656_VA_output_2[c(signa_3656_VA_output_2_pval3),]
dat_3656_out[dat_3656_out$i==659,]
# with FDR correction
signa_3656_VA_output_2$pval <- 10^(-signa_3656_VA_output_2$`-log10pvalue`)
signa_3656_VA_output_2$pval_FDR <- round(p.adjust(signa_3656_VA_output_2$pval, "fdr"), 3)
signa_3656_VA_output_2_pvalFDR <- which(signa_3656_VA_output_2$pval_FDR < 0.05)


signa_3656_VA_output_3 # with Ne=200,200
signa_3656_VA_output_3_pval1 <- which(signa_3656_VA_output_3$`-log10pvalue`> 1.3,) # pval 0.05
signa_3656_VA_output_3_pval2 <- which(signa_3656_VA_output_3$`-log10pvalue`> 2,) # pval 0.01
signa_3656_VA_output_3_pval3 <- which(signa_3656_VA_output_3$`-log10pvalue`> 3,) # pval 0.001
signa_3656_VA_output_3[c(signa_3656_VA_output_3_pval3),]
dat_3656_out[dat_3656_out$i==3589,]
# with FDR correction
signa_3656_VA_output_3$pval <- 10^(-signa_3656_VA_output_3$`-log10pvalue`)
signa_3656_VA_output_3$pval_FDR <- round(p.adjust(signa_3656_VA_output_3$pval, "fdr"), 3)
signa_3656_VA_output_3_pvalFDR <- which(signa_3656_VA_output_3$pval_FDR < 0.05)


x <- NULL
for(vec in c(unique(dat_3656_out$i))){
  print(vec)
  df_to_look <- dat_3656_out[dat_3656_out$i==vec,]
  x[vec] <- df_to_look[1,1]
  signa_3656_VA_output_3$name_snp[vec] <- df_to_look$locus[1]
  if (vec %in% signa_3656_VA_output_3_pvalFDR){
    signa_3656_VA_output_3$FDRoutlier[vec] <- "Yes"
  }
  else{
    signa_3656_VA_output_3$FDRoutlier[vec] <- "No"
  }
}
write.table(signa_3656_VA_output_2, "results_signasel_3656_VA_100Ne")
write.table(signa_3656_VA_output_3, "results_signasel_3656_VA_200Ne")



# SK 288 SNPs
# Ne 100
hist(signa_288_Skjern4pop_Ne100$`-log10pvalue`)
signa_288_Skjern4pop_Ne100_pval1 <- which(signa_288_Skjern4pop_Ne100$`-log10pvalue`> 1.3,) # pval 0.05
signa_288_Skjern4pop_Ne100_pval2 <- which(signa_288_Skjern4pop_Ne100$`-log10pvalue`> 2,) # pval 0.01
signa_288_Skjern4pop_Ne100_pval3 <- which(signa_288_Skjern4pop_Ne100$`-log10pvalue`> 3,) # pval 0.000001
signa_288_Skjern4pop_Ne100[c(signa_288_Skjern4pop_Ne100_pval1),]
dat_out_288[dat_out_288$i==138,]
# with FDR correction
signa_288_Skjern4pop_Ne100$pval <- 10^(-signa_288_Skjern4pop_Ne100$`-log10pvalue`)
signa_288_Skjern4pop_Ne100$pval_FDR <- round(p.adjust(signa_288_Skjern4pop_Ne100$pval, "fdr"), 3)
signa_288_Skjern4pop_Ne100_pvalFDR <- which(signa_288_Skjern4pop_Ne100$pval_FDR < 0.05)

# Ne 200
hist(signa_288_Skjern4pop_Ne200$`-log10pvalue`)
signa_288_Skjern4pop_Ne200_pval1 <- which(signa_288_Skjern4pop_Ne200$`-log10pvalue`> 1.3,) # pval 0.05
signa_288_Skjern4pop_Ne200_pval2 <- which(signa_288_Skjern4pop_Ne200$`-log10pvalue`> 2,) # pval 0.01
signa_288_Skjern4pop_Ne200_pval3 <- which(signa_288_Skjern4pop_Ne200$`-log10pvalue`> 3,) # pval 0.000001
signa_288_Skjern4pop_Ne200[c(signa_288_Skjern4pop_Ne200_pval1),]
dat_out_288[dat_out_288$i==277,]
# with FDR correction
signa_288_Skjern4pop_Ne200$pval <- 10^(-signa_288_Skjern4pop_Ne200$`-log10pvalue`)
signa_288_Skjern4pop_Ne200$pval_FDR <- round(p.adjust(signa_288_Skjern4pop_Ne200$pval, "fdr"), 3)
signa_288_Skjern4pop_Ne200_pvalFDR <- which(signa_288_Skjern4pop_Ne200$pval_FDR < 0.05)
dat_out_288[dat_out_288$i==250,]


x <- NULL
for(vec in c(unique(dat_out_288$i))){
  print(vec)
  df_to_look <- dat_out_288[dat_out_288$i==vec,]
  x[vec] <- df_to_look[1,1]
  signa_288_Skjern4pop_Ne100$name_snp[vec] <- df_to_look$locus[1]
  if (vec %in% signa_288_Skjern4pop_Ne100_pvalFDR){
    signa_288_Skjern4pop_Ne100$FDRoutlier[vec] <- "Yes"
  }
  else{
    signa_288_Skjern4pop_Ne100$FDRoutlier[vec] <- "No"
  }
}
write.table(signa_288_Skjern4pop_Ne100, "results_signasel_288_SK_100Ne")
write.table(signa_288_Skjern4pop_Ne200, "results_signasel_288_SK_200Ne")




# Function to plot all the allele freq trajectories:
plot_allele_freq <- function(freq, titlepdf, outliers){
  pdf(paste(titlepdf, ".pdf"))
  par(mar=c(2,2,2,2))
  par(mfrow=c(3,2))
  for (snpindex in c(1:287)){
    print(snpindex)
    outlier <- freq[freq$i==snpindex,]
    #outlier <- outlier[c(4,7,9,10),]
    outlier <- outlier[c(9,11,13,14),]
    outlier_list <- outliers
    if(nrow(outlier)>0){
      if(snpindex %in% outlier_list){
    plot(c(1:4), t(outlier[,3]),col="black", type="b",lty = 3, lwd = 1, pch = 19,xaxt = 'n',ylim=c(0,1),
         main=paste("Outlier: SNP:",outlier[1,1]),col.main="red",
         xlab="Temporal sample", ylab= "Minor Allele Frequency (MAF)")
    axis(1, at=1:4,cex.axis=0.6, #labels=c("1913", "1930\n-1935", "1936\n-1939", "1940s", "1955","1990", "1999", "2008", "2015"))
         #labels=c("2001", "2011"))
         #labels=c("2001", "2007","2011", "2019"))
         labels=c("2004", "2010","2012", "2016"))
      } else {
    plot(c(1:4), t(outlier[,3]),col="black", type="b",lty = 3, lwd = 1, pch = 19,xaxt = 'n',ylim=c(0,1),
         main=paste("Outlier: SNP:",outlier[1,1]),
         xlab="Temporal sample", ylab= "Minor Allele Frequency (MAF)")
    axis(1, at=1:4,cex.axis=0.6, #labels=c("1913", "1930\n-1935", "1936\n-1939", "1940s", "1955","1990", "1999", "2008", "2015"))
         #labels=c("2001", "2011"))
         #labels=c("2001", "2007","2011", "2019")) 
         labels=c("2004", "2010","2012", "2016"))
      }
      
    }else{"No of rows is not > 0"}
  }
  dev.off()
}


plot_allele_freq(dat_out, "SK_365SNPs")
plot_allele_freq(dat_out, "VA_287SNPs",c(signa_287_VA_output_2_pval1, signa_287_VA_output_3_pval1))
plot_allele_freq(dat_3656_out, "VA_3656SNPs")
plot_allele_freq(dat_3656_out, "SK_3656SNPs")
plot_allele_freq(dat_out, "SK_288SNPs",c(signa_288_Skjern4pop_Ne100_pval1, signa_288_Skjern4pop_Ne200_pval1))

# Function to plot all the allele freqs at each time point and location:
barplot_allele_freq <- function(freq, titlepdf, outliers){
  pdf(paste(titlepdf, ".pdf"))
  par(mar=c(2,2,2,2))
  par(mfrow=c(3,2))
  outlier_list <- outliers
  
  for (snpindex in c(1:288)){
      print(snpindex)
      outlier <- freq[freq$i==snpindex,]
      #outlier <- outlier[c(1,2,7:14),]
      datapoints <- length(outlier$pop)
      if(nrow(outlier)>0){
      if(snpindex %in% outlier_list){
      plot(c(1:datapoints), t(outlier[,3]),ylim=c(0,1),xaxt = 'n',
         main=paste("Outlier: SNP:",outlier[1,1]),col="#69b3a2",
         las=1,xlab="Temporal sample", ylab= "Minor Allele Frequency (MAF)",col.main="red")
      axis(1, at=1:datapoints,cex.axis=0.5, 
         labels=c("HAT", "SK01\nHAT", "SK01\nHIGHAD", "SK01\nLOWAD", "SK07\nHAT","SK07\nHIGHAD", "SK07\nLOWAD", "SK10\nHIGHAD", "SK10\nLOWAD", "SK19\nLOWAD"))
         #labels=c("HAT97", "HAT12\nHAT", "VA04\nHAT", "VA04\nHIGHAD", "VA04\nLOWAD","VA10\nHIGHAD", "VA10\nLOWAD", "VA12\nHIGHAD", "VA12\nLOWAD", "VA16\nLOWAD"))
         #labels=c("2001", "2011"))
         #labels=rownames(table(dat_out$pop)))
      segments(c(1:datapoints),0, c(1:datapoints),t(outlier[,3]))
      } else {
        plot(c(1:datapoints), t(outlier[,3]),ylim=c(0,1),xaxt = 'n',
           main=paste("Outlier: SNP:",outlier[1,1]),col="#69b3a2",
           las=1,xlab="Temporal sample", ylab= "Minor Allele Frequency (MAF)")
        axis(1, at=1:datapoints,cex.axis=0.5, 
           labels=c("HAT", "SK01\nHAT", "SK01\nHIGHAD", "SK01\nLOWAD", "SK07\nHAT","SK07\nHIGHAD", "SK07\nLOWAD", "SK10\nHIGHAD", "SK10\nLOWAD", "SK19\nLOWAD"))
           #labels=c("HAT97", "HAT12\nHAT", "VA04\nHAT", "VA04\nHIGHAD", "VA04\nLOWAD","VA10\nHIGHAD", "VA10\nLOWAD", "VA12\nHIGHAD", "VA12\nLOWAD", "VA16\nLOWAD"))
      #labels=c("2001", "2011"))
      #labels=rownames(table(dat_out$pop)))
        segments(c(1:datapoints),0, c(1:datapoints),t(outlier[,3]))}}}
  dev.off()}
      

freq <- dat_out
head(freq,20)

barplot_allele_freq(dat_out, "288SNPchip_allPops", c(signa_288_Skjern4pop_Ne100_pval1, signa_288_Skjern4pop_Ne200_pval1))
barplot_allele_freq(dat_out, "287SNPchip_allPops", c(signa_287_VA_output_2_pval1, signa_287_VA_output_3_pval1))
dev.off()


# Function to plot only the allele freqs at each time point and location PLUS the Hatchery:
# Only for the outliers
barplot_allele_freq_2 <- function(freq, titlepdf, outliers){

  #count = 1
  outlier_list <- outliers
  for (snpindex in c(1:3656)){
      print(snpindex)
      outlier <- freq[freq$i==snpindex,]
      #outlier <- outlier[c(1,4,7,9,10),] # SKJERN 288
      outlier <- outlier[c(2,9,11,13,14)] # VARDE 287
      #outlier <- outlier[c(1,2,5,7)] # SKJERN 3287
      #outlier <- outlier[c(1,2,11,13,15)] # VARDE 3287
      datapoints <- length(outlier$pop)
      if(snpindex %in% outlier_list){
        if (outlier[1,3] < 0.5){
              lines(c(1:datapoints), t(outlier[,4]),ylim=c(0,1),xaxt = 'n',type = "b",pch = 18,
               main=paste("Outlier: SNP:",outlier[1,1]),col="#69b3a2",
               las=2,xlab="Temporal sample", ylab= "Minor Allele Frequency (MAF)",col.main="red",lt=2)
            #axis(1, at=1:datapoints,cex.axis=0.5, 
               #labels=c("HAT", "SK2001", "SK2007", "SK2010", "SK2019"))
            #segments(c(1:datapoints),0, c(1:datapoints),t(outlier[,4]))
        } else {
            lines(c(1:datapoints), t(outlier[,3]),ylim=c(0,1),xaxt = 'n',type = "b",pch = 18,
               main=paste("Outlier: SNP:",outlier[1,1]),col="#69b3a2",
               las=2,xlab="Temporal sample", ylab= "Minor Allele Frequency (MAF)",lt=2)
          #axis(1, at=1:datapoints,cex.axis=0.5, 
               #labels=c("HAT", "SK2001", "SK2007", "SK2010", "SK2019"))
          #labels=c("HAT97", "HAT12\nHAT", "VA04\nHAT", "VA04\nHIGHAD", "VA04\nLOWAD","VA10\nHIGHAD", "VA10\nLOWAD", "VA12\nHIGHAD", "VA12\nLOWAD", "VA16\nLOWAD"))
          #labels=c("2001", "2011"))
          #labels=rownames(table(dat_out$pop)))
          #segments(c(1:datapoints),0, c(1:datapoints),t(outlier[,3]))
        }
        }
  }
  }




#pdf(paste("test_2021.pdf"))
#par(mar=c(2,2,2,2))
#par(mfrow=c(3,2))

## Plot for SK 288
plot(c(1:datapoints), c(rep(0, datapoints)),ylim=c(0,1),xaxt = 'n',
     main=paste("Outlier: SNP:",outlier[1,1]),col="#69b3a2",
     las=1,xlab="Temporal sample", ylab= "Allele Frequency (AF)",col.main="red")
axis(1, at=1:datapoints,cex.axis=0.5, 
     labels=c("HAT", "SK2001", "SK2007", "SK2010", "SK2019"))
barplot_allele_freq_2(dat_out_288, "test_2021", signa_288_Skjern4pop_Ne200_pvalFDR)
#which snp names?
# the index are: signa_288_Skjern4pop_Ne200_pvalFDR
dat_out_288[dat_out_288$i==55,]
#dev.off()

## Plot for VA 277
plot(c(1:datapoints), c(rep(0, datapoints)),ylim=c(0,1),xaxt = 'n',
     main=paste("Outlier: SNP:",outlier[1,1]),col="#69b3a2",
     las=1,xlab="Temporal sample", ylab= "Allele Frequency (AF)",col.main="red")
axis(1, at=1:datapoints,cex.axis=0.5, 
     labels=c("HAT1997", "VA2004", "VA2010", "VA2012", "VA2016"))
barplot_allele_freq_2(dat_out_277, "test_2021", signa_287_VA_output_3_pvalFDR)
barplot_allele_freq_2(dat_out_277, "test_2021", signa_287_VA_output_2_pvalFDR) # same result
dat_out_277[dat_out_277$i==133,]

# Plot for SK3656
plot(c(1:datapoints), c(rep(0, datapoints)),ylim=c(0,1),xaxt = 'n',
     main=paste("Outlier: SNP:",outlier[1,1]),col="#69b3a2",
     las=1,xlab="Temporal sample", ylab= "Allele Frequency (AF)",col.main="red")
axis(1, at=1:datapoints,cex.axis=0.5, 
     labels=c("HAT1997", "HAT2012", "SK2001", "SK2011"))
barplot_allele_freq_2(dat_3656_out, "test_2021", signa_3656_SK_output_3_pvalFDR)

# Plot for SK3656
plot(c(1:datapoints), c(rep(0, datapoints)),ylim=c(0,1),xaxt = 'n',
     main=paste("Outlier: SNP:",outlier[1,1]),col="#69b3a2",
     las=1,xlab="Temporal sample", ylab= "Allele Frequency (AF)",col.main="red")
axis(1, at=1:datapoints,cex.axis=0.5, 
     labels=c("HAT1997", "HAT2012", "VA2004", "VA2010", "VA2012"))
barplot_allele_freq_2(dat_3656_out, "test_2021", signa_3656_VA_output_3_pvalFDR)
