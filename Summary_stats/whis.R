#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
suppressMessages(library('PopGenome', quietly = TRUE))

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

vcf <- args[1]
sl <- 50000000

command <- paste("zgrep -v '^#' ", vcf, " | wc -l", sep = "")
nsnp <- as.numeric(system(command, intern = TRUE))

command <- paste("zcat ", vcf, " | awk 'NR==2' ", sep = "")
pline <- system(command, intern = TRUE)
priors <- as.numeric(strsplit(pline, split = " ")[[1]][-1])

msg.out <- capture.output(suppressMessages(GENOME.class <- readVCF(vcf, frompos = 1, topos = nsnp, tid = "1", numcols = 10000)))
#GENOME.class <- readVCF(vcf, frompos = 1, topos = nsnp, tid = "1", numcols = 10000)

me <- paste("tsk_", 0:13, sep = "")
mc <- paste("tsk_", 14:19, sep = "")
da <- paste("tsk_", 20:30, sep = "")
dc <- paste("tsk_", 31:44, sep = "")
de <- paste("tsk_", 45:60, sep = "")

msg.out <- capture.output(suppressMessages(GENOME.class <- set.populations(GENOME.class,list(me,mc,da,dc,de), diploid=TRUE)))
#GENOME.class <- set.populations(GENOME.class,list(me,mc,da,dc,de), diploid=TRUE)
#get.sum.data(GENOME.class)
#get.individuals(GENOME.class)

# Neutrality
msg.out <- capture.output(suppressMessages(GENOME.class <- neutrality.stats(GENOME.class)))
TD_me <- get.neutrality(GENOME.class)[[1]][1] # Tajimas D
TD_mc <- get.neutrality(GENOME.class)[[2]][1]
TD_da <- get.neutrality(GENOME.class)[[3]][1]
TD_dc <- get.neutrality(GENOME.class)[[4]][1]
TD_de <- get.neutrality(GENOME.class)[[5]][1]
segsites_me <- get.neutrality(GENOME.class)[[1]][2]/sl # Number of segregating sites per number of simulated bp
segsites_mc <- get.neutrality(GENOME.class)[[2]][2]/sl
segsites_da <- get.neutrality(GENOME.class)[[3]][2]/sl
segsites_dc <- get.neutrality(GENOME.class)[[4]][2]/sl
segsites_de <- get.neutrality(GENOME.class)[[5]][2]/sl
neut <- c(TD_me, TD_mc, TD_da, TD_dc, TD_de, segsites_me, segsites_mc, segsites_da, segsites_dc, segsites_de)

# Diversity within Populations
msg.out <- capture.output(suppressMessages(GENOME.class <- diversity.stats(GENOME.class)))
nuc.diversity.within_me <- get.diversity(GENOME.class)[[1]][1]/sl
hap.diversity.within_me <- get.diversity(GENOME.class)[[1]][2]/sl
nuc.diversity.within_mc <- get.diversity(GENOME.class)[[2]][1]/sl
hap.diversity.within_mc <- get.diversity(GENOME.class)[[2]][2]/sl
nuc.diversity.within_da <- get.diversity(GENOME.class)[[3]][1]/sl
hap.diversity.within_da <- get.diversity(GENOME.class)[[3]][2]/sl
nuc.diversity.within_dc <- get.diversity(GENOME.class)[[4]][1]/sl
hap.diversity.within_dc <- get.diversity(GENOME.class)[[4]][2]/sl
nuc.diversity.within_de <- get.diversity(GENOME.class)[[5]][1]/sl
hap.diversity.within_de <- get.diversity(GENOME.class)[[5]][2]/sl


# Diversity between Populations
msg.out <- capture.output(suppressMessages(GENOME.class <- diversity.stats.between(GENOME.class)))
nuc.diversity.between <- GENOME.class@nuc.diversity.between/sl
nuc.diversity.between_memc <- nuc.diversity.between[1]
nuc.diversity.between_meda <- nuc.diversity.between[2]
nuc.diversity.between_medc <- nuc.diversity.between[3]
nuc.diversity.between_mede <- nuc.diversity.between[4]
nuc.diversity.between_mcda <- nuc.diversity.between[5]
nuc.diversity.between_mcdc <- nuc.diversity.between[6]
nuc.diversity.between_mcde <- nuc.diversity.between[7]
nuc.diversity.between_dadc <- nuc.diversity.between[8]
nuc.diversity.between_dade <- nuc.diversity.between[9]
nuc.diversity.between_dcde <- nuc.diversity.between[10]

# Calculate fixed shared sites
msg.out <- capture.output(suppressMessages(GENOME.class <- calc.fixed.shared(GENOME.class)))
p.fixed.sites <- GENOME.class@n.fixed.sites/sl
p.fixed.sites_memc <- p.fixed.sites[1]
p.fixed.sites_meda <- p.fixed.sites[2]
p.fixed.sites_medc <- p.fixed.sites[3]
p.fixed.sites_mede <- p.fixed.sites[4]
p.fixed.sites_mcda <- p.fixed.sites[5]
p.fixed.sites_mcdc <- p.fixed.sites[6]
p.fixed.sites_mcde <- p.fixed.sites[7]
p.fixed.sites_dadc <- p.fixed.sites[8]
p.fixed.sites_dade <- p.fixed.sites[9]
p.fixed.sites_dcde <- p.fixed.sites[10]
p.shared.sites <- GENOME.class@n.shared.sites/sl
p.shared.sites_memc <- p.shared.sites[1]
p.shared.sites_meda <- p.shared.sites[2]
p.shared.sites_medc <- p.shared.sites[3]
p.shared.sites_mede <- p.shared.sites[4]
p.shared.sites_mcda <- p.shared.sites[5]
p.shared.sites_mcdc <- p.shared.sites[6]
p.shared.sites_mcde <- p.shared.sites[7]
p.shared.sites_dadc <- p.shared.sites[8]
p.shared.sites_dade <- p.shared.sites[9]
p.shared.sites_dcde <- p.shared.sites[10]
p.monomorphic.sites <- GENOME.class@n.monomorphic.sites/sl
p.monomorphic.sites_memc <- p.monomorphic.sites[1]
p.monomorphic.sites_meda <- p.monomorphic.sites[2]
p.monomorphic.sites_medc <- p.monomorphic.sites[3]
p.monomorphic.sites_mede <- p.monomorphic.sites[4]
p.monomorphic.sites_mcda <- p.monomorphic.sites[5]
p.monomorphic.sites_mcdc <- p.monomorphic.sites[6]
p.monomorphic.sites_mcde <- p.monomorphic.sites[7]
p.monomorphic.sites_dadc <- p.monomorphic.sites[8]
p.monomorphic.sites_dade <- p.monomorphic.sites[9]
p.monomorphic.sites_dcde <- p.monomorphic.sites[10]

div <- c(nuc.diversity.within_me, hap.diversity.within_me, nuc.diversity.within_mc, hap.diversity.within_mc, nuc.diversity.within_da, hap.diversity.within_da, nuc.diversity.within_dc, hap.diversity.within_dc, nuc.diversity.within_de, hap.diversity.within_de, nuc.diversity.between_memc, nuc.diversity.between_meda, nuc.diversity.between_medc, nuc.diversity.between_mede, nuc.diversity.between_mcda, nuc.diversity.between_mcdc, nuc.diversity.between_mcde, nuc.diversity.between_dadc, nuc.diversity.between_dade, nuc.diversity.between_dcde, p.fixed.sites_memc, p.fixed.sites_meda, p.fixed.sites_medc, p.fixed.sites_mede, p.fixed.sites_mcda, p.fixed.sites_mcdc, p.fixed.sites_mcde, p.fixed.sites_dadc, p.fixed.sites_dade, p.fixed.sites_dcde, p.shared.sites_memc, p.shared.sites_meda, p.shared.sites_medc, p.shared.sites_mede, p.shared.sites_mcda, p.shared.sites_mcdc, p.shared.sites_mcde, p.shared.sites_dadc, p.shared.sites_dade, p.shared.sites_dcde, p.monomorphic.sites_memc, p.monomorphic.sites_meda, p.monomorphic.sites_medc, p.monomorphic.sites_mede, p.monomorphic.sites_mcda, p.monomorphic.sites_mcdc, p.monomorphic.sites_mcde, p.monomorphic.sites_dadc, p.monomorphic.sites_dade, p.monomorphic.sites_dcde )
#Fst
msg.out <- capture.output(suppressMessages(GENOME.class <- F_ST.stats(GENOME.class)))
nuc.F_ST.pairwise <- get.F_ST(GENOME.class, pairwise = TRUE)[1]
nuc.F_ST_memc <- nuc.F_ST.pairwise[[1]][1]
nuc.F_ST_meda <- nuc.F_ST.pairwise[[1]][2]
nuc.F_ST_medc <- nuc.F_ST.pairwise[[1]][3]
nuc.F_ST_mede <- nuc.F_ST.pairwise[[1]][4]
nuc.F_ST_mcda <- nuc.F_ST.pairwise[[1]][5]
nuc.F_ST_mcdc <- nuc.F_ST.pairwise[[1]][6]
nuc.F_ST_mcde <- nuc.F_ST.pairwise[[1]][7]
nuc.F_ST_dadc <- nuc.F_ST.pairwise[[1]][8]
nuc.F_ST_dade <- nuc.F_ST.pairwise[[1]][9]
nuc.F_ST_dcde <- nuc.F_ST.pairwise[[1]][10]

hap.F_ST.pairwise <- get.F_ST(GENOME.class, pairwise = TRUE)[2]
hap.F_ST_memc <- hap.F_ST.pairwise[[1]][1]
hap.F_ST_meda <- hap.F_ST.pairwise[[1]][2]
hap.F_ST_medc <- hap.F_ST.pairwise[[1]][3]
hap.F_ST_mede <- hap.F_ST.pairwise[[1]][4]
hap.F_ST_mcda <- hap.F_ST.pairwise[[1]][5]
hap.F_ST_mcdc <- hap.F_ST.pairwise[[1]][6]
hap.F_ST_mcde <- hap.F_ST.pairwise[[1]][7]
hap.F_ST_dadc <- hap.F_ST.pairwise[[1]][8]
hap.F_ST_dade <- hap.F_ST.pairwise[[1]][9]
hap.F_ST_dcde <- hap.F_ST.pairwise[[1]][10]

Nei.G_ST.pairwise <- get.F_ST(GENOME.class, pairwise = TRUE)[3]
Nei.G_ST_memc <- Nei.G_ST.pairwise[[1]][1]
Nei.G_ST_meda <- Nei.G_ST.pairwise[[1]][2]
Nei.G_ST_medc <- Nei.G_ST.pairwise[[1]][3]
Nei.G_ST_mede <- Nei.G_ST.pairwise[[1]][4]
Nei.G_ST_mcda <- Nei.G_ST.pairwise[[1]][5]
Nei.G_ST_mcdc <- Nei.G_ST.pairwise[[1]][6]
Nei.G_ST_mcde <- Nei.G_ST.pairwise[[1]][7]
Nei.G_ST_dadc <- Nei.G_ST.pairwise[[1]][8]
Nei.G_ST_dade <- Nei.G_ST.pairwise[[1]][9]
Nei.G_ST_dcde <- Nei.G_ST.pairwise[[1]][10]


fst <- c(nuc.F_ST_memc, nuc.F_ST_meda, nuc.F_ST_medc, nuc.F_ST_mede, nuc.F_ST_mcda, nuc.F_ST_mcdc, nuc.F_ST_mcde, nuc.F_ST_dadc, nuc.F_ST_dade, nuc.F_ST_dcde, hap.F_ST_memc, hap.F_ST_meda, hap.F_ST_medc, hap.F_ST_mede, hap.F_ST_mcda, hap.F_ST_mcdc, hap.F_ST_mcde, hap.F_ST_dadc, hap.F_ST_dade, hap.F_ST_dcde, Nei.G_ST_memc, Nei.G_ST_meda, Nei.G_ST_medc, Nei.G_ST_mede, Nei.G_ST_mcda, Nei.G_ST_mcdc, Nei.G_ST_mcde, Nei.G_ST_dadc, Nei.G_ST_dade, Nei.G_ST_dcde)

# 3D-SFS's of the two mys pops with each one of the dav pops
msg.out <- capture.output(suppressMessages(GENOME.class <- detail.stats(GENOME.class,site.spectrum=TRUE,site.FST=TRUE)))
results <- get.detail(GENOME.class, biallelic.structure=TRUE) 
allele_Freqs <- GENOME.class@region.stats@minor.allele.freqs[[1]] 

df <- expand.grid(0:length(me), 0:length(mc), 0:length(da))
df$freq <- rep(0,nrow(df))
for (i in 1:ncol(allele_Freqs)) {
  freq_pop1 <- allele_Freqs[1,i] * length(me) * 2
  freq_pop2 <- allele_Freqs[2,i] * length(mc) * 2
  freq_pop3 <- allele_Freqs[3,i] * length(da) * 2
  idx <- intersect(intersect(which(df$Var1 == freq_pop1), which(df$Var2 == freq_pop2)), which(df$Var3 == freq_pop3))
  df$freq[idx] <- df$freq[idx] + 1
}
sfs3D_memcda <- df$freq / sl

df <- expand.grid(0:length(me), 0:length(mc), 0:length(dc))
df$freq <- rep(0,nrow(df))
for (i in 1:ncol(allele_Freqs)) {
  freq_pop1 <- allele_Freqs[1,i] * length(me) * 2
  freq_pop2 <- allele_Freqs[2,i] * length(mc) * 2
  freq_pop3 <- allele_Freqs[4,i] * length(dc) * 2
  idx <- intersect(intersect(which(df$Var1 == freq_pop1), which(df$Var2 == freq_pop2)), which(df$Var3 == freq_pop3))
  df$freq[idx] <- df$freq[idx] + 1
}
sfs3D_memcdc <- df$freq / sl

df <- expand.grid(0:length(me), 0:length(mc), 0:length(de))
df$freq <- rep(0,nrow(df))
for (i in 1:ncol(allele_Freqs)) {
  freq_pop1 <- allele_Freqs[1,i] * length(me) * 2
  freq_pop2 <- allele_Freqs[2,i] * length(mc) * 2
  freq_pop3 <- allele_Freqs[5,i] * length(de) * 2
  idx <- intersect(intersect(which(df$Var1 == freq_pop1), which(df$Var2 == freq_pop2)), which(df$Var3 == freq_pop3))
  df$freq[idx] <- df$freq[idx] + 1
}
sfs3D_memcde <- df$freq / sl

# Output
out <- c(priors, neut, div, fst, sfs3D_memcda, sfs3D_memcdc, sfs3D_memcde)
out <- paste(out, collapse = "\t")
write(out, stdout())
##################
## Higher dimensional SFS from minor allele frequencies: ##
# Here is a piece of code that calculates the folded 2D-sfs:
#
#df <- expand.grid(0:length(pop1), 0:length(pop2))
#df$freq <- rep(0,nrow(df))
#for (i in 1:ncol(allele_Freqs)) {
#  freq_pop1 <- allele_Freqs[1,i] * length(pop1) * 2
#  freq_pop2 <- allele_Freqs[2,i] * length(pop2) * 2
#  idx <- intersect(which(df$Var1 == freq_pop1), which(df$Var2 == freq_pop2))
#  df$freq[idx] <- df$freq[idx] + 1
#}
#sfs2D <- df$freq
#
# And here the equivalent for the 3D-sfs:
#
#df <- expand.grid(0:length(pop1), 0:length(pop2), 0:length(pop3))
#df$freq <- rep(0,nrow(df))
#for (i in 1:ncol(allele_Freqs)) {
#  freq_pop1 <- allele_Freqs[1,i] * length(pop1) * 2
#  freq_pop2 <- allele_Freqs[2,i] * length(pop2) * 2
#  freq_pop3 <- allele_Freqs[3,i] * length(pop3) * 2
#  idx <- intersect(intersect(which(df$Var1 == freq_pop1), which(df$Var2 == freq_pop2)), which(df$Var3 == freq_pop3))
#  df$freq[idx] <- df$freq[idx] + 1
#}
#sfs3D <- df$freq


