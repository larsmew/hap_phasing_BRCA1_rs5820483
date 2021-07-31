# Title: Haplotype phasing of rs455055 and BRCA1 mutations
#
# This scripts tried to identify whether a given BRCA1 mutation is in cis or 
# trans with the imputed SNP rs455055 (that is in strong linkage disequilibrium 
# (LD) with SNP rs5820483) by haplotype phasing a set of genotyped SNPs.

library(tidyverse)
library(proxy)

setwd("~/Documents/PhD/haplotype-project/")

# Read data

# Should contain samples as rows and snps as columns,
# and genotypes (0,1,2) as values
geno_data <- read.table("geno.txt", header = T) 


pheno_data <- read.table("pheno.txt", header = T) # 
snp_coords <- read.table("snp_coords_chr17.txt", header = T) # 

# BRCA mutation coordinates
brca1_mutation_coords <- read.table("input/139_ouh_june_2017/BRCA1_mut_position.txt", header = T, sep = "\t")

# geno_data <- read.table("pujana/38_Pujana_OD7-BRCA1/38_brca1_plink_format.txt", header = T)
# pheno_data <- read.table("pujana/38_Pujana_OD7-BRCA1/iCOGS_pheno_export_B1_261112.txt", header = T, sep = "\t", quote = "", fill = T, stringsAsFactors = F)
# snp_coords = read.table("pujana/38_Pujana_OD7-BRCA1/38_pujana_coords.txt", header = T) %>% arrange(position_b37)
# dosage = read.csv("pujana/38_Pujana_OD7-BRCA1/brca1_chr17_41165k_41355k_dosages_filtered.csv", header = T)
# dosage$chr17.41247122.I


getBRCA1info <- function(mut){
    gene <<- "BRCA1"
    mut_start <<- brca1_mutation_coords[which(brca1_mutation_coords$Mut1HGVS == mut), 2]
    mut_stop <<- brca1_mutation_coords[which(brca1_mutation_coords$Mut1HGVS == mut), 3]
    mut_middle <<- (brca_stop - brca_start)/2 + brca_start
    brca_chr <<- 17
}

extractSamples <- function(pheno, geno, chr_coords){

    # Extract samples with given brca mutation
    index <- match(pheno$Onc_ID, geno$SNP)
    matched <- geno[index, ]
    
    # Make sure the snps (columns) are ordered according to physical position
    # And remove snps with no coordinates in snp_coords
    snp_col <- which(names(matched)=="SNP")
    matched <- matched[, c(snp_col, na.omit(match(snp_coords$SNP, colnames(matched))))]
    
    return(matched)
}

firstBreakDist <- function(geno, pheno, snp_coords){
    
    getSnpDist <- function(snp_coords, geno){
        snp_dist = filter(snp_coords, SNP %in% names(geno)) %>%
            arrange(position_b37) %>%
            mutate(brca_distance = as.integer(as.character(position_b37)) - mut_middle)
        return(snp_dist)
    }
    
    distanceSimilarity <- function(x,y){
        # Find homozygote break
        condition = (abs(x-y) == 2 | x == -1 | y == -1)
        # Split vector according to positive or negative distance to brca gene
        cond_pos = condition[(snp_split+1):nrow(snp_dist)]
        cond_neg = condition[1:snp_split.left]
        # get length (in Kb) to brca gene for first homozygous break on either side
        min_index=min(which(cond_pos==T), length(cond_pos))
        pos.len = snp_dist[snp_split+min_index,"brca_distance"] / 1000000
        max_index=max(which(cond_neg==T), 1)
        neg.len = abs(snp_dist[max_index,"brca_distance"]) / 1000000
        # Compute similarity based on length with no homozygous break
        similarity = 1/(pos.len+neg.len)
        return(similarity)
    }
    
    # Define snp distances to BRCA mutation
    snp_dist <<- getSnpDist(snp_coords, geno)
    # Exclude SNPs in mutation, likely not any
    if (removeSNPsInGene){
        snp_dist.dist <<- filter(snp_dist, position_b37 < mut_start | position_b37 > mut_stop)
    }
    # Define position splitting snps on either side of BRCA mutation
    snp.split.left <<- nrow(filter(snp_dist, brca_distance <= 0))
    snp.split <<- nrow(filter(snp_dist, brca_distance < 0))
        
    dst = proxy::dist(geno[,2:ncol(geno)], method=distanceSimilarity)
    
    return(dst)
}


# Genotype values represents the following
#   0: homozygous REF
#   1: heterozygous
#   2: homozygous ALT
genotype_description = c("2"="homref","1"="het","0"="homalt")

brca1_muts <- count(pheno_individual, Mut1HGVS) %>% arrange(desc(n))

height_cutoff = 7


result <- map_dfr(.x = brca1_muts$Mut1HGVS, .f = function(mut){
    getBRCA1info(mut)
    
    # subset pheno info
    pheno_subset <<- filter(pheno_data, Mut1HGVS %in% mut)
    # subset geno data
    geno_subset = extractSamples(pheno_subset, geno_data, chr_coords)
    n_samples = nrow(geno_subset)
    
    if (n_samples < 3){
        pheno_subset$cluster_groups = 1
        k = 1
    } else {
        dir.create("distMatrices", showWarnings = F)
        dst <- firstBreakDist(fam=fam, matched_custom = geno_subset,
                              filename = paste0("distMatrices/dist_",mut_name,".RData"))
        hc <- hclust(dst,"ward.D2")
        k = sum(hc$height>height_cutoff)+1
        if (k>n_samples){
            cluster_groups <<- cutree(hc, n_samples)
            pheno_subset$cluster_groups <- cluster_groups
            k = n_samples
        } else {
            cluster_groups <<- cutree(hc, k)
            pheno_subset$cluster_groups <- cluster_groups
        }
        pheno_subset <<- pheno_subset
    }
    
    result = map_dfr(1:k, function(founder_group){
        geno_founder = geno_subset[geno_subset$SNP %in% subset(pheno_subset, cluster_groups == founder_group)$Onc_ID,]
        
        # Note: that for rs455055 the ALT allele (base A) is the most common, 
        # thus the data is encoded as 0 = homozygote ALT and 2 = homozygote REF
        homref = sum(geno_founder$rs455055 == 2)
        het = sum(geno_founder$rs455055 == 1)
        homalt = sum(geno_founder$rs455055 == 0)
        
        res = data.frame("ICOGS_Person_ID"=geno_founder$SNP, "Mutation"=mut, 
                         "Founder"=founder_group,"Homref"=homref,
                         "Het"=het,"Homalt"=homalt)
        
        # We define the SNP rs455055 to be in trans with the BRCA1 mutation in 
        # a founder group if at least 5 times the number of individuals have 
        # homalt than homref.
        if (homalt>=(5*homref) && homalt > 0){
            # Set samples to be in trans
            res$Splicing = paste("trans", genotype_description[as.character(geno_founder$rs455055)])
            # Individuals in the founder groups not having ALT allele in the 
            # rs455055 SNP cannot de be determined in trans - either they do not
            # belong to the founder group or the ALT allele have been mutated 
            # or recombined away.
            res$Splicing[which(geno_founder$rs455055 == 2)] = "undetermined"
        
        # We define the SNP rs455055 to be in cis with the BRCA1 mutation in 
        # a founder group if at least 5 times the number of individuals have 
        # homref than homalt.
        } else if (homref>=(5*homalt) && homref > 0){
            # Set samples to be in cis
            res$Splicing = paste("cis", genotype_description[as.character(geno_founder$rs455055)])
            # Individuals in the founder groups not having REF allele in the 
            # rs455055 SNP cannot de be determined in cis - either they do not
            # belong to the founder group or the REF allele have been mutated 
            # or recombined away.
            res$Splicing[which(geno_founder$rs455055 == 0)] = "undetermined"
        
        # Otherwise, unable to determine the phasing for the founder group
        } else {
            res$Splicing = "undetermined"
        }
        return(res)
    })
    return(result)
})
summary_table = result %>% count(Mutation,Founder,Splicing,Homref,Het,Homalt, name = "Num_samples")
write.table(result[,c(1:3,7)], file = paste0("output/Splicing_brca1_rs455055_withFounderAnalysis.txt"), quote = F, row.names = F, sep = "\t")
write.table(summary_table, file = paste0("output/Summary_splicing_brca1_rs455055_withFounderAnalysis.txt"), quote = F, row.names = F, sep = "\t")


