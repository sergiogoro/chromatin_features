#
# GENES
#
# Dos datasets:
#   #A) EXON.SHORTINTRONS.n
#   B) GENE.SHORTINTRONS.n
#
# Filtros:
#   B) Genes
#       B.1)
#           data <- na.omit(data) #para el modelo mixto
#

library(car)
library(lattice)
library(ppcor)

rm(list = ls())
mypath <- "/home/sergio/chromatin/analysis/R/chromatin_features/"
setwd(mypath)

# Load datasets
data <- read.table(file="GENE.SHORTINTRONS.n",header=TRUE,sep="\t")
#nrow(data)
data <- subset(data, chromosome != "X")                     #AUTOSOMES
data <- na.omit(data) #para el modelo mixto

### PopGen Summary Statistics ###

pi4f    <- data[["pi_4f"]] / data[["mdmel_4f"]]
pi0f    <- data[["pi_0f"]] / data[["mdmel_0f"]]
piins   <- data[["pi_ins"]] / data[["mdmel_ins"]]
t4f     <- data[["seg_4f"]] / data[["mdmel_4f"]]
t0f     <- data[["seg_0f"]] / data[["mdmel_0f"]]
tins    <- data[["seg_ins"]] / data[["mdmel_ins"]]
k4f     <- ( data[["div_4f"]] + 1 ) / ( data[["mdyak_4f"]] + 1 )
k0f     <- ( data[["div_0f"]] + 1 ) / ( data[["mdyak_0f"]] + 1 )
kins    <- ( ( data$div_ins + 1 ) / ( data[["mdyak_ins"]] + 1 ) )

constraint_4f   <- t0f/t4f
constraint_ins  <- t0f/tins
omega_4f        <- (k0f/k4f)            #Para boxplots
omega_ins       <- (k0f/kins)           #Para boxplots
m               <- data[["mdmel_0f"]] + data[["mdmel_4f"]] + data[["mdmel_2f"]]


### Summary ###
		
#### GENE ###
    #If exists previous file, move it to .bak
    if (file.exists("OUT_GENES-autosomes")) {
        file.rename("OUT_GENES-autosomes", "OUT_GENES-autosomes.bak")
        file.remove("OUT_GENES-autosomes")
    }
    #Remove previous file
    cat("###GENES-autosomes\n\n", file="OUT_GENES-autosomes", append=T)

    #Feature 1
	### Messenger Complexity ###
    cat("##Feature1: Messenger Complexity\n", file="OUT_GENES-autosomes", append=T)

	mcomp <- (data$FBtrs_per_gene/data$num_cds)
    #summary(mcomp)
	mcomp2 = cut(mcomp,c(0,.25,0.4,0.5))
	mcomp5 = cut(mcomp,quantile(mcomp,(0:5)/5))
    #summary(mcomp5)
    cat("#tapply sum m~mcomp5", file="OUT_GENES-autosomes", append=T)
	write.table( tapply(m, mcomp5, sum), file="OUT_GENES-autosomes", quote=T, row.names=T, append=T )

    png( "GENES-autosomes-FEAT1MessengerComplexity_Histogram-mcomp.png", width = 1920, height = 1080 )
	hist(mcomp,breaks=c(50))
    dev.off()
    png( "GENES-autosomes-FEAT1MessengerComplexity_Histogram-contextDistance.png", width = 1920, height = 1080 )
	hist(data$context_distance,breaks=c(100))
    dev.off()

    cat( quantile( mcomp, ( 0:5 ) / 5 )[1], labels="Quantile 0%", file="OUT_GENES-autosomes", fill=T, append=T )
    cat( quantile( mcomp, ( 0:5 ) / 5 )[2], labels="Quantile 20%", file="OUT_GENES-autosomes", fill=T, append=T )
    cat( quantile( mcomp, ( 0:5 ) / 5 )[3], labels="Quantile 40%", file="OUT_GENES-autosomes", fill=T, append=T )
    cat( quantile( mcomp, ( 0:5 ) / 5 )[4], labels="Quantile 60%", file="OUT_GENES-autosomes", fill=T, append=T )
    cat( quantile( mcomp, ( 0:5 ) / 5 )[5], labels="Quantile 80%", file="OUT_GENES-autosomes", fill=T, append=T )
    cat( quantile( mcomp, ( 0:5 ) / 5 )[6], labels="Quantile 100%", file="OUT_GENES-autosomes", fill=T, append=T )
    #quantile(mcomp,(0:5)/5)

    cat(median( mcomp ), file="OUT_GENES-autosomes", fill=T, labels="Median", append=T)
    cat(mean( mcomp ), file="OUT_GENES-autosomes", fill=T, labels="Mean", append=T)
    cat(sd( mcomp ), file="OUT_GENES-autosomes", fill=T, labels="SD", append=T)
	

    png( "GENES-autosomes-FEAT1MessengerComplexity_omega4f-mcomp5.png", width = 1920, height = 1080 )
	boxplot(omega_4f~mcomp5,outline=F,xlab="Number of Transcripts/Number of Exons",ylab="Ka/Ks") #aquellos genes con mas transcritos q exones tienen menor omega q aquellos con menos transcritos q exones.
	abline(h=median(omega_4f),col="black")
    dev.off()
    png( "GENES-autosomes-FEAT1MessengerComplexity_omegains-mcomp5.png", width = 1920, height = 1080 )
	boxplot(omega_ins~mcomp5,outline=F,xlab="Number of Transcripts/Number of Exons",ylab="Ka/Kins") #aquellos genes con mas transcritos q exones tienen menor omega q aquellos con menos transcritos q exones.
	abline(h=median(omega_ins),col="black")	
	        kruskal.test(omega_4f~mcomp5)
	        kruskal.test(omega_ins~mcomp5)
    dev.off()
	
	### chromosome and chromatin ###
    cat("\n", file="OUT_GENES-autosomes", append=T)
    cat("##Feature2: Chromosome and chromatin\n", file="OUT_GENES-autosomes", append=T)

    cat("#tapply sum m~chromosome", file="OUT_GENES-autosomes", append=T)
	write.table( tapply(m, data$chromosome, sum), file="OUT_GENES-autosomes", quote=T, row.names=T, append=T )
    cat("#tapply sum m~chromatin", file="OUT_GENES-autosomes", append=T)
	write.table( tapply(m, data$chromatin, sum), file="OUT_GENES-autosomes", quote=T, row.names=T, append=T )


    cat( kruskal.test( omega_4f ~ data[["chromosome"]] )$statistic, file = "OUT_GENES-autosomes", fill=T, labels="Omega_4f K-W chi-squared", append=T  )
    cat( kruskal.test( omega_4f ~ data[["chromosome"]] )$p.value, file = "OUT_GENES-autosomes", fill=T, labels="Omega_4f K-W p-value", append=T  )
    cat( kruskal.test( omega_ins ~ data[["chromosome"]] )$statistic, file = "OUT_GENES-autosomes", fill=T, labels="Omega_ins K-W chi-squared", append=T  )
    cat( kruskal.test( omega_ins ~ data[["chromosome"]] )$p.value, file = "OUT_GENES-autosomes", fill=T, labels="Omega_ins K-W p-value", append=T  )

    cat( kruskal.test( omega_4f ~ data[["chromatin"]] )$statistic, file = "OUT_GENES-autosomes", fill=T, labels="Omega_4f K-W chi-squared", append=T  )
    cat( kruskal.test( omega_4f ~ data[["chromatin"]] )$p.value, file = "OUT_GENES-autosomes", fill=T, labels="Omega_4f K-W p-value", append=T  )
    cat( kruskal.test( omega_ins ~ data[["chromatin"]] )$statistic, file = "OUT_GENES-autosomes", fill=T, labels="Omega_ins K-W chi-squared", append=T  )
    cat( kruskal.test( omega_ins ~ data[["chromatin"]] )$p.value, file = "OUT_GENES-autosomes", fill=T, labels="Omega_ins K-W p-value", append=T  )


	### Transcripts ###
    cat( "\n", file="OUT_GENES-autosomes", append=T )
    cat( "##Feature3: Transcripts\n", file="OUT_GENES-autosomes", append=T )

    cat("#tapply sum m~FBtrs_per_gene", file="OUT_GENES-autosomes", append=T)
	write.table( tapply( m, data$FBtrs_per_gene, sum ), file="OUT_GENES-autosomes", quote=T, row.names=T, append=T )

	data <-subset(data, FBtrs_per_gene < 11)
	FBtrs = cut(data[["FBtrs_per_gene"]],c(0,1,2,5,75))
    #summary(FBtrs)
    cat("#tapply sum m~FBtrs", file="OUT_GENES-autosomes", append=T)
	write.table( tapply( m, FBtrs, sum ), file="OUT_GENES-autosomes", quote=T, row.names=T, append=T )
	
    png( "GENES-autosomes-FEAT3Transcripts_omega4f-FBtrs.png", width = 1920, height = 1080 )
	boxplot(omega_4f ~ FBtrs, outline=F,xlab="Number of Transcripts/Gene",ylab="Ka/Ks")
	abline(h=median(omega_4f),col="black")
    dev.off()

    png( "GENES-autosomes-FEAT3Transcripts_omegains-FBtrs.png", width = 1920, height = 1080 )
	boxplot(omega_ins ~ FBtrs, outline=F,xlab="Number of Transcripts/Gene",ylab="Ka/Kins")
	abline(h=median(omega_4f),col="black")
    dev.off()

    png( "GENES-autosomes-FEAT3Transcripts_pi4f-FBtrs.png", width = 1920, height = 1080 )
	boxplot(pi4f ~ FBtrs, outline=F, xlab="Number of Transcripts/Gene", ylab="Ka/Ks")
    dev.off()

    png( "GENES-autosomes-FEAT3Transcripts_piins-FBtrs.png", width = 1920, height = 1080 )
	boxplot(piins ~ FBtrs, outline=F, xlab="Number of Transcripts/Gene", ylab="Ka/Kins")
    dev.off()


    cat( cor.test(omega_4f, FBtrs, method="spearman")$estimate, file = "OUT_GENES-autosomes",fill=T, labels="Omega_4f Spearman Rho", append=T)
    cat( cor.test(omega_4f, FBtrs, method="spearman")$p.value, file = "OUT_GENES-autosomes", fill=T, labels="Omega_4f Spearman p-value", append=T)
    cat( cor.test(omega_ins, FBtrs, method="spearman")$estimate, file = "OUT_GENES-autosomes" , fill=T, labels="Omega_ins Spearman Rho", append=T)
    cat( cor.test(omega_ins, FBtrs, method="spearman")$p.value, file = "OUT_GENES-autosomes" , fill=T, labels="Omega_ins Spearman p-value", append=T)

    cat( cor.test(pi4f, FBtrs, method="spearman")$estimate, file = "OUT_GENES-autosomes",fill=T, labels="Omega_4f Spearman Rho", append=T)          #-0.069
    cat( cor.test(pi4f, FBtrs, method="spearman")$p.value, file = "OUT_GENES-autosomes", fill=T, labels="Omega_4f Spearman p-value", append=T)      #-0.069
    cat( cor.test(piins, FBtrs, method="spearman")$estimate, file = "OUT_GENES-autosomes" , fill=T, labels="Omega_ins Spearman Rho", append=T)      #ns
    cat( cor.test(piins, FBtrs, method="spearman")$p.value, file = "OUT_GENES-autosomes" , fill=T, labels="Omega_ins Spearman p-value", append=T)   #ns

    cat( kruskal.test( omega_4f ~ FBtrs )$statistic, file = "OUT_GENES-autosomes", fill=T, labels="Omega_4f K-W chi-squared", append=T  )
    cat( kruskal.test( omega_4f ~ FBtrs )$p.value, file = "OUT_GENES-autosomes", fill=T, labels="Omega_4f K-W p-value", append=T  )
    cat( kruskal.test( omega_ins ~ FBtrs )$statistic, file = "OUT_GENES-autosomes", fill=T, labels="Omega_ins K-W chi-squared", append=T  )
    cat( kruskal.test( omega_ins ~ FBtrs )$p.value, file = "OUT_GENES-autosomes", fill=T, labels="Omega_ins K-W p-value", append=T  )


	### Number of Exons ###
    cat( "\n", file="OUT_GENES-autosomes", append=T )
    cat( "##Feature4: Number of Exons\n", file="OUT_GENES-autosomes", append=T )

    cat("#tapply sum m~num_CDS", file="OUT_GENES-autosomes", append=T)
	write.table( tapply( m, data$num_cds, sum ), file="OUT_GENES-autosomes", quote=T, row.names=T, append=T )

	data <-subset(data, num_cds < 11)
	exons = cut(data[["num_cds"]],c(0,1,2,3,4,5,6,8,11,16,114))     #Revisar factores
    #summary(exons)
    cat("#tapply sum m~exons", file="OUT_GENES-autosomes", append=T)
	write.table( tapply( m, exons, sum ), file="OUT_GENES-autosomes", quote=T, row.names=T, append=T )

    png( "GENES-autosomes-FEAT4NumberOfExons_omega4f-exons.png", width = 1920, height = 1080 )
	boxplot(omega_4f~exons,outline=F,xlab="Number of Exons/Gene",ylab="Ka/Ks")
	abline(h=median(omega_4f),col="black")
    dev.off()

    png( "GENES-autosomes-FEAT4NumberOfExons_omegains-numCDS.png", width = 1920, height = 1080 )
	boxplot(omega_ins~data[["num_cds"]],outline=F,xlab="Number of Exons/Gene",ylab="Ka/Kins")
    dev.off()

    png( "GENES-autosomes-FEAT4NumberOfExons_pi4f-numCDS.png", width = 1920, height = 1080 )
	boxplot(pi4f~data[["num_cds"]],outline=F)
    dev.off()

    png( "GENES-autosomes-FEAT4NumberOfExons_piins-numCDS.png", width = 1920, height = 1080 )
	boxplot(piins~data[["num_cds"]],outline=F)
    dev.off()
	
    cat( cor.test(omega_4f, data$num_cds, method="spearman")$estimate, file = "OUT_GENES-autosomes",fill=T, labels="Omega_4f Spearman Rho", append=T)           # -0.2219826
    cat( cor.test(omega_4f, data$num_cds, method="spearman")$p.value, file = "OUT_GENES-autosomes", fill=T, labels="Omega_4f Spearman p-value", append=T)       # -0.2219826
    cat( cor.test(omega_ins, data$num_cds, method="spearman")$estimate, file = "OUT_GENES-autosomes" , fill=T, labels="Omega_ins Spearman Rho", append=T)       #
    cat( cor.test(omega_ins, data$num_cds, method="spearman")$p.value, file = "OUT_GENES-autosomes" , fill=T, labels="Omega_ins Spearman p-value", append=T)    #
    #cor.test(omega_4f,data$num_cds,method="spearman") # -0.2219826 
    #cor.test(omega_ins,data$num_cds,method="spearman") # 

    cat( cor.test(pi4f, data$num_cds, method="spearman")$estimate, file = "OUT_GENES-autosomes",fill=T, labels="Omega_4f Spearman Rho", append=T)          # +0.037
    cat( cor.test(pi4f, data$num_cds, method="spearman")$p.value, file = "OUT_GENES-autosomes", fill=T, labels="Omega_4f Spearman p-value", append=T)      # +0.037
    cat( cor.test(piins, data$num_cds, method="spearman")$estimate, file = "OUT_GENES-autosomes" , fill=T, labels="Omega_ins Spearman Rho", append=T)      # +0.093
    cat( cor.test(piins, data$num_cds, method="spearman")$p.value, file = "OUT_GENES-autosomes" , fill=T, labels="Omega_ins Spearman p-value", append=T)   # +0.093
    #cor.test(pi4f,data$num_cds,method="spearman") # +0.037
    #cor.test(piins,data$num_cds,method="spearman") # +0.093

    cat( kruskal.test( omega_4f ~ exons )$statistic, file = "OUT_GENES-autosomes", fill=T, labels="Omega_4f K-W chi-squared", append=T  )
    cat( kruskal.test( omega_4f ~ exons )$p.value, file = "OUT_GENES-autosomes", fill=T, labels="Omega_4f K-W p-value", append=T  )
    cat( kruskal.test( omega_ins ~ data[["num_cds"]] )$statistic, file = "OUT_GENES-autosomes", fill=T, labels="Omega_ins K-W chi-squared", append=T  )
    cat( kruskal.test( omega_ins ~ data[["num_cds"]] )$p.value, file = "OUT_GENES-autosomes", fill=T, labels="Omega_ins K-W p-value", append=T  )
    #kruskal.test(omega_4f~exons)
    #kruskal.test(omega_ins~data[["num_cds"]])


#	### Mean Distance between Exons ###
#	plot(log(data$context_distance),log(omega_4f))
#		
#		cor.test(data$context_distance,omega_4f,method="spearman") # -0.24
#		cor.test(data$context_distance,omega_ins,method="spearman") # -0.26
#		
#		cor.test(data$context_distance,data[["num_cds"]],method="spearman") # -0.15
#		cor.test(data$context_distance,m,method="spearman") # -0.05
#
#		# si hay mas distancia entre exones (q es diferente a mas exones) omega se reduce drasticamente. Esto puede ser debido a que fijan menos
#		# mutaciones ligeramente deletereas cuando se dan arrastres en exones colindantes
#
#	xyplot(pi4f~data[["context_distance"]])
#		cor.test(pi4f,data$context_distance,method="spearman") # +0.023
#		cor.test(piins,data$context_distance,method="spearman") # +0.075
#	
#	data$context_distance = cut(data$context_distance,quantile(data$context_distance,(0:5)/5))
#	summary(data[["context_distance"]])
#	with(data,tapply(m,data[["context_distance"]],sum))
#	
#	
#	quantile(data$context_distance,(0:5)/5)
#	median(data$context_distance)
#	mean(data$context_distance)
#	sd(data$context_distance)
#	
#	quantile(data[["stdev_distance"]],(0:5)/5)
#	median(data$stdev_distance)
#	mean(data$stdev_distance)
#	sd(data$stdev_distance)
#	
#	boxplot(omega_4f~data[["context_distance"]],outline=F,xlab="Mean Distance Between Exons",ylab="Ka/Ks")
#	abline(h=median(omega_4f),col="black")
#	
#	boxplot(omega_ins~data[["context_distance"]],outline=F,xlab="Mean Distance Between Exons",ylab="Ka/Kins")
#	kruskal.test(omega_4f~data[["context_distance"]])
#	kruskal.test(omega_ins~data[["context_distance"]])
#	
#	boxplot(m~data[["context_distance"]],outline=F,xlab="Mean Distance Between Exons",ylab="m")
#
#	
#	### Protein/Exon Length ###
#		
#	cor.test(m,omega_4f,method="spearman") #  gene |  exon
#	cor.test(m,omega_ins,method="spearman") #  gene |  exon
#
#	quantile(m,(0:5)/5)
#	median(m)
#	mean(m)
#	sd(m)
#	
#	m5 = cut(m,quantile(m,(0:5)/5))
#	m3 = cut(m,c(0,1533,2372,52503))
#	
#	with(data,tapply(m,m5,sum))
#	summary(m5)
#	
#	boxplot(omega_4f~m3,outline=F,xlab="Protein Length",ylab="Ka/Ks") 
#	abline(h=median(omega_4f),col="black")
#	
#	boxplot(omega_ins~m,outline=F,xlab="Protein Length",ylab="Ka/Kins")
#	kruskal.test(omega_4f~m)
#	kruskal.test(omega_ins~m)	
#				
#	### Expression ###
#	
#	data[["bias_dev_rpkm"]] = cut(data[["bias_dev_rpkm"]],quantile(data[["bias_dev_rpkm"]],(0:5)/5))
#	data[["bias_tis_rpkm"]] = cut(data[["bias_tis_rpkm"]],quantile(data[["bias_tis_rpkm"]],(0:5)/5))
#	data[["bias_str_rpkm"]] = cut(data[["bias_str_rpkm"]],quantile(data[["bias_str_rpkm"]],(0:5)/5))
#	data[["max_dev_rpkm"]] = cut(data[["max_dev_rpkm"]],quantile(data[["max_dev_rpkm"]],(0:5)/5))
#	data[["max_tissue_rpkm"]] = cut(data[["max_tissue_rpkm"]],quantile(data[["max_tissue_rpkm"]],(0:5)/5))
#	data[["max_str_rpkm"]] = cut(data[["max_str_rpkm"]],quantile(data[["max_str_rpkm"]],(0:5)/5))
#
#	with(data,tapply(m,data[["bias_dev_rpkm"]],sum))
#	summary(data[["bias_dev_rpkm"]])
#	with(data,tapply(m,data[["bias_tis_rpkm"]],sum))
#	summary(data[["bias_tis_rpkm"]])	
#	with(data,tapply(m,data[["bias_str_rpkm"]],sum))
#	summary(data[["bias_str_rpkm"]])
#	with(data,tapply(m,data[["max_dev_rpkm"]],sum))
#	summary(data[["max_dev_rpkm"]])
#	with(data,tapply(m,data[["max_tissue_rpkm"]],sum))
#	summary(data[["max_tissue_rpkm"]])	
#	with(data,tapply(m,data[["max_str_rpkm"]],sum))
#	summary(data[["max_str_rpkm"]])
#	
#	
#	quantile(data[["bias_dev_rpkm"]],(0:5)/5)
#	median(data[["bias_dev_rpkm"]])
#	mean(data[["bias_dev_rpkm"]])
#	sd(data[["bias_dev_rpkm"]])	
#
#	quantile(data[["bias_tis_rpkm"]],(0:5)/5)
#	median(data[["bias_tis_rpkm"]])
#	mean(data[["bias_tis_rpkm"]])
#	sd(data[["bias_tis_rpkm"]])		
#
#	quantile(data[["bias_str_rpkm"]],(0:5)/5)
#	median(data[["bias_str_rpkm"]])
#	mean(data[["bias_str_rpkm"]])
#	sd(data[["bias_str_rpkm"]])		
#
#	quantile(data[["max_dev_rpkm"]],(0:5)/5)
#	median(data[["max_dev_rpkm"]])
#	mean(data[["max_dev_rpkm"]])
#	sd(data[["max_dev_rpkm"]])	
#
#	quantile(data[["max_tissue_rpkm"]],(0:5)/5)
#	median(data[["max_tissue_rpkm"]])
#	mean(data[["max_tissue_rpkm"]])
#	sd(data[["max_tissue_rpkm"]])		
#
#	quantile(data[["max_str_rpkm"]],(0:5)/5)
#	median(data[["max_str_rpkm"]])
#	mean(data[["max_str_rpkm"]])
#	sd(data[["max_str_rpkm"]])		
#
#	boxplot(omega_4f~data[["bias_dev_rpkm"]],outline=F,xlab="Developmental Expression Bias",ylab="Ka/Ks")
#	abline(h=median(omega_4f),col="black")
#	boxplot(omega_4f~data[["bias_tis_rpkm"]],outline=F,xlab="Tissue Expression Bias",ylab="Ka/Ks")
#	abline(h=median(omega_4f),col="black")
#	boxplot(omega_4f~data[["bias_str_rpkm"]],outline=F,xlab="Stress Expression Bias",ylab="Ka/Ks")
#	abline(h=median(omega_4f),col="black")
#	
#	boxplot(omega_4f~data[["max_dev_rpkm"]],outline=F,xlab="Log(Maximum Developmental Expression)",ylab="Ka/Ks")
#	abline(h=median(omega_4f),col="black")
#	boxplot(omega_4f~data[["max_tissue_rpkm"]],outline=F,xlab="Log(Maximum Tissue Expression)",ylab="Ka/Ks")
#	abline(h=median(omega_4f),col="black")
#	boxplot(omega_4f~data[["max_str_rpkm"]],outline=F,xlab="Log(Maximum Stress Expression)",ylab="Ka/Ks")
#	abline(h=median(omega_4f),col="black")
#	
#	kruskal.test(omega_4f~data[["bias_dev_rpkm"]])
#	kruskal.test(omega_4f~data[["bias_tis_rpkm"]])
#	kruskal.test(omega_4f~data[["bias_str_rpkm"]])
#	kruskal.test(omega_4f~data[["max_dev_rpkm"]])
#	kruskal.test(omega_4f~data[["max_tissue_rpkm"]])
#	kruskal.test(omega_4f~data[["max_str_rpkm"]])
#	
#	kruskal.test(omega_ins~data[["bias_dev_rpkm"]])
#	kruskal.test(omega_ins~data[["bias_tis_rpkm"]])
#	kruskal.test(omega_ins~data[["bias_str_rpkm"]])
#	kruskal.test(omega_ins~data[["max_dev_rpkm"]])
#	kruskal.test(omega_ins~data[["max_tissue_rpkm"]])
#	kruskal.test(omega_ins~data[["max_str_rpkm"]])
#	
#	cor.test(data[["bias_dev"]],data[["max_dev"]],method="spearman") #
#	cor.test(data[["bias_tis"]],data[["max_tissue"]],method="spearman") #
#	cor.test(data[["bias_str"]],data[["max_str"]],method="spearman") #
#	
#	cor.test(data[["bias_dev_rpkm"]],data[["max_dev_rpkm"]],method="spearman") #
#	cor.test(data[["bias_tis_rpkm"]],data[["max_tissue_rpkm"]],method="spearman") #   
#	cor.test(data[["bias_str_rpkm"]],data[["max_str_rpkm"]],method="spearman") #
#
#	cor.test(data[["bias_dev_rpkm"]],data[["bias_tis_rpkm"]],method="spearman") #
#	cor.test(data[["bias_dev_rpkm"]],data[["bias_str_rpkm"]],method="spearman") #   
#	cor.test(data[["bias_tis_rpkm"]],data[["bias_str_rpkm"]],method="spearman") #	
#
#	cor.test(data[["max_dev_rpkm"]],data[["max_tissue_rpkm"]],method="spearman") #
#	cor.test(data[["max_dev_rpkm"]],data[["max_str_rpkm"]],method="spearman") #   
#	cor.test(data[["max_tissue_rpkm"]],data[["max_str_rpkm"]],method="spearman") #	
#	
#	cor.test(omega_4f,data$bias_dev,method="spearman") # 
#	cor.test(omega_4f,data$max_dev,method="spearman") #
#	cor.test(omega_4f,data$bias_dev_rpkm,method="spearman") #
#	cor.test(omega_4f,data[["max_dev_rpkm"]],method="spearman") # 
#	
#	cor.test(omega_4f,data$bias_tis,method="spearman") # 
#	cor.test(omega_4f,data$max_tissue,method="spearman") # 
#	cor.test(omega_4f,data$bias_tis_rpkm,method="spearman") #
#	cor.test(omega_4f,data[["max_tissue_rpkm"]],method="spearman") #	
#
#	cor.test(omega_4f,data$bias_str,method="spearman") # 
#	cor.test(omega_4f,data[["max_str"]],method="spearman") # 
#	cor.test(omega_4f,data$bias_str_rpkm,method="spearman") # 
#	cor.test(omega_4f,data[["max_str_rpkm"]],method="spearman") # 	
#
#	cor.test(omega_ins,data$bias_dev_rpkm,method="spearman") #
#	cor.test(omega_ins,data[["max_dev_rpkm"]],method="spearman") # 
#	cor.test(omega_ins,data$bias_tis_rpkm,method="spearman") #
#	cor.test(omega_ins,data[["max_tissue_rpkm"]],method="spearman") #
#	cor.test(omega_ins,data$bias_str_rpkm,method="spearman") # 
#	cor.test(omega_ins,data[["max_str_rpkm"]],method="spearman") # 	
#	
#		### Expression Stage ###
#
#		data$stage = factor(data$stage, levels=c("em0-2hr","em2-4hr","em4-6hr","em6-8hr","em8-10hr","em10-12hr","em12-14hr","em14-16hr","em16-18hr","em18-20hr","em20-22hr","em22-24hr","L1","L2","L3_12hr","L3_PS1-2","L3_PS3-6","L3_PS7-9","WPP","WPP_12hr","WPP_24hr","WPP_2days","WPP_3days","WPP_4days","AdF_Ecl_1day","AdM_Ecl_1day","AdF_Ecl_5days","AdM_Ecl_5days","AdF_Ecl_30days","AdM_Ecl_30days"))
#		boxplot(pi4f~data[["stage"]],outline=F)
#		boxplot(piins~data[["stage"]],outline=F)
#
#		data<-subset(data, bias_dev > 0.8)
#		boxplot(omega_4f~data[["stage"]],outline=F)
#		abline(h=median(omega_4f),col="black")
#		boxplot(data$bias_dev~data[["stage"]],outline=F)
#		
#		boxplot(omega_ins~data[["stage"]],outline=F)
#
#		boxplot(data$Dbias_norm~data[["stage"]],outline=F)
#		embryo <- subset(data, stage == "em0-2hr" | stage == "em2-4hr" | stage =="em4-6hr"| stage =="em6-8hr"| stage =="em8-10hr"| stage =="em10-12hr"| stage =="em12-14hr"| stage =="em14-16hr"| stage =="em16-18hr"| stage =="em18-20hr"| stage =="em20-22hr"| stage =="em22-24hr" )
#		larvae <- subset(data, stage == "L1"| stage =="L2"| stage =="L3_12hr"| stage =="L3_PS1-2"| stage =="L3_PS3-6"| stage =="L3_PS7-9")
#		pupae <- subset(data, stage == "WPP"| stage =="WPP_12hr"| stage =="WPP_24hr"| stage =="WPP_2days"| stage =="WPP_3days"| stage =="WPP_4days")
#		adult_male <- subset(data, stage =="AdM_Ecl_1day" | stage =="AdM_Ecl_5days"| stage =="AdM_Ecl_30days")
#		adult_female <- subset(data, stage =="AdF_Ecl_1day" | stage =="AdF_Ecl_5days"| stage =="AdF_Ecl_30days")
#		boxplot((embryo[["div_0f"]]/embryo[["mdmel_0f"]])/(embryo[["div_4f"]]/embryo[["mdmel_4f"]]),(larvae[["div_0f"]]/larvae[["mdmel_0f"]])/(larvae[["div_4f"]]/larvae[["mdmel_4f"]]),(pupae[["div_0f"]]/pupae[["mdmel_0f"]])/(pupae[["div_4f"]]/pupae[["mdmel_4f"]]),(adult_male[["div_0f"]]/adult_male[["mdmel_0f"]])/(adult_male[["div_4f"]]/adult_male[["mdmel_4f"]]),(adult_female[["div_0f"]]/adult_female[["mdmel_0f"]])/(adult_female[["div_4f"]]/adult_female[["mdmel_4f"]]),outline=F)
#
#
#	### Recombination ###
#	
#	data[["mcvean_10kb"]] = cut(data[["mcvean_10kb"]],quantile(data[["mcvean_10kb"]],(0:5)/5))
#	data[["mcvean_100kb"]] = cut(data[["mcvean_100kb"]],quantile(data[["mcvean_100kb"]],(0:5)/5))
#	data[["mcvean_1Mb"]] = cut(data[["mcvean_1Mb"]],quantile(data[["mcvean_1Mb"]],(0:5)/5))
#	data[["comeron_100kb"]] = cut(data[["comeron_100kb"]],quantile(data[["comeron_100kb"]],(0:5)/5))
#	
#	with(data,tapply(m,data[["mcvean_10kb"]],sum))
#	summary(data[["mcvean_10kb"]])
#	with(data,tapply(m,data[["mcvean_100kb"]],sum))
#	summary(data[["mcvean_100kb"]])
#	with(data,tapply(m,data[["comeron_100kb"]],sum))
#	summary(data[["comeron_100kb"]])	
#	
#	boxplot(omega_4f~data[["mcvean_10kb"]],outline=F,xlab="rho 10kb",ylab="Ka/Ks")
#	abline(h=median(omega_4f),col="black")
#	boxplot(omega_4f~data[["mcvean_100kb"]],outline=F,xlab="rho 100kb",ylab="Ka/Ks")
#	abline(h=median(omega_4f),col="black")
#	boxplot(omega_4f~data[["mcvean_1Mb"]],outline=F,xlab="rho 1Mb",ylab="Ka/Ks")
#	abline(h=median(omega_4f),col="black")	
#	boxplot(omega_4f~data[["comeron_100kb"]],outline=F,xlab="cM/Mb 100kb",ylab="Ka/Ks")
#	abline(h=median(omega_4f),col="black")	
#	
#	quantile(data[["mcvean_10kb"]],(0:5)/5)
#	median(data[["mcvean_10kb"]])
#	mean(data[["mcvean_10kb"]])
#	sd(data[["mcvean_10kb"]])
#	quantile(data[["mcvean_100kb"]],(0:5)/5)
#	median(data[["mcvean_100kb"]])
#	mean(data[["mcvean_100kb"]])
#	sd(data[["mcvean_100kb"]])
#	quantile(data[["mcvean_1Mb"]],(0:5)/5)
#	median(data[["mcvean_1Mb"]])
#	mean(data[["mcvean_1Mb"]])
#	sd(data[["mcvean_1Mb"]])
#	quantile(data[["comeron_100kb"]],(0:5)/5)
#	median(data[["comeron_100kb"]])
#	mean(data[["comeron_100kb"]])
#	sd(data[["comeron_100kb"]])	
#	
#		cor.test(data[["mcvean_10kb"]],omega_4f,method="spearman")
#		cor.test(data[["mcvean_10kb"]],omega_ins,method="spearman")
#		cor.test(data[["mcvean_100kb"]],omega_4f,method="spearman")
#		cor.test(data[["mcvean_100kb"]],omega_ins,method="spearman")
#		cor.test(data[["mcvean_1Mb"]],omega_4f,method="spearman")
#		cor.test(data[["mcvean_1Mb"]],omega_ins,method="spearman")
#
#		cor.test(data[["comeron_100kb"]],omega_4f,method="spearman")
#		cor.test(data[["comeron_100kb"]],omega_ins,method="spearman")
#
#	kruskal.test(omega_4f~data[["mcvean_10kb"]])
#	kruskal.test(omega_4f~data[["mcvean_100kb"]])
#	kruskal.test(omega_4f~data[["mcvean_1Mb"]])
#	kruskal.test(omega_4f~data[["comeron_100kb"]])
#
#	kruskal.test(omega_ins~data[["mcvean_10kb"]])
#	kruskal.test(omega_ins~data[["mcvean_100kb"]])
#	kruskal.test(omega_ins~data[["mcvean_1Mb"]])
#	kruskal.test(omega_ins~data[["comeron_100kb"]])	
#	
#
#				### Expression Vs Recombination ###
#				cor.test(data[["bias_dev_rpkm"]],pi4f,method="spearman")
#				cor.test(data[["bias_dev_rpkm"]],t4f,method="spearman")
#				cor.test(data[["max_dev_rpkm"]],pi4f,method="spearman")
#				cor.test(data[["max_dev_rpkm"]],t4f,method="spearman")
#				
#				cor.test(data[["mcvean_10kb"]],data[["bias_dev_rpkm"]],method="spearman")
#				cor.test(data[["mcvean_100kb"]],data[["bias_dev_rpkm"]],method="spearman")
#				cor.test(data[["comeron_100kb"]],data[["bias_dev_rpkm"]],method="spearman")
#				cor.test(data[["mcvean_1Mb"]],data[["bias_dev_rpkm"]],method="spearman")
#				
#				cor.test(data[["mcvean_10kb"]],data[["bias_tis_rpkm"]],method="spearman")
#				cor.test(data[["mcvean_100kb"]],data[["bias_tis_rpkm"]],method="spearman")
#				cor.test(data[["comeron_100kb"]],data[["bias_tis_rpkm"]],method="spearman")
#				cor.test(data[["mcvean_1Mb"]],data[["bias_tis_rpkm"]],method="spearman")
#
#				cor.test(data[["mcvean_10kb"]],data[["bias_str_rpkm"]],method="spearman")
#				cor.test(data[["mcvean_100kb"]],data[["bias_str_rpkm"]],method="spearman")
#				cor.test(data[["comeron_100kb"]],data[["bias_str_rpkm"]],method="spearman")
#				cor.test(data[["mcvean_1Mb"]],data[["bias_str_rpkm"]],method="spearman")				
#				
#				cor.test(data[["mcvean_10kb"]],data[["max_dev_rpkm"]],method="spearman")
#				cor.test(data[["mcvean_100kb"]],data[["max_dev_rpkm"]],method="spearman")
#				cor.test(data[["comeron_100kb"]],data[["max_dev_rpkm"]],method="spearman")
#				cor.test(data[["mcvean_1Mb"]],data[["max_dev_rpkm"]],method="spearman")
#
#				cor.test(data[["mcvean_10kb"]],data[["max_tissue_rpkm"]],method="spearman")
#				cor.test(data[["mcvean_100kb"]],data[["max_tissue_rpkm"]],method="spearman")
#				cor.test(data[["comeron_100kb"]],data[["max_tissue_rpkm"]],method="spearman")
#				cor.test(data[["mcvean_1Mb"]],data[["max_tissue_rpkm"]],method="spearman")
#				
#				cor.test(data[["mcvean_10kb"]],data[["max_str_rpkm"]],method="spearman")
#				cor.test(data[["mcvean_100kb"]],data[["max_str_rpkm"]],method="spearman")
#				cor.test(data[["comeron_100kb"]],data[["max_str_rpkm"]],method="spearman")
#				cor.test(data[["mcvean_1Mb"]],data[["max_str_rpkm"]],method="spearman")
#				
#				boxplot(data[["bias_dev_rpkm"]]~data[["mcvean_10kb"]],outline=F,xlab="rho 10kb",ylab="Developmental Expression Bias")
#				boxplot(pi4f~data[["mcvean_10kb"]],outline=F,xlab="rho 10kb",ylab="Synonymous Pi")
#				boxplot(pi4f~data[["bias_dev_rpkm"]],outline=F,xlab="Developmental Expression Bias",ylab="Synonymous Pi")
#			    
#			    gene <- data.frame(
#                                DevBias = (data$bias_dev_rpkm),        
#                                PiSyn=(pi4f),
#                                rho10kb=(data[["mcvean_10kb"]])
#                                
#                          )
#			    gene <- na.omit(gene)
#			    spcor(gene,method=c("spearman"))
#				# a 10 kb hay correlacion positiva entre Dbias y recombincion, y negativa entre Smax y recombinacion. Los genes housekeeping estan en regiones de baja recomb
#				# o donde estan los genes housekeeping hay subestimas de recombinacion? demasiado bonito para ser verdad, seguramete sera la segunda opcion.... :S
#
#	### Chromosomes ###
#
#	#data <- subset(data,chro_fraction > 0.5)
#	#nrow(data)
#	boxplot(omega_4f~data$chromosome,outline=F,xlab="Chromosomes",ylab="Ka/Ks")
#	abline(h=median(omega_4f),col="black")
#	kruskal.test(omega_ins~data[["chromosome"]])
#
#	### Chromatin ###
#
#	#data <- subset(data,chro_fraction > 0.5)
#	#nrow(data)
#	boxplot(omega_4f~data$chromatin,outline=F,xlab="Chromatin States",ylab="Ka/Ks")
#	abline(h=median(omega_4f),col="black")	
#	kruskal.test(omega_ins~data[["chromatin"]])
#	
#	mA <- ((data[["a_s2_0f"]]+data[["a_bg_0f"]]+ data[["a_s2_4f"]]+data[["a_bg_4f"]]+ data[["a_s2_2f"]]+data[["a_bg_2f"]])/2)
#	mB <- ((data[["b_s2_0f"]]+data[["b_bg_0f"]]+ data[["b_s2_4f"]]+data[["b_bg_4f"]]+ data[["b_s2_2f"]]+data[["b_bg_2f"]])/2)
#	mC <- ((data[["c_s2_0f"]]+data[["c_bg_0f"]]+ data[["c_s2_4f"]]+data[["c_bg_4f"]]+ data[["c_s2_2f"]]+data[["c_bg_2f"]])/2)
#	mD <- ((data[["d_s2_0f"]]+data[["d_bg_0f"]]+ data[["d_s2_4f"]]+data[["d_bg_4f"]]+ data[["d_s2_2f"]]+data[["d_bg_2f"]])/2)
#	mE <- ((data[["e_s2_0f"]]+data[["e_bg_0f"]]+ data[["e_s2_4f"]]+data[["e_bg_4f"]]+ data[["e_s2_2f"]]+data[["e_bg_2f"]])/2)
#	mF <- ((data[["f_s2_0f"]]+data[["f_bg_0f"]]+ data[["f_s2_4f"]]+data[["f_bg_4f"]]+ data[["f_s2_2f"]]+data[["f_bg_2f"]])/2)
#	mG <- ((data[["g_s2_0f"]]+data[["g_bg_0f"]]+ data[["g_s2_4f"]]+data[["g_bg_4f"]]+ data[["g_s2_2f"]]+data[["g_bg_2f"]])/2)
#	mH <- ((data[["h_s2_0f"]]+data[["h_bg_0f"]]+ data[["h_s2_4f"]]+data[["h_bg_4f"]]+ data[["h_s2_2f"]]+data[["h_bg_2f"]])/2)
#	mI <- ((data[["i_s2_0f"]]+data[["i_bg_0f"]]+ data[["i_s2_4f"]]+data[["i_bg_4f"]]+ data[["i_s2_2f"]]+data[["i_bg_2f"]])/2)
#	sum(mA)
#	sum(mB)
#	sum(mC)
#	sum(mD)
#	sum(mE)
#	sum(mF)
#	sum(mG)
#	sum(mH)
#	sum(mI)
#
#	mAs2 <- data[["a_s2_0f"]]+ data[["a_s2_4f"]]+ data[["a_s2_2f"]]
#	mBs2 <-data[["b_s2_0f"]]+ data[["b_s2_4f"]]+ data[["b_s2_2f"]]
#	mCs2 <- data[["c_s2_0f"]]+ data[["c_s2_4f"]]+ data[["c_s2_2f"]]
#	mDs2 <- data[["d_s2_0f"]]+ data[["d_s2_4f"]]+ data[["d_s2_2f"]]
#	mEs2 <- data[["e_s2_0f"]]+ data[["e_s2_4f"]]+ data[["e_s2_2f"]]
#	mFs2 <- data[["f_s2_0f"]]+ data[["f_s2_4f"]]+ data[["f_s2_2f"]]
#	mGs2 <- data[["g_s2_0f"]]+ data[["g_s2_4f"]]+ data[["g_s2_2f"]]
#	mHs2 <- data[["h_s2_0f"]]+ data[["h_s2_4f"]]+ data[["h_s2_2f"]]
#	mIs2 <- data[["i_s2_0f"]]+ data[["i_s2_4f"]]+ data[["i_s2_2f"]]	
#	sum(mAs2)
#	sum(mBs2)
#	sum(mCs2)
#	sum(mDs2)
#	sum(mEs2)
#	sum(mFs2)
#	sum(mGs2)
#	sum(mHs2)
#	sum(mIs2)
#	
#	mAbg <- data[["a_bg_0f"]]+ data[["a_bg_4f"]]+ data[["a_bg_2f"]]
#	mBbg <- data[["b_bg_0f"]]+ data[["b_bg_4f"]]+ data[["b_bg_2f"]]
#	mCbg <- data[["c_bg_0f"]]+ data[["c_bg_4f"]]+ data[["c_bg_2f"]]
#	mDbg <- data[["d_bg_0f"]]+ data[["d_bg_4f"]]+ data[["d_bg_2f"]]
#	mEbg <- data[["e_bg_0f"]]+ data[["e_bg_4f"]]+ data[["e_bg_2f"]]
#	mFbg <- data[["f_bg_0f"]]+ data[["f_bg_4f"]]+ data[["f_bg_2f"]]
#	mGbg <- data[["g_bg_0f"]]+ data[["g_bg_4f"]]+ data[["g_bg_2f"]]
#	mHbg <- data[["h_bg_0f"]]+ data[["h_bg_4f"]]+ data[["h_bg_2f"]]
#	mIbg <- data[["i_bg_0f"]]+ data[["i_bg_4f"]]+ data[["i_bg_2f"]]	
#	sum(mAbg)
#	sum(mBbg)
#	sum(mCbg)
#	sum(mDbg)
#	sum(mEbg)
#	sum(mFbg)
#	sum(mGbg)
#	sum(mHbg)
#	sum(mIbg)	
#	
#	####
#	
#	m <- data[["mdmel_0f"]]+data[["mdmel_4f"]]+data[["mdmel_2f"]]
#	mA <- ((data[["a_s2_0f"]]+data[["a_bg_0f"]]+ data[["a_s2_4f"]]+data[["a_bg_4f"]]+ data[["a_s2_2f"]]+data[["a_bg_2f"]])/2)/(data[["mdmel_0f"]]+data[["mdmel_4f"]]+data[["mdmel_2f"]])
#	mB <- ((data[["b_s2_0f"]]+data[["b_bg_0f"]]+ data[["b_s2_4f"]]+data[["b_bg_4f"]]+ data[["b_s2_2f"]]+data[["b_bg_2f"]])/2)/(data[["mdmel_0f"]]+data[["mdmel_4f"]]+data[["mdmel_2f"]])
#	mC <- ((data[["c_s2_0f"]]+data[["c_bg_0f"]]+ data[["c_s2_4f"]]+data[["c_bg_4f"]]+ data[["c_s2_2f"]]+data[["c_bg_2f"]])/2)/(data[["mdmel_0f"]]+data[["mdmel_4f"]]+data[["mdmel_2f"]])
#	mD <- ((data[["d_s2_0f"]]+data[["d_bg_0f"]]+ data[["d_s2_4f"]]+data[["d_bg_4f"]]+ data[["d_s2_2f"]]+data[["d_bg_2f"]])/2)/(data[["mdmel_0f"]]+data[["mdmel_4f"]]+data[["mdmel_2f"]])
#	mE <- ((data[["e_s2_0f"]]+data[["e_bg_0f"]]+ data[["e_s2_4f"]]+data[["e_bg_4f"]]+ data[["e_s2_2f"]]+data[["e_bg_2f"]])/2)/(data[["mdmel_0f"]]+data[["mdmel_4f"]]+data[["mdmel_2f"]])
#	mF <- ((data[["f_s2_0f"]]+data[["f_bg_0f"]]+ data[["f_s2_4f"]]+data[["f_bg_4f"]]+ data[["f_s2_2f"]]+data[["f_bg_2f"]])/2)/(data[["mdmel_0f"]]+data[["mdmel_4f"]]+data[["mdmel_2f"]])
#	mG <- ((data[["g_s2_0f"]]+data[["g_bg_0f"]]+ data[["g_s2_4f"]]+data[["g_bg_4f"]]+ data[["g_s2_2f"]]+data[["g_bg_2f"]])/2)/(data[["mdmel_0f"]]+data[["mdmel_4f"]]+data[["mdmel_2f"]])
#	mH <- ((data[["h_s2_0f"]]+data[["h_bg_0f"]]+ data[["h_s2_4f"]]+data[["h_bg_4f"]]+ data[["h_s2_2f"]]+data[["h_bg_2f"]])/2)/(data[["mdmel_0f"]]+data[["mdmel_4f"]]+data[["mdmel_2f"]])
#	mI <- ((data[["i_s2_0f"]]+data[["i_bg_0f"]]+ data[["i_s2_4f"]]+data[["i_bg_4f"]]+ data[["i_s2_2f"]]+data[["i_bg_2f"]])/2)/(data[["mdmel_0f"]]+data[["mdmel_4f"]]+data[["mdmel_2f"]])
#









#		### CORRELATIONS ###
#
#		cor.test(pi4f,mA,method="spearman")
#		cor.test(pi4f,mB,method="spearman")
#		cor.test(pi4f,mC,method="spearman")
#		cor.test(pi4f,mD,method="spearman")
#		cor.test(pi4f,mE,method="spearman")
#		cor.test(pi4f,mF,method="spearman")
#		cor.test(pi4f,mG,method="spearman")
#		cor.test(pi4f,mH,method="spearman")
#		cor.test(pi4f,mI,method="spearman")
#
#		cor.test(piins,mA,method="spearman") #se hace mas peqno q el 4f
#		cor.test(piins,mB,method="spearman") #se hace mas grande q el 4f
#		cor.test(piins,mC,method="spearman")
#		cor.test(piins,mD,method="spearman")
#		cor.test(piins,mE,method="spearman")
#		cor.test(piins,mF,method="spearman")
#		cor.test(piins,mG,method="spearman")
#		cor.test(piins,mH,method="spearman")
#		cor.test(piins,mI,method="spearman")
#
#		cor.test(k4f,mA,method="spearman") # cor positiva -> mayor mu? || mayor linked selection -> peor uso de codon
#		cor.test(k4f,mB,method="spearman") # cor negativa -> menor mu? || menor linked selection -> mejor uso de codon
#		cor.test(k4f,mC,method="spearman") # cor negativa -> menor mu? || menor linked selection -> mejor uso de codon
#		cor.test(k4f,mD,method="spearman") # cor negativa -> menor mu? || menor linked selection -> mejor uso de codon
#		cor.test(k4f,mE,method="spearman") # cor positiva -> mayor mu? || mayor linked selection -> peor uso de codon
#		cor.test(k4f,mF,method="spearman") # cor negativa -> menor mu? || menor linked selection -> mejor uso de codon
#		cor.test(k4f,mG,method="spearman") # cor positiva -> mayor mu? || mayor linked selection -> peor uso de codon
#		cor.test(k4f,mH,method="spearman") # cor ns
#		cor.test(k4f,mI,method="spearman") # cor ns
#
#		cor.test(kins,mA,method="spearman") # cor negativa -> menor mu + mayor linked selection -> peor uso de codon
#		cor.test(kins,mB,method="spearman") # cor ns
#		cor.test(kins,mC,method="spearman") # cor ns
#		cor.test(kins,mD,method="spearman") # cor ns
#		cor.test(kins,mE,method="spearman") # cor positiva -> mayor mu + mayor linked selection -> peor uso de codon
#		cor.test(kins,mF,method="spearman") # cor positiva -> mayor mu + menor linked selection -> mejor uso de codon
#		cor.test(kins,mG,method="spearman") # cor positiva -> mayor mu + mayor linked selection -> peor uso de codon
#		cor.test(kins,mH,method="spearman") # cor positiva -> mayor mu
#		cor.test(kins,mI,method="spearman") # cor positiva -> mayor mu
#
#		cor.test(omega_4f,mA,method="spearman")
#		cor.test(omega_4f,mB,method="spearman")
#		cor.test(omega_4f,mC,method="spearman")
#		cor.test(omega_4f,mD,method="spearman")
#		cor.test(omega_4f,mE,method="spearman")
#		cor.test(omega_4f,mF,method="spearman")
#		cor.test(omega_4f,mG,method="spearman")
#		cor.test(omega_4f,mH,method="spearman")
#		cor.test(omega_4f,mI,method="spearman")
#
#		cor.test(omega_ins,mA,method="spearman")
#		cor.test(omega_ins,mB,method="spearman")
#		cor.test(omega_ins,mC,method="spearman")
#		cor.test(omega_ins,mD,method="spearman")
#		cor.test(omega_ins,mE,method="spearman")
#		cor.test(omega_ins,mF,method="spearman")
#		cor.test(omega_ins,mG,method="spearman")
#		cor.test(omega_ins,mH,method="spearman")
#		cor.test(omega_ins,mI,method="spearman")
#
#		cor.test(data$Dbias,mA,method="spearman")
#		cor.test(data$Dbias,mB,method="spearman")
#		cor.test(data$Dbias,mC,method="spearman")
#		cor.test(data$Dbias,mD,method="spearman")
#		cor.test(data$Dbias,mE,method="spearman")
#		cor.test(data$Dbias,mF,method="spearman")
#		cor.test(data$Dbias,mG,method="spearman")
#		cor.test(data$Dbias,mH,method="spearman")
#		cor.test(data$Dbias,mI,method="spearman")
#
#		cor.test(data$Smax,mA,method="spearman")
#		cor.test(data$Smax,mB,method="spearman")
#		cor.test(data$Smax,mC,method="spearman")
#		cor.test(data$Smax,mD,method="spearman")
#		cor.test(data$Smax,mE,method="spearman")
#		cor.test(data$Smax,mF,method="spearman")
#		cor.test(data$Smax,mG,method="spearman")
#		cor.test(data$Smax,mH,method="spearman")
#		cor.test(data$Smax,mI,method="spearman")
#
#
#
#### PLOTS ###
#
#xyplot(omega_ins~comeron_100kb + mcvean_10kb + mcvean_100kb + mcvean_1Mb + stage_rpkm + bias_dev_rpkm + max_dev_rpkm + tissue_rpkm + bias_tis_rpkm + max_tissue_rpkm + stress_rpkm + bias_str_rpkm + max_str_rpkm + context_distance + num_cds + FBtrs_per_gene + m + mcomp + mA + mB + mC + mD + mE + mF + mG + mH + mI + chromosome,data)
#boxplot(omega_4f~data$chromatin,outline=F)
#
#par(mfrow=c(3,1))
#hist(log10(omega_ins))  # parece mejor
#hist(sqrt(omega_ins))
#hist(omega_4f) 
#par(mfrow=c(1,1))
#
#par(mfrow=c(1,2))
#hist(log10(data$comeron_100kb+.1))
#hist(sqrt(data$comeron_100kb)) # parece mejor
#par(mfrow=c(1,1))
#
#par(mfrow=c(1,2))
#hist(log10(data$mcvean_10kb+.1)) # parece mejor
#hist(sqrt(data$mcvean_10kb))
#par(mfrow=c(1,1))
#
#par(mfrow=c(1,2))
#hist(log10(data$mcvean_100kb+.1)) # parece mejor
#hist(sqrt(data$mcvean_100kb))
#par(mfrow=c(1,1))
#
#par(mfrow=c(1,2))
#hist(log10(data$mcvean_1Mb+.1)) 
#hist(sqrt(data$mcvean_1Mb)) # parece mejor
#par(mfrow=c(1,1))
#
#par(mfrow=c(3,1))
#hist(log10(data$bias_dev_rpkm))
#hist(sqrt(data$bias_dev_rpkm)) # parece mejor
#hist(data$bias_dev_rpkm) 
#par(mfrow=c(1,1))
#
#hist(data$max_dev_rpkm) 
#hist(data$max_tissue_rpkm) 
#hist(data$max_str_rpkm) 
#
#par(mfrow=c(3,1))
#hist(log10(data$bias_tis_rpkm))
#hist(sqrt(data$bias_tis_rpkm)) # parece mejor
#hist(data$bias_tis_rpkm) 
#par(mfrow=c(1,1))
#summary(data$bias_tis_rpkm)
#summary(data$bias_dev_rpkm)
#summary(data$bias_str_rpkm)
#
#
#par(mfrow=c(3,1))
#hist(log10(data$bias_str_rpkm))
#hist(sqrt(data$bias_str_rpkm)) # parece mejor
#hist(data$bias_str_rpkm) 
#par(mfrow=c(1,1))
#
#par(mfrow=c(3,1))
#hist(log10(data$context_distance),breaks=c(100)) #este
#hist(sqrt(data$context_distance)) 
#hist(data$context_distance) 
#par(mfrow=c(1,1))
#
#
#par(mfrow=c(3,1))
#hist(log10(data[["num_cds"]]),breaks=c(20)) 
#hist(sqrt((data$num_cds)),breaks=c(20)) #este 
#hist(data$num_cds,breaks=c(100)) 
#par(mfrow=c(1,1))
#
#par(mfrow=c(3,1))
#hist(log10(data[["FBtrs_per_gene"]]),breaks=c(20)) 
#hist(sqrt((data$FBtrs_per_gene)),breaks=c(20)) 
#hist(data$FBtrs_per_gene,breaks=c(100)) #este
#par(mfrow=c(1,1))
#
#par(mfrow=c(3,1))
#hist(log10(m),breaks=c(20)) #este
#hist(sqrt((m)),breaks=c(20)) 
#hist(m,breaks=c(100)) 
#par(mfrow=c(1,1))
#
#par(mfrow=c(3,1))
#hist(log10(mcomp),breaks=c(20)) #este
#hist(sqrt((mcomp)),breaks=c(20)) 
#hist(mcomp,breaks=c(20)) 
#par(mfrow=c(1,1))
#
#par(mfrow=c(3,1))
#hist(log10(mA+1),breaks=c(20)) #este
#hist(sqrt((mA)),breaks=c(20)) 
#hist(mA,breaks=c(100)) 
#par(mfrow=c(1,1))
#
#
#xyplot(log(omega_ins)~log10(mcvean_10kb),data)
#xyplot(log(omega_ins)~log10(mcvean_100kb),data)
#xyplot(log(omega_ins)~log10(mcvean_1Mb),data)
#xyplot(log(omega_ins)~log10(comeron_100kb),data)
#xyplot(log10(omega_ins)~bias_dev_rpkm,data)
#xyplot(log10(omega_ins)~max_dev_rpkm,data)
#xyplot(log10(omega_ins)~log10(m),data)
#xyplot(log10(omega_ins)~log10(context_distance),data)
#xyplot((omega_ins)~(num_cds),data)
#xyplot((omega_ins)~(FBtrs_per_gene),data)
#xyplot((omega_ins)~(mcomp),data)
#
#xyplot(log10(omega_4f)~log10(mA+1),data)
#cor.test(omega_4f,log10(mA+1),method="spearman") #-0.2100069 
#
#summary(data$max_dev_rpkm)
#xyplot(omega_ins~comeron_100kb,data)#+ stage_rpkm + bias_dev_rpkm + max_dev_rpkm + tissue_rpkm + bias_tis_rpkm + max_tissue_rpkm + stress_rpkm + bias_str_rpkm + max_str_rpkm + context_distance + num_cds + FBtrs_per_gene + m + mcomp + mA + mB + mC + mD + mE + mF + mG + mH + mI + chromosome,data)
#
#
#cor.test(omega_ins,data$bias_dev_rpkm,method="spearman") #0.3383352
#cor.test(omega_ins,data$max_dev_rpkm,method="spearman") #-0.1754348
#cor.test(data$bias_dev_rpkm,data$max_dev_rpkm,method="spearman") #-0.230016
#
#cor.test(omega_ins,data$mcvean_10kb,method="spearman") #ns
#cor.test(omega_ins,m,method="spearman") #-0.076
#cor.test(omega_ins,m,method="spearman") #-0.076
#
#cor.test(omega_ins,data$num_cds,method="spearman") #-0.2698525
#cor.test(data$mcvean_10kb,data$context_distance,method="spearman") #0.157
#cor.test(data$bias_dev_rpkm,data$context_distance,method="spearman") #0.127
#cor.test(data$max_dev_rpkm,data$context_distance,method="spearman") #0.042
#cor.test(data$num_cds,data$bias_dev_rpkm,method="spearman") #-0.15
#cor.test(data$num_cds,data$max_dev_rpkm,method="spearman") #-0.072
