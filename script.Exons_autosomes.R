#
# EXONS
#
# Dos datasets:
#   A) EXON.SHORTINTRONS.n
#   #B) GENE.SHORTINTRONS.n
#
# Filtros:
#   A) Exones
#       A.1)
#           Erep <-(data[["num_FBtrs"]]/data[["FBtrs_per_gene"]]) # Exones donde anotaciones 5.49 y 5.50 son distintos!!    #Solo para data exon.shortintrons.n
#           data <- subset(data,Erep <= 1)
#       A.2)
#           data <- na.omit(data) #para el modelo mixto
#

library(car)
library(lattice)
library(ppcor)

rm(list = ls())
mypath <- "/home/sergio/chromatin/analysis/R/chromatin_features/"
setwd(mypath)

# Load datasets
data <- read.table(file="EXON.SHORTINTRONS.n",header=TRUE,sep="\t")

#nrow(data)
data <- subset(data, chromosome != "X")                     #AUTOSOMES
Erep <- (data[["num_FBtrs"]] / data[["FBtrs_per_gene"]])  # Exones donde anotaciones 5.49 y 5.50 son distintos!!    #Solo para data exon.shortintrons.n
data <- subset(data, Erep <= 1)                              #Filtramos los que num_FBtrs/FBtrs_per_gene sean mayores que 1 (Casos que ocurren x usar rev5.49 y rev5.50)

### PopGen Summary Statistics ###

data <- na.omit(data) #para el modelo mixto

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

### EXON ###
    #If exists previous file, move it to .bak
    if (file.exists("OUT_EXONS-autosomes")) {
        file.rename("OUT_EXONS-autosomes", "OUT_EXONS-autosomes.bak")
        file.remove("OUT_EXONS-autosomes") #Remove previous file
    }
    cat("###EXONS-autosomes\n\n", file="OUT_EXONS-autosomes", append=T)

    #Feature 1
	### Transcripts (mRNAs per exon) ###
		cat("##Feature1: mRNAs per exon\n", file="OUT_EXONS-autosomes", append=T)

        #Suma
		cat("#tapply sum m~data$num_FBtrs", file="OUT_EXONS-autosomes", append=T)
        write.table(tapply(m, data$num_FBtrs, sum), file="OUT_EXONS-autosomes", row.names=T, quote=T, append=T)

		num_FBtrs_factor = cut( data[["num_FBtrs"]], c( 0, 1, 2, 75 ) )
		
        png( "EXONS-autosomes-FEAT1Transcripts_omega4f-numFBtrs.png", width = 1920, height = 1080 )
		boxplot( omega_4f~num_FBtrs_factor, outline=F, xlab="Number of Transcripts/Exon", ylab="Ka/Ks", main="Omega 4f" )
		abline( h = median( omega_4f ), col="black" )
        dev.off()

        png( "EXONS-autosomes-FEAT1Transcripts_omegains-numFBtrs.png", width=1920, height=1080 )
		boxplot( omega_ins ~ num_FBtrs_factor, outline=F, xlab="Number of Transcripts/Exon", ylab="Ka/Kins", main="Omega intrones pequeños" )
        abline( h = median( omega_ins ), col="black" )
        dev.off()

        #cor.test(omega_4f, data$num_FBtrs, method="spearman")
        cat( cor.test(omega_4f, data$num_FBtrs, method="spearman")$estimate, labels="Omega_4f Spearman Rho", file = "OUT_EXONS-autosomes", fill=T, append=T)
        cat( cor.test(omega_4f, data$num_FBtrs, method="spearman")$p.value, labels="Omega_4f Spearman p-value", file = "OUT_EXONS-autosomes", fill=T, append=T  )
		cat( cor.test(omega_ins, data$num_FBtrs, method="spearman")$estimate, labels="Omega_ins Spearman Rho", file = "OUT_EXONS-autosomes" , fill=T, append=T  )
		cat( cor.test(omega_ins, data$num_FBtrs, method="spearman")$p.value, labels="Omega_ins Spearman p-value", file = "OUT_EXONS-autosomes" , fill=T, append=T  )
		
        #kruskal.test(omega_4f~num_FBtrs_factor)
        cat( kruskal.test( omega_4f ~ num_FBtrs_factor )$statistic, labels="Omega_4f K-W chi-squared", file = "OUT_EXONS-autosomes", fill=T, append=T  )
		cat( kruskal.test( omega_4f ~ num_FBtrs_factor )$p.value, labels="Omega_4f K-W p-value", file = "OUT_EXONS-autosomes", fill=T, append=T  )
		cat( kruskal.test( omega_ins ~ num_FBtrs_factor )$statistic, labels="Omega_ins K-W chi-squared", file = "OUT_EXONS-autosomes", fill=T, append=T  )
		cat( kruskal.test( omega_ins ~ num_FBtrs_factor )$p.value, labels="Omega_ins K-W p-value", file = "OUT_EXONS-autosomes", fill=T, append=T  )

	
    #Feature 2
	### Exon Representativity ### 
		cat("\n", file="OUT_EXONS-autosomes", append=T)
		cat("##Feature2: Exon Representativity\n", file="OUT_EXONS-autosomes", append=T)

		Erep <- ( data[["num_FBtrs"]] / data[["FBtrs_per_gene"]] )
        #summary( Erep )
	
        #Suma
        alternative <- cut( Erep, c( 0, 0.25, 0.5, 0.99, 1 ) )
		cat("#tapply sum m~alternative", file="OUT_EXONS-autosomes", append=T)
        write.table( tapply( m, alternative, sum ), file="OUT_EXONS-autosomes", quote=T, row.names=T, append=T )

        png( "EXONS-autosomes-FEAT2Erep_omega4f-exonRepresentativity.png", width = 1920, height = 1080 )
		boxplot( omega_4f ~ alternative, outline=F, xlab="Exon Representativity", ylab="Ka/Ks" )
		abline( h = median( omega_4f ), col="black" )
        dev.off()

        png( "EXONS-autosomes-FEAT2Erep_Histogram.png", width = 1920, height = 1080 )
		hist( Erep, breaks = c( 100 ) )
        dev.off()

        cat( quantile( Erep, ( 0:2 ) / 2 )[1], labels="Quantile 0%", file="OUT_EXONS-autosomes", fill=T, append=T )
        cat( quantile( Erep, ( 0:2 ) / 2 )[2], labels="Quantile 50%", file="OUT_EXONS-autosomes", fill=T, append=T )
        cat( quantile( Erep, ( 0:2 ) / 2 )[3], labels="Quantile 100%", file="OUT_EXONS-autosomes", fill=T, append=T )
		
        cat(median( Erep ), file="OUT_EXONS-autosomes", fill=T, labels="Median", append=T)
		cat(mean( Erep ), file="OUT_EXONS-autosomes", fill=T, labels="Mean", append=T)
        cat(sd( Erep ), file="OUT_EXONS-autosomes", fill=T, labels="SD", append=T)

        cat( cor.test(omega_4f, Erep, method="spearman")$estimate, file = "OUT_EXONS-autosomes",fill=T, labels="Omega_4f Spearman Rho", append=T)
        cat( cor.test(omega_4f, Erep, method="spearman")$p.value, file = "OUT_EXONS-autosomes", fill=T, labels="Omega_4f Spearman p-value", append=T)
		cat( cor.test(omega_ins, Erep, method="spearman")$estimate, file = "OUT_EXONS-autosomes" , fill=T, labels="Omega_ins Spearman Rho", append=T)
		cat( cor.test(omega_ins, Erep, method="spearman")$p.value, file = "OUT_EXONS-autosomes" , fill=T, labels="Omega_ins Spearman p-value", append=T)

        cat( kruskal.test( omega_4f ~ alternative )$statistic, file = "OUT_EXONS-autosomes", fill=T, labels="Omega_4f K-W chi-squared", append=T  )
		cat( kruskal.test( omega_4f ~ alternative )$p.value, file = "OUT_EXONS-autosomes", fill=T, labels="Omega_4f K-W p-value", append=T  )
		cat( kruskal.test( omega_ins ~ alternative )$statistic, file = "OUT_EXONS-autosomes", fill=T, labels="Omega_ins K-W chi-squared", append=T  )
		cat( kruskal.test( omega_ins ~ alternative )$p.value, file = "OUT_EXONS-autosomes", fill=T, labels="Omega_ins K-W p-value", append=T  )
		
		
    #Feature 3
	### Order ###
		cat("\n", file="OUT_EXONS-autosomes", append=T)
		cat("##Feature3: Order\n", file="OUT_EXONS-autosomes", append=T)

		cat("#tapply sum m~data$order", file="OUT_EXONS-autosomes", append=T)
        write.table( tapply(m, data$order, sum), file="OUT_EXONS-autosomes", quote=T, row.names=T, append=T)
		
		order_intervals = cut(data$order,c(0,1,2,4,109))    #Revisar factores
        #summary(order_intervals)
		cat("#tapply sum m~order_intervals", file="OUT_EXONS-autosomes", append=T)
        write.table( tapply(m, order_intervals, sum), file="OUT_EXONS-autosomes", quote=T, row.names=T, append=T)

        png( "EXONS-autosomes-FEAT3Order_omega4f-orderIntervals.png", width = 1920, height = 1080 )
		boxplot(omega_4f~order_intervals,outline=F,xlab="Exon Order",ylab="Ka/Ks") 
		abline(h=median(omega_4f),col="black")
        dev.off()
		
        png( "EXONS-autosomes-FEAT3Order_m-orderIntervals.png", width = 1920, height = 1080 )
        boxplot(m~order_intervals,outline=F,xlab="Exon Order",ylab="m") 
        dev.off()

        cat( cor.test(data$order, omega_4f, method="spearman")$estimate, file = "OUT_EXONS-autosomes",fill=T, labels="Omega_4f Spearman Rho", append=T)
        cat( cor.test(data$order, omega_4f, method="spearman")$p.value, file = "OUT_EXONS-autosomes", fill=T, labels="Omega_4f Spearman p-value", append=T)
		cat( cor.test(data$order, omega_ins, method="spearman")$estimate, file = "OUT_EXONS-autosomes" , fill=T, labels="Omega_ins Spearman Rho", append=T)
		cat( cor.test(data$order, omega_ins, method="spearman")$p.value, file = "OUT_EXONS-autosomes" , fill=T, labels="Omega_ins Spearman p-value", append=T)

        cat( kruskal.test( omega_4f ~ order_intervals )$statistic, file = "OUT_EXONS-autosomes", fill=T, labels="Omega_4f K-W chi-squared", append=T  )
		cat( kruskal.test( omega_4f ~ order_intervals )$p.value, file = "OUT_EXONS-autosomes", fill=T, labels="Omega_4f K-W p-value", append=T  )
		cat( kruskal.test( omega_ins ~ order_intervals )$statistic, file = "OUT_EXONS-autosomes", fill=T, labels="Omega_ins K-W chi-squared", append=T  )
		cat( kruskal.test( omega_ins ~ order_intervals )$p.value, file = "OUT_EXONS-autosomes", fill=T, labels="Omega_ins K-W p-value", append=T  )
		

    #Feature 4
	### Exon Length ###
	    cat("\n", file="OUT_EXONS-autosomes", append=T)
	    cat("##Feature4: Exon length\n", file="OUT_EXONS-autosomes", append=T)

	    m5 = cut(m,quantile(m,(0:5)/5))
	    m3 = cut(m,c(0,405,733,12903))
		cat("#tapply sum m~m5", file="OUT_EXONS-autosomes", append=T)
	    write.table( tapply(m,m5,sum), file="OUT_EXONS-autosomes", quote=T, row.names=T, append=T)
        #summary(m5)

        cat( quantile( m, ( 0:5 ) / 5 )[1], labels="Quantile 0%", file="OUT_EXONS-autosomes", fill=T, append=T )
        cat( quantile( m, ( 0:5 ) / 5 )[2], labels="Quantile 20%", file="OUT_EXONS-autosomes", fill=T, append=T )
        cat( quantile( m, ( 0:5 ) / 5 )[3], labels="Quantile 40%", file="OUT_EXONS-autosomes", fill=T, append=T )
        cat( quantile( m, ( 0:5 ) / 5 )[4], labels="Quantile 60%", file="OUT_EXONS-autosomes", fill=T, append=T )
        cat( quantile( m, ( 0:5 ) / 5 )[5], labels="Quantile 80%", file="OUT_EXONS-autosomes", fill=T, append=T )
        cat( quantile( m, ( 0:5 ) / 5 )[6], labels="Quantile 100%", file="OUT_EXONS-autosomes", fill=T, append=T )

        cat(median( m ), file="OUT_EXONS-autosomes", fill=T, labels="Median", append=T)
		cat(mean( m ), file="OUT_EXONS-autosomes", fill=T, labels="Mean", append=T)
        cat(sd( m ), file="OUT_EXONS-autosomes", fill=T, labels="SD", append=T)

        png( "EXONS-autosomes-FEAT4Elength_omega4f-m3.png", width = 1920, height = 1080 )
	    boxplot( omega_4f~m3, outline=F, xlab="Exon Length", ylab="Ka/Ks" )
	    abline( h = median( omega_4f ), col="black" )
        dev.off()

        png( "EXONS-autosomes-FEAT4Elength_omegains-m3.png", width = 1920, height = 1080 )
	    boxplot( omega_ins~m3, outline=F, xlab="Exon Length", ylab="Ka/Kins")
	    dev.off()

        cat( cor.test(m, omega_4f, method="spearman")$estimate, file = "OUT_EXONS-autosomes",fill=T, labels="Omega_4f Spearman Rho", append=T)
        cat( cor.test(m, omega_4f, method="spearman")$p.value, file = "OUT_EXONS-autosomes", fill=T, labels="Omega_4f Spearman p-value", append=T)
		cat( cor.test(m, omega_ins, method="spearman")$estimate, file = "OUT_EXONS-autosomes" , fill=T, labels="Omega_ins Spearman Rho", append=T)
		cat( cor.test(m, omega_ins, method="spearman")$p.value, file = "OUT_EXONS-autosomes" , fill=T, labels="Omega_ins Spearman p-value", append=T)

        cat( kruskal.test( omega_4f ~ m5 )$statistic, file = "OUT_EXONS-autosomes", fill=T, labels="Omega_4f K-W chi-squared", append=T  )
		cat( kruskal.test( omega_4f ~ m5 )$p.value, file = "OUT_EXONS-autosomes", fill=T, labels="Omega_4f K-W p-value", append=T  )
		cat( kruskal.test( omega_ins ~ m5 )$statistic, file = "OUT_EXONS-autosomes", fill=T, labels="Omega_ins K-W chi-squared", append=T  )
		cat( kruskal.test( omega_ins ~ m5 )$p.value, file = "OUT_EXONS-autosomes", fill=T, labels="Omega_ins K-W p-value", append=T  )
		
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
