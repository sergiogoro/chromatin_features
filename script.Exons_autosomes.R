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
setwd("/home/sergio/chromatin/analysis/R/chromatin_features/")

# Load datasets
data<-read.table(file="EXON.SHORTINTRONS.n",header=TRUE,sep="\t")

#nrow(data)
data <- subset(data, chromosome != "X")                     #AUTOSOMES
Erep <- ( data[["num_FBtrs"]] / data[["FBtrs_per_gene"]] )  # Exones donde anotaciones 5.49 y 5.50 son distintos!!    #Solo para data exon.shortintrons.n
data <- subset(data,Erep <= 1)                              #Filtramos los que num_FBtrs/FBtrs_per_gene sean mayores que 1 (Casos que ocurren x usar rev5.49 y rev5.50)

### PopGen Summary Statistics ###

data <- na.omit(data) #para el modelo mixto

pi4f    <- data[["pi_4f"]] / data[["mdmel_4f"]]
pi0f    <- data[["pi_0f"]] / data[["mdmel_0f"]]
piins   <- data[["pi_ins"]] / data[["mdmel_ins"]]
t4f     <- data[["seg_4f"]] / data[["mdmel_4f"]]
t0f     <- data[["seg_0f"]] / data[["mdmel_0f"]]
tins    <- data[["seg_ins"]] / data[["mdmel_ins"]]
k4f     <- ( data[["div_4f"]] +1 ) / ( data[["mdyak_4f"]] +1 )
k0f     <- ( data[["div_0f"]] +1 ) / ( data[["mdyak_0f"]] +1 )
kins    <- ( ( data$div_ins +1 ) / ( data[["mdyak_ins"]] +1 ) )

constraint_4f   <- t0f/t4f
constraint_ins  <- t0f/tins
omega_4f        <- ( k0f / k4f )            #Para boxplots
omega_ins       <- ( k0f / kins )           #Para boxplots
m               <- data[["mdmel_0f"]] + data[["mdmel_4f"]] + data[["mdmel_2f"]]


### Summary ###

### EXON ###

    #Feature 1
	### Transcripts (mRNAs per exon) ###
	    cat( with( data, tapply( m, num_FBtrs, sum ) ), file = "OUT_Exons-autosomes" )
		num_FBtrs_factor = cut( data[["num_FBtrs"]], c(0, 1, 2, 75) )
		
        png( "EXONS-autosomes_omega4f-numFBtrs.png", width=1920, height=1080 )
		boxplot( omega_4f~num_FBtrs_factor, outline=F, xlab="Number of Transcripts/Exon", ylab="Ka/Ks", main="Omega 4f" )
		abline( h = median( omega_4f ), col="black" )
        #dev.off()
        png( "EXONS-autosomes_omegains-numFBtrs.png", width=1920, height=1080 )
		boxplot( omega_ins ~ num_FBtrs_factor, outline=F, xlab="Number of Transcripts/Exon", ylab="Ka/Kins", main="Omega intrones pequeños" )
        abline( h = median( omega_ins ), col="black" )
        #dev.off()

        #cor.test(omega_4f, data$num_FBtrs_factor, method="spearman")
            #cor.test(omega_4f,data$num_FBtrs_factor, method="spearman")[3]  # Spearman rho
            #cor.test(omega_4f,data$num_FBtrs_factor, method="spearman")[4]  # Spearman p-value
        write.table(cor.test(omega_4f,data$num_FBtrs_factor, method="spearman")[3], file = "OUT_Exons-autosomes", quote=T, row.names=F, append=T)
        write.table( cor.test(omega_4f, data$num_FBtrs_factor, method="spearman")[4], file = "OUT_Exons-autosomes", quote=T, row.names=F, append=T  )
		write.table( cor.test(omega_ins, data$num_FBtrs_factor, method="spearman")[3], file = "OUT_Exons-autosomes" , quote=T, row.names=F, append=T  )
		write.table( cor.test(omega_ins, data$num_FBtrs_factor, method="spearman")[4], file = "OUT_Exons-autosomes" , quote=T, row.names=F, append=T  )
		
        #kruskal.test(omega_4f~num_FBtrs_factor)
            #kruskal.test(omega_4f~num_FBtrs_factor)[1] # K-W chi-squared
            #kruskal.test(omega_4f~num_FBtrs_factor)[3] # p-value
		write.table( kruskal.test( omega_4f ~ num_FBtrs_factor )[1], file = "OUT_Exons-autosomes" , quote=T, row.names=F, append=T  )
		write.table( kruskal.test( omega_4f ~ num_FBtrs_factor )[3], file = "OUT_Exons-autosomes" , quote=T, row.names=F, append=T  )
		write.table( kruskal.test( omega_ins ~ num_FBtrs_factor )[1], file = "OUT_Exons-autosomes" , quote=T, row.names=F, append=T  )
		write.table( kruskal.test( omega_ins ~ num_FBtrs_factor )[3], file = "OUT_Exons-autosomes" , quote=T, row.names=F, append=T  )

	
#    #Feature 2
#	### Exon Representativity ### 
#		Erep <-(data[["num_FBtrs"]]/data[["FBtrs_per_gene"]])
#		summary(Erep)
#		
#		alternative = cut(Erep,c(0,0.25,0.5,0.99,1))
#		boxplot(omega_4f~alternative,outline=F,xlab="Exon Representativity",ylab="Ka/Ks")
#		abline(h=median(omega_4f),col="black")
#		with(data,tapply(m,alternative,sum))
#
#		hist(Erep,breaks=c(100))
#		quantile(Erep,(0:2)/2)
#		median(Erep)
#		mean(Erep)
#		sd(Erep)
#		cor.test(omega_4f,Erep,method="spearman") #
#		cor.test(omega_ins,Erep,method="spearman")
#		kruskal.test(omega_4f~alternative)
#		kruskal.test(omega_ins~alternative)
#		
#	### Order ###
#		with(data,tapply(m,order,sum))	
#		cor.test(data$order,omega_4f,method="spearman") 
#		cor.test(data$order,omega_ins,method="spearman") 
#		
#		order_intervals = cut(data$order,c(0,1,2,4,109))    #Revisar factores
#		summary(order_intervals)
#		with(data,tapply(m,order_intervals,sum))	
#
#		boxplot(omega_4f~order_intervals,outline=F,xlab="Exon Order",ylab="Ka/Ks") 
#		abline(h=median(omega_4f),col="black")
#		kruskal.test(omega_4f~order_intervals)
#		kruskal.test(omega_ins~order_intervals)
#		
#		boxplot(m~order_intervals,outline=F,xlab="Exon Order",ylab="m") 
#
#	### Exon Length ###
#		
#	cor.test(m,omega_4f,method="spearman") #  gene |  exon
#	cor.test(m,omega_ins,method="spearman") #  gene |  exon
#
#	m5 = cut(m,quantile(m,(0:5)/5))
#	m3 = cut(m,c(0,405,733,12903))
#	
#	with(data,tapply(m,m5,sum))
#	summary(m5)
#	quantile(m,(0:5)/5)
#	median(m)
#	mean(m)
#	sd(m)
#	
#	boxplot(omega_4f~m3,outline=F,xlab="Exon Length",ylab="Ka/Ks") 
#	abline(h=median(omega_4f),col="black")
#	
#	boxplot(omega_ins~m3,outline=F,xlab="Exon Length",ylab="Ka/Kins")
#	kruskal.test(omega_4f~m5)
#	kruskal.test(omega_ins~m5)
#		
#### GENE ###
#
#	### Messenger Complexity ###
#	mcomp <- (data$FBtrs_per_gene/data$num_cds)
#	summary(mcomp)
#	hist(mcomp,breaks=c(50))
#	hist(data$context_distance,breaks=c(100))
#	
#	quantile(mcomp,(0:5)/5)
#	median(mcomp)
#	mean(mcomp)
#	sd(mcomp)
#	
#	mcomp2 = cut(mcomp,c(0,.25,0.4,0.5))
#	mcomp5 = cut(mcomp,quantile(mcomp,(0:5)/5))
#	boxplot(omega_4f~mcomp5,outline=F,xlab="Number of Transcripts/Number of Exons",ylab="Ka/Ks") #aquellos genes con mas transcritos q exones tienen menor omega q aquellos con menos transcritos q exones.
#	abline(h=median(omega_4f),col="black")
#	boxplot(omega_ins~mcomp5,outline=F,xlab="Number of Transcripts/Number of Exons",ylab="Ka/Kins") #aquellos genes con mas transcritos q exones tienen menor omega q aquellos con menos transcritos q exones.
#	abline(h=median(omega_ins),col="black")	
#	        kruskal.test(omega_4f~mcomp5)
#	        kruskal.test(omega_ins~mcomp5)
#
#	with(data,tapply(m,mcomp5,sum))
#	summary(mcomp5)
#	
#	### chromosome and chromatin ###
#	with(data,tapply(m,chromosome,sum))        
#	with(data,tapply(m,chromatin,sum))        
#	kruskal.test(omega_4f~data[["chromosome"]])
#	kruskal.test(omega_ins~data[["chromosome"]])   
#	kruskal.test(omega_4f~data[["chromatin"]])
#	kruskal.test(omega_ins~data[["chromatin"]]) 	
#	
#
#	### Transcripts ###
#	with(data,tapply(m,FBtrs_per_gene,sum))
#	data <-subset(data, FBtrs_per_gene < 11)
#	FBtrs = cut(data[["FBtrs_per_gene"]],c(0,1,2,5,75))
#	summary(FBtrs)
#	with(data,tapply(m,FBtrs,sum))
#	
#		boxplot(omega_4f~FBtrs,outline=F,xlab="Number of Transcripts/Gene",ylab="Ka/Ks")
#		abline(h=median(omega_4f),col="black")
#		boxplot(omega_ins~data[["FBtrs_per_gene"]],outline=F,xlab="Number of Transcripts/Gene",ylab="Ka/Kins")
#		cor.test(omega_4f,data$FBtrs_per_gene,method="spearman")
#		cor.test(omega_ins,data$FBtrs_per_gene,method="spearman")
#		kruskal.test(omega_4f~FBtrs)
#		kruskal.test(omega_ins~data[["FBtrs_per_gene"]])
#
#		boxplot(pi4f~data[["FBtrs_per_gene"]],outline=F)
#		boxplot(piins~data[["FBtrs_per_gene"]],outline=F)
#		cor.test(pi4f,data$FBtrs_per_gene,method="spearman") #-0.069
#		cor.test(piins,data$FBtrs_per_gene,method="spearman") #ns
#
#	### Number of Exons ###
#	with(data,tapply(m,num_cds,sum))
#	data <-subset(data, num_cds < 11)
#	exons = cut(data[["num_cds"]],c(0,1,2,3,4,5,6,8,11,16,114))     #Revisar factores
#	summary(exons)
#	with(data,tapply(m,exons,sum))
#
#		boxplot(omega_4f~exons,outline=F,xlab="Number of Exons/Gene",ylab="Ka/Ks")
#		abline(h=median(omega_4f),col="black")
#		
#		boxplot(omega_ins~data[["num_cds"]],outline=F,xlab="Number of Exons/Gene",ylab="Ka/Kins")
#		
#		
#		cor.test(omega_4f,data$num_cds,method="spearman") # -0.2219826 
#		cor.test(omega_ins,data$num_cds,method="spearman") # 
#		kruskal.test(omega_4f~exons)
#		kruskal.test(omega_ins~data[["num_cds"]])
#
#		boxplot(pi4f~data[["num_cds"]],outline=F)
#		boxplot(piins~data[["num_cds"]],outline=F)
#		cor.test(pi4f,data$num_cds,method="spearman") # +0.037
#		cor.test(piins,data$num_cds,method="spearman") # +0.093
#
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
