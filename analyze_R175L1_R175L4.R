
source("/Users/dlundberg/Documents/abt6/scripts/R/functions/microbiome_custom_functions.R")
library(vegan)
library(reshape)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(gplots)

date=format(Sys.Date(), format="%Y%m%d")

LANE="L4_REDO"

#16S
if(LANE=="L4"){   #
	otutab="/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/Load_HiSeq175_lane4_20200428.txt"
	tax="/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/HiSeq0175_all_trimmed_lane4.tax"
	metadataT="/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/metadata_HiSeq175.txt"
}
if(LANE=="L1_REDO"){
	otutab="/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175_REDO/R175L1_zOTUtab_20200517.txt"
	tax="/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175_REDO/R175L1_taxonomy_20200517.tax"  
	metadataT="/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175_REDO/REDO_metadata.txt"
}
if(LANE=="L4_REDO"){
	otutab="/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175_REDO/HiSeq0175_lane4_20200517_OTUtabv2.txt"
	tax="/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175_REDO/HiSeq0175_lane4_20200517_allv2.tax"
	metadataT="/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175_REDO/REDO_metadata.txt"
}
if(LANE=="L1meta"){
	otutab="/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/Load_HiSeq175_lane1meta_20200430.txt"
	tax="/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/HiSeq0175_all_trimmed_lane1meta.tax"
	metadataT="/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/metadata_HiSeq175.txt"
}
if(LANE=="L1other"){
	otutab="/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/Load_HiSeq175_lane1other_20200430.txt"
	tax="/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/HiSeq0175_all_trimmed_lane1other.tax"
	metadataT="/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/metadata_HiSeq175.txt"
}




#load OTU table
	read.countTable(file=otutab)->otureads
	
#load metadata
	as.matrix(read.table(metadataT, sep="\t"))->metadata
	metadata[1,]->colnames(metadata)
	metadata[2:nrow(metadata),]->metadata

#fix sample names based on metadata
if(LANE!="L1_REDO" & LANE!="L4_REDO"){
	colnames(otureads)=metadata[match(colnames(otureads), metadata[,"name"]),"sample"]
}
if(LANE=="L1_REDO" | LANE=="L4_REDO"){
	print("REDO")
	colnames(otureads)=metadata[match(colnames(otureads), metadata[,1]),"OTUtableID"]
}


#load taxonomy
	as.matrix(read.table(file=tax, sep="\t"))->taxonomy
	
#remove contaminants
	otureads=remove_contaminant_taxa(otureads, taxonomy, keywords=c("Chloroplast", "Mitochondria"))
			
#keep samples with more than 1000 total reads
	otureads=otureads[,which(colSums(otureads)>=1000)]
	
#find low abundance OTUs to later remove from whole dataset.
	#minimum_readcount=25
	#minimum_samples=5
	#quantitative_OTUs=rownames(otureads)[which(apply(otureads, 1, max)>=minimum_readcount)]
	#quantitative_OTUs=c("HOST", quantitative_OTUs)

#normalize
	#otureads=normalize100(otureads)

#remove low abundance rows 0.5% used 20200229 for DISPLAY ONLY. 0.05 good for 16S comparison
	#otureads=low_abundance(otureads, percent=0.05)
	
#order full table by rowsums
	order(rowSums(otureads), decreasing=FALSE)->ordering_vector
		otureads[ordering_vector,]->otureads

#next set specific order
	toporder="low_abundance"
	bottomorder=c("Otu3", "HOSTArabidopsisthalianaGiganteaGI502")		
	otureads=topOrder(otureads, toporder, bottomorder)->currentreads




#>
	#>
		#>
			#> #pepper infiltration LANE4
		#>
	#>
#>
#first, load growth curve data
pepperCFU<-read.table("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/pepper_growth_curve/pepper_CFU_counts.txt")
as.matrix(pepperCFU[1,])->names(pepperCFU)
pepperCFU[2:nrow(pepperCFU),]->pepperCFU
pepperCFU=as.data.frame(pepperCFU)
pepperCFU=pepperCFU[,c("group", "sample", "experiment", "CFU_cm2")]
pepperCFU$groupsample=paste(pepperCFU$group, pepperCFU$sample, sep="_")
pepperCFU$CFU_cm2=as.numeric(as.matrix(pepperCFU$CFU_cm2))
infiltrationCFU=pepperCFU[pepperCFU$experiment=="infiltration",]

infiltrationCFU[which(infiltrationCFU$sample==1),]

CFU_4=infiltrationCFU$CFU_cm2[which(infiltrationCFU$group=="10exp4")]
CFU_5=infiltrationCFU$CFU_cm2[which(infiltrationCFU$group=="10exp5")]
CFU_6=infiltrationCFU$CFU_cm2[which(infiltrationCFU$group=="10exp6")]
CFU_7=infiltrationCFU$CFU_cm2[which(infiltrationCFU$group=="10exp7")]
CFU_8=infiltrationCFU$CFU_cm2[which(infiltrationCFU$group=="10exp8")]
CFU_infiltration=c(CFU_4, CFU_5, CFU_6, CFU_7, CFU_8)

#next load qPCR data
inqPCR<-read.table("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/qPCR_pepper/qPCR_inf.txt")
qPCR_4=inqPCR[,2][which(inqPCR[,3]=="in10e4")]
qPCR_5=inqPCR[,2][which(inqPCR[,3]=="in10e5")]
qPCR_6=inqPCR[,2][which(inqPCR[,3]=="in10e6")]
qPCR_7=inqPCR[,2][which(inqPCR[,3]=="in10e7")]
qPCR_8=inqPCR[,2][which(inqPCR[,3]=="in10e8")]
qPCR_in=c(qPCR_4, qPCR_5, qPCR_6, qPCR_7, qPCR_8)


#pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq168/LOAD_HiSeq_", current, "_untrimmed", date, ".pdf", sep="", collapse=""), width = 3, height = 2, useDingbats=FALSE)
#ggplot(histDL, aes(x=histDL$group, y=histDL$CFU_cm2, group=histDL$group, color=histDL$group)) +
#  geom_line() + geom_boxplot(color="#444444") + geom_jitter(size=2, shape=16, position=position_jitter(0.2)) +
#  theme_minimal() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
#  theme(legend.position="none") + 
#  scale_y_continuous(name = "log10 CFU / cm2", breaks = 0:9, labels = 0:9, limits = c(0,9)) + 
#  scale_color_manual(values=rep("#000000", 20)) 


InfiltRep1=grep("S61xaltV4CaGIxin10e4.1|S62xaltV4CaGIxin10e5.1|S63xaltV4CaGIxin10e6.1|S64xaltV4CaGIxin10e7.1|S65xaltV4CaGIxin10e8.1", colnames(otureads), value=TRUE)
InfiltRep2=grep("S66xaltV4CaGIxin10e4.2|S67xaltV4CaGIxin10e5.2|S68xaltV4CaGIxin10e6.2|S69xaltV4CaGIxin10e7.2|S70xaltV4CaGIxin10e8.2", colnames(otureads), value=TRUE)
InfiltRep3=grep("S71xaltV4CaGIxin10e4.3|S72xaltV4CaGIxin10e5.3|S73xaltV4CaGIxin10e6.3|S74xaltV4CaGIxin10e7.3|S75xaltV4CaGIxin10e8.3", colnames(otureads), value=TRUE)
InfiltRep4=grep("S76xaltV4CaGIxin10e4.4|S77xaltV4CaGIxin10e5.4|S78xaltV4CaGIxin10e6.4|S79xaltV4CaGIxin10e7.4|S80xaltV4CaGIxin10e8.4", colnames(otureads), value=TRUE)
InfiltLevels=c("in10e4","in10e5","in10e6","in10e7","in10e8")


InfiltReps=otureads[,c("S61xaltV4CaGIxin10e4.1","S66xaltV4CaGIxin10e4.2","S71xaltV4CaGIxin10e4.3","S76xaltV4CaGIxin10e4.4","S62xaltV4CaGIxin10e5.1","S67xaltV4CaGIxin10e5.2","S72xaltV4CaGIxin10e5.3","S77xaltV4CaGIxin10e5.4","S63xaltV4CaGIxin10e6.1","S68xaltV4CaGIxin10e6.2","S73xaltV4CaGIxin10e6.3","S78xaltV4CaGIxin10e6.4","S64xaltV4CaGIxin10e7.1","S69xaltV4CaGIxin10e7.2","S74xaltV4CaGIxin10e7.3","S79xaltV4CaGIxin10e7.4","S65xaltV4CaGIxin10e8.1","S70xaltV4CaGIxin10e8.2","S75xaltV4CaGIxin10e8.3","S80xaltV4CaGIxin10e8.4")]
InfiltReps=InfiltReps[which(rowSums(InfiltReps)>0),]

InfiltReps=normalize100(InfiltReps)

#FIND KEY OTUS
if(LANE=="L4"){
	HOST=colSums(matrix_format(InfiltReps[grep("Otu4$|Otu12890$", rownames(InfiltReps)),]))
	InfiltReps=InfiltReps[setdiff(1:nrow(InfiltReps), grep("Otu4$|Otu12890$", rownames(InfiltReps))),]
	Xe8510=InfiltReps["Otu1",]
	InfiltReps=InfiltReps[setdiff(rownames(InfiltReps), c("Otu1")),]
}

if(LANE=="L4_REDO"){
	HOST=colSums(matrix_format(InfiltReps[grep("Otu4$|Otu12886$", rownames(InfiltReps)),]))
	InfiltReps=InfiltReps[setdiff(1:nrow(InfiltReps), grep("Otu4$|Otu12886$", rownames(InfiltReps))),]
	Xe8510=InfiltReps["Otu1",]
	InfiltReps=InfiltReps[setdiff(rownames(InfiltReps), c("Otu1")),]
}

#convert to family, add host afterwards
InfiltReps<-filter_by_taxonomyic_level(countTable=InfiltReps, taxonomy=taxonomy, keywords=c("c"), last_filter=TRUE)
InfiltReps=InfiltReps[which(rowSums(InfiltReps)>0),]
InfiltReps=InfiltReps[apply(InfiltReps, 1, median)>=0.1,]

#add back host
InfiltReps=rbind(InfiltReps, Xe8510, HOST)


order(rowSums(InfiltReps), decreasing=FALSE)->ordering_vector
		InfiltReps[ordering_vector,]->InfiltReps

#next set specific order
	toporder=c("empty", "low_abundance")
	#bottomorder=c("Comamonadaceae", "Geodermatophilaceae", "Oxalobacteraceae", "Sphingomonadaceae", "Pseudomonadaceae", "Streptomycetaceae", "Nocardioidaceae", "Chitinophagaceae", "Xanthomonadaceae", "Xe8510", "HOST")		
	bottomorder=c("Actinobacteria", "Gammaproteobacteria", "Betaproteobacteria", "Alphaproteobacteria", "HOST","Xe8510")
	InfiltReps=topOrder(InfiltReps, toporder, bottomorder)
	InfiltReps=InfiltReps[nrow(InfiltReps):1,]

melt(InfiltReps)->histDL
c("organism", "sample_name", "abundance")->names(histDL)
gsub(".*xin", "in", histDL$sample_name)->histDL$sample_category
	gsub("\\..*", "", histDL$sample_category)->histDL$sample_category
gsub(".*\\.", "", histDL$sample_name)->histDL$replicate
paste(histDL$replicate, histDL$organism, sep="_")->histDL$org_rep
histDL$sample_category=factor(histDL$sample_category, levels=InfiltLevels)
	
	unique(as.matrix(histDL$organism))->taxa_color_pairs
	COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
	#ALPHA_COLUMN_IN_COLORSLIST=7
	source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
	"black"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Hpa"),2]
	"blue"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Xe8510"),2]
	"green"->taxa_color_pairs[which(taxa_color_pairs[,1]=="HOST"),2]
	"darkgreen"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Nocardioidaceae"),2]

	taxa_color_pairs[match(levels(histDL$organism), taxa_color_pairs[,1]),2]->histDL$colors
	
	
#relative abundance
#pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq168/LOAD_HiSeq_Titration_ReAb_", date, ".pdf", sep="", collapse=""), width = 5, height = 3, useDingbats=FALSE)
ggplot(histDL, aes(x=histDL$sample_category, y=histDL$abundance, group=histDL$org_rep, color=histDL$organism)) +
  geom_line() + geom_point(size=2) + scale_color_manual(values=histDL$colors) + 
  theme_minimal() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
  theme(legend.position="none") 
#dev.off() 

pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/PepperInfiltBact_logReAb_", date, ".pdf", sep="", collapse=""), width = 7, height = 4, useDingbats=FALSE)
ggplot(histDL, aes(x=sample_category, y=log10(abundance), group=org_rep, color=organism)) +
  geom_line() + geom_point(size=2) + scale_color_manual(values=histDL$colors) + 
  theme_classic() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
  theme(legend.position="none") +
  theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold")) + 
    theme(axis.text.x = element_text(angle=90, hjust=0.99, vjust=0.9, color="black")) + 
	theme(axis.ticks = element_line(colour = 'black')) + 
	   	scale_x_discrete(name = "injection") + 
  scale_y_continuous(name = "RA (%)", breaks=c(-3:2), labels=10^c(-3:2), limits = c(-3,2)) #+ 
 dev.off()
 	 

InfiltReps["HOST", grep("in10e8",colnames(InfiltReps))]


InfiltReps_load=InfiltReps

divide_by_host=function(tab){
	tab["HOST", which(tab["HOST",]==0)]<-0.00000001
	for(c in 1:ncol(tab)){
		tab[,c]/tab["HOST",c]->tab[,c]
		if( max(tab[,c])>100 ){
			tab[,c]=200*tab[,c]/sum(tab[,c])
		}
	}
	#tab[which(tab>100)]<-100
	return(tab)
}
	
InfiltReps_load=divide_by_host(InfiltReps_load)

InfiltReps_load=rbind(InfiltReps_load, CFU_infiltration, qPCR_in)


melt(InfiltReps_load)->histDL
c("organism", "sample_name", "abundance")->names(histDL)
gsub(".*xin", "in", histDL$sample_name)->histDL$sample_category
	gsub("\\..*", "", histDL$sample_category)->histDL$sample_category
gsub(".*\\.", "", histDL$sample_name)->histDL$replicate
paste(histDL$replicate, histDL$organism, sep="_")->histDL$org_rep
histDL$sample_category=factor(histDL$sample_category, levels=InfiltLevels)
histDL[histDL$organism!="HOST",]->histDL
#histDL[histDL$organism!="CFU_infiltration",]->histDL
histDL$organism=factor(histDL$organism, levels=as.character(unique(histDL$organism)))
	levels(histDL$organism)->taxa_color_pairs
	COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
	#ALPHA_COLUMN_IN_COLORSLIST=7
	source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
	"blue"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Xe8510"),2]
	"green"->taxa_color_pairs[which(taxa_color_pairs[,1]=="HOST"),2]
	"darkgreen"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Nocardioidaceae"),2]
	"#7FC97F"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Actinobacteria"),2]
	"#BEAED4"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Gammaproteobacteria"),2]
	"#FDC086"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Betaproteobacteria"),2]
	"#FFFF99"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Alphaproteobacteria"),2]
	taxa_color_pairs[match(levels(histDL$organism), taxa_color_pairs[,1]),2]->histDL$colors
	
#boxplot just 16S data
#histDL2=histDL[grep("Oxalobacteraceae|Sphingomonadaceae|Pseudomonadaceae|Streptomycetaceae|Nocardioidaceae|Chitinophagaceae|Xanthomonadaceae|Xe8510", histDL$organism, invert=FALSE),]
histDL2=histDL[grep("Actinobacteria|Gammaproteobacteria|Betaproteobacteria|Alphaproteobacteria|Xe8510", histDL$organism, invert=FALSE),]
histDL2$organism=factor(histDL2$organism, levels=unique(as.character(histDL2$organism)))
log10(histDL2$abundance+.00000000000001)->histDL2$abundance
date=format(Sys.Date(), format="%Y%m%d")
pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/PepperInfiltBact_", date, ".pdf", sep="", collapse=""), width = 7, height = 4, useDingbats=FALSE)
ggplot(histDL2, aes(fill=organism, y=abundance, x=sample_category)) + 
   	 	geom_boxplot(outlier.size=0, outlier.shape=NA, alpha=0.5, width=.7, position=position_dodge(.9)) + 
   	 	geom_point(aes(fill=organism), size = 1.5, stroke=0.1, shape = 21, position = position_jitterdodge(jitter.width=0.1, dodge=.9)) +
   	 	theme_classic() + scale_color_manual(values=histDL2$colors) + scale_fill_manual(values=histDL2$colors) + 
    	theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=0.99, vjust=0.9, color="black")) + 
    	theme(axis.ticks = element_line(colour = 'black')) +
    	scale_x_discrete(name = "injection") + 
 		scale_y_continuous(name = "log10 abundance", breaks=c(-3:3), labels=10^c(-3:3), limits = c(-3,3)) 
dev.off()



#boxplot ONLY CFU
melt(InfiltReps_load)->histDL
c("organism", "sample_name", "abundance")->names(histDL)
gsub(".*xin", "in", histDL$sample_name)->histDL$sample_category
	gsub("\\..*", "", histDL$sample_category)->histDL$sample_category
gsub(".*\\.", "", histDL$sample_name)->histDL$replicate
paste(histDL$replicate, histDL$organism, sep="_")->histDL$org_rep
histDL$sample_category=factor(histDL$sample_category, levels=InfiltLevels)
histDL[histDL$organism!="HOST",]->histDL
histDL$organism=factor(histDL$organism, levels=as.character(unique(histDL$organism)))
	levels(histDL$organism)->taxa_color_pairs
	COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
	#ALPHA_COLUMN_IN_COLORSLIST=7
	source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
	"blue"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Xe8510"),2]
	"lightblue"->taxa_color_pairs[which(taxa_color_pairs[,1]=="HOST"),2]
	"darkgreen"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Nocardioidaceae"),2]
	"lightblue"->taxa_color_pairs[which(taxa_color_pairs[,1]=="CFU_infiltration"),2]
	taxa_color_pairs[match(levels(histDL$organism), taxa_color_pairs[,1]),2]->histDL$colors
	
	
histDL3=histDL[grep("CFU_infiltration|Xe8510|qPCR_in", histDL$organism),]
histDL3$organism=factor(histDL3$organism, levels=as.character(c("CFU_infiltration", "qPCR_in","Xe8510")))

hamPCR_median=median(histDL3$abundance[which(histDL3$organism=="Xe8510" & histDL3$sample_category=="in10e6")])
qPCR_median=median(histDL3$abundance[which(histDL3$organism=="qPCR_in" & histDL3$sample_category=="in10e6")])
CFU_median=median(histDL3$abundance[which(histDL3$organism=="CFU_infiltration" & histDL3$sample_category=="in10e6")])

qPCR_scale_factor=CFU_median/qPCR_median
hamPCR_scale_factor=CFU_median/hamPCR_median
histDL3$abundance[which(histDL3$organism=="Xe8510")]=hamPCR_scale_factor*histDL3$abundance[which(histDL3$organism=="Xe8510")]
histDL3$abundance[which(histDL3$organism=="qPCR_in")]=qPCR_scale_factor*histDL3$abundance[which(histDL3$organism=="qPCR_in")]

log10(histDL3$abundance)->histDL3$abundance

#ALIGN hamPCR to the QPCR.
hamPCR_values=histDL3$abundance[which(histDL3$organism=="Xe8510")]
qPCR_values=histDL3$abundance[which(histDL3$organism=="qPCR_in")]
CFU_values=histDL3$abundance[which(histDL3$organism=="CFU_infiltration")]
ham_v_qpcr=lm(hamPCR_values~qPCR_values)[[1]]
slope=ham_v_qpcr[2]
intercept=ham_v_qpcr[1]
#abline(intercept, slope)
par(mfrow=c(1,2))
plot(qPCR_values, hamPCR_values)
plot(slope*qPCR_values + intercept, hamPCR_values)
histDL3$abundance[which(histDL3$organism=="qPCR_in")]<-(slope*histDL3$abundance[which(histDL3$organism=="qPCR_in")]) + intercept


#ALIGN hamPCR and QPCR to the CFU, but only the middle concentrations
middle_range=which(match(histDL3$sample_category, c("in10e4","in10e5","in10e6", "in10e7", "in10e8"), nomatch=0)>0)
hamPCR_values=histDL3$abundance[intersect(which(histDL3$organism=="Xe8510"), middle_range)]
qPCR_values=histDL3$abundance[intersect(which(histDL3$organism=="qPCR_in"), middle_range)]
CFU_values=histDL3$abundance[intersect(which(histDL3$organism=="CFU_infiltration"), middle_range)]
ham_v_CFU=lm(CFU_values~hamPCR_values)[[1]]
slope=ham_v_CFU[2]
intercept=ham_v_CFU[1]
abline(intercept, slope)
par(mfrow=c(1,2))
plot(hamPCR_values, CFU_values)
plot(slope*hamPCR_values + intercept, qPCR_values)
histDL3$abundance[which(histDL3$organism=="Xe8510")]<-(slope*histDL3$abundance[which(histDL3$organism=="Xe8510")]) + intercept
histDL3$abundance[which(histDL3$organism=="qPCR_in")]<-(slope*histDL3$abundance[which(histDL3$organism=="qPCR_in")]) + intercept



date=format(Sys.Date(), format="%Y%m%d")
pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/PepperInfiltCFU_", date, ".pdf", sep="", collapse=""), width = 4, height = 4, useDingbats=FALSE)
ggplot(histDL3, aes(fill=organism, y=abundance, x=sample_category)) + 
   	 	geom_boxplot(outlier.size=0, outlier.shape=NA, alpha=0.75, width=.7, position=position_dodge(.9)) + 
   	 	geom_point(aes(fill=organism), size = 1.5, stroke=0.1, shape = 21, position = position_jitterdodge(jitter.width=0.1, dodge=.9)) +
   	 	theme_classic() + scale_color_manual(values=c("lightblue", "orange", "blue")) + scale_fill_manual(values=c("lightblue", "orange", "blue")) + 
    	theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=0.99, vjust=0.9, color="black")) + 
    	theme(axis.ticks = element_line(colour = 'black')) +
    	scale_x_discrete(name = "injection") + 
 		scale_y_continuous(name = "log10 abundance", breaks=c(0:9), labels=c(0:9), limits = c(0,9)) 
dev.off()

#line plots only
#melt(InfiltReps_load)->histDL
#c("organism", "sample_name", "abundance")->names(histDL)
#gsub(".*xin", "in", histDL$sample_name)->histDL$sample_category
#	gsub("\\..*", "", histDL$sample_category)->histDL$sample_category
#gsub(".*\\.", "", histDL$sample_name)->histDL$replicate
#paste(histDL$replicate, histDL$organism, sep="_")->histDL$org_rep
#histDL$sample_category=factor(histDL$sample_category, levels=InfiltLevels)
#histDL[histDL$organism!="HOST",]->histDL
#histDL$organism=factor(histDL$organism, levels=as.character(c("CFU_infiltration", "Xanthomonadaceae","Pseudomonadaceae","Chitinophagaceae","Sphingomonadaceae","Enterobacteriaceae")))
#log10(histDL$abundance)->histDL$abundance
#	rownames(InfiltReps_load)->taxa_color_pairs
#	COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
#	#ALPHA_COLUMN_IN_COLORSLIST=7
#	source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
#	taxa_color_pairs[match(histDL$organism, taxa_color_pairs[,1]),2]->histDL$colors
#
#
##pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/MiSeq_1/Pepper_infiltration_lines_", date, ".pdf", sep="", collapse=""), width = 5, height = 3, useDingbats=FALSE)
#ggplot(histDL, aes(x=histDL$sample_category, y=histDL$abundance, group=histDL$org_rep, color=histDL$organism)) +
#  geom_line() + geom_point(size=0) + scale_color_manual(values=histDL$colors) + 
#  theme_minimal() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
#  theme(legend.position="none")  +
#   scale_x_discrete(name = "injection") + 
# 	scale_y_continuous(name = "log10 abundance", breaks=c(0:7), limits = c(0,7)) +
#	theme(panel.grid.minor = element_line(size = 0.15), panel.grid.major = element_line(size = .15)) + 
#	theme(panel.grid.minor = element_line(colour = "black"), panel.grid.major = element_line(colour = "black")) + 
#	theme(panel.grid.minor = element_blank()) 
##dev.off()
#
#

#correlate growth curve vs. 16S data  
melt(InfiltReps_load)->histDL
c("organism", "sample_name", "abundance")->names(histDL)
gsub(".*xin", "in", histDL$sample_name)->histDL$sample_category
	gsub("\\..*", "", histDL$sample_category)->histDL$sample_category
gsub(".*\\.", "", histDL$sample_name)->histDL$replicate
paste(histDL$replicate, histDL$organism, sep="_")->histDL$org_rep
histDL$sample_category=factor(histDL$sample_category, levels=InfiltLevels)
histDL[histDL$organism!="HOST",]->histDL
histDL$organism=factor(histDL$organism, levels=as.character(unique(histDL$organism)))

histDL3=histDL[grep("CFU_infiltration|Xe8510", histDL$organism),]
histDL3$organism=factor(histDL3$organism, levels=as.character(c("Xe8510","CFU_infiltration")))
correlation=cor(log10(InfiltReps_load["CFU_infiltration",]), log10(InfiltReps_load["Xe8510",]))^2
histDL3$X=log10(histDL3$abundance[which(histDL3$organism=="Xe8510")])
histDL3$Y=log10(histDL3$abundance[which(histDL3$organism=="CFU_infiltration")])
	ggplot(histDL3, aes(x=X, y=Y)) + geom_point(size=1.5, alpha=0.2) + 
			theme_minimal() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
			theme(legend.position="none") + 
			coord_fixed(ratio=1/1) +  
				theme(panel.grid.minor = element_line(size = 0.15), panel.grid.major = element_line(size = .15)) + 
			theme(panel.grid.minor = element_line(colour = "black"), panel.grid.major = element_line(colour = "black")) + 
			theme(panel.grid.minor = element_blank()) +
			geom_smooth(method='lm', size=0.5, se = FALSE) +
			scale_x_continuous(name = "a + b", breaks = c(-3:3), labels = c(-3:3), limits = c(-3,3)) + 
			scale_y_continuous(name = "c + d", breaks = c(2:7), labels = c(2:7), limits = c(2,7)) 
	# 
			#annotate(geom="text", size=3, x = 16^(1/exponent), y = 36^(1/exponent), label = deparse(bquote(R^2~"="~.(correlation))), parse=T) +
			#annotate(geom="text", size=3, x = 16^(1/exponent), y = 24^(1/exponent), label = deparse(paste("K-S = ", round(ks[[2]],3), sep="")), parse=T) 


e4=as.vector(InfiltReps_load[grep("Actinobacteria|Gammaproteobacteria|Betaproteobacteria|Alphaproteobacteria", rownames(InfiltReps_load)),grep(".*10e4.*",colnames(InfiltReps_load))])
e5=as.vector(InfiltReps_load[grep("Actinobacteria|Gammaproteobacteria|Betaproteobacteria|Alphaproteobacteria", rownames(InfiltReps_load)),grep(".*10e5.*",colnames(InfiltReps_load))])
e6=as.vector(InfiltReps_load[grep("Actinobacteria|Gammaproteobacteria|Betaproteobacteria|Alphaproteobacteria", rownames(InfiltReps_load)),grep(".*10e6.*",colnames(InfiltReps_load))])
e7=as.vector(InfiltReps_load[grep("Actinobacteria|Gammaproteobacteria|Betaproteobacteria|Alphaproteobacteria", rownames(InfiltReps_load)),grep(".*10e7.*",colnames(InfiltReps_load))])
e8=as.vector(InfiltReps_load[grep("Actinobacteria|Gammaproteobacteria|Betaproteobacteria|Alphaproteobacteria", rownames(InfiltReps_load)),grep(".*10e8.*",colnames(InfiltReps_load))])

all=melt(rbind(e4, e5, e6, e7, e8))
names(all)=c("category", "b", "abundance")

pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/PepperINF_increaseCommensal_", date, ".pdf", sep="", collapse=""), width = 4, height = 4, useDingbats=FALSE)
ggplot(all, aes(y=log10(abundance), x=category)) + 
   	 	geom_boxplot(outlier.size=0, outlier.shape=NA, alpha=0.5, width=.7, position=position_dodge(.9)) + 
   	 	geom_point(aes(fill=category), size = 2.5, color="black", alpha=0.3, stroke=0.1, shape = 21, position = position_jitterdodge(jitter.width=0.9, dodge=.9)) +
   	 	theme_classic() + scale_color_manual(values=rep("black", 5)) + scale_fill_manual(values=rep("black", 5)) + 
    	theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=0.99, vjust=0.9, color="black")) + 
    	theme(axis.ticks = element_line(colour = 'black')) +
    	scale_x_discrete(name = "injection") + 
  		scale_y_continuous(name = "log10 abundance", breaks=c(-4:0), labels=c(-4:0), limits = c(-4,.5)) 
dev.off()

wilcox.test(e8, e4, "greater")





#>
	#>
		#>
			#> #pepper growth curve Lane 4
		#>
	#>
#>
#first, load growth curve data
pepperCFU<-read.table("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/pepper_growth_curve/pepper_CFU_counts.txt")
as.matrix(pepperCFU[1,])->names(pepperCFU)
pepperCFU[2:nrow(pepperCFU),]->pepperCFU
pepperCFU=as.data.frame(pepperCFU)
pepperCFU=pepperCFU[,c("group", "sample", "experiment", "CFU_cm2")]
pepperCFU$groupsample=paste(pepperCFU$group, pepperCFU$sample, sep="_")
pepperCFU$CFU_cm2=as.numeric(as.matrix(pepperCFU$CFU_cm2))
gcCFU=pepperCFU[pepperCFU$experiment=="growth_curve",]

CFU_Rep1=gcCFU$CFU_cm2[which(gcCFU$sample==1)]
CFU_Rep2=gcCFU$CFU_cm2[which(gcCFU$sample==2)]
CFU_Rep3=gcCFU$CFU_cm2[which(gcCFU$sample==3)]
CFU_Rep4=gcCFU$CFU_cm2[which(gcCFU$sample==4)]
CFU_Rep5=gcCFU$CFU_cm2[which(gcCFU$sample==5)]
CFU_Rep6=gcCFU$CFU_cm2[which(gcCFU$sample==6)]
CFU_gc=c(CFU_Rep1, CFU_Rep2, CFU_Rep3, CFU_Rep4, CFU_Rep5, CFU_Rep6)


#next load qPCR data
gcqPCR<-read.table("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/qPCR_pepper/qPCR_growthcurve.txt")
qPCR_Rep1=gcqPCR[,2][which(gcqPCR[,4]==1)]
qPCR_Rep2=gcqPCR[,2][which(gcqPCR[,4]==2)]
qPCR_Rep3=gcqPCR[,2][which(gcqPCR[,4]==3)]
qPCR_Rep4=gcqPCR[,2][which(gcqPCR[,4]==4)]
qPCR_Rep5=gcqPCR[,2][which(gcqPCR[,4]==5)]
qPCR_Rep6=gcqPCR[,2][which(gcqPCR[,4]==6)]
qPCR_gc=c(qPCR_Rep1, qPCR_Rep2, qPCR_Rep3, qPCR_Rep4, qPCR_Rep5, qPCR_Rep6)


gcRep1=grep("S25xaltV4CaGIxgc11dpi.1|S31xaltV4CaGIxgc9dpi.1|S37xaltV4CaGIxgc7dpi.1|S43xaltV4CaGIxgc4dpi.1|S49xaltV4CaGIxgc2dpi.1|S55xaltV4CaGIxgc0dpi.1", colnames(otureads), value=TRUE)
gcRep2=grep("S26xaltV4CaGIxgc11dpi.2|S32xaltV4CaGIxgc9dpi.2|S38xaltV4CaGIxgc7dpi.2|S44xaltV4CaGIxgc4dpi.2|S50xaltV4CaGIxgc2dpi.2|S56xaltV4CaGIxgc0dpi.2", colnames(otureads), value=TRUE)
gcRep3=grep("S27xaltV4CaGIxgc11dpi.3|S33xaltV4CaGIxgc9dpi.3|S39xaltV4CaGIxgc7dpi.3|S45xaltV4CaGIxgc4dpi.3|S51xaltV4CaGIxgc2dpi.3|S57xaltV4CaGIxgc0dpi.3", colnames(otureads), value=TRUE)
gcRep4=grep("S28xaltV4CaGIxgc11dpi.4|S34xaltV4CaGIxgc9dpi.4|S40xaltV4CaGIxgc7dpi.4|S46xaltV4CaGIxgc4dpi.4|S52xaltV4CaGIxgc2dpi.4|S58xaltV4CaGIxgc0dpi.4", colnames(otureads), value=TRUE)
gcRep5=grep("S29xaltV4CaGIxgc11dpi.5|S35xaltV4CaGIxgc9dpi.5|S41xaltV4CaGIxgc7dpi.5|S47xaltV4CaGIxgc4dpi.5|S53xaltV4CaGIxgc2dpi.5|S59xaltV4CaGIxgc0dpi.5", colnames(otureads), value=TRUE)
gcRep6=grep("S30xaltV4CaGIxgc11dpi.6|S36xaltV4CaGIxgc9dpi.6|S42xaltV4CaGIxgc7dpi.6|S48xaltV4CaGIxgc4dpi.6|S54xaltV4CaGIxgc2dpi.6|S60xaltV4CaGIxgc0dpi.6", colnames(otureads), value=TRUE)
gcLevels=c("gc0dpi","gc2dpi","gc4dpi","gc7dpi","gc9dpi","gc11dpi")


gcReps=otureads[,c("S55xaltV4CaGIxgc0dpi.1","S49xaltV4CaGIxgc2dpi.1","S43xaltV4CaGIxgc4dpi.1","S37xaltV4CaGIxgc7dpi.1","S31xaltV4CaGIxgc9dpi.1","S25xaltV4CaGIxgc11dpi.1","S56xaltV4CaGIxgc0dpi.2","S50xaltV4CaGIxgc2dpi.2","S44xaltV4CaGIxgc4dpi.2","S38xaltV4CaGIxgc7dpi.2","S32xaltV4CaGIxgc9dpi.2","S26xaltV4CaGIxgc11dpi.2","S57xaltV4CaGIxgc0dpi.3","S51xaltV4CaGIxgc2dpi.3","S45xaltV4CaGIxgc4dpi.3","S39xaltV4CaGIxgc7dpi.3","S33xaltV4CaGIxgc9dpi.3","S27xaltV4CaGIxgc11dpi.3","S58xaltV4CaGIxgc0dpi.4","S52xaltV4CaGIxgc2dpi.4","S46xaltV4CaGIxgc4dpi.4","S40xaltV4CaGIxgc7dpi.4","S34xaltV4CaGIxgc9dpi.4","S28xaltV4CaGIxgc11dpi.4","S59xaltV4CaGIxgc0dpi.5","S53xaltV4CaGIxgc2dpi.5","S47xaltV4CaGIxgc4dpi.5","S41xaltV4CaGIxgc7dpi.5","S35xaltV4CaGIxgc9dpi.5","S29xaltV4CaGIxgc11dpi.5","S60xaltV4CaGIxgc0dpi.6","S54xaltV4CaGIxgc2dpi.6","S48xaltV4CaGIxgc4dpi.6","S42xaltV4CaGIxgc7dpi.6","S36xaltV4CaGIxgc9dpi.6","S30xaltV4CaGIxgc11dpi.6")]
#gcReps=otureads[,c(gcRep1, gcRep2, gcRep3, gcRep4, gcRep5, gcRep6)]
gcReps=gcReps[which(rowSums(gcReps)>0),]
gcReps=normalize100(gcReps)

#FIND KEY OTUS
if(LANE=="L4"){
	HOST=colSums(matrix_format(gcReps[grep("Otu4$|Otu12890$", rownames(gcReps)),]))
	gcReps=gcReps[setdiff(1:nrow(gcReps), grep("Otu4$|Otu12890$", rownames(gcReps))),]
	Xe8510=gcReps["Otu1",]
	gcReps=gcReps[setdiff(rownames(gcReps), c("Otu1")),]
}

if(LANE=="L4_REDO"){
	HOST=colSums(matrix_format(gcReps[grep("Otu4$|Otu12886$", rownames(gcReps)),]))
	gcReps=gcReps[setdiff(1:nrow(gcReps), grep("Otu4$|Otu12886$", rownames(gcReps))),]
	Xe8510=gcReps["Otu1",]
	gcReps=gcReps[setdiff(rownames(gcReps), c("Otu1")),]
}

#convert to family, add host afterwards
gcReps<-filter_by_taxonomyic_level(countTable=gcReps, taxonomy=taxonomy, keywords=c("c"), last_filter=TRUE)



gcReps=gcReps[which(rowSums(gcReps)>0),]
gcReps=gcReps[apply(gcReps, 1, min)>0.0,]

#add back host
gcReps=rbind(gcReps, Xe8510, HOST)

order(rowSums(gcReps), decreasing=FALSE)->ordering_vector
		gcReps[ordering_vector,]->gcReps

#next set specific order
	toporder=c("empty", "low_abundance")
	#bottomorder=c("Pseudonocardiaceae","Nocardioidaceae","Comamonadaceae","Enterobacteriaceae","Rhodobacteraceae","Streptomycetaceae","Moraxellaceae","Oxalobacteraceae","Geodermatophilaceae","Micrococcaceae","Sphingomonadaceae","Pseudomonadaceae","HOST","Xe8510")
	bottomorder=c("Actinobacteria", "Gammaproteobacteria", "Betaproteobacteria", "Alphaproteobacteria", "HOST","Xe8510")
	gcReps=topOrder(gcReps, toporder, bottomorder)
	gcReps=gcReps[nrow(gcReps):1,]



melt(gcReps)->histDL
c("organism", "sample_name", "abundance")->names(histDL)
gsub(".*xgc", "gc", histDL$sample_name)->histDL$sample_category
	gsub("\\..*", "", histDL$sample_category)->histDL$sample_category
gsub(".*\\.", "", histDL$sample_name)->histDL$replicate
paste(histDL$replicate, histDL$organism, sep="_")->histDL$org_rep
histDL$sample_category=factor(histDL$sample_category, levels=gcLevels)
	
	unique(as.matrix(histDL$organism))->taxa_color_pairs
	COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
	#ALPHA_COLUMN_IN_COLORSLIST=7
	source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
	"black"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Hpa"),2]
	"blue"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Xe8510"),2]
	"green"->taxa_color_pairs[which(taxa_color_pairs[,1]=="HOST"),2]
	taxa_color_pairs[match(levels(histDL$organism), taxa_color_pairs[,1]),2]->histDL$colors
	
	
#relative abundance
#pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq168/LOAD_HiSeq_Titration_ReAb_", date, ".pdf", sep="", collapse=""), width = 5, height = 3, useDingbats=FALSE)
ggplot(histDL, aes(x=sample_category, y=abundance, group=org_rep, color=organism)) +
  geom_line() + geom_point(size=2) + scale_color_manual(values=histDL$colors) + 
  theme_minimal() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
  theme(legend.position="none") 
#dev.off()  

pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/PepperGCBact_logReAb_", date, ".pdf", sep="", collapse=""), width = 7, height = 4, useDingbats=FALSE)
ggplot(histDL, aes(x=sample_category, y=log10(abundance), group=org_rep, color=organism)) +
  geom_line() + geom_point(size=2) + scale_color_manual(values=histDL$colors) + 
  theme_classic() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
  theme(legend.position="none") +
  theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold")) + 
    theme(axis.text.x = element_text(angle=90, hjust=0.99, vjust=0.9, color="black")) + 
	theme(axis.ticks = element_line(colour = 'black')) + 
	   	scale_x_discrete(name = "injection") + 
  scale_y_continuous(name = "RA (%)", breaks=c(-3:2), labels=10^c(-3:2), limits = c(-3,2)) #+ 
 dev.off()
 	
 		
#dev.off()  
gcReps["HOST", grep("11dpi",colnames(gcReps))]



gcReps_load=gcReps

divide_by_host=function(tab){
	tab["HOST", which(tab["HOST",]==0)]<-0.00000001
	for(c in 1:ncol(tab)){
		tab[,c]/tab["HOST",c]->tab[,c]
		if( max(tab[,c])>100 ){
			tab[,c]=200*tab[,c]/sum(tab[,c])
		}
	}
	#tab[which(tab>100)]<-100
	return(tab)
}

#skip this step for a relative abundance comparison	
gcReps_load=divide_by_host(gcReps_load)

gcReps_load=rbind(gcReps_load, CFU_gc, qPCR_gc)


melt(gcReps_load)->histDL
c("organism", "sample_name", "abundance")->names(histDL)
gsub(".*xgc", "gc", histDL$sample_name)->histDL$sample_category
	gsub("\\..*", "", histDL$sample_category)->histDL$sample_category
	histDL$sample_category=factor(histDL$sample_category, levels=gcLevels)
gsub(".*\\.", "", histDL$sample_name)->histDL$replicate
paste(histDL$replicate, histDL$organism, sep="_")->histDL$org_rep
histDL[histDL$organism!="HOST",]->histDL
#histDL[histDL$organism!="CFU_gc",]->histDL
histDL$organism=factor(histDL$organism, levels=unique(as.character(histDL$organism)))
	levels(histDL$organism)->taxa_color_pairs
	COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
	#ALPHA_COLUMN_IN_COLORSLIST=7
	source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
	"blue"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Xe8510"),2]
	"lightblue"->taxa_color_pairs[which(taxa_color_pairs[,1]=="CFU_gc"),2]
	"darkgreen"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Nocardioidaceae"),2]
	"#7FC97F"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Actinobacteria"),2]
	"#BEAED4"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Gammaproteobacteria"),2]
	"#FDC086"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Betaproteobacteria"),2]
	"#FFFF99"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Alphaproteobacteria"),2]
	taxa_color_pairs[match(levels(histDL$organism), taxa_color_pairs[,1]),2]->histDL$colors


#boxplot just 16S data
#histDL2=histDL[grep("Streptomycetaceae|Moraxellaceae|Oxalobacteraceae|Geodermatophilaceae|Micrococcaceae|Sphingomonadaceae|Pseudomonadaceae|Xe8510", histDL$organism, invert=FALSE),]
histDL2=histDL[grep("Actinobacteria|Gammaproteobacteria|Betaproteobacteria|Alphaproteobacteria|Xe8510", histDL$organism, invert=FALSE),]
histDL2$organism=factor(histDL2$organism, levels=unique(as.character(histDL2$organism)))
log10(histDL2$abundance+.00000000000001)->histDL2$abundance
date=format(Sys.Date(), format="%Y%m%d")
pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/PepperGCBact_", date, ".pdf", sep="", collapse=""), width = 7, height = 4, useDingbats=FALSE)
ggplot(histDL2, aes(fill=organism, y=abundance, x=sample_category)) + 
   	 	geom_boxplot(outlier.size=0, outlier.shape=NA, alpha=0.5, width=.7, position=position_dodge(.9)) + 
   	 	geom_point(aes(fill=organism), size = 1.5, stroke=0.1, shape = 21, position = position_jitterdodge(jitter.width=0.1, dodge=.9)) +
   	 	theme_classic() + scale_color_manual(values=histDL2$colors) + scale_fill_manual(values=histDL2$colors) + 
    	theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=0.99, vjust=0.9, color="black")) + 
    	theme(axis.ticks = element_line(colour = 'black')) +
    	scale_x_discrete(name = "injection") + 
 		scale_y_continuous(name = "log10 abundance", breaks=c(-3:3), labels=10^c(-3:3), limits = c(-3,3)) #+ 
 		#geom_hline(yintercept=log10(30), linetype="dashed", 
         #       color = "red", size=2)
dev.off()


#boxplot ONLY CFU
histDL3=histDL[grep("CFU_gc|qPCR_gc|Xe8510", histDL$organism),]
histDL3$organism=factor(histDL3$organism, levels=as.character(c("CFU_gc", "qPCR_gc", "Xe8510")))
#################################

hamPCR_median=median(histDL3$abundance[which(histDL3$organism=="Xe8510" & histDL3$sample_category=="gc4dpi")])
qPCR_median=median(histDL3$abundance[which(histDL3$organism=="qPCR_gc" & histDL3$sample_category=="gc4dpi")])
CFU_median=median(histDL3$abundance[which(histDL3$organism=="CFU_gc" & histDL3$sample_category=="gc4dpi")])

qPCR_scale_factor=CFU_median/qPCR_median
hamPCR_scale_factor=CFU_median/hamPCR_median
histDL3$abundance[which(histDL3$organism=="Xe8510")]=hamPCR_scale_factor*histDL3$abundance[which(histDL3$organism=="Xe8510")]
histDL3$abundance[which(histDL3$organism=="qPCR_gc")]=qPCR_scale_factor*histDL3$abundance[which(histDL3$organism=="qPCR_gc")]

log10(histDL3$abundance)->histDL3$abundance

#ALIGN hamPCR to the QPCR.
hamPCR_values=histDL3$abundance[which(histDL3$organism=="Xe8510")]
qPCR_values=histDL3$abundance[which(histDL3$organism=="qPCR_gc")]
CFU_values=histDL3$abundance[which(histDL3$organism=="CFU_gc")]
ham_v_qpcr=lm(qPCR_values~hamPCR_values)[[1]]
slope=ham_v_qpcr[2]
intercept=ham_v_qpcr[1]
#abline(intercept, slope)
#par(mfrow=c(1,2))
#plot(hamPCR_values, qPCR_values)
#plot(slope*hamPCR_values + intercept, qPCR_values)
histDL3$abundance[which(histDL3$organism=="Xe8510")]<-(slope*histDL3$abundance[which(histDL3$organism=="Xe8510")]) + intercept


#ALIGN hamPCR and QPCR to the CFU, but only the middle concentrations
middle_range=which(match(histDL3$sample_category, c("gc0dpi","gc2dpi","gc4dpi", "gc7dpi", "gc9dpi", "gc11dpi"), nomatch=0)>0)
hamPCR_values=histDL3$abundance[intersect(which(histDL3$organism=="Xe8510"), middle_range)]
qPCR_values=histDL3$abundance[intersect(which(histDL3$organism=="qPCR_gc"), middle_range)]
CFU_values=histDL3$abundance[intersect(which(histDL3$organism=="CFU_gc"), middle_range)]
ham_v_CFU=lm(CFU_values~hamPCR_values)[[1]]
slope=ham_v_CFU[2]
intercept=ham_v_CFU[1]
#abline(intercept, slope)
#par(mfrow=c(1,2))
#plot(hamPCR_values, CFU_values)
#plot(slope*hamPCR_values + intercept, qPCR_values)
histDL3$abundance[which(histDL3$organism=="Xe8510")]<-(slope*histDL3$abundance[which(histDL3$organism=="Xe8510")]) + intercept
histDL3$abundance[which(histDL3$organism=="qPCR_gc")]<-(slope*histDL3$abundance[which(histDL3$organism=="qPCR_gc")]) + intercept
#####################################

date=format(Sys.Date(), format="%Y%m%d")
pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/PepperGC_CFUcomp_", date, ".pdf", sep="", collapse=""), width = 4, height = 4, useDingbats=FALSE)
ggplot(histDL3, aes(fill=organism, y=abundance, x=sample_category)) + 
   	 	geom_boxplot(outlier.size=0, outlier.shape=NA, alpha=0.75, width=.7, position=position_dodge(.9)) + 
   	 	geom_point(aes(fill=organism), size = 1.0, stroke=0.1, shape = 21, position = position_jitterdodge(jitter.width=0.1, dodge=.9)) +
   	 	theme_classic() + scale_color_manual(values=histDL3$colors) + scale_fill_manual(values=c("lightblue", "orange", "blue")) + 
    	theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=0.99, vjust=0.9, color="black")) + 
    	theme(axis.ticks = element_line(colour = 'black')) +
    	scale_x_discrete(name = "injection") + 
 		scale_y_continuous(name = "log10 abundance", breaks=c(0:9), labels=c(0:9), limits = c(0,9)) 
dev.off()


X=histDL3$abundance[which(histDL3$organism=="CFU_gc")]
Y=histDL3$abundance[which(histDL3$organism=="Xe8510")]

dpi11=as.vector(gcReps_load[grep("Actinobacteria|Gammaproteobacteria|Betaproteobacteria|Alphaproteobacteria", rownames(gcReps_load)),grep(".*11dpi.*",colnames(gcReps_load))])
dpi9=as.vector(gcReps_load[grep("Actinobacteria|Gammaproteobacteria|Betaproteobacteria|Alphaproteobacteria", rownames(gcReps_load)),grep(".*9dpi.*",colnames(gcReps_load))])
dpi7=as.vector(gcReps_load[grep("Actinobacteria|Gammaproteobacteria|Betaproteobacteria|Alphaproteobacteria", rownames(gcReps_load)),grep(".*7dpi.*",colnames(gcReps_load))])
dpi4=as.vector(gcReps_load[grep("Actinobacteria|Gammaproteobacteria|Betaproteobacteria|Alphaproteobacteria", rownames(gcReps_load)),grep(".*4dpi.*",colnames(gcReps_load))])
dpi2=as.vector(gcReps_load[grep("Actinobacteria|Gammaproteobacteria|Betaproteobacteria|Alphaproteobacteria", rownames(gcReps_load)),grep(".*2dpi.*",colnames(gcReps_load))])
dpi0=as.vector(gcReps_load[grep("Actinobacteria|Gammaproteobacteria|Betaproteobacteria|Alphaproteobacteria", rownames(gcReps_load)),grep(".*0dpi.*",colnames(gcReps_load))])

wilcox.test(dpi7, dpi0, "greater")

all=melt(rbind(dpi0, dpi2, dpi4))
names(all)=c("category", "b", "abundance")

pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/PepperGC_increaseCommensal_", date, ".pdf", sep="", collapse=""), width = 3, height = 4, useDingbats=FALSE)
ggplot(all, aes(y=log10(abundance), x=category)) + 
   	 	geom_boxplot(outlier.size=0, outlier.shape=NA, alpha=0.5, width=.7, position=position_dodge(.9)) + 
   	 	geom_point(aes(fill=category), size = 2.5, color="black", alpha=0.3, stroke=0.1, shape = 21, position = position_jitterdodge(jitter.width=0.9, dodge=.9)) +
   	 	theme_classic() + scale_color_manual(values=rep("black", 5)) + scale_fill_manual(values=rep("black", 5)) + 
    	theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=0.99, vjust=0.9, color="black")) + 
    	theme(axis.ticks = element_line(colour = 'black')) +
    	scale_x_discrete(name = "injection") + 
  		scale_y_continuous(name = "log10 abundance", breaks=c(-4:0), labels=c(-4:0), limits = c(-4,.5)) 
dev.off()




#>
	#>
		#>
			#> #worm titration LANE 4
		#>
	#>
#>  	

WormRep1=grep("S49xV5V6V7csq1xWorm_0D7C|S50xV5V6V7csq1xWorm_1D6C|S51xV5V6V7csq1xWorm_2D5C|S52xV5V6V7csq1xWorm_3D4C|S53xV5V6V7csq1xWorm_4D3C|S54xV5V6V7csq1xWorm_5D2C|S55xV5V6V7csq1xWorm_6D1C|S56xV5V6V7csq1xWorm_7D0C", colnames(otureads), value=TRUE)
WormRep2=grep("S57xV5V6V7csq1xWorm_0D7C|S58xV5V6V7csq1xWorm_1D6C|S59xV5V6V7csq1xWorm_2D5C|S60xV5V6V7csq1xWorm_3D4C|S61xV5V6V7csq1xWorm_4D3C|S62xV5V6V7csq1xWorm_5D2C|S63xV5V6V7csq1xWorm_6D1C|S64xV5V6V7csq1xWorm_7D0C" , colnames(otureads), value=TRUE)
WormRep3=grep("S65xV5V6V7csq1xWorm_0D7C|S66xV5V6V7csq1xWorm_1D6C|S67xV5V6V7csq1xWorm_2D5C|S68xV5V6V7csq1xWorm_3D4C|S69xV5V6V7csq1xWorm_4D3C|S70xV5V6V7csq1xWorm_5D2C|S71xV5V6V7csq1xWorm_6D1C|S72xV5V6V7csq1xWorm_7D0C", colnames(otureads), value=TRUE)

WormTit=otureads[,c(WormRep1, WormRep2, WormRep3)]

WormTit=normalize100(WormTit)

if(LANE=="L4"){
	HOST=WormTit[c("Otu2"),]
	WormTit=WormTit[setdiff(rownames(WormTit), c("Otu2")),]
	Lysinibacillus=WormTit[c("Otu18"),]
	WormTit=WormTit[setdiff(rownames(WormTit), c("Otu18")),]
	Pseudomonas=WormTit[c("Otu45"),]
	WormTit=WormTit[setdiff(rownames(WormTit), c("Otu45")),]
	Ecoli=WormTit[c("Otu7"),]
	WormTit=WormTit[setdiff(rownames(WormTit), c("Otu7")),]
}
if(LANE=="L4_REDO"){
	HOST=WormTit[c("Otu2"),]
	WormTit=WormTit[setdiff(rownames(WormTit), c("Otu2")),]
	Lysinibacillus=WormTit[c("Otu16"),]
	WormTit=WormTit[setdiff(rownames(WormTit), c("Otu16")),]
	Pseudomonas=WormTit[c("Otu44"),]
	WormTit=WormTit[setdiff(rownames(WormTit), c("Otu44")),]
	Ecoli=WormTit[c("Otu6"),]
	WormTit=WormTit[setdiff(rownames(WormTit), c("Otu6")),]
}


WormTit=rbind(WormTit, Lysinibacillus, Pseudomonas, Ecoli)
SUMS=colSums(WormTit)
WormTit=rbind(WormTit, SUMS, HOST)


WormTit=WormTit[which(rowSums(WormTit)>0),]
WormTit=WormTit[apply(WormTit, 1, max)>.2,]
sort(rowSums(WormTit))


melt(WormTit)->histDL
c("organism", "sample_name", "abundance")->names(histDL)
gsub(".*csq1x", "", histDL$sample_name)->histDL$sample_category
rep(0, nrow(histDL))->histDL$replicate
	histDL$replicate[grep(paste("S",49:56, "x", sep="",collapse="|"), histDL$sample_name)]=1
	histDL$replicate[grep(paste("S",57:64, "x", sep="",collapse="|"), histDL$sample_name)]=2
	histDL$replicate[grep(paste("S",65:72, "x", sep="",collapse="|"), histDL$sample_name)]=3
paste(histDL$replicate, histDL$organism, sep="_")->histDL$org_rep
histDL$sample_category=as.matrix(histDL$sample_category)
paste(histDL$sample_category, histDL$replicate, sep="_")->histDL$sample_rep
histDL$org_rep=factor(histDL$org_rep, levels=unique(histDL$org_rep))
histDL$sample_category=factor(histDL$sample_category, levels=c("Worm_0D7C","Worm_1D6C","Worm_2D5C","Worm_3D4C","Worm_4D3C","Worm_5D2C","Worm_6D1C","Worm_7D0C"))
histDL[histDL$organism!="SUMS",]->histDL


	levels(histDL$organism)->taxa_color_pairs
	COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
	source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
	"green"->taxa_color_pairs[which(taxa_color_pairs[,1]=="HOST"),2]
	"#fdba17"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Lysinibacillus"),2]
	"#3e78bd"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Pseudomonas"),2]
	"#e9427e"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Ecoli"),2]
	
	taxa_color_pairs[match(histDL$organism, taxa_color_pairs[,1]),2]->histDL$colors
	
#relative abundance
histDL2=histDL[grep("Lysinibacillus|Pseudomonas|Ecoli|HOST|SUMS", histDL$organism),]
histDL2$organism=factor(histDL2$organism, levels=unique(as.character(histDL2$organism)))
date=format(Sys.Date(), format="%Y%m%d")
pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/WormTit_", "ReAb", "_", date, ".pdf", sep="", collapse=""), width = 2.5, height = 4, useDingbats=FALSE)
ggplot(histDL2, aes(x=sample_category, y=abundance, group=org_rep, color=organism)) +
  geom_line() + geom_point(size=2) + scale_color_manual(values=histDL2$colors) + 
  theme_minimal() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
  theme(legend.position="none") +
  theme(panel.grid.minor = element_line(size = 0.15), panel.grid.major = element_line(size = .15)) + 
			theme(panel.grid.minor = element_line(colour = "black"), panel.grid.major = element_line(colour = "black")) + 
			theme(panel.grid.minor = element_blank()) +
			scale_y_continuous(limits=c(0, 100))	
dev.off()  


histDL3=histDL2[grep("Lysinibacillus|Pseudomonas|Ecoli", histDL2$organism),]
histDL3$organism=factor(histDL3$organism, levels=unique(as.character(histDL3$organism)))
pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/WormTit_", "ReAbHist", "_", date, ".pdf", sep="", collapse=""), width = 2.5, height = 4, useDingbats=FALSE)
ggplot(histDL3,aes(x=sample_rep,y=abundance, fill=organism)) + 
 	 geom_bar(position="fill", stat="identity", width=1.0) + 
 	 scale_fill_manual(values=histDL3$colors) + theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
 	 xlab("") + ylab("") + guides(fill=FALSE) +
 	 theme_classic() + theme(legend.position="none") + 
		scale_alpha_continuous(range = c(0.3, 1)) #+
		#scale_y_continuous(limits=c(0, Ylim))	
dev.off()




WormTit_load=WormTit
#Yaxis="linear"
Yaxis="sqrt"

divide_by_host=function(tab){
	tab["HOST", which(tab["HOST",]==0)]<-0.00000001
	for(c in 1:ncol(tab)){
		tab[,c]/tab["HOST",c]->tab[,c]
		if( max(tab[,c])>100 ){
			tab[,c]=tab[,c]/sum(tab[,c])
		}
	}
	#tab[which(tab>100)]<-100
	return(tab)
}
	
WormTit_load=divide_by_host(WormTit_load)
#WormTit_load=add.empty.taxa(WormTit_load)

colSums(WormTit_load)



melt(WormTit_load)->histDL
c("organism", "sample_name", "abundance")->names(histDL)
gsub(".*csq1x", "", histDL$sample_name)->histDL$sample_category
rep(0, nrow(histDL))->histDL$replicate
	histDL$replicate[grep(paste("S",49:56, "x", sep="",collapse="|"), histDL$sample_name)]=1
	histDL$replicate[grep(paste("S",57:64, "x", sep="",collapse="|"), histDL$sample_name)]=2
	histDL$replicate[grep(paste("S",65:72, "x", sep="",collapse="|"), histDL$sample_name)]=3
paste(histDL$replicate, histDL$organism, sep="_")->histDL$org_rep
histDL$sample_category=as.matrix(histDL$sample_category)
paste(histDL$sample_category, histDL$replicate, sep="_")->histDL$sample_rep
histDL$org_rep=factor(histDL$org_rep, levels=unique(histDL$org_rep))
histDL$sample_category=factor(histDL$sample_category, levels=c("Worm_0D7C","Worm_1D6C","Worm_2D5C","Worm_3D4C","Worm_4D3C","Worm_5D2C","Worm_6D1C","Worm_7D0C"))


histDL[histDL$organism!="HOST",]->histDL
histDL$organism=factor(histDL$organism, levels=as.character(unique(histDL$organism)))

if(Yaxis=="sqrt"){ 
	print("HI")
	sqrt(sqrt(histDL$abundance))->histDL$abundance
	Yaxmin=0
	Yaxmax=sqrt(sqrt(.3))
	breaks=c(0, .2, .4, .6, .8)
	breaks=sqrt(sqrt(c(0, .1, .2, .3)))
	breaklabels=c(0, .2, .4, .6, .8)
	breaklabels=c(0, .1, .2, .3)
	
}
if(Yaxis=="linear"){ 
	Yaxmin=0
	Yaxmax=0.3
	breaks=c(0, .1, .2, .3)
	breaklabels=c(0, .1, .2, .3)
}

	levels(histDL$organism)->taxa_color_pairs
	COLUMN_IN_COLORSLIST=10    # 3 is colors, 4 is black / white / green / purple
	source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
	"black"->taxa_color_pairs[which(taxa_color_pairs[,1]=="SUMS"),2]
	"#fdba17"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Lysinibacillus"),2]
	"#3e78bd"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Pseudomonas"),2]
	"#e9427e"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Ecoli"),2]

	
	
	taxa_color_pairs[match(levels(histDL$organism), taxa_color_pairs[,1]),2]->histDL$colors
	
histDL2=histDL[grep("Lysinibacillus|Pseudomonas|Ecoli|SUMS", histDL$organism),]
histDL2$organism=factor(histDL2$organism, levels=unique(as.character(histDL2$organism)))
#histDL2=histDL
date=format(Sys.Date(), format="%Y%m%d")
pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/WormTit_", Yaxis, "_", date, ".pdf", sep="", collapse=""), width = 2.5, height = 4, useDingbats=FALSE)
ggplot(histDL2, aes(x=sample_category, y=abundance, group=org_rep, color=organism)) +
  geom_line() + geom_point(size=2) + scale_color_manual(values=histDL2$colors) + 
  theme_minimal() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
   theme(legend.position="none") +
  theme(panel.grid.minor = element_line(size = 0.15), panel.grid.major = element_line(size = .15)) + 
			theme(panel.grid.minor = element_line(colour = "black"), panel.grid.major = element_line(colour = "black")) + 
			theme(panel.grid.minor = element_blank()) +
	scale_y_continuous(name = "abundance", breaks = breaks, labels = breaklabels, limits = c(Yaxmin,Yaxmax))

dev.off()

histDL3=histDL2[grep("Lysinibacillus|Pseudomonas|Ecoli", histDL2$organism),]
histDL3$organism=factor(histDL3$organism, levels=unique(as.character(histDL3$organism)))
pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/WormTit_", "LOADHist", "_", date, ".pdf", sep="", collapse=""), width = 2.5, height = 4, useDingbats=FALSE)
ggplot(histDL3,aes(x=sample_rep,y=abundance, fill=organism)) + 
 	 geom_bar(position="stack", stat="identity", width=1.0) + 
 	 scale_fill_manual(values=histDL3$colors) + theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
 	 xlab("") + ylab("") + guides(fill=FALSE) +
 	 theme_classic() + theme(legend.position="none") + 
		scale_y_continuous(name = "abundance", breaks = breaks, labels = breaklabels, limits = c(Yaxmin,Yaxmax))
dev.off()


#>
	#>
		#>
			#> GAUTAM TITRATION LANE 4
		#>
	#>
#>  	

ActinRep1=grep("S25xV4GI502actinxHPA_0D7C|S26xV4GI502actinxHPA_1D6C|S27xV4GI502actinxHPA_2D5C|S28xV4GI502actinxHPA_3D4C|S29xV4GI502actinxHPA_4D3C|S30xV4GI502actinxHPA_5D2C|S31xV4GI502actinxHPA_6D1C|S32xV4GI502actinxHPA_7D0C", colnames(otureads), value=TRUE)
ActinRep2=grep("S33xV4GI502actinxHPA_0D7C|S34xV4GI502actinxHPA_1D6C|S35xV4GI502actinxHPA_2D5C|S36xV4GI502actinxHPA_3D4C|S37xV4GI502actinxHPA_4D3C|S38xV4GI502actinxHPA_5D2C|S39xV4GI502actinxHPA_6D1C|S40xV4GI502actinxHPA_7D0C", colnames(otureads), value=TRUE)
ActinRep3=grep("S41xV4GI502actinxHPA_0D7C|S42xV4GI502actinxHPA_1D6C|S43xV4GI502actinxHPA_2D5C|S44xV4GI502actinxHPA_3D4C|S45xV4GI502actinxHPA_4D3C|S46xV4GI502actinxHPA_5D2C|S47xV4GI502actinxHPA_6D1C|S48xV4GI502actinxHPA_7D0C", colnames(otureads), value=TRUE)

ActinTit=otureads[,c(ActinRep1, ActinRep2, ActinRep3)]

ActinTit=normalize100(ActinTit)

if(LANE=="L4"){
	HOST=ActinTit[c("Otu5"),]
	ActinTit=ActinTit[setdiff(rownames(ActinTit), c("Otu5")),]
	Hpa=ActinTit[c("Otu10"),]
	ActinTit=ActinTit[setdiff(rownames(ActinTit), c("Otu10")),]
	DC3000=colSums(ActinTit[c("Otu11", "Otu12"),])
	ActinTit=ActinTit[setdiff(rownames(ActinTit), c("Otu11", "Otu12")),]
}
if(LANE=="L4_REDO"){
	HOST=ActinTit[c("Otu5"),]
	ActinTit=ActinTit[setdiff(rownames(ActinTit), c("Otu5")),]
	Hpa=ActinTit[c("Otu9"),]
	ActinTit=ActinTit[setdiff(rownames(ActinTit), c("Otu9")),]
	SUMS=colSums(ActinTit)
	DC3000=colSums(ActinTit[c("Otu24", "Otu10"),])
	ActinTit=ActinTit[setdiff(rownames(ActinTit), c("Otu24", "Otu10")),]
}

#max = 66, min = 6

#convert to family, add host afterwards
ActinTit<-filter_by_taxonomyic_level(countTable=ActinTit, taxonomy=taxonomy, keywords=c("f"), last_filter=TRUE)
	SUMS=colSums(ActinTit)

ActinTit=ActinTit[which(rowSums(ActinTit)>0),]
ActinTit=ActinTit[apply(ActinTit, 1, min)>.02,]

#add back host and HPA
ActinTit=rbind(ActinTit, HOST, DC3000, Hpa, SUMS)

#ActinTit=add.empty.taxa(ActinTit)


order(rowSums(ActinTit), decreasing=FALSE)->ordering_vector
		ActinTit[ordering_vector,]->ActinTit

#next set specific order
	toporder=c("empty", "low_abundance")
	bottomorder=c("Caulobacteraceae", "Methylobacteriaceae", "Microbacteriaceae", "Comamonadaceae", "Flavobacteriaceae", "Oxalobacteraceae", "Sphingobacteriaceae", "Cytophagaceae", "Sphingomonadaceae", "Pseudomonadaceae", "HOST","DC3000","Hpa")		
	ActinTit=topOrder(ActinTit, toporder, bottomorder)


melt(ActinTit)->histDL
c("organism", "sample_name", "abundance")->names(histDL)
gsub(".*actinx", "", histDL$sample_name)->histDL$sample_category
rep(0, nrow(histDL))->histDL$replicate
	histDL$replicate[grep(paste("S",25:32, "x", sep="",collapse="|"), histDL$sample_name)]=1
	histDL$replicate[grep(paste("S",33:40, "x", sep="",collapse="|"), histDL$sample_name)]=2
	histDL$replicate[grep(paste("S",41:48, "x", sep="",collapse="|"), histDL$sample_name)]=3
paste(histDL$replicate, histDL$organism, sep="_")->histDL$org_rep
histDL$sample_category=as.matrix(histDL$sample_category)
paste(histDL$sample_category, histDL$replicate, sep="_")->histDL$sample_rep
histDL$org_rep=factor(histDL$org_rep, levels=unique(histDL$org_rep))
histDL$sample_category=factor(histDL$sample_category, levels=c("HPA_0D7C","HPA_1D6C","HPA_2D5C","HPA_3D4C","HPA_4D3C","HPA_5D2C","HPA_6D1C","HPA_7D0C"))
histDL$organism=factor(histDL$organism, levels=as.character(unique(histDL$organism)))


	levels(histDL$organism)->taxa_color_pairs
	COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
	#ALPHA_COLUMN_IN_COLORSLIST=7
	source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
	"yellow"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Hpa"),2]
	"#5cffc0"->taxa_color_pairs[which(taxa_color_pairs[,1]=="DC3000"),2]
	"green"->taxa_color_pairs[which(taxa_color_pairs[,1]=="HOST"),2]
	"black"->taxa_color_pairs[which(taxa_color_pairs[,1]=="SUMS"),2]
	"grey"->taxa_color_pairs[which(taxa_color_pairs[,1]=="low_abundance"),2]
	"#FFBC00"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Caulobacteraceae"),2]
	"#f43df0"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Flavobacteriaceae"),2]
	"#c7f214"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Sphingobacteriaceae"),2]
	taxa_color_pairs[match(levels(histDL$organism), taxa_color_pairs[,1]),2]->histDL$colors
	
#relative abundance
date=format(Sys.Date(), format="%Y%m%d")
histDL=histDL[grep("SUMS", histDL$organism, invert=TRUE),]
histDL$organism=factor(histDL$organism, levels=unique(as.character(histDL$organism)))
pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/GautamTit_", "ReAb", "_", date, ".pdf", sep="", collapse=""), width = 2.5, height = 4, useDingbats=FALSE)
ggplot(histDL, aes(x=histDL$sample_category, y=histDL$abundance, group=histDL$org_rep, color=histDL$organism)) +
  geom_line() + geom_point(size=2) + scale_color_manual(values=histDL$colors) + 
  theme_minimal() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
  theme(legend.position="none") +
  theme(panel.grid.minor = element_line(size = 0.15), panel.grid.major = element_line(size = .15)) + 
			theme(panel.grid.minor = element_line(colour = "black"), panel.grid.major = element_line(colour = "black")) + 
			theme(panel.grid.minor = element_blank()) +
				scale_y_continuous(name = "abundance", breaks = c(0, 35, 70), labels = c(0, 35, 70), limits = c(0,70))
dev.off()  


histDL3=histDL[grep("HOST|Hpa|SUMS", histDL$organism, invert=TRUE),]
histDL3$organism=factor(histDL3$organism, levels=unique(as.character(histDL3$organism)))
#pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/GautamTit_", "ReAbHist", "_", date, ".pdf", sep="", collapse=""), width = 2.5, height = 4, useDingbats=FALSE)
ggplot(histDL3,aes(x=sample_rep,y=abundance, fill=organism)) + 
 	 geom_bar(position="fill", stat="identity", width=1.0) + 
 	 scale_fill_manual(values=histDL3$colors) + theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
 	 xlab("") + ylab("") + guides(fill=FALSE) +
 	 theme_classic() + theme(legend.position="none") + 
		scale_alpha_continuous(range = c(0.3, 1)) +
		scale_y_continuous(name = "abundance", breaks = c(0, .5, 1), labels = c(0, 50, 100))
#dev.off()


ActinTit_load=ActinTit

Yaxis="linear"
#Yaxis="sqrt"

divide_by_host=function(tab){
	tab["HOST", which(tab["HOST",]==0)]<-0.00000001
	for(c in 1:ncol(tab)){
		tab[,c]/tab["HOST",c]->tab[,c]
		if( max(tab[,c])>100 ){
			tab[,c]=tab[,c]/sum(tab[,c])
		}
	}
	#tab[which(tab>100)]<-100
	return(tab)
}
	
ActinTit_load=divide_by_host(ActinTit_load)
#ActinTit_load=add.empty.taxa(ActinTit_load)

toporder=c("empty", "low_abundance")
	bottomorder=c("Caulobacteraceae", "Methylobacteriaceae", "Microbacteriaceae", "Comamonadaceae", "Flavobacteriaceae", "Oxalobacteraceae", "Sphingobacteriaceae", "Cytophagaceae", "Sphingomonadaceae", "Pseudomonadaceae", "HOST","DC3000","Hpa")		
	ActinTit_load=topOrder(ActinTit_load, toporder, bottomorder)


melt(ActinTit_load)->histDL
c("organism", "sample_name", "abundance")->names(histDL)
gsub(".*actinx", "", histDL$sample_name)->histDL$sample_category
rep(0, nrow(histDL))->histDL$replicate
	histDL$replicate[grep(paste("S",25:32, "x", sep="",collapse="|"), histDL$sample_name)]=1
	histDL$replicate[grep(paste("S",33:40, "x", sep="",collapse="|"), histDL$sample_name)]=2
	histDL$replicate[grep(paste("S",41:48, "x", sep="",collapse="|"), histDL$sample_name)]=3
paste(histDL$replicate, histDL$organism, sep="_")->histDL$org_rep
histDL$sample_category=as.matrix(histDL$sample_category)
paste(histDL$sample_category, histDL$replicate, sep="_")->histDL$sample_rep
histDL$org_rep=factor(histDL$org_rep, levels=unique(histDL$org_rep))
histDL$sample_category=factor(histDL$sample_category, levels=c("HPA_0D7C","HPA_1D6C","HPA_2D5C","HPA_3D4C","HPA_4D3C","HPA_5D2C","HPA_6D1C","HPA_7D0C"))
histDL$organism=factor(histDL$organism, levels=as.character(unique(histDL$organism)))


if(Yaxis=="sqrt"){ 
	print("HI")
	sqrt(sqrt(histDL$abundance))->histDL$abundance
	Yaxmin=0
	Yaxmax=sqrt(sqrt(7))
	#breaks=c(0, .8, 1.6)
	#breaklabels=c(0, .8, 1.6)
	breaks=sqrt(sqrt(c(0, 3.5, 7)))
	breaklabels=c(0, 3.5, 7)
	
}
if(Yaxis=="linear"){ 
	Yaxmin=0
	Yaxmax=7
	breaks=c(0, 3.5, 7)
	breaklabels=c(0, 3.5, 7)
}

 
	levels(histDL$organism)->taxa_color_pairs
	COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
	#ALPHA_COLUMN_IN_COLORSLIST=7
	source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
	"yellow"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Hpa"),2]
	"#5cffc0"->taxa_color_pairs[which(taxa_color_pairs[,1]=="DC3000"),2]
	"green"->taxa_color_pairs[which(taxa_color_pairs[,1]=="HOST"),2]
	"black"->taxa_color_pairs[which(taxa_color_pairs[,1]=="SUMS"),2]
	"grey"->taxa_color_pairs[which(taxa_color_pairs[,1]=="low_abundance"),2]
	"#FFBC00"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Caulobacteraceae"),2]
	"#f43df0"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Flavobacteriaceae"),2]
	"#c7f214"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Sphingobacteriaceae"),2]
	taxa_color_pairs[match(levels(histDL$organism), taxa_color_pairs[,1]),2]->histDL$colors
	
	
	
histDL2=histDL[grep("HOST", histDL$organism, invert=TRUE),]
histDL2$organism=factor(histDL2$organism, levels=unique(as.character(histDL2$organism)))
date=format(Sys.Date(), format="%Y%m%d")
#pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/GautamTit_", Yaxis, "_", date, ".pdf", sep="", collapse=""), width = 2.5, height = 4, useDingbats=FALSE)
ggplot(histDL2, aes(x=sample_category, y=abundance, group=org_rep, color=organism)) +
  geom_line() + geom_point(size=2) + scale_color_manual(values=histDL2$colors) + 
  theme_minimal() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
  theme(legend.position="none") +
  theme(panel.grid.minor = element_line(size = 0.15), panel.grid.major = element_line(size = .15)) + 
			theme(panel.grid.minor = element_line(colour = "black"), panel.grid.major = element_line(colour = "black")) + 
			theme(panel.grid.minor = element_blank()) +
	scale_y_continuous(name = "abundance", breaks = breaks, labels = breaklabels , limits = c(Yaxmin,Yaxmax)) 	
dev.off()


histDL3=histDL[grep("HOST|SUMS", histDL$organism, invert=TRUE),]
histDL3$organism=factor(histDL3$organism, levels=unique(as.character(histDL3$organism)))
#pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/GautamTit_", "LOADHist", "_", date, ".pdf", sep="", collapse=""), width = 2.5, height = 4, useDingbats=FALSE)
ggplot(histDL3,aes(x=sample_rep,y=abundance, fill=organism)) + 
 	 geom_bar(position="stack", stat="identity", width=1.0) + 
 	 scale_fill_manual(values=histDL3$colors) + theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
 	 xlab("") + ylab("") + guides(fill=FALSE) +
 	 theme_classic() + theme(legend.position="none") + 
		scale_alpha_continuous(range = c(0.3, 1)) +
		scale_y_continuous(name = "abundance", breaks = c(0, 3.5, 7), labels = c(0, 3.5, 7))
dev.off()

#>
	#>
		#>
			#> GAUTAM EXPERIMENT LANE1
		#>
	#>
#>  	
#note uses different OTU set than the titration :(



Col_hpa1_dc=grep("S1xV4GI502actinxCol_hpa1_dc|S2xV4GI502actinxCol_hpa1_dc|S3xV4GI502actinxCol_hpa1_dc|S4xV4GI502actinxCol_hpa1_dc", colnames(otureads), value=TRUE)
eds1_hpa1_dc=grep("S5xV4GI502actinxeds1_hpa1_dc|S6xV4GI502actinxeds1_hpa1_dc|S7xV4GI502actinxeds1_hpa1_dc|S8xV4GI502actinxeds1_hpa1_dc", colnames(otureads), value=TRUE)
Col_hpa1_mg=grep("S9xV4GI502actinxCol_hpa1_mg|S10xV4GI502actinxCol_hpa1_mg|S11xV4GI502actinxCol_hpa1_mg|S12xV4GI502actinxCol_hpa1_mg", colnames(otureads), value=TRUE)
eds1_hpa1_mg=grep("S13xV4GI502actinxeds1_hpa1_mg|S14xV4GI502actinxeds1_hpa1_mg|S15xV4GI502actinxeds1_hpa1_mg|S16xV4GI502actinxeds1_hpa1_mg", colnames(otureads), value=TRUE)
Col_mg_mg=grep("S17xV4GI502actinxCol_mg_mg|S18xV4GI502actinxCol_mg_mg|S19xV4GI502actinxCol_mg_mg|S20xV4GI502actinxCol_mg_mg", colnames(otureads), value=TRUE)
eds1_mg_mg=grep("S21xV4GI502actinxeds1_mg_mg|S22xV4GI502actinxeds1_mg_mg|S23xV4GI502actinxeds1_mg_mg|S24xV4GI502actinxeds1_mg_mg", colnames(otureads), value=TRUE)
GautamLevels=c("Col_hpa1_dc","eds1_hpa1_dc","Col_hpa1_mg","eds1_hpa1_mg", "Col_mg_mg","eds1_mg_mg")
   
#Gautam=otureads[,c(Col_hpa1_dc,eds1_hpa1_dc,Col_hpa1_mg,eds1_hpa1_mg, Col_mg_mg,eds1_mg_mg)]

Gautam=otureads[,c("S1xV4GI502actinxCol_hpa1_dc","S2xV4GI502actinxCol_hpa1_dc","S3xV4GI502actinxCol_hpa1_dc","S4xV4GI502actinxCol_hpa1_dc","S5xV4GI502actinxeds1_hpa1_dc","S6xV4GI502actinxeds1_hpa1_dc","S7xV4GI502actinxeds1_hpa1_dc","S8xV4GI502actinxeds1_hpa1_dc","S9xV4GI502actinxCol_hpa1_mg","S10xV4GI502actinxCol_hpa1_mg","S11xV4GI502actinxCol_hpa1_mg","S12xV4GI502actinxCol_hpa1_mg","S13xV4GI502actinxeds1_hpa1_mg","S14xV4GI502actinxeds1_hpa1_mg","S15xV4GI502actinxeds1_hpa1_mg","S16xV4GI502actinxeds1_hpa1_mg","S17xV4GI502actinxCol_mg_mg","S18xV4GI502actinxCol_mg_mg","S19xV4GI502actinxCol_mg_mg","S20xV4GI502actinxCol_mg_mg","S21xV4GI502actinxeds1_mg_mg","S22xV4GI502actinxeds1_mg_mg","S23xV4GI502actinxeds1_mg_mg","S24xV4GI502actinxeds1_mg_mg")]

Gautam=normalize100(Gautam)
Gautam=Gautam[which(rowSums(Gautam)>0),]

if(LANE!="L1_REDO"){
	HOST=Gautam[c("Otu1"),]
	Gautam=Gautam[setdiff(rownames(Gautam), c("Otu1")),]
	Hpa=Gautam[c("Otu14"),]
	Gautam=Gautam[setdiff(rownames(Gautam), c("Otu14")),]
	DC3000=colSums(Gautam[c("Otu2", "Otu9"),])
	Gautam=Gautam[setdiff(rownames(Gautam), c("Otu2", "Otu9")),]
}
if(LANE=="L1_REDO"){
	HOST=Gautam[c("Otu1"),]
	Gautam=Gautam[setdiff(rownames(Gautam), c("Otu1")),]
	Hpa=Gautam[c("Otu19"),]
	Gautam=Gautam[setdiff(rownames(Gautam), c("Otu19")),]
	DC3000=colSums(Gautam[c("Otu2", "Otu5"),])
	Gautam=Gautam[setdiff(rownames(Gautam), c("Otu2", "Otu5")),]
}


#convert to family, add host afterwards
Gautam<-filter_by_taxonomyic_level(countTable=Gautam, taxonomy=taxonomy, keywords=c("f"), last_filter=TRUE)
Gautam=Gautam[which(rowSums(Gautam)>0),]
Gautam=low_abundance(Gautam, percent=0.05)

#add back host
Gautam=rbind(Gautam, DC3000, Hpa, HOST)

Gautam_load=Gautam


divide_by_host=function(tab){
	tab["HOST", which(tab["HOST",]==0)]<-0.00000001
	for(c in 1:ncol(tab)){
		tab[,c]/tab["HOST",c]->tab[,c]
		if( max(tab[,c])>100 ){
			tab[,c]=tab[,c]/sum(tab[,c])
		}
	}
	#tab[which(tab>100)]<-100
	return(tab)
}
	
Gautam_load=divide_by_host(Gautam_load)

#remove host amplicon
Gautam_load=Gautam_load[which(rownames(Gautam_load)!="HOST"),]

norm="ReAb"
if(norm=="LOAD"){
	Gautam_load=add.empty.taxa(Gautam_load, rescale=FALSE)
	Ylim=80
}
if(norm=="ReAb"){
	Gautam_load=Gautam
	Gautam_load=Gautam_load[which(rownames(Gautam_load)!="HOST"),]
	Gautam_load=Gautam_load[which(rownames(Gautam_load)!="Hpa"),]
	Gautam_load=normalize100(Gautam_load)
	Ylim=101
}

order_tab=Gautam_load[,grep("mg_mg",colnames(Gautam_load), invert=TRUE)]
order(rowSums(order_tab), decreasing=FALSE)->ordering_vector
		Gautam_load[ordering_vector,]->Gautam_load

#next set specific order
toporder=c("empty", "low_abundance")
	bottomorder=c("Caulobacteraceae", "Methylobacteriaceae", "Microbacteriaceae", "Comamonadaceae", "Flavobacteriaceae", "Oxalobacteraceae", "Sphingobacteriaceae", "Cytophagaceae", "Sphingomonadaceae", "Pseudomonadaceae", "HOST","DC3000","Hpa")		
	Gautam_load=topOrder(Gautam_load, toporder, bottomorder)

####### LOAD INCREASER
Gautam_load["Hpa",]=4*Gautam_load["Hpa",]


melt(Gautam_load)->histDL
c("organism", "sample_name", "abundance")->names(histDL)
gsub(".*actinx", "", histDL$sample_name)->histDL$sample_category
histDL$sample_category=factor(histDL$sample_category, levels=GautamLevels)
L=length(unique(histDL$organism))
histDL$organism=factor(histDL$organism, levels=unique(histDL$organism))
histDL$sample_name=factor(histDL$sample_name, levels=unique(histDL$sample_name))
	
	levels(histDL$organism)->taxa_color_pairs
	COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
	#ALPHA_COLUMN_IN_COLORSLIST=7
	source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
	"yellow"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Hpa"),2]
	"#5cffc0"->taxa_color_pairs[which(taxa_color_pairs[,1]=="DC3000"),2]
	"green"->taxa_color_pairs[which(taxa_color_pairs[,1]=="HOST"),2]
	"black"->taxa_color_pairs[which(taxa_color_pairs[,1]=="SUMS"),2]
	"grey"->taxa_color_pairs[which(taxa_color_pairs[,1]=="low_abundance"),2]
	"#FFBC00"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Caulobacteraceae"),2]
	"#f43df0"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Flavobacteriaceae"),2]
	"#c7f214"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Sphingobacteriaceae"),2]
	taxa_color_pairs[match(levels(histDL$organism), taxa_color_pairs[,1]),2]->histDL$colors

	
#histDL$abundance=sqrt(histDL$abundance)
date=format(Sys.Date(), format="%Y%m%d")
#pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/GautamExp_", norm, "_", date, ".pdf", sep="", collapse=""), width = 5, height = 3, useDingbats=FALSE)
	ggplot(histDL,aes(x=sample_name,y=abundance, fill=organism)) + 
 	 geom_bar(position="stack", stat="identity", width=1.0) + 
 	 scale_fill_manual(values=histDL$colors) + theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
 	 xlab("") + ylab("") + guides(fill=FALSE) +
 	 theme_classic() + theme(legend.position="none") + 
		scale_alpha_continuous(range = c(0.3, 1)) 	+
		scale_y_continuous(name = "abundance", breaks = c(0, Ylim), labels = c(0, Ylim), limits=c(0, Ylim))
#dev.off()



bcount=histDL[grep("Col_hpa1_dc|eds1_hpa1_dc", histDL$sample_category),]
bcount=bcount[grep("DC3000", bcount$organism),]
bcount=rbind(bcount, bcount, bcount)
bcount$abundance[1:8]=c(4.730171,2.063581,1.546129,1.522790,11.652902,23.703262,9.868589,6.891590)
bcount$abundance[9:16]=c(52105.26316, 15000, 4615.384615, 7500, 24705.88235, 97500, 24000, 62068.96552)
bcount$abundance[17:24]=c(39.67922,38.30553,21.02404,24.44420,53.86480,35.43408,39.20907,39.43859)
bcount$sample_category=factor(
			c(rep("Col_hpa1_dc", 4), rep("eds1_hpa1_dc", 4),
			  rep("Col_hpa1_dc_CFU", 4), rep("eds1_hpa1_dc_CFU", 4),
			  rep("Col_hpa1_dc_ReAb", 4), rep("eds1_hpa1_dc_ReAb", 4)),
			levels=c("Col_hpa1_dc_CFU","eds1_hpa1_dc_CFU","Col_hpa1_dc","eds1_hpa1_dc","Col_hpa1_dc_ReAb","eds1_hpa1_dc_ReAb"))
scale_factor1=median(bcount$abundance[9:12])/median(bcount$abundance[1:4])
bcount$abundance[1:8]=bcount$abundance[1:8]*scale_factor1
scale_factor2=median(bcount$abundance[9:12])/median(bcount$abundance[17:20])
bcount$abundance[17:24]=bcount$abundance[17:24]*scale_factor2
bcount$abundance[1:8]=bcount$abundance[1:8]/median(bcount$abundance[1:4])
bcount$abundance[9:16]=bcount$abundance[9:16]/median(bcount$abundance[9:12])
bcount$abundance[17:24]=bcount$abundance[17:24]/median(bcount$abundance[17:20])


date=format(Sys.Date(), format="%Y%m%d")
pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/GautamExp_CFUplot_", date, ".pdf", sep="", collapse=""), width = 4, height = 5, useDingbats=FALSE)
ggplot(bcount,aes(x=sample_category,y=abundance, fill=sample_category)) + 
  	 	geom_boxplot(aes(color=sample_category), outlier.size=0, outlier.shape=NA, alpha=0, width=.7, position=position_dodge(.7)) + 
   	 	geom_point(aes(fill=sample_category), size = 3, stroke=0.1, shape = 21, position = position_jitterdodge(jitter.width=0.4)) +
   	 	theme_classic() + scale_color_manual(values=rep(c("darkgrey", "black"), 3)) + scale_fill_manual(values=rep(c("darkgrey", "black"), 3)) + 
    	theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=0.99, vjust=0.9, color="black")) + 
    	theme(axis.ticks = element_line(colour = 'black')) +
    	scale_x_discrete(name = "injection") + 
 		scale_y_continuous(limits = c(0,15)) 
 dev.off()
 	
 	



bcount2=Gautam_load[grep("empty|low_abundance|Hpa", rownames(Gautam_load), invert=TRUE),]
bcount2=colSums(bcount2)
#for mg comparisons

melt(bcount2)->bcount2
rownames(bcount2)->bcount2$names
names(bcount2)[1]<-"abundance"
class(bcount2$abundance)
bcount2=bcount2[grep("Col_hpa1_mg|Col_mg_mg|eds1_hpa1_mg|eds1_mg_mg", bcount2$names),]
bcount2$sample_category=factor(c(rep("Col_hpa1_mg", 4), rep("eds1_hpa1_mg", 4), rep("Col_mg_mg", 4), rep("eds1_mg_mg", 4)), levels=c("Col_hpa1_mg", "Col_mg_mg", "eds1_hpa1_mg", "eds1_mg_mg"))
date=format(Sys.Date(), format="%Y%m%d")
pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/GautamExp_HpaMg_", date, ".pdf", sep="", collapse=""), width = 3.4, height = 5, useDingbats=FALSE)
ggplot(bcount2,aes(x=sample_category,y=abundance, fill=sample_category)) + 
  	 	geom_boxplot(aes(color=sample_category), outlier.size=0, outlier.shape=NA, alpha=0, width=.7, position=position_dodge(.7)) + 
   	 	geom_point(aes(fill=sample_category), size = 3, stroke=0.1, shape = 21, position = position_jitterdodge(jitter.width=0.4)) +
   	 	theme_classic() + scale_color_manual(values=rep(c("black", "black","black", "black"), 3)) + scale_fill_manual(values=rep(c("black", "black","black", "black"), 3)) + 
    	theme(axis.text=element_text(size=15, color="black"), axis.title=element_text(size=14,face="bold")) + 
    	theme(axis.text.x = element_text(angle=90, hjust=0.99, vjust=0.9, color="black")) + 
    	theme(axis.ticks = element_line(colour = 'black')) +
    	scale_x_discrete(name = "injection") + 
 		scale_y_continuous(limits = c(0,10)) 
dev.off()


log10(2) + log10(2)
log10(4)


#>
	#>
		#>
			#> #Athaliana metagenome (don't normalize to 100 first) lane1
		#>
	#>
#>
#first load metagenome load data
chromCounts2="/Users/dlundberg/Documents/abt6/Regalado2017/count_tables/chromCountsNew.txt"
read.countTable(chromCounts2)->chromCounts2
chromCounts2[grep("S", rownames(chromCounts2)),]->chromCounts2
chromCounts2->chromCounts

	metatab_B="/Users/dlundberg/Documents/abt6/Regalado2017/count_tables/BacteriaFamilyRaw.txt"
	read.MetacountTable(metatab_B)->metareads_B
		metareads_B=metareads_B[setdiff(1:nrow(metareads_B),match(c("Bacteria", "TotalSeq"), rownames(metareads_B))),]
		chromCounts[match(names(chromCounts), colnames(metareads_B))]->HOST
	metareads_B=rbind(HOST, metareads_B)

#80	S176
#63	S159
#71	S167
#36	S132
#83	S179
#61	S157
#19	S19
#15	S15
		
metareads_B=metareads_B[,match(c("S176","S159","S167","S132","S179","S157","S19","S15"), colnames(metareads_B))]
colnames(metareads_B)=c("S80","S63","S71","S36","S83","S61","S19","S15")



AtV3V4rep1=match(c("S25xV3V4GI502xmeta_S80", "S26xV3V4GI502xmeta_S63", "S27xV3V4GI502xmeta_S71", "S28xV3V4GI502xmeta_S36", "S29xV3V4GI502xmeta_S83", "S30xV3V4GI502xmeta_S61", "S31xV3V4GI502xmeta_S19", "S32xV3V4GI502xmeta_S15"),colnames(otureads))
AtV3V4rep2=match(c("S33xV3V4GI502xmeta_S80", "S34xV3V4GI502xmeta_S63", "S35xV3V4GI502xmeta_S71", "S36xV3V4GI502xmeta_S36", "S37xV3V4GI502xmeta_S83", "S38xV3V4GI502xmeta_S61", "S39xV3V4GI502xmeta_S19", "S40xV3V4GI502xmeta_S15"),colnames(otureads))
AtV3V4rep3=match(c("S41xV3V4GI502xmeta_S80", "S42xV3V4GI502xmeta_S63", "S43xV3V4GI502xmeta_S71", "S44xV3V4GI502xmeta_S36", "S45xV3V4GI502xmeta_S83", "S46xV3V4GI502xmeta_S61", "S47xV3V4GI502xmeta_S19", "S48xV3V4GI502xmeta_S15"),colnames(otureads))

AtV4rep1=match(c("S49xV4r806GI502xmeta_S80", "S50xV4r806GI502xmeta_S63", "S51xV4r806GI502xmeta_S71", "S52xV4r806GI502xmeta_S36", "S53xV4r806GI502xmeta_S83", "S54xV4r806GI502xmeta_S61", "S55xV4r806GI502xmeta_S19", "S56xV4r806GI502xmeta_S15"),colnames(otureads))
AtV4rep2=match(c("S57xV4r806GI502xmeta_S80", "S58xV4r806GI502xmeta_S63", "S59xV4r806GI502xmeta_S71", "S60xV4r806GI502xmeta_S36", "S61xV4r806GI502xmeta_S83", "S62xV4r806GI502xmeta_S61", "S63xV4r806GI502xmeta_S19", "S64xV4r806GI502xmeta_S15"),colnames(otureads))
AtV4rep3=match(c("S65xV4r806GI502xmeta_S80", "S66xV4r806GI502xmeta_S63", "S67xV4r806GI502xmeta_S71", "S68xV4r806GI502xmeta_S36", "S69xV4r806GI502xmeta_S83", "S70xV4r806GI502xmeta_S61", "S71xV4r806GI502xmeta_S19", "S72xV4r806GI502xmeta_S15"),colnames(otureads))

AtV5V6V7rep1=match(c("S1xV5V6V7GI466xmeta_S80", "S2xV5V6V7GI466xmeta_S63", "S3xV5V6V7GI466xmeta_S71", "S4xV5V6V7GI466xmeta_S36", "S5xV5V6V7GI466xmeta_S83", "S6xV5V6V7GI466xmeta_S61", "S7xV5V6V7GI466xmeta_S19", "S8xV5V6V7GI466xmeta_S15"),colnames(otureads))
AtV5V6V7rep2=match(c("S9xV5V6V7GI466xmeta_S80", "S10xV5V6V7GI466xmeta_S63", "S11xV5V6V7GI466xmeta_S71", "S12xV5V6V7GI466xmeta_S36", "S13xV5V6V7GI466xmeta_S83", "S14xV5V6V7GI466xmeta_S61", "S15xV5V6V7GI466xmeta_S19", "S16xV5V6V7GI466xmeta_S15"),colnames(otureads))
AtV5V6V7rep3=match(c("S17xV5V6V7GI466xmeta_S80", "S18xV5V6V7GI466xmeta_S63", "S19xV5V6V7GI466xmeta_S71", "S20xV5V6V7GI466xmeta_S36", "S21xV5V6V7GI466xmeta_S83", "S22xV5V6V7GI466xmeta_S61", "S23xV5V6V7GI466xmeta_S19", "S24xV5V6V7GI466xmeta_S15"),colnames(otureads))


dataset="V4"

if(LANE!="L1_REDO"){
if(dataset=="V3V4"){
	AtMeta=otureads[,c(AtV3V4rep1, AtV3V4rep2, AtV3V4rep3)]
	AtMeta=normalize100(AtMeta[,which(colSums(AtMeta)>=1500)])
	HOST=colSums(AtMeta[c("Otu1", "Otu14372"),])
	AtMeta=AtMeta[setdiff(rownames(AtMeta), c("Otu1", "Otu14372")),]
	}
if(dataset=="V4"){
	AtMeta=otureads[,c(AtV4rep1, AtV4rep2, AtV4rep3)]
	AtMeta=normalize100(AtMeta[,which(colSums(AtMeta)>=1500)])
	HOST=colSums(AtMeta[c("Otu1", "Otu14372"),])
	AtMeta=AtMeta[setdiff(rownames(AtMeta), c("Otu1", "Otu14372")),]
	}
if(dataset=="V5V6V7"){
	AtMeta=otureads[,c(AtV5V6V7rep1, AtV5V6V7rep2, AtV5V6V7rep3)]
	AtMeta=normalize100(AtMeta[,which(colSums(AtMeta)>=1500)])
	HOST=colSums(AtMeta[c("Otu1", "Otu14372"),])
	AtMeta=AtMeta[setdiff(rownames(AtMeta), c("Otu1", "Otu14372")),]
	}
if(dataset=="all"){
	AtMeta=otureads[,c(AtV3V4rep1, AtV3V4rep2, AtV3V4rep3, AtV4rep1, AtV4rep2, AtV4rep3, AtV5V6V7rep1, AtV5V6V7rep2, AtV5V6V7rep3)]
	AtMeta=normalize100(AtMeta[,which(colSums(AtMeta)>=1500)])
	HOST=colSums(AtMeta[c("Otu1", "Otu14372"),])
	AtMeta=AtMeta[setdiff(rownames(AtMeta), c("Otu1", "Otu14372")),]
	}
if(dataset=="meta"){
	AtMeta=metareads_B
	AtMeta=normalize100(AtMeta[,which(colSums(AtMeta)>=15000)])
	HOST=AtMeta["HOST",]
	AtMeta=AtMeta[setdiff(rownames(AtMeta), "HOST"),]
	}
}
if(LANE=="L1_REDO"){
if(dataset=="V3V4"){
	AtMeta=otureads[,c(AtV3V4rep1, AtV3V4rep2, AtV3V4rep3)]
	AtMeta=normalize100(AtMeta[,which(colSums(AtMeta)>=1500)])
	HOST=colSums(AtMeta[c("Otu1", "Otu17444"),])
	AtMeta=AtMeta[setdiff(rownames(AtMeta), c("Otu1", "Otu17444")),]
	}
if(dataset=="V4"){
	AtMeta=otureads[,c(AtV4rep1, AtV4rep2, AtV4rep3)]
	AtMeta=normalize100(AtMeta[,which(colSums(AtMeta)>=1500)])
	HOST=colSums(AtMeta[c("Otu1", "Otu17444"),])
	AtMeta=AtMeta[setdiff(rownames(AtMeta), c("Otu1", "Otu17444")),]
	}
if(dataset=="V5V6V7"){
	AtMeta=otureads[,c(AtV5V6V7rep1, AtV5V6V7rep2, AtV5V6V7rep3)]
	AtMeta=normalize100(AtMeta[,which(colSums(AtMeta)>=1500)])
	HOST=colSums(AtMeta[c("Otu1", "Otu17444"),])
	AtMeta=AtMeta[setdiff(rownames(AtMeta), c("Otu1", "Otu17444")),]
	}
if(dataset=="all"){
	AtMeta=otureads[,c(AtV3V4rep1, AtV3V4rep2, AtV3V4rep3, AtV4rep1, AtV4rep2, AtV4rep3, AtV5V6V7rep1, AtV5V6V7rep2, AtV5V6V7rep3)]
	AtMeta=normalize100(AtMeta[,which(colSums(AtMeta)>=1500)])
	HOST=colSums(AtMeta[c("Otu1", "Otu17444"),])
	AtMeta=AtMeta[setdiff(rownames(AtMeta), c("Otu1", "Otu17444")),]
	}
if(dataset=="meta"){
	AtMeta=metareads_B
	AtMeta=normalize100(AtMeta[,which(colSums(AtMeta)>=15000)])
	HOST=AtMeta["HOST",]
	AtMeta=AtMeta[setdiff(rownames(AtMeta), "HOST"),]
	}
}


if(dataset!="meta"){
	AtMeta<-filter_by_taxonomyic_level(countTable=AtMeta, taxonomy=taxonomy, keywords=c("f"), last_filter=TRUE)
}

AtMeta=AtMeta[which(rowSums(AtMeta)>0),]
#AtMeta=AtMeta[apply(AtMeta, 1, max)>0.01,]
AtMeta=low_abundance(AtMeta, percent=0.01)
	
	
#add back host
AtMeta=rbind(AtMeta, HOST)

melt(AtMeta)->histDL
c("organism", "sample_name", "abundance")->names(histDL)
gsub(".*xmeta_", "", histDL$sample_name)->histDL$sample_category
histDL$sample_category=factor(histDL$sample_category, levels=c("S80","S63","S71","S36","S83","S61","S19","S15"))
rep(0, nrow(histDL))->histDL$replicate
if(dataset=="V3V4"){
	histDL$replicate[grep(paste("S",25:32, "x", sep="",collapse="|"), histDL$sample_name)]=1
	histDL$replicate[grep(paste("S",33:40, "x", sep="",collapse="|"), histDL$sample_name)]=2
	histDL$replicate[grep(paste("S",41:48, "x", sep="",collapse="|"), histDL$sample_name)]=3
}
if(dataset=="V4"){
	histDL$replicate[grep(paste("S",49:56, "x", sep="",collapse="|"), histDL$sample_name)]=1
	histDL$replicate[grep(paste("S",57:64, "x", sep="",collapse="|"), histDL$sample_name)]=2
	histDL$replicate[grep(paste("S",65:72, "x", sep="",collapse="|"), histDL$sample_name)]=3
}
if(dataset=="V5V6V7"){
	histDL$replicate[grep(paste("S",1:8, "x", sep="",collapse="|"), histDL$sample_name)]=1
	histDL$replicate[grep(paste("S",9:16, "x", sep="",collapse="|"), histDL$sample_name)]=2
	histDL$replicate[grep(paste("S",17:24, "x", sep="",collapse="|"), histDL$sample_name)]=3
}
if(dataset=="meta"){
	histDL$replicate=rep(1, nrow(histDL))
}
paste(histDL$replicate, histDL$organism, sep="_")->histDL$org_rep
histDL$org_rep=factor(histDL$org_rep, levels=unique(histDL$org_rep))
	
	rownames(AtMeta)->taxa_color_pairs
	COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
	#ALPHA_COLUMN_IN_COLORSLIST=7
	source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
	"green"->taxa_color_pairs[which(taxa_color_pairs[,1]=="HOST"),2]
	taxa_color_pairs[match(levels(histDL$organism), taxa_color_pairs[,1]),2]->histDL$colors
	
	
#relative abundance
#pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq168/LOAD_HiSeq_Titration_ReAb_", date, ".pdf", sep="", collapse=""), width = 5, height = 3, useDingbats=FALSE)
ggplot(histDL, aes(x=histDL$sample_category, y=histDL$abundance, group=histDL$org_rep, color=histDL$organism)) +
  geom_line() + geom_point(size=2) + scale_color_manual(values=histDL$colors) + 
  theme_minimal() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
  theme(legend.position="none") 
#dev.off()  


norm="ReAbavg"
norm="LOADavg"
AtMeta_load=AtMeta

divide_by_host=function(tab){
	tab["HOST", which(tab["HOST",]==0)]<-0.00000001
	for(c in 1:ncol(tab)){
		tab[,c]/tab["HOST",c]->tab[,c]
		if( max(tab[,c])>100 ){
			tab[,c]=200*tab[,c]/sum(tab[,c])
		}
	}
	#tab[which(tab>100)]<-100
	return(tab)
}
	
AtMeta_load=divide_by_host(AtMeta_load)
AtMeta_load=AtMeta_load[order(rowSums(AtMeta_load)),]


#remove host amplicon
AtMeta_load=AtMeta_load[which(rownames(AtMeta_load)!="HOST"),]


#AtMeta_load=AtMeta_load[grep("Pseudomonadaceae", rownames(AtMeta_load), invert=TRUE),]
if(norm=="LOAD"){ AtMeta_load=add.empty.taxa(AtMeta_load, rescale=FALSE) }
if(norm=="LOADavg"){ 
	if(dataset!="meta"){
		AtMeta_load=(AtMeta_load[,1:8]+AtMeta_load[,9:16]+AtMeta_load[,17:24])/3
	}
		AtMeta_load=add.empty.taxa(AtMeta_load, rescale=FALSE)
}
if(norm=="ReAb"){ AtMeta_load=normalize100(AtMeta_load) }
if(norm=="ReAbavg"){ 
	AtMeta_load=normalize100(AtMeta_load)
	if(dataset!="meta"){
		AtMeta_load=(AtMeta_load[,1:8]+AtMeta_load[,9:16]+AtMeta_load[,17:24])/3
	}
	AtMeta_load=normalize100(AtMeta_load)
}

#next set specific order
	toporder=c("empty", "low_abundance")
	bottomorder=c("Nakamurellaceae", "Armatimonadaceae", "Nocardiaceae", "Polyangiaceae", "Methylophilaceae", "Burkholderiales_incertae_sedis", "Mycobacteriaceae", "Rhodobacteraceae", "Alcaligenaceae", "Xanthomonadaceae", "Verrucomicrobiaceae", "Pseudonocardiaceae", "Moraxellaceae", "Bradyrhizobiaceae", "Kineosporiaceae", "Micrococcaceae", "Acetobacteraceae", "Geodermatophilaceae", "Nocardioidaceae", "Micromonosporaceae", "Aurantimonadaceae", "Chitinophagaceae", "Hyphomicrobiaceae", "Rhizobiaceae", "Caulobacteraceae", "Methylobacteriaceae", "Microbacteriaceae", "Comamonadaceae", "Flavobacteriaceae", "Oxalobacteraceae", "Sphingobacteriaceae", "Cytophagaceae", "Sphingomonadaceae", "Pseudomonadaceae", "HOST")	
	AtMeta_load=topOrder(AtMeta_load, toporder, bottomorder)

abundant_families=grep("low_abundance", rownames(AtMeta_load[which(apply(AtMeta_load, 1, min)>-1),]), invert=TRUE, value=TRUE)
abundant_families=abundant_families[(length(abundant_families)-9):length(abundant_families)]
abundant_families=c("Caulobacteraceae", "Methylobacteriaceae", "Microbacteriaceae", "Comamonadaceae", "Flavobacteriaceae", "Oxalobacteraceae", "Sphingobacteriaceae", "Cytophagaceae", "Sphingomonadaceae", "Pseudomonadaceae")

melt(AtMeta_load)->histDL
c("organism", "sample_name", "abundance")->names(histDL)
gsub(".*xmeta_", "", histDL$sample_name)->histDL$sample_category
histDL$sample_category=factor(histDL$sample_category, levels=c("S80","S63","S71","S36","S83","S61","S19","S15"))
rep(0, nrow(histDL))->histDL$replicate
if(dataset=="V3V4"){
	histDL$replicate[grep(paste("S",25:32, "x", sep="",collapse="|"), histDL$sample_name)]=1
	histDL$replicate[grep(paste("S",33:40, "x", sep="",collapse="|"), histDL$sample_name)]=2
	histDL$replicate[grep(paste("S",41:48, "x", sep="",collapse="|"), histDL$sample_name)]=3
	if(norm=="LOADavg"){histDL$replicate=rep(1, nrow(histDL))}
	breaks=c(0, .5, 1, 1.5, 2)
	breaklabels=c(0, .5, 1, 1.5, 2)
	Yaxmin=0
	Yaxmax=2
	histloadmax=8
}
if(dataset=="V4"){
	histDL$replicate[grep(paste("S",49:56, "x", sep="",collapse="|"), histDL$sample_name)]=1
	histDL$replicate[grep(paste("S",57:64, "x", sep="",collapse="|"), histDL$sample_name)]=2
	histDL$replicate[grep(paste("S",65:72, "x", sep="",collapse="|"), histDL$sample_name)]=3
	if(norm=="LOADavg"){histDL$replicate=rep(1, nrow(histDL))}
	breaks=c(0, .5, 1, 1.5, 2)
	breaklabels=c(0, .5, 1, 1.5, 2)
	Yaxmin=0
	Yaxmax=2
	histloadmax=20
}
if(dataset=="V5V6V7"){
	histDL$replicate[grep(paste("S",1:8, "x", sep="",collapse="|"), histDL$sample_name)]=1
	histDL$replicate[grep(paste("S",9:16, "x", sep="",collapse="|"), histDL$sample_name)]=2
	histDL$replicate[grep(paste("S",17:24, "x", sep="",collapse="|"), histDL$sample_name)]=3
	if(norm=="LOADavg"){histDL$replicate=rep(1, nrow(histDL))}
	breaks=c(0, .5, 1, 1.5, 2)
	breaklabels=c(0, .5, 1, 1.5, 2)
	Yaxmin=0
	Yaxmax=2
	histloadmax=3
}
if(dataset=="meta"){
	histDL$replicate=rep(1, nrow(histDL))
	breaks=c(0, .5, 1)
	breaklabels=c(0, .5, 1)
	Yaxmin=0
	Yaxmax=1
	histloadmax=.8
}


paste(histDL$replicate, histDL$organism, sep="_")->histDL$org_rep
paste(histDL$sample_category, histDL$replicate, sep="_")->histDL$sample_rep
histDL$org_rep=factor(histDL$org_rep, levels=unique(histDL$org_rep))
histDL$sample_rep=factor(histDL$sample_rep, levels=c("S80_1","S80_2","S80_3","S63_1","S63_2","S63_3","S71_1","S71_2","S71_3","S36_1","S36_2","S36_3","S83_1","S83_2","S83_3","S61_1","S61_2","S61_3","S19_1","S19_2","S19_3","S15_1","S15_2","S15_3"))
histDL$organism=factor(histDL$organism, levels=unique(as.character(histDL$organism)))


	rownames(AtMeta)->taxa_color_pairs
	COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
	#ALPHA_COLUMN_IN_COLORSLIST=7
	source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
	"green"->taxa_color_pairs[which(taxa_color_pairs[,1]=="HOST"),2]
	"gray"->taxa_color_pairs[which(taxa_color_pairs[,1]=="low_abundance"),2]
	"#FFBC00"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Caulobacteraceae"),2]
	"#f43df0"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Flavobacteriaceae"),2]
	"#c7f214"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Sphingobacteriaceae"),2]
	taxa_color_pairs[match(levels(histDL$organism), taxa_color_pairs[,1]),2]->histDL$colors
	histDL$alphas=rep(0.3, nrow(histDL))
	histDL$alphas[match(histDL$organism, abundant_families, nomatch=0)>0]=1

if("a"=="b"){		
histDL2=histDL
histDL2$abundance=sqrt(sqrt(histDL2$abundance))
date=format(Sys.Date(), format="%Y%m%d")
#pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq168/LOAD_HiSeq_Titration_Corrected_4ROOT", date, ".pdf", sep="", collapse=""), width = 5, height = 3, useDingbats=FALSE)
ggplot(histDL2, aes(x=sample_category, y=abundance, group=org_rep, color=organism)) +
  geom_line() + geom_point(size=2) + scale_color_manual(values=histDL$colors) + 
  theme_minimal() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
  theme(legend.position="none") # +
	#scale_y_continuous(limits = c(0,2))
}

#Histogram with triplicates
#if(dataset=="V3V4"){Ylim=10}
#if(dataset=="V4"){Ylim=25}
#if(dataset=="V5V6V7"){Ylim=25}
#if(norm=="ReAb"){ Ylim=101 }
date=format(Sys.Date(), format="%Y%m%d")
#pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/AtMeta_", dataset, "_", norm, "_", date, ".pdf", sep="", collapse=""), width = 5, height = 2, useDingbats=FALSE)
ggplot(histDL,aes(x=sample_rep,y=abundance, fill=organism, alpha=histDL$alphas)) + 
 	 geom_bar(position="stack", stat="identity", width=1.0) + 
 	 scale_fill_manual(values=histDL$colors) + theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
 	 xlab("") + ylab("") + guides(fill=FALSE) +
 	 theme_classic() + theme(legend.position="none") + 
		scale_alpha_continuous(range = c(0.3, 1)) +
		scale_y_continuous(limits=c(0, histloadmax))	
#dev.off()


#LINES
histDL2=histDL[grep("HOST", histDL$organism, invert=TRUE),]
histDL2=histDL2[match(histDL2$organism, abundant_families, nomatch=0)>0,]
histDL2$organism=factor(histDL2$organism, levels=unique(as.character(histDL2$organism)))

date=format(Sys.Date(), format="%Y%m%d")
pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/AtMeta_", dataset, "_lines_", norm, "_", date, ".pdf", sep="", collapse=""), width = 5, height = 4, useDingbats=FALSE)
ggplot(histDL2, aes(x=sample_category, y=sqrt(sqrt(abundance)), group=org_rep, color=organism)) +
  geom_line() + geom_point(size=2) + scale_color_manual(values=histDL2$colors) + 
  theme_minimal() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
  theme(legend.position="none") +
  theme(panel.grid.minor = element_line(size = 0.15), panel.grid.major = element_line(size = .15)) + 
			theme(panel.grid.minor = element_line(colour = "black"), panel.grid.major = element_line(colour = "black")) + 
			theme(panel.grid.minor = element_blank()) +
	scale_y_continuous(name = "abundance", breaks = breaks, labels = breaklabels , limits = c(Yaxmin,Yaxmax)) 	
dev.off()


#V5V6V7abundant_families=abundant_families


######>
	######> METAGENOME COMPARISON
######> 
	
metaB=normalize100(metareads_B[,which(colSums(metareads_B)>=15000)])
HOST=metaB["HOST",]
metaB=metaB[setdiff(rownames(metaB), "HOST"),]
metaB=metaB[which(rowSums(metaB)>0),]
metaB=low_abundance(metaB, percent=0.01)
#add back host
metaB=rbind(metaB, HOST)
metaB=divide_by_host(metaB)

currentreadsX=metaB
currentreadsY=AtMeta_load
currentreadsX=(135/5)*metaB
currentreadsY=(2/5)*AtMeta_load

intersect(rownames(currentreadsX), rownames(currentreadsY))->common
common=setdiff(common, "low_abundance")
#common=c("Pseudomonadaceae", "Sphingomonadaceae" ,    "Sphingobacteriaceae" ,  "Oxalobacteraceae")
currentreadsX=currentreadsX[common,]
currentreadsY=currentreadsY[common,]

rownames(currentreadsX)->taxa_color_pairs

melt(currentreadsX)->histDL1
histDL1[,c(2, 1, 3)]->histDL1
c("sample1", "taxa1", "abundance1")->names(histDL1)
factor(histDL1$taxa1, levels=rownames(currentreadsX))->histDL1$taxa1

melt(currentreadsY)->histDL2
histDL2[,c(2, 1, 3)]->histDL2
c("sample2", "taxa2", "abundance2")->names(histDL2)
factor(histDL2$taxa2, levels=rownames(currentreadsY))->histDL2$taxa2

cbind(histDL1, histDL2)->histDL3
zeros=intersect(which(histDL3$abundance1==0), which(histDL3$abundance2==0))
histDL3=histDL3[setdiff(1:nrow(histDL3), zeros),]

COLUMN_IN_COLORSLIST=3    # 3 is colors, 4 is black / white / green / purple
ALPHA_COLUMN_IN_COLORSLIST=7
source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
	"gray"->taxa_color_pairs[which(taxa_color_pairs[,1]=="low_abundance"),2]
	"#FFBC00"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Caulobacteraceae"),2]
	"#f43df0"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Flavobacteriaceae"),2]
	"#c7f214"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Sphingobacteriaceae"),2]
	taxa_color_pairs[match(histDL3$taxa1, taxa_color_pairs[,1]),2]->histDL3$colors
	taxa_color_pairs[match(histDL3$taxa1, taxa_color_pairs[,1]),3]->histDL3$alphas
	as.numeric(histDL3$alphas)->histDL3$alphas
	histDL3$alphas=rep(0.3, nrow(histDL3))
	histDL3$alphas[match(histDL3$taxa1, abundant_families, nomatch=0)>0]=1
	histDL3$taxa1=factor(histDL3$taxa1, levels=unique(histDL3$taxa1))
	histDL3$taxa2=factor(histDL3$taxa2, levels=unique(histDL3$taxa2))

correlation=round(cor(sqrt(sqrt(histDL3$abundance1)), sqrt(sqrt(histDL3$abundance2)), method="pearson"), 4)^2

#SCATTER
date=format(Sys.Date(), format="%Y%m%d")
#pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/AtMeta_Meta_vs_", dataset, "_SCATTER_", norm, "_", date, ".pdf", sep="", collapse=""), width = 5, height = 4, useDingbats=FALSE)
ggplot(histDL3, aes(x=abs(sqrt(sqrt(histDL3$abundance1))), y=abs(sqrt(sqrt(histDL3$abundance2))), color=histDL3$taxa1, alpha=histDL3$alphas)) + 
	geom_point(shape=16, size=3) + 
	theme_minimal() + theme(axis.text.x = element_text(size=4),axis.text.y = element_text(size=4)) +  
	scale_color_manual(values=histDL3$colors) + theme(legend.position="none") +  
	scale_y_continuous(name = "", breaks = c(0, .5, 1, 1.5, 2), labels = c(0, .5, 1, 1.5, 2)^4, limits = c(0,2)) + 
	scale_x_continuous(name = "", breaks = c(0, .5, 1), labels = c(0, .5, 1)^4, limits = c(0,1)) + 
	annotate(geom="text", size=3, x = .7, y = .8, label = deparse(bquote(R^2~"="~.(correlation))), parse=T) + 
	theme(panel.grid.minor = element_line(size = 0.15), panel.grid.major = element_line(size = .15)) + 
	theme(panel.grid.minor = element_line(colour = "black"), panel.grid.major = element_line(colour = "black")) + 
	theme(panel.grid.minor = element_blank())
#dev.off()

ggplot(histDL3, aes(x=abs(sqrt(sqrt(histDL3$abundance1))), y=abs(sqrt(sqrt(histDL3$abundance2))), color=histDL3$taxa1, alpha=histDL3$alphas)) + 
	geom_point(shape=16, size=3) + 
	theme_minimal() + theme(axis.text.x = element_text(size=4),axis.text.y = element_text(size=4)) +  
	scale_color_manual(values=histDL3$colors) + theme(legend.position="none") +  
	scale_y_continuous(name = "", breaks = c(0, .5, 1, 1.5, 2), labels = c(0, .5, 1, 1.5, 2)^4, limits = c(0,2)) + 
	scale_x_continuous(name = "", breaks = c(0, .5, 1, 1.5, 2), labels = c(0, .5, 1, 1.5, 2)^4, limits = c(0,2)) + 
	annotate(geom="text", size=3, x = .7, y = .8, label = deparse(bquote(R^2~"="~.(correlation))), parse=T) + 
	theme(panel.grid.minor = element_line(size = 0.15), panel.grid.major = element_line(size = .15)) + 
	theme(panel.grid.minor = element_line(colour = "black"), panel.grid.major = element_line(colour = "black")) + 
	theme(panel.grid.minor = element_blank())


#>
	#>
		#>
			#> #Compare wildDNA concentration
		#>
	#>
#>


WildConc=grep("S73xV4GI502xwildDNA5ng|S74xV4GI502xwildDNA10ng|S75xV4GI502xwildDNA25ng|S76xV4GI502xwildDNA50ng|S77xV4GI502xwildDNA100ng|S78xV4GI502xwildDNA200ng|S79xV4GI502xwildDNA300ng|S80xV4GI502xwildDNA500ng", colnames(otureads), value=TRUE)

WildConc=otureads[,c("S73xV4GI502xwildDNA5ng","S74xV4GI502xwildDNA10ng","S75xV4GI502xwildDNA25ng","S76xV4GI502xwildDNA50ng","S77xV4GI502xwildDNA100ng","S78xV4GI502xwildDNA200ng","S79xV4GI502xwildDNA300ng","S80xV4GI502xwildDNA500ng")]


#compare16S=otureads[,c(Og16S, OgHam, Rm16S, RmHam, Og16S_7_10, OgHam_7_10,plasmid_2,plasmid_7_10)]
#compare16S=compare16S[which(rowSums(compare16S)>0),]


#rename HOST and add it back
if(LANE!="L1_REDO"){
	HOST=colSums(WildConc[c("Otu1", "Otu5748"),])
	WildConc=WildConc[setdiff(rownames(WildConc), c("Otu1", "Otu5748")),]
}
if(LANE=="L1_REDO"){
	HOST=colSums(WildConc[c("Otu1", "Otu17444"),])
	WildConc=WildConc[setdiff(rownames(WildConc), c("Otu1", "Otu17444")),]
}
WildConc=rbind(WildConc, HOST)

WildConc=normalize100(WildConc)
#AtMeta<-filter_by_taxonomyic_level(countTable=AtMeta, taxonomy=taxonomy, keywords=c("f"), last_filter=TRUE)
WildConc=WildConc[which(rowSums(WildConc)>0),]
WildConc=low_abundance(WildConc, percent=0.01)
	
#get rid of low abundance
#compare16S=compare16S[grep("low_abundance", rownames(compare16S), invert=TRUE),]


#order by abundance
order(rowSums(WildConc), decreasing=FALSE)->ordering_vector
		WildConc[ordering_vector,]->WildConc

#next set specific order
	toporder=c("empty", "low_abundance")
	bottomorder=c("HOST")		
	WildConc=topOrder(WildConc, toporder, bottomorder)


melt(WildConc)->histDL
c("organism", "sample_name", "abundance")->names(histDL)
gsub(".*actinx", "", histDL$sample_name)->histDL$sample_category
L=length(unique(histDL$organism))
histDL$organism=factor(histDL$organism, levels=unique(histDL$organism))
histDL$sample_name=factor(histDL$sample_name, levels=unique(histDL$sample_name))
		
	unique(as.matrix(histDL$organism))->taxa_color_pairs
	COLUMN_IN_COLORSLIST=10    # 3 is colors, 4 is black / white / green / purple
	#ALPHA_COLUMN_IN_COLORSLIST=7
	source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
	"black"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Hpa"),2]
	"yellow"->taxa_color_pairs[which(taxa_color_pairs[,1]=="DC3000"),2]
	"green"->taxa_color_pairs[which(taxa_color_pairs[,1]=="HOST"),2]
	
	taxa_color_pairs[match(levels(histDL$organism), taxa_color_pairs[,1]),2]->histDL$colors
	
#histDL$abundance=sqrt(histDL$abundance)
#relative abundance histogram
date=format(Sys.Date(), format="%Y%m%d")
pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/WildTotalConcentrations_", date, ".pdf", sep="", collapse=""), width = 5, height = 4, useDingbats=FALSE)
	ggplot(histDL,aes(x=sample_name,y=abundance, fill=organism)) + 
 	 geom_bar(position="fill", stat="identity", width=1.0) + 
 	 scale_fill_manual(values=histDL$colors) + theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
 	 xlab("") + ylab("") + guides(fill=FALSE) +
 	 theme_classic() + theme(legend.position="none") + 
		scale_alpha_continuous(range = c(0.3, 1)) 	
dev.off()


#histDL$abundance=sqrt(histDL$abundance)
#relative abundance histogram
date=format(Sys.Date(), format="%Y%m%d")
histDL3=histDL[grep("HOST", histDL$organism, invert=TRUE),]
#histDL3$organism=factor(histDL3$organism, levels=unique(as.character(histDL3$organism)))
pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/WildTotalConcentration_ReAb_", date, ".pdf", sep="", collapse=""), width = 5, height = 4, useDingbats=FALSE)
	ggplot(histDL3,aes(x=sample_name,y=abundance, fill=organism)) + 
 	 geom_bar(position="fill", stat="identity", width=1.0) + 
 	 scale_fill_manual(values=histDL3$colors) + theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
 	 xlab("") + ylab("") + guides(fill=FALSE) +
 	 theme_classic() + theme(legend.position="none") + 
		scale_alpha_continuous(range = c(0.3, 1)) 	
dev.off()

compare16S=sqrt(sqrt(compare16S))

breaklabels=sqrt(sqrt(c(0, 0.05, 1, 4, 16, 36, 64, 100)))
exponent=4
axmax=2.5

#16Sog vs. 16Sog
	X=c(compare16S[,Og16S[1]], compare16S[,Og16S[2]])
	Y=c(compare16S[,Og16S[3]], compare16S[,Og16S[4]])
	#threshold rare OTUs
	thresholded=which(apply(cbind(X, Y), 1, min)>=0.05)
	X=X[thresholded]
	Y=Y[thresholded]
	ks=ks.test(X, Y)
		as.data.frame(cbind(X, Y))->histDL3
		c("X", "Y")->names(histDL3)
		rownames(histDL3)->histDL3$unique_comparison
		correlation=cor(histDL3$X, histDL3$Y)^2
		rep("black", nrow(histDL3))->histDL3$color
		#0, 0.02, 0.1, 1, 2, 4, 8, 16
		breaklabel=c(0, 0.3760603, 0.5623413, 1, 1.189207, 1.414214, 1.681793, 2, 2.213364)
		#pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq168/Correlation_Wild_16S_v_16S_", date, ".pdf", sep="", collapse=""), width = 3, height = 3, useDingbats=FALSE)
		ggplot(histDL3, aes(x=X, y=Y)) + geom_point(size=1.5, color=histDL3$color, alpha=0.2) + 
			theme_minimal() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
			theme(legend.position="none") + 
			coord_fixed(ratio=1/1) +  
			scale_x_continuous(name = "a + b", breaks = breaklabels, labels = breaklabels^exponent, limits = c(0,axmax)) + 
			scale_y_continuous(name = "c + d", breaks = breaklabels, labels = breaklabels^exponent, limits = c(0,axmax)) +
			theme(panel.grid.minor = element_line(size = 0.15), panel.grid.major = element_line(size = .15)) + 
			theme(panel.grid.minor = element_line(colour = "black"), panel.grid.major = element_line(colour = "black")) + 
			theme(panel.grid.minor = element_blank()) +
			geom_smooth(method='lm', size=0.5, se = FALSE) +
			annotate(geom="text", size=3, x = 16^(1/exponent), y = 36^(1/exponent), label = deparse(bquote(R^2~"="~.(correlation))), parse=T) +
			annotate(geom="text", size=3, x = 16^(1/exponent), y = 24^(1/exponent), label = deparse(paste("K-S = ", round(ks[[2]],3), sep="")), parse=T) 
		#dev.off()	


#>
	#>
		#>
			#> #Compare plasmid DNA concentration
		#>
	#>
#>

PlasmidConc=otureads[,c("S81xV4GI502ITSoxplasmid0.48pg","S82xV4GI502ITSoxplasmid0.96pg","S83xV4GI502ITSoxplasmid1.92pg","S84xV4GI502ITSoxplasmid2.4pg","S85xV4GI502ITSoxplasmid4.8pg","S86xV4GI502ITSoxplasmid9.6pg","S87xV4GI502ITSoxplasmid19.2pg","S88xV4GI502ITSoxplasmid24pg","S89xV4GI502ITSoxplasmid48pg","S90xV4GI502ITSoxplasmid96pg","S91xV4GI502ITSoxplasmid192pg","S92xV4GI502ITSoxplasmid240pg","S93xV4GI502ITSoxplasmid480pg","S94xV4GI502ITSoxplasmid960pg","S95xV4GI502ITSoxplasmid1920pg","S96xV4GI502ITSoxplasmid2400pg")]


#rename HOST and ITS add it back
HOST=PlasmidConc[c("Otu1"),]
PlasmidConc=PlasmidConc[setdiff(rownames(PlasmidConc), c("Otu1")),]
ITS=PlasmidConc[c("Otu4"),]
PlasmidConc=PlasmidConc[setdiff(rownames(PlasmidConc), c("Otu4")),]
DC3000=PlasmidConc[c("Otu2"),]
PlasmidConc=PlasmidConc[setdiff(rownames(PlasmidConc), c("Otu2")),]


PlasmidConc=rbind(HOST, ITS, DC3000)

PlasmidConc=normalize100(PlasmidConc)

melt(PlasmidConc)->histDL
c("organism", "sample_name", "abundance")->names(histDL)
gsub(".*actinx", "", histDL$sample_name)->histDL$sample_category
histDL$sample_category=factor(histDL$sample_category, levels=GautamLevels)
L=length(unique(histDL$organism))
histDL$organism=factor(histDL$organism, levels=unique(histDL$organism))
histDL$sample_name=factor(histDL$sample_name, levels=unique(histDL$sample_name))
		
	unique(as.matrix(histDL$organism))->taxa_color_pairs
	COLUMN_IN_COLORSLIST=10    # 3 is colors, 4 is black / white / green / purple
	#ALPHA_COLUMN_IN_COLORSLIST=7
	source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
	"black"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Hpa"),2]
	"yellow"->taxa_color_pairs[which(taxa_color_pairs[,1]=="DC3000"),2]
	"green"->taxa_color_pairs[which(taxa_color_pairs[,1]=="HOST"),2]
	
	taxa_color_pairs[match(levels(histDL$organism), taxa_color_pairs[,1]),2]->histDL$colors
	
#histDL$abundance=sqrt(histDL$abundance)
#relative abundance histogram
	ggplot(histDL,aes(x=sample_name,y=abundance, fill=organism)) + 
 	 geom_bar(position="fill", stat="identity", width=1.0) + 
 	 scale_fill_manual(values=histDL$colors) + theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
 	 xlab("") + ylab("") + guides(fill=FALSE) +
 	 theme_classic() + theme(legend.position="none") + 
		scale_alpha_continuous(range = c(0.3, 1)) 	



 #PLOT HISTOGRAM 

#pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq168/LOAD_HiSeq_", current, "_untrimmed", date, ".pdf", sep="", collapse=""), width = 3, height = 2, useDingbats=FALSE)
	ggplot(histDL,aes(x=histDL$sample,y=histDL$abundance, fill=histDL$taxa)) + 
 	 geom_bar(position="fill", stat="identity", width=1.0) + 
 	 scale_fill_manual(values=histDL$colors) + theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
 	 xlab("") + ylab("") + guides(fill=FALSE) +
 	 theme_classic() + theme(legend.position="none") + 
		scale_alpha_continuous(range = c(0.3, 1)) 	
#dev.off()
	


#>
	#>
		#>
			#> #effect of primer concentration on product LANE 4
		#>
	#>
#>


AlterPrim=otureads[,c("S9xV4GI502xwild_GIto16S_10to1", "S10xV4GI502xwild_GIto16S_6to1", "S11xV4GI502xwild_GIto16S_3to1", "S12xV4GI502xwild_GIto16S_1to1", "S13xV4GI502xwild_GIto16S_1to1", "S14xV4GI502xwild_GIto16S_1to3", "S15xV4GI502xwild_GIto16S_1to6", "S16xV4GI502xwild_GIto16S_1to10")]

AlterPrim=normalize100(AlterPrim)

#rename HOST and add it back
HOST=AlterPrim[c("Otu5"),]
AlterPrim=AlterPrim[setdiff(rownames(AlterPrim), c("Otu5")),]


AlterPrim<-filter_by_taxonomyic_level(countTable=AlterPrim, taxonomy=taxonomy, keywords=c("p"), last_filter=TRUE)
AlterPrim=AlterPrim[which(rowSums(AlterPrim)>0),]
AlterPrim=low_abundance(AlterPrim, percent=0.01)
	
AlterPrim=rbind(AlterPrim, HOST)

#order by abundance
order(rowSums(AlterPrim), decreasing=FALSE)->ordering_vector
		AlterPrim[ordering_vector,]->AlterPrim

#next set specific order
	toporder=c("empty", "low_abundance")
	bottomorder=c("HOST")		
	AlterPrim=topOrder(AlterPrim, toporder, bottomorder)

melt(AlterPrim)->histDL

c("organism", "sample_name", "abundance")->names(histDL)
#gsub(".*actinx", "", histDL$sample_name)->histDL$sample_category
L=length(unique(histDL$organism))
histDL$organism=factor(histDL$organism, levels=unique(histDL$organism))
histDL$sample_name=factor(histDL$sample_name, levels=unique(histDL$sample_name))
		
	unique(as.matrix(histDL$organism))->taxa_color_pairs
	COLUMN_IN_COLORSLIST=10    # 3 is colors, 4 is black / white / green / purple
	#ALPHA_COLUMN_IN_COLORSLIST=7
	source("/Users/dlundberg/Documents/abt6/Regalado2017/family_colorscheme.R")
	"black"->taxa_color_pairs[which(taxa_color_pairs[,1]=="Hpa"),2]
	"yellow"->taxa_color_pairs[which(taxa_color_pairs[,1]=="DC3000"),2]
	"green"->taxa_color_pairs[which(taxa_color_pairs[,1]=="HOST"),2]
	
	taxa_color_pairs[match(levels(histDL$organism), taxa_color_pairs[,1]),2]->histDL$colors
	
#histDL$abundance=sqrt(histDL$abundance)
#relative abundance histogram
date=format(Sys.Date(), format="%Y%m%d")
pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/AlteredPrimerConcentrations_", date, ".pdf", sep="", collapse=""), width = 5, height = 4, useDingbats=FALSE)
	ggplot(histDL,aes(x=sample_name,y=abundance, fill=organism)) + 
 	 geom_bar(position="fill", stat="identity", width=1.0) + 
 	 scale_fill_manual(values=histDL$colors) + theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
 	 xlab("") + ylab("") + guides(fill=FALSE) +
 	 theme_classic() + theme(legend.position="none") + 
		scale_alpha_continuous(range = c(0.3, 1)) 	
dev.off()


#calculate GI / 16S 
product_ratio=vector(length=8)
for(i in 1:8){
	Hsum=AlterPrim[grep("HOST", rownames(AlterPrim)),i]
	Msum=sum(AlterPrim[grep("HOST", rownames(AlterPrim), invert=TRUE),i])
	product_ratio[i]=(Msum/Hsum)
	
}

primer_ratio=c(1/10, 1/6, 1/3, 1/1, 1/1, 3/1, 6/1, 10/1)

pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/AlteredPrimerConcentrations_scatter1_", date, ".pdf", sep="", collapse=""), width = 4, height = 4, useDingbats=FALSE)
	plot(primer_ratio, product_ratio, pch=16, cex=1, xlim=c(0, 10), ylim=c(0, 300))
dev.off()


plot(primer_ratio, sqrt(product_ratio), pch=16, cex=1, xlim=c(0, 10), ylim=c(0, 20))

#histDL$abundance=sqrt(histDL$abundance)
#relative abundance histogram
date=format(Sys.Date(), format="%Y%m%d")
histDL3=histDL[grep("HOST", histDL$organism, invert=TRUE),]
#histDL3$organism=factor(histDL3$organism, levels=unique(as.character(histDL3$organism)))
pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq175/WildTotalConcentration_ReAb_", date, ".pdf", sep="", collapse=""), width = 5, height = 4, useDingbats=FALSE)
	ggplot(histDL3,aes(x=sample_name,y=abundance, fill=organism)) + 
 	 geom_bar(position="fill", stat="identity", width=1.0) + 
 	 scale_fill_manual(values=histDL3$colors) + theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
 	 xlab("") + ylab("") + guides(fill=FALSE) +
 	 theme_classic() + theme(legend.position="none") + 
		scale_alpha_continuous(range = c(0.3, 1)) 	
dev.off()

compare16S=sqrt(sqrt(compare16S))

breaklabels=sqrt(sqrt(c(0, 0.05, 1, 4, 16, 36, 64, 100)))
exponent=4
axmax=2.5

#16Sog vs. 16Sog
	X=c(compare16S[,Og16S[1]], compare16S[,Og16S[2]])
	Y=c(compare16S[,Og16S[3]], compare16S[,Og16S[4]])
	#threshold rare OTUs
	thresholded=which(apply(cbind(X, Y), 1, min)>=0.05)
	X=X[thresholded]
	Y=Y[thresholded]
	ks=ks.test(X, Y)
		as.data.frame(cbind(X, Y))->histDL3
		c("X", "Y")->names(histDL3)
		rownames(histDL3)->histDL3$unique_comparison
		correlation=cor(histDL3$X, histDL3$Y)^2
		rep("black", nrow(histDL3))->histDL3$color
		#0, 0.02, 0.1, 1, 2, 4, 8, 16
		breaklabel=c(0, 0.3760603, 0.5623413, 1, 1.189207, 1.414214, 1.681793, 2, 2.213364)
		#pdf(paste("/Users/dlundberg/Documents/abt6/Pratchaya/PCR_Load_Quantification/HiSeq168/Correlation_Wild_16S_v_16S_", date, ".pdf", sep="", collapse=""), width = 3, height = 3, useDingbats=FALSE)
		ggplot(histDL3, aes(x=X, y=Y)) + geom_point(size=1.5, color=histDL3$color, alpha=0.2) + 
			theme_minimal() + theme(axis.text.x = element_text(size=5),axis.text.y = element_text(size=5)) +  
			theme(legend.position="none") + 
			coord_fixed(ratio=1/1) +  
			scale_x_continuous(name = "a + b", breaks = breaklabels, labels = breaklabels^exponent, limits = c(0,axmax)) + 
			scale_y_continuous(name = "c + d", breaks = breaklabels, labels = breaklabels^exponent, limits = c(0,axmax)) +
			theme(panel.grid.minor = element_line(size = 0.15), panel.grid.major = element_line(size = .15)) + 
			theme(panel.grid.minor = element_line(colour = "black"), panel.grid.major = element_line(colour = "black")) + 
			theme(panel.grid.minor = element_blank()) +
			geom_smooth(method='lm', size=0.5, se = FALSE) +
			annotate(geom="text", size=3, x = 16^(1/exponent), y = 36^(1/exponent), label = deparse(bquote(R^2~"="~.(correlation))), parse=T) +
			annotate(geom="text", size=3, x = 16^(1/exponent), y = 24^(1/exponent), label = deparse(paste("K-S = ", round(ks[[2]],3), sep="")), parse=T) 
		#dev.off()	


