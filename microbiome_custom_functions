#matrix_format(matrix)

#low_abundance(readtable, percent=5)
	#takes taxa that have abundance below a certain threshold and puts them in a low abundance category
	#requires readtable to have rownames and colnames

#read.otutable(filename)
	#reads the table and takes the first column as otus and the first row as samples, converts the rest of the table to a numeric matrix

#remove_contaminant_taxa(countTable, taxonomy, keywords=c("Chloroplast", "Mitochondria", "Archaea"))
	#needs countTable with row and column names, a USEARCH taxonomy table, and some keywords to look for and eliminate. spits out only the reduced count table
	
#filter_by_taxonomyic_level()

#find_difference_taxa_tables<-function(matrix1, matrix2)
	#returns the difference of subtracting the second matrix from the first matrix
	





options(scipen=999)
	
read.countTable<-function(filename){
	reads=read.table(filename, comment.char = "", sep="\t", quote="")	#load OTU table
	matrix_format(reads)->reads
	reads[2:nrow(reads), 1]->otu_names
	matrix_format(reads[1, 2:ncol(reads)])->otu_samples
	matrix_format(reads[2:nrow(reads), 2:ncol(reads)])->reads
	if(nrow(reads)==1){
		t(reads)->reads
	}
	apply(reads,2,as.numeric)->reads
	otu_names->rownames(reads)
	otu_samples->colnames(reads)
	return(reads)
}

	
add.other.taxa<-function(matrix, subtract=TRUE){
	c("other", rownames(matrix))->newrows
	rbind(rep(0, ncol(matrix)), matrix)->matrix
	newrows->rownames(matrix)
	if(subtract==TRUE){
		max(colSums(matrix))->max_microbe_abundance
	}
	if(subtract!=TRUE){
		100->max_microbe_abundance
	}
	for(c in 1:ncol(matrix)){
		max_microbe_abundance-sum(matrix[,c])->matrix[1,c]
		100*matrix[,c]/sum(matrix[,c])->matrix[,c]
	}
	return(matrix)
}
add.empty.taxa<-function(matrix, subtract=TRUE, rescale=TRUE){
	c("empty", rownames(matrix))->newrows
	rbind(rep(0, ncol(matrix)), matrix)->matrix
	newrows->rownames(matrix)
	if(subtract==TRUE){
		max(colSums(matrix))->max_microbe_abundance
	}
	if(subtract!=TRUE){
		100->max_microbe_abundance
	}
	for(c in 1:ncol(matrix)){
		max_microbe_abundance-sum(matrix[,c])->matrix[1,c]
		if(rescale==TRUE){
			100*matrix[,c]/sum(matrix[,c])->matrix[,c]
		}
		if(rescale!=TRUE){
			max_microbe_abundance*matrix[,c]/sum(matrix[,c])->matrix[,c]
		}
	}
	return(matrix)
}

add.plant<-function(matrix){
	c("empty", "plant", rownames(matrix))->newrows
	rbind(rep(0, ncol(matrix)), rep(0, ncol(matrix)), matrix)->matrix
	newrows->rownames(matrix)
	max(colSums(matrix))->max_microbe_abundance
	for(c in 1:ncol(matrix)){
		max_microbe_abundance-sum(matrix[,c])->matrix[1,c]
		100->matrix[2,c]
		100*matrix[,c]/sum(matrix[,c])->matrix[,c]
	}
	return(matrix)
}

date.2<-function(){
	return(format(Sys.Date(), format="%Y%m%d"))
}


getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

plant_color_list<-matrix(c(
	'Athaliana', "#AACF37",
	'Arabdidopsis', "#AACF37",
	
    'Moss', "#2078B4", 
    
    'FuzzyRosemary', "#F59999",
    
    'Soil', "#9F4222", 
    
    'Draba', "#3DB549", 
    
    'SpringKraut', "#DA4699", 
	'Impatiens', "#DA4699", 
    
    'Cardamine', "#0C6D38", 
    
    'Plantago', "#FF0000",
    
    'Thistle', "#FF8B00", 
    "Sonchus",  "#FF8B00", 
    
    'BigClover', "#727272",
    'Clover', "#57C8E8",
    'SmallClover', "#57C8E8", 
    "Trifolium", "#57C8E8",
    
    'Grass', "#5A00EE",
    
    'Dandelion', "#fbe955",
    'Taraxacum', "#fbe955",
    
    "?_S4.2_Pseudobulk", "red",  
    
    'NeckarKing', "#000000",
    "NK", "#000000",
    
    "ShortTree", "#22F9B8",
	
	"blank1_blank", "#CAB2D6",
	"blank1_blank", "#000000"
), ncol=2, byrow=TRUE)




read.MetacountTable<-function(filename){
	reads=t(read.table(filename, sep="\t", comment.char = ""))	
	matrix_format(reads)->reads
	reads[2:nrow(reads), 1]->otu_names
	reads[1, 2:ncol(reads)]->otu_samples
	reads[2:nrow(reads), 2:ncol(reads)]->reads
	apply(reads,2,as.numeric)->reads
	otu_names->rownames(reads)
	otu_samples->colnames(reads)
	return(reads)
}
	
	
		
remove_contaminant_taxa<-function(countTable, taxonomy, keywords=c("Chloroplast", "Mitochondria", "Archaea"), invert=FALSE){
	taxonomy[match(rownames(countTable), taxonomy[,1], nomatch=0),4]->table_taxonomy  
	#nomatch=0 removed 20200227... messed everything up!!!
	#taxonomy[match(rownames(countTable), taxonomy[,1], nomatch=0),4]->table_taxonomy 
	vector(length=0)->contaminants
	for(word in keywords){
		c(contaminants, grep(word, table_taxonomy))->contaminants
	}
	unique(contaminants)->contaminants
	if(invert==FALSE){
		setdiff(1:length(table_taxonomy), contaminants)->keep_taxa
	}
	if(invert==TRUE){ contaminants->keep_taxa }
	table_taxonomy[keep_taxa]->table_taxonomy
	countTable[keep_taxa,]->countTable
	return(countTable)
}


filter_by_taxonomyic_level<-function(countTable, taxonomy, keywords=c("p"), last_filter=TRUE, clean_last_filter=TRUE){
	taxa_order=c("d:", "p:", "c:", "o:", "f:", "g:") 
	grep(keywords, taxa_order)->taxa_pos
	#nomatch=0 removed 20200227... messed everything up!!!
	#taxonomy[match(rownames(countTable), taxonomy[,1], nomatch=0),4]->table_taxonomy
	taxonomy[match(rownames(countTable), taxonomy[,1]),4]->table_taxonomy
	grep(paste(keywords,":", sep=""), table_taxonomy)->taxa
	countTable[taxa,]->countTable
	table_taxonomy[taxa]->table_taxonomy
	countTable->taxatable
	#clean up everything higher than the taxa named, only if last_filter=TRUE
	if(last_filter==TRUE){
		gsub(paste(".*", keywords, ":", sep=""),"", table_taxonomy)->table_taxonomy
		#clean up the crap after the chosen taxa
		if((taxa_pos)!=length(taxa_order)){
			for(t in (taxa_pos+1):length(taxa_order)){
				print(t)
				if(clean_last_filter==TRUE){
					gsub(paste(",",taxa_order[t],".*",sep=""), "", table_taxonomy)->table_taxonomy
					gsub(paste(taxa_order[t],".*",sep=""), "", table_taxonomy)->table_taxonomy  #check
				}
			}
		}
		taxa=unique(table_taxonomy)
		taxatable=matrix(ncol=ncol(countTable), nrow=length(taxa))
		for(t in taxa){
			current_row=which(t==taxa)
			apply(matrix_format(countTable[which(table_taxonomy==t),]), 2, sum)->taxatable[current_row,]
		}
		taxa->rownames(taxatable)
		colnames(countTable)->colnames(taxatable)
	}	
	return(taxatable)
}
filter_by_taxonomic_level<-filter_by_taxonomyic_level

low_abundance<-function(readtable, percent=5, return_low_abundance_rows=FALSE){

	#make a category for low abundance for phyla that have less than 5%
	rownames(readtable)->phyla
	colnames(readtable)->samples
	c(phyla, "low_abundance")->phyla
	#puts an extra row at the bottom of readtable for the low abundance
	rbind(readtable, vector(length=ncol(readtable)))->readtable
	
	#order(rowSums(readtable), decreasing=TRUE)->ordering_vector

	#phyla[ordering_vector]->phyla
	#readtable[ordering_vector,]->readtable

		#helper function
		lowAbund<-function(vector){
			which(vector<=percent)->LA
			sum(vector[LA])->theSum
			vector[LA]<-0
			theSum->vector[length(vector)]
			return(vector)
		}

	apply(readtable, 2, lowAbund)->readtable

	#remove phyla rows that are now "low abundance"
	which(rowSums(readtable)!=0)->save_rows
	which(rowSums(readtable)==0)->low_abundance_rows
	save_rows[length(save_rows):1]->save_rows
	readtable[save_rows,]->readtable
	phyla[save_rows]->phyla

	phyla->rownames(readtable)
	samples->colnames(readtable)
	if(return_low_abundance_rows==FALSE){
		return(readtable)
	}
	if(return_low_abundance_rows==TRUE){
		return(low_abundance_rows)
	}
	
}


topOrder<-function(countTable, topOrder, bottomOrder){
	match(topOrder,rownames(countTable), nomatch=0)->top
	setdiff(1:nrow(countTable), top)->bottom
	countTable[c(top, bottom),]->countTable
	
	match(bottomOrder,rownames(countTable), nomatch=0)->footer
	setdiff(1:nrow(countTable), footer)->header
	countTable[c(header, footer),]->countTable
}

	

#MINREADSxMINSAMPLES threshold
readsXsampleTHRES<-function(countTable, minreads=1, minsamples=5){
	countTable[rowSums(countTable)>(minreads*minsamples),]->countTable
	keep=vector(length=nrow(countTable))
	for(r in 1:nrow(countTable)){
		length(which(countTable[r,]>minreads))->keep[r]
	}
	countTable[which(keep>minsamples),]->countTable
	return(countTable)
}


matrix_format<-function(matrix, type=""){
	as.matrix(matrix)->matrix
	if(type=="" | type=="H"){
		if(ncol(matrix)==1){
			t(matrix)->matrix
		}
	}
	if(type=="V"){}
	return(matrix)
}	



find_difference_taxa_tables<-function(matrix1, matrix2){	

	#returns the difference of subtracting the second matrix from the first matrix
	
	#first make sure rows are the same, so direct subtraction of one table from the other is possble
	unique(c(rownames(matrix1),rownames(matrix2)))->all_rownames
	updated_matrix1=matrix(nrow=length(all_rownames), ncol=ncol(matrix1), data=0)
	updated_matrix2=matrix(nrow=length(all_rownames), ncol=ncol(matrix2), data=0)
	rownames(updated_matrix1)=rownames(updated_matrix2)=all_rownames
	colnames(updated_matrix1)=colnames(updated_matrix2)=colnames(matrix1)
	
	setdiff(all_rownames, rownames(matrix1))
	rownames(matrix1[match(all_rownames, rownames(matrix1), nomatch=0),])->destination_rownames
	matrix1[match(all_rownames, rownames(matrix1), nomatch=0),]->updated_matrix1[destination_rownames,]
	
	setdiff(all_rownames, rownames(matrix2))
	rownames(matrix2[match(all_rownames, rownames(matrix2), nomatch=0),])->destination_rownames
	matrix2[match(all_rownames, rownames(matrix2), nomatch=0),]->updated_matrix2[destination_rownames,]

	#Next subtract matrix 2 from matrix 1
	updated_matrix1-updated_matrix2->Mat1minusMat2
		
	#add a category for "blank space" for whatever doesn't add up to 100%, so that when graphed 
	#the difference result isn't scaled to a new 100% and renormalized
	positive_blank_space=rep(0, ncol(Mat1minusMat2))
	negative_blank_space=rep(0, ncol(Mat1minusMat2))
	rbind(positive_blank_space, negative_blank_space, Mat1minusMat2)->Mat1minusMat2
	
	for(c in 1:ncol(Mat1minusMat2)){
		100-sum(Mat1minusMat2[which(Mat1minusMat2[,c]>0),c])->pos
		-100-sum(Mat1minusMat2[which(Mat1minusMat2[,c]<0),c])->neg
		Mat1minusMat2["positive_blank_space", c]<-pos
		Mat1minusMat2["negative_blank_space", c]<-neg
	}
		
	return(Mat1minusMat2)
}		


normalize100<-function(matrix){
	for(c in 1:ncol(matrix)){
		100*matrix[,c]/sum(matrix[,c])->matrix[,c]
	}
	return(matrix)
}
	


combine_taxa_tables<-function(matrix1, matrix2, return="both"){	

	#returns the union of two taxa tables.. adds matrix2 to the right of matrix1, adding taxa and 0's as necessary
	#if return=1, adjusts the first table so that it could be combined with the second. But does not append second and returns 1st.
	#if return=2, same but returns second table
	
	#first make sure rows are the same, so direct subtraction of one table from the other is possble
	unique(c(rownames(matrix1),rownames(matrix2)))->all_rownames
	updated_matrix1=matrix(nrow=length(all_rownames), ncol=ncol(matrix1), data=0)
	updated_matrix2=matrix(nrow=length(all_rownames), ncol=ncol(matrix2), data=0)
	rownames(updated_matrix1)=rownames(updated_matrix2)=all_rownames
	colnames(updated_matrix1)=colnames(matrix1)
	colnames(updated_matrix2)=colnames(matrix2)
		
	setdiff(all_rownames, rownames(matrix1))
	rownames(matrix1[match(all_rownames, rownames(matrix1), nomatch=0),])->destination_rownames
	matrix1[match(all_rownames, rownames(matrix1), nomatch=0),]->updated_matrix1[destination_rownames,]
	
	setdiff(all_rownames, rownames(matrix2))
	rownames(matrix2[match(all_rownames, rownames(matrix2), nomatch=0),])->destination_rownames
	matrix2[match(all_rownames, rownames(matrix2), nomatch=0),]->updated_matrix2[destination_rownames,]

	#Next bind matrix 1 onto matrix 2
	if(return=="both"){
		cbind(updated_matrix1,updated_matrix2)->Mat1combMat2
		return(Mat1combMat2)
	}
	if(return==1){
		return(updated_matrix1)
	}
	if(return==2){
		return(updated_matrix2)
	}
}		




#rarefy
#the rarefaction function from the VEGAN library
rrarefy<-function (x, sample) 
{
    x <- as.matrix(x)
    if (ncol(x) == 1) 
        x <- t(x)
    if (length(sample) > 1 && length(sample) != nrow(x)) 
        stop(gettextf("length of 'sample' and number of rows of 'x' do not match"))
    sample <- rep(sample, length = nrow(x))
    colnames(x) <- colnames(x, do.NULL = FALSE)
    nm <- colnames(x)
    for (i in 1:nrow(x)) {
        row <- sample(rep(nm, times = x[i, ]), sample[i])
        row <- table(row)
        ind <- names(row)
        x[i, ] <- 0
        x[i, ind] <- row
    }
    x
}





ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
   # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
   # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.

   if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
   if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
   else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
   else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
   m[tri] <- t(m)[tri]
   return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
   # Returns a matrix (M) of distances between geographic points.
   # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
   # (df.geopoints$lat[j], df.geopoints$lon[j]).
   # The row and column names are given by df.geopoints$name.

   GeoDistanceInMetres <- function(g1, g2){
      # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
      # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
      # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
      # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
      # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
      DistM <- function(g1, g2){
         require("Imap")
         return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
      }
      return(mapply(DistM, g1, g2))
   }

   n.geopoints <- nrow(df.geopoints)

   # The index column is used to ensure we only do calculations for the upper triangle of points
   df.geopoints$index <- 1:n.geopoints

   # Create a list of lists
   list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})

   # Get a matrix of distances (in metres)
   mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")

   # Set the row and column names
   rownames(mat.distances) <- df.geopoints$name
   colnames(mat.distances) <- df.geopoints$name

   return(mat.distances)
}





#reverse complement
rc<-function(sequence, reverse=TRUE, complement=TRUE){
	as.matrix(sequence)->sequence
	if(complement==TRUE){
		tolower(sequence)->sequence
		gsub("a","T", sequence)->sequence
		gsub("c","G", sequence)->sequence
		gsub("g","C", sequence)->sequence
		gsub("t","A", sequence)->sequence
		gsub("u","U", sequence)->sequence
		gsub("r","Y", sequence)->sequence
		gsub("y","R", sequence)->sequence
		gsub("m","K", sequence)->sequence
		gsub("k","M", sequence)->sequence
		gsub("w","W", sequence)->sequence
		gsub("s","S", sequence)->sequence
		gsub("b","V", sequence)->sequence
		gsub("d","H", sequence)->sequence
		gsub("h","D", sequence)->sequence
		gsub("v","B", sequence)->sequence
		gsub("n","N", sequence)->sequence
	}
	if(reverse==TRUE){
		reverser<-function(sequence){
			strsplit(sequence,split="")[[1]]->sequence
			sequence[length(sequence):1]->sequence
			paste(sequence, collapse="", sep="")->sequence
			return(sequence)
		}
		apply(sequence, 1, reverser)->sequence
	}
	as.matrix(sequence)->sequence
	return(sequence)
}



split_FASTA<-function(FASTA, add_integer_column=TRUE, leave_carrot=TRUE, reverse=FALSE){
	as.matrix(FASTA)->FASTA
	
	if(reverse==FALSE){
		#Finds carrots
		grep(">",FASTA)->name_locations
		
		#Makes a space to hold the sequences, one for each carrot
		vector(length=length(name_locations))->sequences
		
		#Makes a new 'name locations.2' that pretends there is one more carrot at the end of the last sequence.
		#Then we assume all lines occuring between carrots are a sequence that should be concatenated.
		c(name_locations, max(dim(FASTA))+1)->name_locations.2
		for(s in 1:length(name_locations)){
			paste(FASTA[(name_locations.2[s]+1):(name_locations.2[s+1]-1),],collapse="")->sequences[s]
		}
		
		#Uses the name locations to pull ou the names
		names=FASTA[name_locations]
		
		#Binds names and sequences into a split matrix
		cbind(names, sequences)->FASTA
		
		#Calculates substitute integer names for each sequence
		if(add_integer_column==TRUE){
			integer_column=1:length(name_locations)
			cbind(FASTA, integer_column)->FASTA
		}
		
		#Remove the carrot from the name field
		if(leave_carrot==FALSE){
			gsub(">", "", FASTA[,1])->FASTA[,1]
		}
		gsub("U", "T", FASTA[,2])->FASTA[,2]
		gsub(" ", "", FASTA[,2])->FASTA[,2]
		return(FASTA)
	}
	if(reverse==TRUE){
	unsplit=vector(length=2*nrow(FASTA))
		if(length(grep(">", FASTA[,1], fixed=TRUE))==0){
			cbind(rep(">",nrow(FASTA)),FASTA[,1])->carrot
			apply(carrot, 1, paste, collapse="")->FASTA[,1]
		}
	FASTA[,1]->unsplit[seq(1,(2*nrow(FASTA)-1),by=2)]
	FASTA[,2]->unsplit[seq(2,(2*nrow(FASTA)), by=2)]
	return(unsplit)
	}
}






fastq_to_fasta<-function(fastq){
	as.matrix(fastq)->fastq
	seq(1, nrow(fastq), by=4)->names
	seq(1, nrow(fastq), by=4)+1->seqs
	fastq[sort(c(names, seqs)),]->fastq
	gsub("@",">",fastq)->fastq
	return(fastq)
}
	
	
	library(ggplot2)
library(gridExtra)

theme_black = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "white",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      #panel.grid.major = element_line(color = "grey35"),  
      #panel.grid.minor = element_line(color = "grey20"),  
      panel.margin = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}

colorscheme<-function(color_key){
	n <- 60
	qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
	col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
	DerekPalette<-c("#000000","#ffff00","#f922b9","#32fff4","#1F78B4","#FB9A99","#FDBF6F","#E31A1C","#6A3D9A","#CAB2D6","#FF7F00","#FFFF99","#A6CEE3")
	#DerekPalette<-c("#ffff00","#000000","#ffff00","#f922b9","#000000","#ffff00","#f922b9","#32fff4","#1F78B4","#FB9A99","#FDBF6F","#E31A1C","#6A3D9A","#CAB2D6","#FF7F00","#FFFF99","#A6CEE3")
	DerekPalette<-c("#f922b9","#32fff4","#1F78B4","#FB9A99","#FDBF6F","#E31A1C","#6A3D9A","#CAB2D6","#FF7F00","#FFFF99","#A6CEE3")
		
	colorslist=read.table(color_key,  sep="\t",comment.char="")
	as.matrix(colorslist)->colorslist
	colorslist[1,]->colnames(colorslist)
	
	c(DerekPalette,DerekPalette,col_vector,col_vector,col_vector)->DerekPalette
	c(col_vector,col_vector,col_vector)->DerekPalette
	rep(DerekPalette, 10)->DerekPalette
	#set standard color scheme
	taxa_color_pairs[length(taxa_color_pairs):1]->taxa_color_pairs
	cbind(taxa_color_pairs, DerekPalette[1:length(taxa_color_pairs)], rep(0.3, length(taxa_color_pairs)))->taxa_color_pairs
	
	
	
	for(r in 1:nrow(colorslist)){
		taxa_color_pairs[which(taxa_color_pairs[,1]==colorslist[r,2]),2]<-colorslist[r,COLUMN_IN_COLORSLIST]
	}
	
	
	if(exists("ALPHA_COLUMN_IN_COLORSLIST")==FALSE){ALPHA_COLUMN_IN_COLORSLIST=8}
	for(r in 1:nrow(colorslist)){
		taxa_color_pairs[which(taxa_color_pairs[,1]==colorslist[r,2]),3]<-colorslist[r,ALPHA_COLUMN_IN_COLORSLIST]
	}
}
	
	
