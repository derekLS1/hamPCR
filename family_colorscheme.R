n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
DerekPalette<-c("#000000","#ffff00","#f922b9","#32fff4","#1F78B4","#FB9A99","#FDBF6F","#E31A1C","#6A3D9A","#CAB2D6","#FF7F00","#FFFF99","#A6CEE3")
#DerekPalette<-c("#ffff00","#000000","#ffff00","#f922b9","#000000","#ffff00","#f922b9","#32fff4","#1F78B4","#FB9A99","#FDBF6F","#E31A1C","#6A3D9A","#CAB2D6","#FF7F00","#FFFF99","#A6CEE3")
DerekPalette<-c("#f922b9","#32fff4","#1F78B4","#FB9A99","#FDBF6F","#E31A1C","#6A3D9A","#CAB2D6","#FF7F00","#FFFF99","#A6CEE3")


#to black out most stuff use these 2 lines
#DerekPalette<-c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000")
#col_vector<-c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000")



if(file.exists("/Users/dlundberg/Documents/abt6/Regalado2017/compare16StoMETA/fungal_colors.txt")==TRUE){
	colorslist=read.table("/Users/dlundberg/Documents/abt6/Regalado2017/compare16StoMETA/fungal_colors.txt", sep="\t", comment.char="")
}
if(file.exists("/Users/dlundberg/Documents/abt6/Regalado2017/compare16StoMETA/bacterial_colors.txt")==TRUE){
	colorslist=read.table("/Users/dlundberg/Documents/abt6/Regalado2017/compare16StoMETA/bacterial_colors.txt", sep="\t",comment.char="")
}
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

