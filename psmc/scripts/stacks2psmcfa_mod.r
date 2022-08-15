# Generating input files for PSMC from tag files and allele files.
# Tag files and allele files are generated through pstacks (Stacks v1.13).
# Jan 27, 2015, by Shenglin Liu.
# Aprl 17, 2022 modified by Kumar

setwd("/nobackup/beegfs/workspace/ks575/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/psmc/data")
#============parameters============#
f.tag<-"ROC31.tags.tsv"	# name of the tag file
f.allele<-"ROC31.alleles.tsv"	# name of the allele file

f.output<-"ROC31.psmcfa"	# name of the output file
f.log<-"psmc_log.txt"	# name of the log file

chromosome<-c("Backbone", "contig", "scaffold")	# keywords for the chromosomes that you want to select
interval<-3000000	# number of lines to be read each time from the input stream
#==================================#

write(sprintf("Started\tparsing %s at %s",f.tag,Sys.time()),f.log,append=T)
#
ex<-data.frame(integer(0),character(0),integer(0))
a<-read.table(f.tag, sep="\t",comment.char="",stringsAsFactors=F, skip = 1)
index<-as.vector(unlist(sapply(chromosome,function(x){grep(x,a[,4])})))
ex<-a[index,3:5]
#
write(sprintf("Ended\tparsing %s at %s\n",f.tag,Sys.time()),f.log,append=T)
#
print(dim(ex))
#
zygosity<-rep("T",nrow(ex))
zygosity[1]<-"K"
ex<-cbind(ex,zygosity)
ex[1,4]<-"T"
rownames(ex)<-ex[,1]

b<-read.table(f.allele,sep="\t",stringsAsFactors=F)

ex[as.character(unique(b[,3])),4]<-"K"
ex<-na.omit(ex)

sink(f.output)
chr<-unique(ex[,2])
for(i.chr in chr)
{
	subex<-ex[ex[,2]==i.chr,]
	subex<-subex[order(subex[,3]),]
	cat(">",i.chr,"\n",sep="")
	cat(as.character(subex[,4]),"\n",sep="")
}
sink()
