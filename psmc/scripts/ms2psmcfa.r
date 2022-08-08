# By Shenglin Liu, Jun 2, 2015.
# This function randomly samples RAD loci from chromosomes simulated by ms (Hudson 2002),
#	and transforms the RAD data into a ".psmcfa" file for PSMC analysis.
# The input file is the simulation result from ms; it must have an extension name of ".ms".
# If the number of repeats in the input file (say n) is larger than n.chrom,
#	the function will generate n%/%n.chrom files,
#	each having n.chrom chromosomes.

# file: a file containing simulated chromosomes created by ms;
#	it must have an extension name of ".ms".
# l.chrom: the length of one chromosome.
# n.chrom: the number of chromosome pairs for the individual.
# l.bin: bin size.
# shrink: shrink the genome by "shrink" times by randomly sampling bins.
# consecutive: the number of consecutive bins taken at each RAD locus.
ms2psmcfa<-function(file,l.chrom,n.chrom=1,l.bin=100,shrink=1,consecutive=1)
{
	a<-scan(file,what="",quiet=T,sep="\n")
	a<-grep("positions:",a,value=T)
	a<-sub("positions: ","",a)
	
	n.samp<-length(a)%/%n.chrom
	
	output<-paste(sub(".ms$","",file),"_%0",nchar(n.samp),"d",sep="")
	output<-sprintf("%s.CHR%02d",output,n.chrom)
	if(shrink>1)
	{
		output<-sprintf("%s.RAD%d",output,shrink)
		if(consecutive-1)
		{
			output<-sprintf("%s.CS%d",output,consecutive)
		}
	}
	if(consecutive<1)
	{
		stop("'consecutive' must be larger than or equal to 1!")
	}
	if(shrink<1)
	{
		stop("'shrink' must be larger than or equal to 1!")
	}
	if(n.chrom<1)
	{
		stop("'n.chrom' must be larger than or equal to 1!")
	}
	
	for(i in 1:n.samp)
	{
		sink(sprintf(paste(output,".psmcfa",sep=""),i))
		for(j in 1:n.chrom)
		{
			cat(sprintf(">seq%02d\n",j))
			write(a[(i-1)*n.chrom+j],"temp")
			c<-scan("temp",quiet=T)
			c<-c*l.chrom
			c<-round(c)
			b<-rep("T",l.chrom%/%l.bin)
			c<-(c-1)%/%l.bin+1
			c<-unique(c)
			c<-c[c<(l.chrom%/%l.bin)]
			b[c]="K"
			c<-sample(1:(l.chrom%/%l.bin),l.chrom%/%l.bin%/%shrink)
			if(consecutive-1)
			{
				temp<-c
				for(k in 1:(consecutive-1))
				{
					temp<-c(temp,c+k)
				}
				c<-temp
			}
			c<-sort(c)
			c<-c[c<(l.chrom%/%l.bin)]
			b<-b[c]
			cat(b,sep="")
			cat("\n")
		}
		sink()
	}
}

