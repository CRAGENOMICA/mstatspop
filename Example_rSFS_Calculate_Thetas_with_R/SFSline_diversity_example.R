#functions to estimate diversity using the rSFS (option -o 92 in mstatspop) in relation to the total population

#structure of the input data: infile:	./chr16.DUROC_100.RAW.ID.tfa.gz	scaffold_name:	16	start_window:	250212	end_window:	322075	missing:	0	iteration:	0	npermutations:	0	seed:	123456	Length:	71864.00	Lengtht:	71864	mh:	0	Ratio_S/V:	2.082	Ratio_Missing:	NA	Variants:	1279	npops:	1	nsam[0]:	100	Eff_length1_pop[0]:	71864.00	Eff_length2_pop[0]:	71864.00	Eff_length3_pop[0]:	71864.00	
#linexfreqy_pop[0]: .. HERE there is a vector divided in line (rows) x freqs (columns): in this case is 100 lines and 50 diff folded frequencies (no fixed included)
#field 16 (nsam[0]:) and field 20 (linexfreqy_pop[0]:) are necessary. to estimate diversity per position the field 18 (Eff_length2_pop[0]:) is also necessary.

#Functions to calculate Theta Watterson (1979) Theta Tajima (1983) and Theta FuLi (1993) from folded SFS using the Achaz approach (2009)
deltak <- function(i,j) {
  if(i==j) return(1)
  return(0)
}

weight.watt.folded <-function(nsam) {
  fold <- 1 #sfs folded
  w <- array(0,dim=c(floor(nsam/(fold+1))))
  for(i in 1:length(w)) {
    w[i] <- nsam/(i*(nsam-i)*(1+deltak(i,nsam-i)))
  }
  w
}

weight.taj.folded <-function(nsam) {
  fold <- 1 #sfs folded
  w <- array(0,dim=c(floor(nsam/(fold+1))))
  for(i in 1:length(w)) {
    w[i] <- nsam/(1+deltak(i,nsam-i))
  }
  w
}

weight.fuli.folded <-function(nsam) {
  fold <- 1 #sfs folded
  w <- array(0,dim=c(floor(nsam/(fold+1))))
  w[1] <- nsam
  w
}

psi.fold.ij <- function(nsam,subnsam) {
  w <- array(0,dim=c(floor(nsam/2)))
  for(i in 1:length(w)) { 
    w[i] <- (1-dhyper(x=0,k=subnsam,m=i,n=nsam-i))
  }
  w
}

Calc.Theta.folded <- function(sfs,w,phi,psi) {
  th <- 0
  for(i in 1:length(sfs)) {
    th <- th + w[i] * sfs[i] * 1/phi[i] * 1/psi[i]
  }
  th <- th/(sum(w))
}

#EXAMPLE HERE:
nsam <- 200
npops <- 100
path_file <- "./"

#ESTIMATE RELATIVE NUCLEOTIDE DIVERSITY PER POP (ind) (that is, IN RELATION TO TOTAL, option -o 92)
#define data frames for sending results

#take the genes and their locations
name.list <- "MICROBIOTE_BILEACID_genes"
name.gene <- read.table(file=sprintf("%s/GCF_000003025.6_Sscrofa11.1_%s_chromosome18.txt",path_file,name.list))
name.gene <- name.gene[!duplicated(name.gene),]

columns.fuli <- c("NAME","nsam","scaffold","Gene","Start","End","Length",sprintf("ThetaRFuLi_%0.f",c(1:npops)))
data.results.Fuli <- data.frame(matrix(nrow=0,ncol=length(columns.fuli)))
colnames(data.results.Fuli) <- columns.fuli

columns.watt <- c("NAME","nsam","scaffold","Gene","Start","End","Length",sprintf("ThetaRWatt_%0.f",c(1:npops)))
data.results.Watt <- data.frame(matrix(nrow=0,ncol=length(columns.watt)))
colnames(data.results.Watt) <- columns.watt

columns.taj  <- c("NAME","nsam","scaffold","Gene","Start","End","Length",sprintf("ThetaRTaj_%0.f",c(1:npops)))
data.results.Taj  <- data.frame(matrix(nrow=0,ncol=length(columns.taj)))
colnames(data.results.Taj) <- columns.taj

#for each chromosome (1:18)
data.freqs <- read.table(file="chr18.DUROC_100.results_mstatspop_o92_MICROBIOTE_BILEACID_genes.txt")

#for each gene extract the sfs per lineage and calculate their theta statistics
init.mat <- grep("rSFS",data.freqs[1,])
num.gene <- 1 #test
for(num.gene in 1:dim(data.freqs)[1]) {
  p.len.seq <- grep("Length:",data.freqs[num.gene,]) + 1
  sfs.line <- matrix(as.numeric(data.freqs[num.gene,((init.mat+1):(init.mat+(npops*floor(nsam/2))))]),byrow=T,nrow=npops,ncol=(floor(nsam/2)))
  
  w.fuli <- weight.fuli.folded(nsam)
  w.watt <- weight.watt.folded(nsam)
  w.taj  <- weight.taj.folded(nsam)
  
  phi <- w.watt
  w.psi <- psi.fold.ij(nsam,subnsam=2)
  
  Theta.FuLi.line <- round(apply(sfs.line,1,Calc.Theta.folded,w=w.fuli,phi=phi,psi=w.psi)/data.freqs[num.gene,p.len.seq],5)
  Theta.Watt.line <- round(apply(sfs.line,1,Calc.Theta.folded,w=w.watt,phi=phi,psi=w.psi)/data.freqs[num.gene,p.len.seq],5)
  Theta.Taji.line <- round(apply(sfs.line,1,Calc.Theta.folded,w=w.taj ,phi=phi,psi=w.psi)/data.freqs[num.gene,p.len.seq],5)
  
  ngene <- name.gene[name.gene[,1]==data.freqs[num.gene,c(4)] & name.gene[,2]==data.freqs[num.gene,c(6)] & name.gene[,3]==data.freqs[num.gene,c(8)],4]
  
  data.results.Fuli <- rbind(data.results.Fuli,data.frame(name.list,nsam,(data.freqs[num.gene,c(4)]),ngene,(data.freqs[num.gene,c(6)]),(data.freqs[num.gene,c(8)]),data.freqs[num.gene,p.len.seq],t(Theta.FuLi.line)))
  data.results.Watt <- rbind(data.results.Watt,data.frame(name.list,nsam,(data.freqs[num.gene,c(4)]),ngene,(data.freqs[num.gene,c(6)]),(data.freqs[num.gene,c(8)]),data.freqs[num.gene,p.len.seq],t(Theta.Watt.line)))
  data.results.Taj  <- rbind(data.results.Taj, data.frame(name.list,nsam,(data.freqs[num.gene,c(4)]),ngene,(data.freqs[num.gene,c(6)]),(data.freqs[num.gene,c(8)]),data.freqs[num.gene,p.len.seq],t(Theta.Taji.line)))
}
colnames(data.results.Fuli) <- columns.fuli
colnames(data.results.Watt) <- columns.watt
colnames(data.results.Taj)  <- columns.taj

write.table(file=sprintf("%s/Results_RDiversity_%s_RThetaFuli.txt",path_file,name.list),x=data.results.Fuli,quote=F,row.names=F)
write.table(file=sprintf("%s/Results_RDiversity_%s_RThetaWatt.txt",path_file,name.list),x=data.results.Watt,quote=F,row.names=F)
write.table(file=sprintf("%s/Results_RDiversity_%s_RThetaTaji.txt",path_file,name.list),x=data.results.Taj, quote=F,row.names=F)
