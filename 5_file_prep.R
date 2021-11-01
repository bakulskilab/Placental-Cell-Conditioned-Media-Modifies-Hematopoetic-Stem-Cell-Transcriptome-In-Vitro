##############################
# packages needed,
# install if don't have
##############################
install.packages('refGenome')
install.packages('png')

install.packages("http://hartleys.github.io/QoRTs/QoRTs_STABLE.tar.gz",
                 repos=NULL, 
                 type="source");

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt", version = "3.8")

##############################
# meta data
##############################
BeWo <- c('minus','plus','plus','minus','plus','plus','minus','plus','plus','minus','plus','plus')
vc <- c('vc','vc','h10','vc','vc','h10','vc','vc','h10','vc','vc','h10')
date <- c(rep('041118',3),rep('041718',3),rep('042418',3),rep('050118',3))
folders <- list.dirs(path='/nfs/turbo/bakulski1/People/johndou/BeWo_RNAseq/Run_2567_Output/featureCounts',recursive=FALSE)
samp.names <- gsub('/nfs/turbo/bakulski1/People/johndou/BeWo_RNAseq/Run_2567_Output/featureCounts/','',folders)

meta <- data.frame(unique.ID=samp.names,group.ID=paste0(vc,'_',BeWo,'BeWo'),day=date)
meta$group.ID <- relevel(meta$group.ID, ref='vc_minusBeWo')

###############################
#QoRTs plots
###############################

#done after running the QoRTs program on STAR output
library(QoRTs)
library(png)

res <- read.qc.results.data('/nfs/turbo/bakulski1/People/johndou/BeWo_RNAseq/Run_2567_Output/QoRTs_QC/', decoder=meta)

makeMultiPlot.all(res, outfile.dir='/nfs/turbo/bakulski1/People/johndou/BeWo_RNAseq/Run_2567_Output/QoRTs_QC/',plot.device.name='pdf')


###############################
#feature counts summaries
###############################

#get all the file paths
folders <- list.dirs(path='/nfs/turbo/bakulski1/People/johndou/BeWo_RNAseq/Run_2567_Output/featureCounts',recursive=FALSE)
summary.paths <- file.path(folders,'feature_counts.summary')

#read first summary in to guide looping
sum.example <- read.delim(summary.paths[1])

#blank dataframe to fill in
fc.summary <- data.frame(matrix(0,nrow=length(folders),ncol=nrow(sum.example)))
rownames(fc.summary) <- sum.example$Status
samp.names <- gsub('/nfs/turbo/bakulski1/People/johndou/BeWo_RNAseq/Run_2567_Output/featureCounts/','',folders)
colnames(fc.summary) <- samp.names

#fill in info for each sample
for(i in 1:length(summary.paths)){
  summary <- read.delim(summary.paths[i])
  fc.summary[,i] <- summary$Aligned.out.sorted
}

#save
write.csv(fc.summary,file='/nfs/turbo/bakulski1/People/johndou/BeWo_RNAseq/Run_2567_Output/featureCounts/feature counts summary.csv')
fc.summary <- read.csv('/nfs/turbo/bakulski1/People/johndou/BeWo_RNAseq/Run_2567_Output/featureCounts/feature counts summary.csv', row.names=1)

###############################
### load feature counts files
###############################

setwd('/nfs/turbo/bakulski1/People/johndou/BeWo_RNAseq/')

#get all the file paths
folders <- list.dirs(path='/nfs/turbo/bakulski1/People/johndou/BeWo_RNAseq/Run_2567_Output/featureCounts',recursive=FALSE)
counts.paths <- file.path(folders,'feature_counts')

#read in
counts <- lapply(counts.paths, function(fn) read.table(file.path(fn), skip=2))

#check to make sure gene ID column all same, if this if false something is wrong
all(sapply(counts, function(a) all(a$V1 == counts[[1]]$V1))) 

#make a matrix for counts
tbl.counts <- sapply(counts, function(x) x$V7)
colnames(tbl.counts) <- samp.names
rownames(tbl.counts) <- counts[[1]]$V1
head(tbl.counts)

############################### 
#gene identifiers:
#in count matrix ensembl ids are used
#using biomart, matrix to match to
#entrez ids and hgnc names is made
###############################

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ens.names <- rownames(tbl.counts)
ens.names <- gsub('\\..*','',ens.names)
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id","entrezgene","hgnc_symbol","chromosome_name","start_position","end_position"),
  values=ens.names,
  mart=mart)
genes.ord <- genes[match(genes$ensembl_gene_id,ens.names),]
genes.ord.uni <- genes.ord[!duplicated(genes.ord$ensembl_gene_id),]

#use feature counts output as base, add info
gene.info <- read.delim("/nfs/turbo/bakulski1/People/johndou/BeWo_RNAseq/Run_2567_Output/featureCounts/Sample_119411/feature_counts",skip=1)
gene.info <- gene.info[,-7]
rownames(gene.info) <- gene.info$Geneid
gene.info$Geneid <- gsub('\\..*','',gene.info$Geneid)
gene.info$Row.names <- rownames(gene.info)
gene.info <- merge(gene.info,genes.ord.uni,by.x='Geneid',by.y='ensembl_gene_id',all.x=T,all.y=F)
rownames(gene.info) <- gene.info$Row.names
gene.info <- gene.info[rownames(DGElist$counts),]
gene.info <- gene.info[,c('Geneid','entrezgene','hgnc_symbol','Chr','Start','End','Strand','Length')]
names(gene.info)[1] <- 'ensembl_gene_id'

head(gene.info)
