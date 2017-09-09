#### R CODE FOR PARAFIT ANALYSIS ####


#### SET WORKING DIRECTORY
setwd("~name")


#### LOAD LIBRARIES 
library(ape)
library(phangorn)
library(parallel)


#### SET ALPHA LEVEL
p_level<-0.055


#### READ IN 'PARASITE' TREE 
p.tree<-read.nexus("name.tre")


#### REMOVE 'PARASITE' OUTGROUPS OR DUPLICATE TAXA 
parasite_tree <- drop.tip(p.tree, c('outgroup1', 'outgroup2', 'duplicate1'))


#### READ IN 'HOST' TREE
host_tree<- read.tree("name.tre") # may also use read.nexus


#### READ IN HOST-PARASITE MATRIX
hpmatrix<-as.matrix(read.table("matrix.txt", header=T))


#### CALCULATE PATRISTIC DISTANCES AND ROTATE HEADINGS TO FIT MATRIX
#### PARASITES
parasite.patristic.dist<-as.matrix(cophenetic(parasite_tree))
parasite.patristic.dist <- parasite.patristic.dist[colnames(hpmatrix), colnames(hpmatrix)]


#### HOSTS
host.patristic.dist<-as.matrix(cophenetic(host_tree))
host.patristic.dist <- host.patristic.dist[rownames(hpmatrix), rownames(hpmatrix)] 


#### RUN PARAFIT (ONE RUN)
runParaFit<- function(i){
  tmp <- parafit(host.patristic.dist, parasite.patristic.dist, hpmatrix, nperm=999, test.links=T, correction='cailliez')
  return (tmp)
  }
  
# alternative: correction='lingoes'


#### RUN PARALLEL PARAFIT ANALYSES
cat("Scheduled",length(host_tree),"jobs to run in parallel ... Running ... Nothing will be printed on screen until all runs are completed ...",sep=" ")


#### LOOP MULTIPLE PARAFIT RUNS (100 IN THIS CASE)
#### IDEAL TO RUN > 1 RUN BECAUSE THE RESULTS VARY SLIGHTLY OVER MULTIPLE RUNS
res<- lapply(1:100, function(i){
  simplify2array(mclapply(1,runParaFit,mc.preschedule=T,mc.cores=10,mc.set.seed=F))
})

res # THIS SHOWS ALL RUNS AND P.GLOBAL VALUES (ALL SEPARATE 'OBJECTS' FROM ONE ANOTHER)

data.class(res) # SHOULD BE 'LIST'


#### EXTRACT P.GLOBAL VALUES FROM 'res'
res_p<-sapply(res, "[[",2) 
# res_p # TO SEE OUTPUT
data.class(res_p) # NUMERIC
as.matrix(res_p) # MAKE MATRIX FOR NEXT STEP


#### NOW, YOU CAN GET THE AVERAGE P.GLOBAL VALUE
p_global_mean<-mean(res_p[[1]])
p_global_mean 


#### NOW, NEED TO COMBINE ALL THE PARAFIT RUNS AND P.GLOBAL VALUES INTO 1 'OBJECT' IN ORDER TO EXTRACT P-VALUES (WHICH THEN CAN BE USED TO GET THE (OVERALL) MEAN P.GLOBAL AS WELL AS THE MEAN P-VALUES FOR EACH LINK OVER MULTIPLE PARAFIT RUNS). ALSO THIS IS USEFUL FOR MAKING SURE THE ORDER THAT IS LISTED IN THE PARAFIT OUTPUT IS THE SAME AS THE ORDER THE P-VALUES APPEAR (BECAUSE R LIKES TO REORDER THESE OBJECTS IN STRANGE WAYS... WHICH CAN LEAD TO THE WRONG LINKS BEING ASSOCIATED WITH 'SIGNIFICANCE')
res_link<-c()
for (i in 1:100){ 
  res_link<-cbind(res_link,res[[i]])
}
# res_link # TO SEE OUTPUT
data.class(res_link) # MATRIX, ALL VALUES FROM 'res' ARE COMBINED INTO 'res_link'


#### CHECK ORDER OF PARAFIT OUTPUT AND P-VALUES (VERY IMPORTANT!)
links_parafit_order=cbind(rownames(hpmatrix)[res_link[,1]$link.table[,1]],colnames(hpmatrix)[res_link[,1]$link.table[,2]])
# links_parafit_order # TO CHECK IF R CORRECTLY ORDERED OUTPUT/LINKS/P-VALUES



#### CREATE MATRIX OF P VALUES; COLUMNS ARE HP LINKS, ROWS ARE DIFFERENT PARAFIT RUNS [CERTAINLY, INSTEAD OF THE LOOP, ONE CAN USE LAPPLY ETC HERE]
p_table=NULL 
for (i in 1:ncol(res_link)){ 
  p.F1=res_link[,i]$link.table[,4] #col4= p.F1, col6=p.F2
  p_table=rbind(p_table,p.F1)
  #cat(p.F1, "***\n\n")
}
p_table


#### ADJUST THOSE P-VALUES; COLUMNS ARE HP LINKS, ROWS ARE DIFFERENT PARAFIT RUNS
p_table_adj<-apply(p_table,2,p.adjust,method="BH") 
# p_table_adj # TO SEE OUTPUT
data.class(p_table_adj) # MATRIX


#### ADJUSTED P-VALUE MEANS FOR EACH INDIVIDUAL LINK
p_means<-colMeans(p_table_adj)
p_means_table<-cbind(links_parafit_order,(p_means))
p_means_table_df<-as.data.frame(p_means_table)
p_means_table_df
write.table(p_means_table_df, file='name.txt',sep='\t', row.names=F, col.names=T,quote=F)


#### WRITING P-VALUE TABLES
ptable<-cbind(links_parafit_order,t(p_table))
colnames(ptable)<-c("Host","Parasite",sprintf("p.F1_tree_%s",seq(1:ncol(res_link))))
ptable
write.table(ptable, file='name.txt', sep='\t', row.names=F, col.names=T,quote=F)


#### WRITING P-VALUE TABLES WITH ADJUSTED VALUES
p_tableadj<-cbind(links_parafit_order,t(p_table_adj))
colnames(p_tableadj)<-c("Host","Parasite",sprintf("p.F1.adj_tree_%s",seq(1:ncol(res_link))))
write.table(p_tableadj, file='name.txt', sep='\t', row.names=F, col.names=T,quote=F)


#### PLOTTING
#### ADJUST COLORS AND LINE TYPES
p_colors=NULL

for (i in 1:ncol(p_table_adj)){ 
  p_min=min(p_table_adj[,i]);p_max=max(p_table_adj[,i])
  color="red3" # MIX OF SIGNIFICANT & NOT-SIGNIFICANT P-VALUES
  if(p_min>p_level && p_max>p_level){color="gray55"} # NOT-SIGNIFICANT
  if(p_min<=p_level && p_max<=p_level){color="green3"} # SIGNIFICANT
  p_colors[i]=color
}

p_type=NULL

for (i in 1:ncol(p_table_adj)){ 
  p_min=min(p_table_adj[,i]);p_max=max(p_table_adj[,i])
  lty="dotdash" # MIX OF SIGNIFICANT & NOT-SIGNIFICANT P-VALUES
  if(p_min>p_level && p_max>p_level){lty="longdash"} # NOT-SIGNIFICANT
  if(p_min<=p_level && p_max<=p_level){lty="solid"} # SIGNIFICANT
  p_type[i]=lty
}

#### TANGLEGRAM TO MINIMIZE LINK CROSSING
library(phytools) # DOWNLOADED FROM GITHUB (NOT CRAN) https://github.com/liamrevell/phytools

ph1 <- host_tree
ph2 <- parasite_tree
assoc_links <- links_parafit_order
cophy <- cophylo(ph1, ph2, assoc=assoc_links, rotate=T)

pdf("name.pdf", height = 30, width = 45) # sets size of printing area 
cophyloplot_min<-plot(cophy, link.lty=p_type, link.col=p_colors, link.lwd=8, fsize=3.5,pts=F)
dev.off()

#### IF LABEL NAMES ARE TOO LONG
# EMMANUEL PARADIS AT http://grokbase.com/t/r/r-sig-phylo/12c3djkx3h/modification-of-figures-produced-using-cophyloplot
# FIRST MAKE A COPY OF THE INTERNAL FUNCTION
# titi<-plotCophylo2
# then do fix(titi) and modify what you want (here look for calls to
# text(... cex = ...), save and close. Do fix(cophyloplot) and change
# "plotCophylo2(..." by "titi(...", save and close.
# fix(titi) # 
# if (show.tip.label) {
#   text(a[1:N.tip.x, ], cex = 0.5, font = font, pos = 4, labels = x$tip.label)
#   text(b2[1:N.tip.y, ], cex = 0.5, font = font, pos = 2, 
#        labels = y$tip.label)
# }
#
# titi
# fix(cophyloplot) # Do fix(cophyloplot) and change "plotCophylo2(..." by "titi(...", save and close.
