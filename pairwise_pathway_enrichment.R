load("C:/Users/saadkhan/kegg2entrez.RData")
#######convert gene to pathway matrix############
kegg_pathway_matrix <- as.matrix(xtabs(~pathway_id+gene_id,data=keggpathway2gene))
autism_gene_list <- read.table("case_control.pairs.large.txt",header = T,sep = "\t")
mygenes <- union(mylarge_file$Entrez1,mylarge_file$Entrez2)
geneid <- as.matrix(unique(keggpathway2gene$gene_id))
SigPairEnrich=function(autism_gene_list,geneid,uni_gene,kegg_pathway_matrix,FDR=0.05,fdrmethod=c('bonferroni')){
#---------------------------------------------------------
###--Sort gene in pathway
###
s=as.matrix(order(uni_gene));
UniG=as.matrix(uni_gene[s,]);
Path=as.matrix(path_info_mtx[,s]);

#----------------------------------
# Profile Gene in pathway
IntG=as.matrix(sort(intersect(geneid,UniG)));#The genes involved in the pathway which contained in geneid
#---------------------------------------
# pathway gene in Profile
a=UniG %in% IntG;
Path=Path[,a];#IntG genes in profile
L=length(IntG);#background genes
m=autism_gene_list[,1] %in% IntG;
n=autism_gene_list[,2] %in% IntG;

k=(m+n)==2;
L_sigP=sum(k);#Gene pairs involved in pathway,(num of interest pairs)
if(L_sigP==0){
    stop("There is no interest pairs involved in pathway");
}
BP=choose(L,2);#background pairs
##-----------------------------------------------------------
#--SigP enrich

LP=length(Path[,1]);#path num
HyperP=matrix(,nrow=LP,ncol=6);#"pathway_index","p_value","adjusted p-value","k","n","m".
HyperP[,3]=L_sigP;
HyperP[,4]=BP;

for (i in 1:LP){
  Temp1=Path[i,];
  TempG=IntG[which(Temp1==1),];
  mytemplist[[i]] <- TempG
  HyperP[i,1]=i;
  if (length(TempG)<2){
    
    HyperP[i,2]=1;
    HyperP[i,5]=0;
    HyperP[i,6]=0;
    next;
  }else{
    BPinPath=choose(length(TempG),2);
    HyperP[i,6]=BPinPath;
    a=autism_gene_list[,1] %in% TempG;
    b=autism_gene_list[,2] %in% TempG;
    c=(a+b)==2;
    if (sum(c)==0){
      HyperP[i,2]=1;
      HyperP[i,5]=0;
      
      next;
    }else{
      L_sigPPath=sum(c);
      HyperP[i,5]=L_sigPPath;
    }
  }
  p1=matrix(1-phyper(L_sigPPath-1,BPinPath,BP-BPinPath,L_sigP));
  HyperP[i,2]=p1;
}
return(HyperP)
} 
