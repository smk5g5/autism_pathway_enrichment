library(dplyr)
library(KEGGREST)
#load("./kegg2entrez.RData")
load("C:/Users/saadkhan/kegg2entrez.RData")
#######convert gene to pathway matrix############
autism_gene_list <- read.table("case_control.pairs.large.txt",header = T,sep = "\t")
autism_gene_list <- autism_gene_list[complete.cases(autism_gene_list),] # remove NA values
autism_gene_list <- autism_gene_list %>% 
    filter(! Entrez1 == Entrez2 ) %>%
    mutate(G1 = pmin(Entrez1, Entrez2), G2 = pmax(Entrez1, Entrez2)) %>%
    select(G1, G2, max.g) %>%
    rename(Entrez1 = G1, Entrez2 = G2) %>%
    group_by(Entrez1, Entrez2) %>%
    summarize(mg = max(max.g)) %>%
    mutate(max.g = mg) %>% select(Entrez1, Entrez2, max.g)
autism_gene_list <- as.data.frame(autism_gene_list)

uni_gene <- unique(keggpathway2gene$gene_id) #uni gene should be all the genes in pathways (universal genes)
geneid <- union(autism_gene_list$Entrez1,autism_gene_list$Entrez2) #gene id should be all unique genes in autism results 
kegg_pathway_matrix <- as.matrix(xtabs(~pathway_id+gene_id,data=keggpathway2gene))
SigPairEnrich=function(autism_gene_list,geneid,uni_gene,kegg_pathway_matrix,FDR=0.05,fdrmethod=c('bonferroni')){
#---------------------------------------------------------
###--Sort gene in pathway
###
s=as.matrix(order(uni_gene)); #ordered indices of gene ids
UniG=as.matrix(uni_gene[s]); # sorted matrix of uni_gene 
Path=as.matrix(kegg_pathway_matrix[,s]);

#----------------------------------
# Profile Gene in pathway
IntG=as.matrix(sort(intersect(geneid,UniG)));#The genes involved in the pathway which contained in geneid #353
#---------------------------------------
# pathway gene in Profile
a=UniG %in% IntG; #IntG is all unique genes from autism results that have KEGG pathway assigned
Path=Path[,a];#IntG genes in profile
L=length(IntG);#background genes
m=autism_gene_list[,1] %in% IntG; #order may be important
## Number of gene that intersect with intg that are there in entrez 1
n=autism_gene_list[,2] %in% IntG; #
## Number of gene that intersect with intg that are there in entrez 2
k=(m+n)==2;
##Gene pairs that are present in the list of pathways at the same time 
L_sigP=sum(k);#Gene pairs involved in pathway,(num of interest pairs)
if(L_sigP==0){
    stop("There is no interest pairs involved in pathway");
}
BP=choose(L,2);#background pairs
#BP is the number of all possible combinations for pairwise background genes
#BP=nC2
##-----------------------------------------------------------
#--SigP enrich

LP=length(Path[,1]);#path num #total number of pathways
multiLP = choose(LP,2) #Total number of pairwise pathway combinations
HyperP=matrix(,nrow=multiLP,ncol=8);#"pathway_index","p_value","adjusted p-value","k","n","m".
#create a hypergeometric distribution matrix with length LP(total number of pathways)
HyperP[,3]=L_sigP; #(Column 3 is gene pairs involved in one of the LP pathways)
HyperP[,4]=BP; #(#Column 4 all possible combination of background genes)

#combinations_path <- rownames(Path)
#pairwise_comb <- as.list(combn(combinations_path,2,simplify = F))
multiLP = length(pairwise_comb)
for(i in 1:multiLP){
   pathnames_vec = unlist(pairwise_comb[i]) 
    Temp1 = Path[pathnames_vec[1],] 
    Temp2=Path[pathnames_vec[2],]
    TempG1=IntG[which(Temp1==1),]; #All the genes in pathway1
    TempG2= IntG[which(Temp2==1),] #All the genes in pathway2
    TempG2_complement <- TempG2[(TempG1!=TempG2)]
    TempG2_complement <- TempG2_complement[is.na(TempG2_complement)]      
    TempG1_complement <- TempG1[(TempG1!=TempG2)] 
    Total_genes = c(TempG1_complement,TempG2_complement) #I can still have all the genes
    HyperP[i,1]=i;
    HyperP[i,7]=pathnames_vec[1];
    HyperP[i,8]=pathnames_vec[2];
  if ((length(TempG1_complement)<1)&(length(TempG2_complement)<1)){
    HyperP[i,2]=1;
    HyperP[i,5]=0;
    HyperP[i,6]=0;
    next;
  }else{
    BPinPath=choose(length(Total_genes),2); #this is x which is all possible combinations of pathways that may contain these genes
    HyperP[i,6]=BPinPath; #
    a=autism_gene_list[,1] %in% TempG1_complement;
    b=autism_gene_list[,2] %in% TempG2_complement;
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
}
pathway_info <- keggList("pathway","hsa") #Getting Kegg pathway names using KEGGREST
selpathnames <- as.character(paste("path",rownames(Path),sep=":"))
pathway_info <- pathway_info[selpathnames]
q1=as.matrix(p.adjust(HyperP[,2],method="BH",length(HyperP[,2])),ncol=1);
SigPath=cbind(HyperP[,1:2],q1,HyperP[,3:6]);
q1fdr=(q1<=FDR);
SigPath_info=pathway_info[q1fdr];
SigPath_info <- unname(SigPath_info)
SigPath=SigPath[q1fdr,];
if(length(SigPath[,1])==0){
    print("Can't find the pathway for statistically significant enrichment under the threshold value FDR!\n");
}
else{
res=as.data.frame(cbind(SigPath_info,SigPath));
colnames(res)=c("pathway_name","pathway_index","p_value","adjusted p-value","k","n","x","m")
return(res);
}
}
