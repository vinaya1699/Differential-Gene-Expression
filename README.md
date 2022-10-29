# Differential-Gene-Expression

# To identify the genes which are differential in tumor vs control samples

logcpm=source("logCPM.rds")

#logcpm contains log of CPM of provided bulk cell data .

![image](https://user-images.githubusercontent.com/110582335/198826384-45ef4420-c025-41b5-8e61-4246804b8de7.png)

# Summary of data used i.e. logcpm 
![image](https://user-images.githubusercontent.com/110582335/198827084-553700f5-ecc4-491d-8716-9b3d9c052f29.png)

# Creating a dataframe
mat=matrix(NA,ncol=4,nrow = nrow(logcpm))
rownames(mat)=rownames(logcpm)
colnames(mat)=c('meanTumor','meanControl','pvalue','log2FC')

# Column 1 to 7 contain tumorous genes whereas 8 to 11 contain control gene samples
for(i in 1:nrow(logcpm)){
  vector1 = as.numeric(logcpm[i, 1:7])
  
  vector2 = as.numeric(logcpm[i, 8:11])
  
  res=t.test(vector1, vector2, paired = F, alternative = "two.sided")
  
  mat[i,1]=res$estimate[[1]]
  
  mat[i,2]=res$estimate[[2]]
  
  mat[i,3]=res$p.value
  
  mat[i,4]=mat[i,1]-mat[i,2]
  
}

mat=as.data.frame(mat)
num=which(is.nan(mat$pvalue))
mat[num,'pvalue']=1
![image](https://user-images.githubusercontent.com/110582335/198826249-3ac30081-2b46-4f7d-ae8f-73d14b617b5b.png)

# Plotting a volcano plot 
library(EnhancedVolcano)
EnhancedVolcano(mat,lab = rownames(mat),x = 'log2FC' ,y ='pvalue')![volcano_plot](https://user-images.githubusercontent.com/110582335/197959094-6d2c9411-f515-419d-a0a0-7308a520442d.jpg)

# Interpretation :- 
1) Highly expressed upregulated genes :-
  ENSG00000139364
  ENSG00000013588
  ENSG00000198223
  ENSG00000121570
2) Highly expressed downregulated genes :-
  ENSG00000102804
  ENSG00000196549
  ENSG00000231298
  ENSG00000111371
  ENSG00000135821
