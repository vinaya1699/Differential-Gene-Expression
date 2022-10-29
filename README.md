# Differential-Gene-Expression

# To identify the genes which are differential in tumor vs control samples

logcpm=source("logCPM.rds")
![image](https://user-images.githubusercontent.com/110582335/198826145-5d8a97b2-a13a-4767-8a0f-695a4eaeeecb.png)
![image](https://user-images.githubusercontent.com/110582335/198826296-e5510a2d-2fd8-4f56-905f-b719d43913d3.png)

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
