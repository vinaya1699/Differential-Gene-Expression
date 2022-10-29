# Differential-Gene-Expression

# To identify the genes which are differential in tumor vs control samples

logcpm=source("logCPM.rds")

#logcpm contains log of CPM of provided bulk cell data .
# #Reading an input cancer data file
# data_file=read.csv("C:/Users/91973/Desktop/3rd_semester/Cancer_analysis/GSE149650_Read_counts.csv",sep=",",header=T,row.names = 1)

# #Create a count per matrix
# cpmatrix=data_file
# for(i in 1:ncol(data_file)){
#   cpmatrix[,i]=(data_file[,i]/sum(data_file[,i]))*1000000
# }

# #Calculate a log of cpm
# logcpm=log2(cpmatrix+1)
![image](https://user-images.githubusercontent.com/110582335/198826384-45ef4420-c025-41b5-8e61-4246804b8de7.png)

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
