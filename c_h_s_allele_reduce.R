#Chance vs history vs selection project
library("vegan")
library("plyr")
library("RColorBrewer")
library("ggplot2")
library("data.table")
library("dplyr")
library("reshape2")

theme_set(theme_bw())

# Save all filenames in directory ending in .txt as a variable temp
temp = list.files(pattern="*.txt")
all_snp <- lapply(temp,read.delim) #create a list with each object in list a data frame of mutation table
# headers we want to keep
headers <- c("aa_new_seq","aa_position","aa_ref_seq","codon_new_seq","codon_number","codon_position","codon_ref_seq","gene_name","gene_position","gene_product","frequency","mutation_category","new_seq","position","seq_id","size","snp_type","type","title")
#take a look at the first data frame in the list
df_snp_01 <- all_snp[[1]]
df_snp_01 <- df_snp_01[,headers] #select only headers we want to keep
head(df_snp_01)
# convert all data frames in the list to one large data frame
all_snp_df <- as.data.frame(data.table::rbindlist(all_snp,fill=T),colClasses = c("character"))
all_snp_df <- all_snp_df[,headers] #keep headers we want
all_snp_df$mutation <- paste(all_snp_df$aa_ref_seq,all_snp_df$aa_position,all_snp_df$aa_new_seq, sep=":") #create mutation column
View(all_snp_df)

all_snp_df$title <- gsub("sample", "breseq",all_snp_df$title)
all_snp_df$sample <- all_snp_df$title

system("iconv -c -f utf-8 -t ascii sample_key.csv > sample_key_noutf.csv")

sample_key <- read.csv("sample_key_noutf.csv",header=T)
rownames(sample_key) <- sample_key$sample
all_snp_df <- (merge(all_snp_df,sample_key,by="sample"))

write.csv(all_snp_df,file="df_snp_noref.csv")



##### reruns ####
system("iconv -c -f utf-8 -t ascii c-h-a_rerun_snp_Breseq_Output.csv > c-h-a_rerun_snp_Breseq_Output_noutf.csv")
system("iconv -c -f utf-8 -t ascii  c-h-a_rerun_mc_Breseq_Output.csv >c-h-a_rerun_mc_Breseq_Output_noutf.csv")
system("iconv -c -f utf-8 -t ascii  c-h-a_rerun_nje_Breseq_Output.csv >c-h-a_rerun_nje_Breseq_Output_noutf.csv")

df_snp <- read.csv("c-h-a_rerun_snp_Breseq_Output_noutf.csv",header=TRUE)
df_mc <- read.csv("c-h-a_rerun_mc_Breseq_Output_noutf.csv",header=TRUE)
df_nje <- read.csv("c-h-a_rerun_nje_Breseq_Output_noutf.csv",header=TRUE)

df_snp$Position <- gsub(",","",df_snp$Position)
df_snp$chr_pos <- paste(df_snp$Seq.ID,df_snp$Position,sep="_")
#Mutation data frames found in common reference = false positives
df_1_10ref_snps <- read.csv("/Users/chrismarshall/Documents/Pitt/Cooper_Lab/Abaum_UO1/experimental/references_breseq_out/final_ref_NZ_CP012004_pAb123/snp_final_ref.csv", header=TRUE)
df_1_10ref_snps$Position <- gsub(",","",df_1_10ref_snps$Position)
df_1_10ref_snps$chr_pos <- paste(df_1_10ref_snps$Seq.ID,df_1_10ref_snps$Position,sep="_")
#now do same with gene level instead of position level (more aggressive, could miss things)
df_snp_noref <- df_snp[ !(df_snp$chr_pos %in% df_1_10ref_snps$chr_pos), ]

write.csv(df_snp_noref,file="c-h-a_rerun_snp_Breseq_Output_noref.csv")


df_nje$gene_annot <- paste(df_nje$Gene,df_nje$Annotation, sep=":") 

df_1_10ref_nje <- read.csv("/Users/chrismarshall/Documents/Pitt/Cooper_Lab/Abaum_UO1/experimental/references_breseq_out/final_ref_NZ_CP012004_pAb123/nje_final_ref_Breseq_Output_noutf.csv", header=TRUE)
df_1_10ref_nje$gene_annot <- paste(df_1_10ref_nje$Gene, df_1_10ref_nje$Annotation,sep=":") 

df_nje_noref <- df_nje[ !(df_nje$gene_annot %in% df_1_10ref_nje$gene_annot), ]
write.csv(df_nje_noref,file="c-h-a_rerun_nje_noref.csv",row.names=F)


 
df_1_10ref_mc <- read.csv("/Users/chrismarshall/Documents/Pitt/Cooper_Lab/Abaum_UO1/experimental/references_breseq_out/final_ref_NZ_CP012004_pAb123/mc_final_ref_Breseq_Output_noutf.csv", header=TRUE)

df_mc_noref <- df_mc[ !(df_mc$Gene%in% df_1_10ref_nje$Gene), ]
View(df_mc_noref)


### rerun 2 ####

#install the packages if necessary
if(!require("tidyverse")) install.packages("tidyverse")
if(!require("fs")) install.packages("fs")
if(!require("readxl")) install.packages("readxl")

#load packages
library(tidyverse)
library(fs)
library(readxl)

snp <- read_excel("c_h_s_rerun2_Breseq_Output.xlsx",1) #imports as tibble
mc <- read_excel("c_h_s_rerun2_Breseq_Output.xlsx",2)
nje <- read_excel("c_h_s_rerun2_Breseq_Output.xlsx",3,trim_ws = TRUE)

sample_key <- read_csv("sample_key_rerun2.csv",col_names = F)
colnames(sample_key) <- c("Sample","names")
sample_key$Sample <- gsub('3','breseq_3',sample_key$Sample)

snp <- snp %>% inner_join(sample_key,by="Sample")

snp$Position <- gsub(",","",snp$Position)
snp$chr_pos <- paste(snp$`Seq ID`,snp$Position,sep="_")
#Mutation data frames found in common reference = false positives
df_1_10ref_snps <- read_csv("~/final_ref_NZ_CP012004_pAb123/snp_final_ref.csv", col_names=TRUE) #this is tibble
df_1_10ref_snps$Position <- gsub(",","",df_1_10ref_snps$Position)
df_1_10ref_snps$chr_pos <- paste(df_1_10ref_snps$`Seq ID`,df_1_10ref_snps$Position,sep="_")
df_snp_noref <- snp[ !(snp$chr_pos %in% df_1_10ref_snps$chr_pos), ]

write_csv(df_snp_noref,path="c-h-a_rerun2_snp_Breseq_Output_noref.csv")
library("writexl")
write_xlsx(df_snp_noref,path="c-h-a_rerun2_snp_Breseq_Output_noref.xlsx")


nje
nje$gene_annot <- paste(nje$Gene,nje$Annotation, sep=":")
#nje$id_pos <- paste(nje$`Seq Id`,nje$Position, sep=":") 
nje$gene_annot <- gsub(" ","", nje$gene_annot,fixed=T) # didnt work

#nje$id_pos <- gsub(" ","", nje$id_pos)


df_1_10ref_nje <- read_csv("~/final_ref_NZ_CP012004_pAb123/nje_final_ref_Breseq_Output_noutf.csv", col_names=TRUE)
#df_1_10ref_nje$gene_annot <- paste(df_1_10ref_nje$Gene, df_1_10ref_nje$Annotation,sep=":") 
#df_1_10ref_nje$id_pos <- paste(df_1_10ref_nje$`Seq Id`, df_1_10ref_nje$Position,sep=":") 

#df_nje_noref <- nje[ !(nje$gene_annot %in% df_1_10ref_nje$gene_annot), ]
df_nje_noref <- nje[ !(nje$Gene %in% df_1_10ref_nje$Gene), ]

write.csv(df_nje_noref,file="c-h-a_rerun2_nje_noref.csv",row.names=F)
write_xlsx(df_nje_noref,path="c-h-a_rerun2_nje_Breseq_Output_noref.xlsx")


df_1_10ref_mc <- read_csv("~/final_ref_NZ_CP012004_pAb123/mc_final_ref_Breseq_Output_noutf.csv", col_names=TRUE)

mc
df_mc_noref <- mc[ !(mc$Gene%in% df_1_10ref_nje$Gene), ]
View(df_mc_noref)
