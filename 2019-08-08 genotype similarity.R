library(reshape2)
library(vegan)

setwd("C:/Users/cbt12/My Documents/History Chance and Adaptation Antibiotics")

mut_list <- read.csv("2019-08-08 mutation frequencies partially melted.csv", header=TRUE )

#set NAs to zero
mut_list[is.na(mut_list)] <- 0



chance <- function(within){
  #returns mean distance between replicate populations
  sum <- 0
  count <- 0
  evolved <- grepl("\\.",row.names(within))
  within <- within[evolved, evolved]
  for(i in 1:(nrow(within)-1)){
    for(j in (i+1):ncol(within)){
      if(rownames(within)[i] != colnames(within)[j] & substr(rownames(within)[i],4,4)==substr(colnames(within)[j],4,4) ){
        sum <- sum+within[i,j]
        count <- count+1
      }
    }
  }
  
  return(sum/count*0.5)
}

history <- function(within){
  #Return mean distance between populations with different ancestors minus Chance
  sum <- 0
  count <- 0
  evolved <- grepl("\\.",row.names(within))
  within <- within[evolved, evolved]
  for(i in 1:(nrow(within)-1)){
    for(j in (i+1):ncol(within)){
      if(substr(rownames(within)[i],4,4)!=substr(colnames(within)[j],4,4) ){
        sum <- sum+within[i,j]
        count <- count+1
      }
    }
  }

  return(sum/count*0.5-chance(within))
}

adaptation <- function(within, comm){
  #returns mean distance between ancestral and evolved populations
  sum <- 0
  count <- 0

  for(i in 1:(nrow(within)-1)){
    for(j in (i+1):ncol(within)){
      if(!grepl("\\.", rownames(within)[i]) & grepl("\\.", colnames(within)[j]) & substr(rownames(within)[i],4,4)==substr(colnames(within)[j],4,4) ){
        sum <- sum+within[i,j]
        count <- count+1
      }
    }
  }

  #make a new matrix with ancestral mutations removed.
  #this allows the history component of only new mutations to be subtracted from the adaptation value
  anc_removed <- comm
  ancestors <- c(1,5,9,13,17,21)
  for(i in 1:ncol(comm)){
    for(j in ancestors){
      if(comm[j,i] <= comm[j+1,i]){
        anc_removed[j+1,i] <- comm[j+1,i]-comm[j,i]
      }
      if(comm[j,i] <= comm[j+2,i]){
        anc_removed[j+2,i] <- comm[j+2,i]-comm[j,i]  
      }
      if(comm[j,i] <= comm[j+3,i]){
        anc_removed[j+3,i] <- comm[j+3,i]-comm[j,i]
      }  
    }
  }
  anc_removed_manhattan <- as.matrix(vegdist(anc_removed, method="manhattan"))
  

  return((sum/count)-chance(within)-history(anc_removed_manhattan))

}


hca <- function(mut_list){
  #reshape the mutation list into a community data matrix
  
  community <- dcast(mut_list, variable~gene_product,value.var="value",sum)
  comm <- community[,-1]
  rownames(comm) <- community[,1]
  
  manhattan <- as.matrix(vegdist(comm, method="manhattan"))
  
  #run the three analyses, history, chance and adaptation
  print("chance:")
  print(chance(manhattan))
  print("history:")
  print(history(manhattan))
  print("adaptation:")
  print(adaptation(manhattan, comm))
  
}




#set up ceftazidime mutations
ceft_bio <- mut_list[mut_list$environment=="ceftazidime" & mut_list$Environment_history=="biofilm",]
melt_ceft_bio <- melt(ceft_bio)
ceft_plank <- mut_list[mut_list$environment=="ceftazidime"& mut_list$Environment_history=="planktonic",]
plank_names <- colnames(mut_list)
plank_names[6:17] <- c("anc4","anc4.1","anc4.2","anc4.3", "anc5", "anc5.1","anc5.2","anc5.3","anc6","anc6.1","anc6.2","anc6.3")
colnames(ceft_plank) <- plank_names
melt_ceft_plank <- melt(ceft_plank)

melt_ceft_all <- rbind(melt_ceft_bio, melt_ceft_plank)
hca(melt_ceft_all)


#set up imipenem mutations
imip_bio <- mut_list[mut_list$environment=="imipenem"& mut_list$Environment_history=="biofilm",]
melt_imip_bio <- melt(imip_bio)
imip_plank <- mut_list[mut_list$environment=="imipenem"& mut_list$Environment_history=="planktonic",]
colnames(imip_plank) <- plank_names
melt_imip_plank <- melt(imip_plank)

melt_imip_all <- rbind(melt_imip_bio, melt_imip_plank)
hca(melt_imip_all)


