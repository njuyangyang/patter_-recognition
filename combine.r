#####################################
##
## pattern progject
## author: YANG YANG
## Data: 04/2015
## 
#####################################

## install packages and set enviroment
install.packages("e1071")
library(e1071)

install.packages("HiDimDA")
library(HiDimDA)

install.packages("class")
libraray(class)

install("plyr")
library(plyr)


setwd("/Users/dajiwu/Desktop/pattern-project")
TrainData<-read.table("Training_Data.txt",header=T,sep='\t')
TestData<-read.table("Testing_Data.txt",header=T,sep='\t')

#######################################
##
##function definition
##
#######################################
##
##support vectors mechine
##
##resubstitution for svm
###################################
resubstitution_error_svm=function(TrainCase,Train_target){
  class_model<-svm(x=TrainCase,y=Train_target,cost=0.1)
  tab=table(predict(class_model),Train_target)
  num1=sum(diag(tab))
  denom1=sum(tab)
  signif(1-num1/denom1,3)
}


#####################################
##cross-validate for svm
#######################
leave_one_error_svm=function(TrainCase,Train_target){
  error=0
  TrainCase<-data.frame(TrainCase)
  for(i in 1:nrow(TrainCase)){
    tempcase=TrainCase[-i,]
    temptarget=Train_target[-i]
    temp_model=svm(x=tempcase,y=temptarget,cost=0.1)
    f=predict(temp_model,TrainCase[i,])
    f=as.numeric(levels(f))[f]
    g=as.numeric(levels(Train_target))[Train_target[i]]
    error=error+(g-f)^2
  }
  error/nrow(TrainCase)
}

#################################
## 
##Diagonal linear discriminant analysis
##
##resubstitution for dlda
###################################
resubstitution_error_dlda=function(TrainCase,Train_target){
  class_model<-Dlda(data=as.matrix(TrainCase),grouping=Train_target,prior=c(0.5,0.5))
  tab=table(predict(class_model,TrainCase)$class,Train_target)
  num1=sum(diag(tab))
  denom1=sum(tab)
  signif(1-num1/denom1,3)
}


#########################
##cross-validate for dlda
#######################
leave_one_error_dlda=function(TrainCase,Train_target){
  error=0
  TrainCase<-data.frame(TrainCase)
  for(i in 1:nrow(TrainCase)){
    tempcase=TrainCase[-i,]
    temptarget=Train_target[-i]
    temp_model=Dlda(data=as.matrix(TrainCase),grouping=Train_target,prior=c(0.5,0.5))
    f=predict(temp_model,TrainCase[i,])$class
    f=as.numeric(levels(f))[f]
    g=as.numeric(levels(Train_target))[Train_target[i]]
    error=error+(g-f)^2
  }
  error/nrow(TrainCase)
}

###################################
##
##k nearest neighbour classification, k=3
##
##resubstitution for knn
###################################
resubstitution_error_knn=function(TrainCase,Train_target){
  TrainCase=as.matrix(TrainCase)
  class_model<-knn(train=TrainCase,test=TrainCase,cl=Train_target,k=3)
  tab=table(class_model,Train_target)
  num1=sum(diag(tab))
  denom1=sum(tab)
  signif(1-num1/denom1,3)
}


#########################
##cross-validate for knn
#######################
leave_one_error_knn=function(TrainCase,Train_target){
  error=0
  TrainCase<-data.frame(TrainCase)
  for(i in 1:nrow(TrainCase)){
    tempcase=TrainCase[-i,]
    temptarget=Train_target[-i]
    temp_model=knn(train=tempcase,test=TrainCase[i,],cl=temptarget,k=3)
    error=error+(as.numeric(levels(Train_target[i]))[Train_target[i]]-as.numeric(levels(temp_model))[temp_model])^2
  }
  error/nrow(TrainCase)
}


##############################
##            
##                
##            balabala      
##
################################


##############################
##
##main procedure
##
##############################

####################
##exhausive search for 1 using resubstitution
####################
svm_resub_1_error<-data.frame(gene_num=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:num_fea)
{
  feature<-i
  TrainCase<-TrainData[,feature]
  TestCase<-TestData[,feature]
  Train_target<-as.factor(TrainData[,ncol(TrainData)])
  Test_target<-as.factor(TestData[,ncol(TrainData)])
  result<-resubstitution_error_svm(TrainCase,Train_target)
  new_item<-list(gene_num=i,result=result)
  svm_resub_1_error<-rbind(svm_resub_1_error,new_item)
}

svm_resub_1_error<-arrange(svm_resub_1_error,result)


#######################
##exhausive search for 2 using resubstitution
#######################

svm_resub_2_error<-data.frame(gene_num_1=numeric(),gene_num_2=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:(num_fea-1)){
  for(j in (i+1):num_fea){
    feature<-c(i,j)
    TrainCase<-TrainData[,feature]
    TestCase<-TestData[,feature]
    Train_target<-as.factor(TrainData[,ncol(TrainData)])
    Test_target<-as.factor(TestData[,ncol(TrainData)])
    result<-resubstitution_error_svm(TrainCase,Train_target)
    new_item<-list(gene_num_1=i,gene_num_2=j,result=result)
    svm_resub_2_error<-rbind(svm_resub_2_error,new_item)
    
  }
}

svm_resub_2_error<-arrange(svm_resub_2_error,result)

#######################
##exhausive search for 3 using resubstitution
#######################

svm_resub_3_error<-data.frame(gene_num_1=numeric(),gene_num_2=numeric(),gene_num_3=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:(num_fea-2)){
  for(j in (i+1):(num_fea-1)){
    for(k in (j+1):num_fea){
      feature<-c(i,j,k)
      TrainCase<-TrainData[,feature]
      TestCase<-TestData[,feature]
      Train_target<-as.factor(TrainData[,ncol(TrainData)])
      Test_target<-as.factor(TestData[,ncol(TrainData)])
      result<-resubstitution_error_svm(TrainCase,Train_target)
      new_item<-list(gene_num_1=i,gene_num_2=j,gene_num_3=k,result=result)
      svm_resub_3_error<-rbind(svm_resub_3_error,new_item)
    }
  }
}

svm_resub_3_error<-arrange(svm_resub_3_error,result)

####################
##exhausive search for 1 using leave one out
####################
svm_leave_1_error<-data.frame(gene_num=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:num_fea)
{
  feature<-i
  TrainCase<-TrainData[,feature]
  TestCase<-TestData[,feature]
  Train_target<-as.factor(TrainData[,ncol(TrainData)])
  Test_target<-as.factor(TestData[,ncol(TrainData)])
  result<-leave_one_error_svm(TrainCase,Train_target)
  new_item<-list(gene_num=i,result=result)
  svm_leave_1_error<-rbind(svm_leave_1_error,new_item)
}

svm_leave_1_error<-arrange(svm_leave_1_error,result)


#######################
##exhausive search for 2 using leave one out
#######################

svm_leave_2_error<-data.frame(gene_num_1=numeric(),gene_num_2=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:(num_fea-1)){
  for(j in (i+1):num_fea){
    feature<-c(i,j)
    TrainCase<-TrainData[,feature]
    TestCase<-TestData[,feature]
    Train_target<-as.factor(TrainData[,ncol(TrainData)])
    Test_target<-as.factor(TestData[,ncol(TrainData)])
    result<-leave_one_error_svm(TrainCase,Train_target)
    new_item<-list(gene_num_1=i,gene_num_2=j,result=result)
    svm_leave_2_error<-rbind(svm_leave_2_error,new_item)
    
  }
}

svm_leave_2_error<-arrange(svm_leave_2_error,result)

#######################
##exhausive search for 3 using leave one out
#######################

svm_leave_3_error<-data.frame(gene_num_1=numeric(),gene_num_2=numeric(),gene_num_3=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:(num_fea-2)){
  for(j in (i+1):(num_fea-1)){
    for(k in (j+1):num_fea){
      feature<-c(i,j,k)
      TrainCase<-TrainData[,feature]
      TestCase<-TestData[,feature]
      Train_target<-as.factor(TrainData[,ncol(TrainData)])
      Test_target<-as.factor(TestData[,ncol(TrainData)])
      result<-leave_one_error_svm(TrainCase,Train_target)
      new_item<-list(gene_num_1=i,gene_num_2=j,gene_num_3=k,result=result)
      svm_leave_3_error<-rbind(svm_leave_3_error,new_item)
    }
  }
}

svm_leave_3_error<-arrange(svm_leave_3_error,result)


##########
##
####################
##exhausive search for 1 using resubstitution
####################
dlda_resub_1_error<-data.frame(gene_num=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:num_fea)
{
  feature<-i
  TrainCase<-TrainData[,feature]
  TestCase<-TestData[,feature]
  Train_target<-as.factor(TrainData[,ncol(TrainData)])
  Test_target<-as.factor(TestData[,ncol(TrainData)])
  result<-resubstitution_error_dlda(TrainCase,Train_target)
  new_item<-list(gene_num=i,result=result)
  dlda_resub_1_error<-rbind(dlda_resub_1_error,new_item)
}

dlda_resub_1_error<-arrange(dlda_resub_1_error,result)


#######################
##exhausive search for 2 using resubstitution
#######################

dlda_resub_2_error<-data.frame(gene_num_1=numeric(),gene_num_2=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:(num_fea-1)){
  for(j in (i+1):num_fea){
    feature<-c(i,j)
    TrainCase<-TrainData[,feature]
    TestCase<-TestData[,feature]
    Train_target<-as.factor(TrainData[,ncol(TrainData)])
    Test_target<-as.factor(TestData[,ncol(TrainData)])
    result<-resubstitution_error_dlda(TrainCase,Train_target)
    new_item<-list(gene_num_1=i,gene_num_2=j,result=result)
    dlda_resub_2_error<-rbind(dlda_resub_2_error,new_item)
    
  }
}

dlda_resub_2_error<-arrange(dlda_resub_2_error,result)

#######################
##exhausive search for 3 using resubstitution
#######################

dlda_resub_3_error<-data.frame(gene_num_1=numeric(),gene_num_2=numeric(),gene_num_3=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:(num_fea-2)){
  for(j in (i+1):(num_fea-1)){
    for(k in (j+1):num_fea){
      feature<-c(i,j,k)
      TrainCase<-TrainData[,feature]
      TestCase<-TestData[,feature]
      Train_target<-as.factor(TrainData[,ncol(TrainData)])
      Test_target<-as.factor(TestData[,ncol(TrainData)])
      result<-resubstitution_error_dlda(TrainCase,Train_target)
      new_item<-list(gene_num_1=i,gene_num_2=j,gene_num_3=k,result=result)
      dlda_resub_3_error<-rbind(dlda_resub_3_error,new_item)
    }
  }
}

dlda_resub_3_error<-arrange(dlda_resub_3_error,result)

####################
##exhausive search for 1 using leave one out
####################
dlda_leave_1_error<-data.frame(gene_num=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:num_fea)
{
  feature<-i
  TrainCase<-TrainData[,feature]
  TestCase<-TestData[,feature]
  Train_target<-as.factor(TrainData[,ncol(TrainData)])
  Test_target<-as.factor(TestData[,ncol(TrainData)])
  result<-leave_one_error_dlda(TrainCase,Train_target)
  new_item<-list(gene_num=i,result=result)
  dlda_leave_1_error<-rbind(dlda_leave_1_error,new_item)
}

dlda_leave_1_error<-arrange(dlda_leave_1_error,result)


#######################
##exhausive search for 2 using leave one out
#######################

dlda_leave_2_error<-data.frame(gene_num_1=numeric(),gene_num_2=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:(num_fea-1)){
  for(j in (i+1):num_fea){
    feature<-c(i,j)
    TrainCase<-TrainData[,feature]
    TestCase<-TestData[,feature]
    Train_target<-as.factor(TrainData[,ncol(TrainData)])
    Test_target<-as.factor(TestData[,ncol(TrainData)])
    result<-leave_one_error_dlda(TrainCase,Train_target)
    new_item<-list(gene_num_1=i,gene_num_2=j,result=result)
    dlda_leave_2_error<-rbind(dlda_leave_2_error,new_item)
    
  }
}

dlda_leave_2_error<-arrange(dlda_leave_2_error,result)

#######################
##exhausive search for 3 using leave one out
#######################

dlda_leave_3_error<-data.frame(gene_num_1=numeric(),gene_num_2=numeric(),gene_num_3=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:(num_fea-2)){
  for(j in (i+1):(num_fea-1)){
    for(k in (j+1):num_fea){
      feature<-c(i,j,k)
      TrainCase<-TrainData[,feature]
      TestCase<-TestData[,feature]
      Train_target<-as.factor(TrainData[,ncol(TrainData)])
      Test_target<-as.factor(TestData[,ncol(TrainData)])
      result<-leave_one_error_dlda(TrainCase,Train_target)
      new_item<-list(gene_num_1=i,gene_num_2=j,gene_num_3=k,result=result)
      dlda_leave_3_error<-rbind(dlda_leave_3_error,new_item)
    }
  }
}

dlda_leave_3_error<-arrange(dlda_leave_3_error,result)


####################
##exhausive search for 1 using resubstitution
####################
knn_resub_1_error<-data.frame(gene_num=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:num_fea)
{
  feature<-i
  TrainCase<-TrainData[,feature]
  TestCase<-TestData[,feature]
  Train_target<-as.factor(TrainData[,ncol(TrainData)])
  Test_target<-as.factor(TestData[,ncol(TrainData)])
  result<-resubstitution_error_knn(TrainCase,Train_target)
  new_item<-list(gene_num=i,result=result)
  knn_resub_1_error<-rbind(knn_resub_1_error,new_item)
}

knn_resub_1_error<-arrange(knn_resub_1_error,result)


#######################
##exhausive search for 2 using resubstitution
#######################

knn_resub_2_error<-data.frame(gene_num_1=numeric(),gene_num_2=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:(num_fea-1)){
  for(j in (i+1):num_fea){
    feature<-c(i,j)
    TrainCase<-TrainData[,feature]
    TestCase<-TestData[,feature]
    Train_target<-as.factor(TrainData[,ncol(TrainData)])
    Test_target<-as.factor(TestData[,ncol(TrainData)])
    result<-resubstitution_error_knn(TrainCase,Train_target)
    new_item<-list(gene_num_1=i,gene_num_2=j,result=result)
    knn_resub_2_error<-rbind(knn_resub_2_error,new_item)
    
  }
}

knn_resub_2_error<-arrange(knn_resub_2_error,result)

#######################
##exhausive search for 3 using resubstitution
#######################

knn_resub_3_error<-data.frame(gene_num_1=numeric(),gene_num_2=numeric(),gene_num_3=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:(num_fea-2)){
  for(j in (i+1):(num_fea-1)){
    for(k in (j+1):num_fea){
      feature<-c(i,j,k)
      TrainCase<-TrainData[,feature]
      TestCase<-TestData[,feature]
      Train_target<-as.factor(TrainData[,ncol(TrainData)])
      Test_target<-as.factor(TestData[,ncol(TrainData)])
      result<-resubstitution_error_knn(TrainCase,Train_target)
      new_item<-list(gene_num_1=i,gene_num_2=j,gene_num_3=k,result=result)
      knn_resub_3_error<-rbind(knn_resub_3_error,new_item)
    }
  }
}

knn_resub_3_error<-arrange(knn_resub_3_error,result)

####################
##exhausive search for 1 using leave one out
####################
knn_leave_1_error<-data.frame(gene_num=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:num_fea)
{
  feature<-i
  TrainCase<-TrainData[,feature]
  TestCase<-TestData[,feature]
  Train_target<-as.factor(TrainData[,ncol(TrainData)])
  Test_target<-as.factor(TestData[,ncol(TrainData)])
  result<-leave_one_error_knn(TrainCase,Train_target)
  new_item<-list(gene_num=i,result=result)
  knn_leave_1_error<-rbind(knn_leave_1_error,new_item)
}

knn_leave_1_error<-arrange(knn_leave_1_error,result)


#######################
##exhausive search for 2 using leave one out
#######################

knn_leave_2_error<-data.frame(gene_num_1=numeric(),gene_num_2=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:(num_fea-1)){
  for(j in (i+1):num_fea){
    feature<-c(i,j)
    TrainCase<-TrainData[,feature]
    TestCase<-TestData[,feature]
    Train_target<-as.factor(TrainData[,ncol(TrainData)])
    Test_target<-as.factor(TestData[,ncol(TrainData)])
    result<-leave_one_error_knn(TrainCase,Train_target)
    new_item<-list(gene_num_1=i,gene_num_2=j,result=result)
    knn_leave_2_error<-rbind(knn_leave_2_error,new_item)
    
  }
}

knn_leave_2_error<-arrange(knn_leave_2_error,result)

#######################
##exhausive search for 3 using leave one out
#######################

knn_leave_3_error<-data.frame(gene_num_1=numeric(),gene_num_2=numeric(),gene_num_3=numeric(),result=numeric())
num_fea<-ncol(TrainData)-1
for(i in 2:(num_fea-2)){
  for(j in (i+1):(num_fea-1)){
    for(k in (j+1):num_fea){
      feature<-c(i,j,k)
      TrainCase<-TrainData[,feature]
      TestCase<-TestData[,feature]
      Train_target<-as.factor(TrainData[,ncol(TrainData)])
      Test_target<-as.factor(TestData[,ncol(TrainData)])
      result<-leave_one_error_knn(TrainCase,Train_target)
      new_item<-list(gene_num_1=i,gene_num_2=j,gene_num_3=k,result=result)
      knn_leave_3_error<-rbind(knn_leave_3_error,new_item)
    }
  }
}

knn_leave_3_error<-arrange(knn_leave_3_error,result)

##############
##output all work
##############
svm_resub_1_error_head<-head(svm_resub_1_error,20)
svm_resub_2_error_head<-head(svm_resub_2_error,20)
svm_resub_3_error_head<-head(svm_resub_3_error,20)
svm_leave_1_error_head<-head(svm_leave_1_error,20)
svm_leave_2_error_head<-head(svm_leave_2_error,20)
svm_leave_3_error_head<-head(svm_leave_3_error,20)
dlda_resub_1_error_head<-head(dlda_resub_1_error,20)
dlda_resub_2_error_head<-head(dlda_resub_2_error,20)
dlda_resub_3_error_head<-head(dlda_resub_3_error,20)
dlda_leave_1_error_head<-head(dlda_leave_1_error,20)
dlda_leave_2_error_head<-head(dlda_leave_2_error,20)
dlda_leave_3_error_head<-head(dlda_leave_3_error,20)
knn_resub_1_error_head<-head(knn_resub_1_error,20)
knn_resub_2_error_head<-head(knn_resub_2_error,20)
knn_resub_3_error_head<-head(knn_resub_3_error,20)
knn_leave_1_error_head<-head(knn_leave_1_error,20)
knn_leave_2_error_head<-head(knn_leave_2_error,20)
knn_leave_3_error_head<-head(knn_leave_3_error,20)



write.csv(svm_resub_1_error_head,"svm_resub_1_error_head.csv")
write.csv(svm_resub_2_error_head,"svm_resub_2_error_head.csv")
write.csv(svm_resub_3_error_head,"svm_resub_3_error_head.csv")
write.csv(svm_leave_1_error_head,"svm_leave_1_error_head.csv")
write.csv(svm_leave_2_error_head,"svm_leave_2_error_head.csv")
write.csv(svm_leave_3_error_head,"svm_leave_3_error_head.csv")

write.csv(dlda_resub_1_error_head,"dlda_resub_1_error_head.csv")
write.csv(dlda_resub_2_error_head,"dlda_resub_2_error_head.csv")
write.csv(dlda_resub_3_error_head,"dlda_resub_3_error_head.csv")
write.csv(dlda_leave_1_error_head,"dlda_leave_1_error_head.csv")
write.csv(dlda_leave_2_error_head,"dlda_leave_2_error_head.csv")
write.csv(dlda_leave_3_error_head,"dlda_leave_3_error_head.csv")

write.csv(svm_resub_1_error_head,"svm_resub_1_error_head.csv")
write.csv(svm_resub_2_error_head,"svm_resub_2_error_head.csv")
write.csv(svm_resub_3_error_head,"svm_resub_3_error_head.csv")
write.csv(svm_leave_1_error_head,"svm_leave_1_error_head.csv")
write.csv(svm_leave_2_error_head,"svm_leave_2_error_head.csv")
write.csv(svm_leave_3_error_head,"svm_leave_3_error_head.csv")








