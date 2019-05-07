d1 <-merge(copy_number,RNA_expression, by="Gene")
mydata <- merge(d1, protein_expression, by="Gene")
clean_data <-mydata[complete.cases(mydata),]
library(dplyr)
#linear regression and generalized regression
#fit <-lm(MYC_Protein~MYC_RNA +MYC_TRUE_CN, data = clean_data)
mod <-glm(MYC_Protein~MYC_RNA +MYC_TRUE_CN, data = clean_data, family = gaussian(link = "identity"))

# calculate and store predicted values
cn_data <-mutate(clean_data,predicted=predict(mod,type = "response"))

#order by copy numer and then RNA then protein
cn_data <-arrange(cn_data,MYC_RNA, MYC_TRUE_CN)
scatter.smooth(x=cn_data$MYC_Protein, y=cn_data$predicted)

#classification and regression tree prediction (train=500,test=379)
need_data <-clean_data[,c(1,3,4,5)]
library(rpart)
library(rpart.plot)
train_indices <-sample(1:nrow(need_data),500, replace = FALSE)
train <- need_data[train_indices,]
test <- need_data[-train_indices,]
mymodel <-rpart(MYC_Protein~ MYC_TRUE_CN+MYC_RNA,data = train,method = "anova")
rpart.plot(mymodel,type = 3,digits = 3,fallen.leaves = TRUE)
predict_data <-predict(mymodel, test)
MAE <-function(actual, predicted){mean(abs(actual-predicted))}
MAE(test$MYC_Protein,predict_data)             
total_error <-sum(abs(test$MYC_Protein-predict_data))          
total_error

#random forecast
library(randomForest)
library(MASS)
mod1 <- randomForest(MYC_Protein~MYC_TRUE_CN+MYC_RNA, need_data,ntree=500)
m <- predict(mod1, need_data)
error2 <-sum(abs(m-need_data$MYC_Protein))
MAE(m,need_data$MYC_Protein)

#gradient boosting
library(xgboost)
library(data.table)
gb <-xgb.DMatrix(data.matrix(need_data[,2:3]), label =as.numeric(need_data[,4]))
fit <-xgb.train(data=gb,
                eta = 0.3,
                max.depth=8,
                nrounds = 50,
                eval_metric = "merror",
                objective ="reg:linear")
n <- predict(fit, as.matrix(need_data[,2:3]))
error3 <-sum(abs(n-need_data$MYC_Protein))
error3
MAE(n,need_data$MYC_Protein)
             