
## Download R packages msgl that we need
library(msgl)


## Import data.golubxin represents the expanded input matrix.leibiaoqian represents the class labels.tezhengqunweizhi represents the grouping order.
x<-read.table("golubxin.txt",header=F)  
x<-as.matrix(x)      
y=read.csv("leibiaoqian.csv",header=F)
y<-t(y)  
y<-as.factor(y)       
tezheng<-read.table("tezhengqunweizhi.txt",header=F)
tezheng<-t(tezheng)   
tezheng<-as.matrix(tezheng)



#### we randomly divide the data set into two parts: training set and test set according to the strategic dividing method. The following is one of random dividing methods.
xtrain<-x[c(1,4,5,7:25,31:38,55:72),]
xtrain<-as.matrix(xtrain) 
ytrain<-y[c(1,4,5,7:25,31:38,55:72)]
ytrain<-as.factor(ytrain) 
xtest<-x[-c(1,4,5,7:25,31:38,55:72),]
xtest<-as.matrix(xtest) 
ytest<-y[-c(1,4,5,7:25,31:38,55:72)]
ytest<-as.factor(ytest) 


##Solve MOGL
cl <- makeCluster(2)
registerDoParallel(cl)
fit.cv <- msgl::cv(xtrain, classes=ytrain,grouping=tezheng, fold =5,alpha=0.5, lambda = 0.05, use_parallel = TRUE)     
stopCluster(cl)
fit.cv
Err(fit.cv,type="rate") 
fit<-msgl::fit(xtrain,classes=ytrain,grouping=tezheng,alpha = 0.5,lambda = 0.05)
fit	
features(fit)[[best_model(fit.cv)]]  
parameters(fit)[[best_model(fit.cv)]] 
coef(fit, best_model(fit.cv))  


## Make prediction   
res <- predict(fit,xtest)
yyuce<-res$classes[,best_model(fit.cv)]   
















