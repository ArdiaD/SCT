#' Script that perform Monte Carlo simulations and store results as CSV files
#' 2024-11-23 - Ardia & Sessinou

rm(list = ls())
library("readxl")
library("pacman")
p_load(parallel,foreach,doParallel,ggplot2,install=FALSE,update=FALSE)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Functions.R")

# Setup in the paper !!! Takes several days !!!
cT = 250; cK = cN = 10; h0 = 0
NB = 500; nb = 100    
seq_cN <- c(2, 10, 50, 100, 400)
seq_cK <- c(2, 10, 50, 100)
seq_garchmodel <- 1:12
seq_denseORsparse <- c(0, 1)
cl <- makeCluster(detectCores() - 1)

# Setup to go fast and test 
cT = 250; cK = cN = 10; h0 = 0
NB = 200; nb = 10  
seq_cN <- c(2, 10)
seq_cK <- c(2, 10)
seq_denseORsparse <- c(0)
seq_garchmodel <- c(1)
cl <- makeCluster(1)

clusterEvalQ(cl, source("Functions.R"))

for(cN in seq_cN) {
      for(denseORsparse in seq_denseORsparse) {
            for(cK in seq_cK) {
                  for (garchmodel in seq_garchmodel){ 
                        
                        resultats<-list()
                        h<-0
                        # ncp=0
                        for (ncp in seq(-0.4,0.4,0.05))
                        {
                              N=cN
                              K=cK
                              ct=cT
                              if(h0==0){
                                    alpha_0=rep(ncp,N)
                                    if(denseORsparse==0)alpha_0[seq(floor(N*0.5))]=0;
                                    Beta_0=rep(1,K)%*%t(rep(1+ncp,N))
                              }
                              if(h0==1){
                                    alpha_0=rep(0,N)
                                    Beta_0=rep(1,K)%*%t(rep(1,N))
                              }
                              
                              Sx=toeplitz((0.8)^seq(0,K-1))
                              Sy=toeplitz((0.5)^seq(0,N-1))
                              
                              Pi=K
                              ids=cut(seq(ct),Pi,labels=F)
                              # NB=1000;nb=100
                              list.rep<-rep(NB/nb,nb)
                              list.rep0<-cumsum(rep(NB/nb,nb))-10
                              ncp<-ncp
                              
                              output <- parallel::clusterApplyLB(
                                    cl       = cl,
                                    x        = seq(list.rep),
                                    fun      = f_many_simulation,
                                    N = N, K = K, ct = ct, alpha_0 = alpha_0, 
                                    Beta_0 = Beta_0, Sx = Sx, Sy = Sy, Pi = Pi, ids = ids, ncp = ncp, 
                                    garchmodel = garchmodel, list.rep = list.rep, list.rep0 = list.rep0
                              )
                              
                              output <- foreach::foreach(a=output,.combine = "c")%do%a
                              
                              outputs<-lapply(seq(output[[1]]), function(i)sapply(seq(NB),function(j){
                                    output[[j]][[i]]
                              }))
                              
                              outputsBis<-c(lapply(seq(length(outputs)-2), function(i){
                                    out<-sapply(c(0.01,0.05,0.1), function(a){
                                          outputs[[i]]<-as.matrix(outputs[[i]])
                                          if(ncol(outputs[[i]])>1) return(rowMeans(outputs[[i]]<a))
                                          else return(mean(outputs[[i]]<a))
                                    })
                                    if(i<4)rownames(out)<-c(sapply(LETTERS[1:5], function(a)paste0(a,c("1/3","1/2","2/3"))))
                                    out
                                    
                              }),
                              
                              list(sapply(c(0.01,0.05,0.1), function(leveltest)
                                    mean(apply(outputs[[length(outputs)-1]],2, function(a){
                                          Reject<-0
                                          accept<-reject<-0
                                          if(a[1]>leveltest) accept<-1
                                          if(a[2]<=leveltest) reject<-1
                                          
                                          if((reject==1)&(accept==1))Reject<-0
                                          else{
                                                if(reject==1) Reject<-1
                                          }
                                          Reject
                                    })))),
                              
                              list(sapply(c(0.01,0.05,0.1), function(leveltest)
                                    mean(apply(outputs[[length(outputs)]],2, function(a){
                                          Reject<-0
                                          accept<-reject<-0
                                          if(a[1]>leveltest) accept<-1
                                          if(a[2]<=leveltest) reject<-1
                                          
                                          if((reject==1)&(accept==1))Reject<-0
                                          else{
                                                if(reject==1) Reject<-1
                                          }
                                          Reject
                                    }))))
                              
                              )
                              
                              resultats[[(h<-h+1)]]<-c(list(ncp=ncp),outputsBis)
                              
                              for( i in seq(outputsBis)){
                                    
                                    print(paste0("K=",K," N=",N, "T=", ct,"D=",i-1,"H0=",h0,"ncp=",ncp,sep=""))
                                    print(outputsBis[[i]])
                              }
                              
                        }
                        
                        resultats[[1]][[6]]<-as.matrix(resultats[[1]][[6]])
                        
                        Ans1<-foreach(j=seq(2,5),.combine = "rbind")%do%sapply(seq(resultats), function(i)resultats[[i]][[j]][seq(1,15,3),2])
                        if(ncol(resultats[[1]][[6]])>1){
                              Ans2<-sapply(seq(resultats), function(i)resultats[[i]][[6]][,2])
                        }else{
                              Ans2<-sapply(seq(resultats), function(i)resultats[[i]][[6]][2])
                        } 
                        
                        rownames(Ans1)<-c(sapply(seq(0,3), function(a)paste0(c("CCTd","CCTda","CCTd/a","CCTa","CCTd/a,a"),"L",a)))
                        
                        Ans3<-sapply(seq(resultats), function(i)resultats[[i]][[7]][2])
                        Ans4<-sapply(seq(resultats), function(i)resultats[[i]][[8]][2])
                        
                        OutputData=t(rbind(Ans1,Ans2,mveGL=Ans3,SpanningGL=Ans4))
                        
                        print("===================================================")
                        if(denseORsparse==1) print(paste0("DenseWyRobStudModel_",garchmodel,"_K_",cK,"_N_",cN,".csv"))
                        if(denseORsparse==0) print(paste0("SparseWyRobStudModel_",garchmodel,"_K_",cK,"_N_",cN,".csv"))
                        print(OutputData)
                        print("===================================================")
                        
                        if(denseORsparse==1)write.csv(OutputData,file = paste0("_csv/DenseWyRobStudModel_",garchmodel,"_K_",cK,"_N_",cN,".csv"))
                        if(denseORsparse==0)write.csv(OutputData,file = paste0("_csv/SparseWyRobStudModel_",garchmodel,"_K_",cK,"_N_",cN,".csv"))
                  }
            }
      }
      
}
parallel::stopCluster(cl)

