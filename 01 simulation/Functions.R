#dirname(rstudioapi::getActiveDocumentContext()$path)

cauchypv<-function(p){
      .5-atan(mean(tan((.5-p)*pi)))/pi
}


cct<-function(x,y)
{
      require(testcorr)
      require(stats)
      cauchypv(c(sapply(seq(ncol(x)), function(i)(unlist(testcorr::cc.test(x[,i],y[,i],max.lag = 50,plot = F,table = F)[["pvqtilde"]])))
      ))
}
ttest<-function(eps,func=colMeans,k=1/1.5)
{
      eps=as.matrix(eps)
      k=floor(nrow(eps)^k)#ceiling(nrow(eps)/30)#
      folds=cut(seq(nrow(eps)),k,labels = F)
      estvar=as.matrix(sapply(seq(k),function(i)func(as.matrix(eps[folds==i,]))))
      if(ncol(estvar)>1)estvar=t(estvar)
      return(list(student=apply(estvar, 2, function(u)t.test(u)$p.value),normal=2*pnorm(-abs(mean(eps)/apply(estvar, 2, function(u)sd(u))))))
}

waldtest<-function(h,b)
{
      cauchypv(apply(h, 2,function(a) {
            require(sandwich)
            require(lmtest)
            reg<-lm(a~b-1)
            NW_VCOV <- NeweyWest(reg, prewhite = T, adjust = T)
            coeftest(reg, vcov = NW_VCOV)[,4]
      }))
}
testbm<-function(eps,...)
{
      n=ncol(eps)
      out<-ttest(eps,...)
      sapply(out,function(u)cauchypv(u))
}

garchm<-function(ct=250,cseed=123,model=garchmodel)
{
      require(fGarch)
      require(stats)
      
      if(model==1) return(rnorm(ct))#normal
      if(model==2) return(rt(ct,5))#st
      if(model==3) return(rsstd(ct,nu=4,xi = 0.9))#sst
      
      if(model==4){#garch normal
            spec = garchSpec(model = list(omega = 0.1, alpha = 0.1, beta = 0.8))
            return(garchSim(spec, n = ct))#runif(1,0.5,1)*
      }
      if(model==5){#garch st
            spec = garchSpec(model = list(omega = 0.1, alpha = 0.1, beta = 0.8),cond.dist="std")
            return(garchSim(spec, n = ct))#runif(1,0.5,1)*
      }
      if(model==6){#garch sst
            spec = garchSpec(model = list(omega = 0.1, alpha = 0.1, beta = 0.8),cond.dist="sstd")
            return(garchSim(spec, n = ct))#runif(1,0.5,1)*
      }
      
      if(model==7){#AR-garch normal
            spec = garchSpec(model = list(omega = 0.1, alpha = 0.1, beta = 0.8,ar = 0.2))
            return(garchSim(spec, n = ct))#runif(1,0.5,1)*
      }
      if(model==8){#AR-garch st
            spec = garchSpec(model = list(omega = 0.1, alpha = 0.1, beta = 0.8,ar = 0.2),cond.dist="std")
            return(garchSim(spec, n = ct))#runif(1,0.5,1)*
      }
      if(model==9){#AR-garch sst
            spec = garchSpec(model = list(omega = 0.1, alpha = 0.1, beta = 0.8,ar = 0.2),cond.dist="sstd")
            return(garchSim(spec, n = ct))#runif(1,0.5,1)*
      }  
      
      if(model == 10) return(arima.sim(n = ct,rand.gen = function(n) rnorm(ct), list(ar = 0.2 )))
      if(model == 11) return(arima.sim(n = ct,rand.gen = function(n) rt(n, df = 5), list(ar = 0.2 )))
      if(model == 12) return(arima.sim(n = ct,rand.gen = function(n) rsstd(ct,nu = 4, xi = 0.9), list(ar = 0.2 )))
      
}

ranklex <- function(x,uu) {
      n  <- length(x)
      vv <- cbind(x,uu) 
      r  <- matrix(0, nr=n, nc=1)
      for (i in 1:(n-1)){
            if (vv[n,1] > vv[i,1]){
                  r[i] <- 1     
            }
            if ((vv[n,1] == vv[i,1]) & (vv[n,2] > vv[i,2])) {
                  r[i] <- 1   
            }
      }
      sum(r)+1
}

HK<-function(R1,R2){
      
      #generate R
      R<-cbind(R1,R2)
      p=ncol(R1)
      p2=ncol(R2)
      N=nrow(R1)
      
      #CCT
      reg<-lm(R2-R1[,1]~cbind(-R1[,1],R1[,-1]-R1[,1]))
      outpv<-apply(sapply(coef(summary(reg)), function(a)a[1:2,4]),1,cauchypv)
      names(outpv)<-c("CCTa","CCTd")
      # View(outpv)
      
      #Commpute F-stats
      mu<-matrix(colMeans(R))
      ones<-matrix(1,p+p2)
      V<-cov(R)
      iV<-solve(V)
      a<-crossprod(mu,iV)%*%mu
      b<-crossprod(mu,iV)%*%ones
      c<-crossprod(ones,iV)%*%ones
      d<-a*c-b*b
      
      V11<-cov(R1)
      iV11<-solve(V11)
      mu1<-matrix(colMeans(R1))
      ones1=matrix(1,p)
      a1<-crossprod(mu1,iV11)%*%mu1
      b1<-crossprod(mu1,iV11)%*%ones1
      c1<-crossprod(ones1,iV11)%*%ones1
      d1<-a1*c1-b1*b1
      
      U<-(c1+d1)/(c+d)
      f<-((N-p-p2)/p2)*(1/sqrt(U)-1)
      f1<-((N-p-p2)/p2)*((a-a1)/(1+a1))
      f2<-((N-p-p2+1)/p2)*( ((c+d)/(c1+d1))*((1+a1)/(1+a))-1 )
      
      pval<-function(f,f1,f2,N,p,p2)
      {
            pv<-pf(f,df1=2*p2,df2=2*(N-p-p2),lower.tail = F)
            pv1<-pf(f1,df1=p2,df2=(N-p-p2),lower.tail = F)
            pv2<-pf(f2,df1=p2,df2=(N-p-p2+1),lower.tail = F)
            return(c(HK=pv,F1=pv1,F2=pv2))
      }
      
      alternative<-function(R,p2,p)
      {
            y<-R[,1]
            x<-R[,1]-R[,-1]
            #GMVP weigths
            GMVPweigth<-function(u)
            {
                  a<-solve(cov(u))
                  matrix(colSums(a)/sum(a))
            }
            GMVPweigth(R)
            reg<-lm(y~x)
            require(car)
            gmvpspanning<-linearHypothesis(reg,hypothesis.matrix =paste0("x",seq(p,p+p2-1) ))  
            x2<-R[,-1]-R[,1]
            y2<-rep(1,N)
            reg2<-lm(y2~x2-1)
            # coef(reg2)/sum(coef(reg2))
            # crossprod(solve(cov(x2)),colMeans(x2))/sum(crossprod(solve(cov(x2)),colMeans(x2)))
            ctangspanning<-linearHypothesis(reg2,hypothesis.matrix =paste0("x2",seq(p,p+p2-1) ))  
            a<-c(gmvpspanning$`Pr(>F)`[-1],ctangspanning$`Pr(>F)`[-1])
            return(a)
      }
      
      
      res1<-c(pval(f,f1,f2,N,p,p2),outpv,CCTad=cauchypv(outpv))
      res4<-alternative(R,p2,p)
      return(c(res1,MK_BJ=res4))
      
}
GungorLungerPesaranY<-function(Y,X,seed=1){
      
      TT<-nrow(X)
      K<-ncol(X)
      N<-ncol(Y)
      totsim <- 500
      ones <- matrix(1,TT,1)
      
      XX   <- matrix(cbind(ones, X),TT,(K+1))
      
      
      
      # OLS parameter estimates
      
      
      Xtemp    <- solve( t(XX) %*% XX)
      
      Bhat1 <- Xtemp %*% t(XX) %*% Y
      
      Ehat1 <- Y - XX %*% Bhat1         # Unrestricted
      SSRu <- t(Ehat1) %*% Ehat1
      SigmaU <- SSRu / TT
      
      
      H         <- matrix(0,1,K+1)
      H[1,1]    <- 1
      C <- matrix(0,1,N)
      
      Bhat0 <- Bhat1 - Xtemp %*% t(H) %*% solve( H %*% Xtemp %*% t(H) ) %*% (H %*% Bhat1 - C )		# Restricted
      Ehat0 <- Y - XX %*% Bhat0
      SSRr <- t(Ehat0) %*% Ehat0
      SigmaR <- SSRr / TT
      
      
      # GRS test
      
      #
      if ((TT - K - N) >= 1 ) {
            ahat      <- Bhat1[1,]
            shat      <- SigmaU
            
            fact     <- t(X)
            mufactor     <- matrix(0,K,1)
            ssqm_temp    <- matrix(0, TT, K)
            
            for (i in 1:K){
                  mufactor[i]     <- (1/TT)*sum(X[,i])
                  ssqm_temp[,i]   <- X[,i]-mufactor[i]
            }
            
            ssqm <-  1/TT*( t(ssqm_temp)%*%ssqm_temp )
            
            GRS       <- ((TT-N-K)/N)*solve(1+( t(mufactor) %*% solve(ssqm) %*% mufactor ))*(t(ahat)%*%solve(shat)%*%ahat)
            pval_GRS  <-  1-pf(GRS,N,(TT-N-K))
            pvalHK<-HK(X,Y)
      }
      
      
      # Pesaran-Yamagata tests
      
      
      v <- TT - K - 1
      
      t2 <- matrix(0, N, 1)
      
      eye <- matrix(0, TT, TT)
      diag(eye) <- 1
      
      MX <- eye - X %*% solve( t(X) %*% X ) %*% t(X)
      
      num <- matrix( (t(ones) %*% MX %*% ones) * v, N, 1)
      
      t2 <- (Bhat1[1,]^2) * num / ( TT * diag(SigmaU) )
      
      pN <- 0.05/(N-1)
      thetaN <- qnorm( 1-pN/2 )^2
      rhobar <- 0
      for ( i in 2 : N ){
            for ( j in 1 : (i-1) ){
                  temp <- SigmaU[i,j] / sqrt( SigmaU[i,i] * SigmaU[j,j] )
                  temp2 <- temp^2
                  if ( v*temp2 >= thetaN){
                        rhobar <- rhobar + temp2
                  }
            }
      }
      rhobar <- rhobar * 2 / ( N * ( N - 1 ) )
      
      Jalpha2 <- sum(t2 - v/(v-2))/sqrt(N)
      
      den <- ( v/ (v-2) ) * sqrt( 2 * (v-1) * (1 + (N-1)*rhobar)  / (v-4) )
      Jalpha2 <- Jalpha2 / den
      
      pval_Jalpha2 <- 1-pnorm(Jalpha2)
      
      ## MC Fmax tests
      temp <-  (diag(SSRr) - diag(SSRu) )/( diag(SSRu)  )            
      Fmax  <- max( temp )  
      Fmax_actual  <- Fmax 	# Gives the Fmax from the actual (not simulated) sample
      
      LMCstats <- matrix(0,totsim,1)
      LMCstats[totsim] <- Fmax
      
      BMCstats <- matrix(0,totsim,1)
      BMCstats[totsim,1] <- Fmax 
      
      Ehat0data <- Ehat0   
      Bhat0data <- Bhat0     
      
      set.seed(123456*seed)
      
      for (isim in 1:(totsim-1)){
            
            esim <- sign(rnorm(TT)) * Ehat0data
            
            Ysim <-  XX %*% Bhat0data +  esim
            
            Y <- Ysim
            
            Bhat1 <- Xtemp %*% t(XX) %*% Y    
            Ehat1 <- Y - XX %*% Bhat1
            SSRu <- t(Ehat1) %*% Ehat1
            
            Bhat0 <- Bhat1 - Xtemp %*% t(H) %*% solve( H %*% Xtemp %*% t(H) ) %*% (H %*% Bhat1 - C )
            Ehat0 <- Y - XX %*% Bhat0
            SSRr <- t(Ehat0) %*% Ehat0
            
            # LMC test		    		  		    
            
            temp <- (diag(SSRr) - diag(SSRu) )/( diag(SSRu)  )                        
            Fmax  <- max( temp )		    		       		    		    	    	    	    
            LMCstats[isim] <- Fmax	            
            
            # BMC test
            
            SSRr <- t(esim) %*% esim
            temp <- (diag(SSRr) - diag(SSRu) )/( diag(SSRu)  )                                                          
            Fmax  <- max( temp )				    		    		                		 		                	    		             		                	    
            BMCstats[isim] <- Fmax
            
      }    	
      
      uu <- runif(totsim)
      
      temp <- ranklex(LMCstats,uu) 
      LMCpvalueFmax <- (totsim - temp + 1)/totsim
      
      temp <- ranklex(BMCstats,uu) 
      BMCpvalueFmax <- (totsim - temp + 1)/totsim    
      
      
      # print('============================================')
      
      # print(c('Period =', datee[B],'to', datee[E], 'Nsize =', N, 'Ksize =', K))
      
      # print('============================================')
      
      # if ((TT - K - N) >= 1 ) {
      # 	print(c('GRS        =', GRS))
      # 	print(c('p-value    =', pval_GRS))
      # }
      
      # print('============================================')
      # 
      # print(c('Jalpha2    =',Jalpha2))
      # print(c('p-value    =', pval_Jalpha2))
      # 
      # print('============================================')
      # 
      # print('F-Max')
      # print(c('------','F-max           =', Fmax_actual))
      # 
      # if (LMCpvalueFmax > 0.05 ) {
      # 	print(c('ACCEPT ----------','LMCp-value       =', LMCpvalueFmax))
      # }
      # 
      # if (BMCpvalueFmax <= 0.05 ) {
      # 	print(c('REJECT ----------','BMC p-value      =', BMCpvalueFmax))
      # }
      # 
      # if ( (LMCpvalueFmax <= 0.05 ) & (BMCpvalueFmax > 0.05 ) ) {
      # 	print(c('INCONCLUSIVE ----','LMC p-value      =', LMCpvalueFmax))
      # 	print(c('INCONCLUSIVE ----','BMC p-value      =', BMCpvalueFmax))
      # }
      if ((TT - K - N) >= 1 ) {
            return(list(c(GRS=pval_GRS,PY=pval_Jalpha2,pvalHK), c(LMCpvalueFmax=LMCpvalueFmax,BMCpvalueFmax=BMCpvalueFmax) )  )
      } else{
            return(list(PY=pval_Jalpha2, c(LMCpvalueFmax=LMCpvalueFmax,BMCpvalueFmax=BMCpvalueFmax)))
      }
      
      
}

SpanningHKGL<-function(Y,X,printRes=0){
      TT<-nrow(X)
      K<-ncol(X)
      N<-ncol(Y)
      totsim <- 500
      ones <- matrix(1,TT,1)
      
      XX   <- matrix(cbind(ones, X),TT,(K+1))
      
      
      
      # OLS parameter estimates
      
      Xtemp    <- solve( t(XX) %*% XX)
      
      Bhat1 <- Xtemp %*% t(XX) %*% Y
      
      
      Ehat1 <- Y - XX %*% Bhat1         # Unrestricted
      SSRu <- t(Ehat1) %*% Ehat1
      SigmaU <- SSRu / TT
      
      
      H         <- matrix(0,2,K+1)
      H[1,1]    <- 1
      H[2,2:(K+1)] <- 1
      C <- matrix(0,2,N)
      C[2,]  <- 1
      
      Bhat0 <- Bhat1 - Xtemp %*% t(H) %*% solve( H %*% Xtemp %*% t(H) ) %*% (H %*% Bhat1 - C )   # Restricted
      Ehat0 <- Y - XX %*% Bhat0
      SSRr <- t(Ehat0) %*% Ehat0
      SigmaR <- SSRr / TT
      
      
      ## HK test
      
      # if (TRUE){
      #   
      #   A  <- 100   # Rescaling scalar: This takes care of the overflow during the numerical computation of the determinant
      #   
      #   if (2*(TT - K - N) >= 1 ) {
      #     UU <- det(SigmaR*A)/det(SigmaU*A)        
      #     HK <- (TT - K - N)*(1/N)* (sqrt(UU) - 1)
      #     pval_HK <- 1- pf(HK, 2*N, 2*(TT - K - N) )
      #   }            
      # }  
      
      
      ## New MC tests
      
      temp <-  (diag(SSRr) - diag(SSRu) )/( diag(SSRu)  )      
      Fmax  <- max( temp ) 
      
      Fmax_actual  <- Fmax # Gives the Fmax stat from the actual (not simulated) sample
      
      LMCstats <- matrix(0,totsim,1)
      LMCstats[totsim] <- Fmax
      
      BMCstats <- matrix(0,totsim,1)
      BMCstats[totsim] <- Fmax 
      
      Ehat0data <- Ehat0   
      Bhat0data <- Bhat0    
      
      set.seed(123456) 
      
      for (isim in 1:(totsim-1)){
            
            esim <- sign(rnorm(TT)) * Ehat0data
            
            Ysim <-  XX %*% Bhat0data +  esim
            
            Y <- Ysim
            
            Bhat1 <- Xtemp %*% t(XX) %*% Y    
            Ehat1 <- Y - XX %*% Bhat1
            SSRu <- t(Ehat1) %*% Ehat1
            
            Bhat0 <- Bhat1 - Xtemp %*% t(H) %*% solve( H %*% Xtemp %*% t(H) ) %*% (H %*% Bhat1 - C )
            Ehat0 <- Y - XX %*% Bhat0
            SSRr <- t(Ehat0) %*% Ehat0
            
            # LMC test	
            
            temp <- (diag(SSRr) - diag(SSRu) )/( diag(SSRu)  )                         
            Fmax  <- max( temp )                   		    	    	    	    
            LMCstats[isim] <- Fmax	            
            
            # BMC test
            
            SSRr <- t(esim) %*% esim		    		    		                		  
            temp <- (diag(SSRr) - diag(SSRu) )/( diag(SSRu)  )                        
            Fmax  <- max( temp )        	                	    		             		                	   
            BMCstats[isim] <- Fmax
            
      }    	
      
      uu <- runif(totsim)
      
      temp <- ranklex(LMCstats,uu) 
      LMCpvalueFmax <- (totsim - temp + 1)/totsim
      
      temp <- ranklex(BMCstats,uu) 
      BMCpvalueFmax <- (totsim - temp + 1)/totsim    
      
      if(printRes==1){
            print('============================================')
            print('============================================')
            
            print('F-Max')
            print(c('------','F-max           =',Fmax_actual))
            
            if (LMCpvalueFmax > 0.05 ) {
                  print(c('ACCEPT ----','LMCp-value       =', LMCpvalueFmax))
            }
            
            if (BMCpvalueFmax <= 0.05 ) {
                  print(c('REJECT ----','BMC p-value      =', BMCpvalueFmax))
            }
            
            if ( (LMCpvalueFmax <= 0.05 ) & (BMCpvalueFmax > 0.05 ) ) {
                  print(c('INCONCLUSIVE ----','LMC p-value      =', LMCpvalueFmax))
                  print(c('INCONCLUSIVE ----','BMC p-value      =', BMCpvalueFmax))
            }    
      }
      
      
      c(SpanLMCpvalueFmax=LMCpvalueFmax,SpanBSpanMCpvalueFmax=BMCpvalueFmax)
      
}

f_many_simulation <- function(nb,N,K,ct,alpha_0,Beta_0,Sx,Sy,Pi,ids,ncp,garchmodel,list.rep,list.rep0) {
      #browser()
      seed <- list.rep[nb] + list.rep0[nb]
      out <- vector("list", list.rep[nb])
      for (i in 1:list.rep[nb]) {
            out[[i]] <- f_one_simulation(i+seed,N,K,ct,alpha_0,Beta_0,Sx,Sy,Pi,ids,ncp,garchmodel)
      }
      #browser()
      return(out) 
}

f_one_simulation <- function(seed,N,K,ct,alpha_0,Beta_0,Sx,Sy,Pi,ids,ncp,garchmodel) {
      #browser()
      set.seed(seed)
      z <- sapply(seq(K), function(i) c(garchm(ct = ct, cseed = seed, model = garchmodel))) %*% chol(Sx) #test set
      y <- alpha_0 + z %*% Beta_0 + sapply(seq(N), function(i)c(garchm(ct = ct, cseed = seed * 100000, model = garchmodel))) %*% chol(Sy) #Benchmark set
      x <- cbind(z[,1], z[,-1] + z[,1])
      
      one = rep(1, nrow(y))
      
      prods<-function(score,k){
            if(k>0){
                  a=sapply(seq(k), function(k)rnorm(nrow(score),1,1))
                  apply(as.matrix(a), 1, prod)
            }else{
                  1
            }
      }
      
      getpv<-function(u,x){
            require(sandwich)
            R1=u
            R2=x
            # y=y[,1]-x[,1]
            y=u-x[,1]
            x=cbind(x[,1],x[,-1]-x[,1])
            
            reg<-lm(y~x)
            regA<-lm(x[,1]~cbind(y,x[,-1]))
            regB<-lm(one~cbind(y,x)-1)
            reg3<-lm(y~x-1)
            regC<-lm(x[,1]~cbind(y,x[,-1])-1)
            # sn=cD
            
            forsn<-function(sn){
                  ########
                  ew=(resid(regA))
                  ew1=(resid(regB))
                  
                  ####################################################
                  res=(resid((reg)))
                  score2=score=as.matrix(cbind(res*ew))
                  
                  ########
                  ew3=(resid(regC))
                  ####################################################
                  res3=(resid(reg3))
                  score3=as.matrix(cbind(res3*ew3))
                  
                  #delta=0
                  scoreb=score*prods(score,sn)
                  # scoreb=(score)
                  scoretest=c(test1=testbm(scoreb,k=1/3)[1],test1=testbm(scoreb,k=1/2)[1],test2=testbm(scoreb,k=2/3)[1]
                  )
                  #delta=alpha=0
                  score=as.matrix(cbind(score2,res*ew1))
                  scoreb=score*prods(score,sn)
                  # scoreb=(score)
                  scoretest2=c(test1=testbm(scoreb,k=1/3)[1],test1=testbm(scoreb,k=1/2)[1],test2=testbm(scoreb,k=2/3)[1]
                  )
                  #delta=0/alpha=0
                  score=as.matrix(score3)
                  scoreb=score*prods(score,sn)
                  # scoreb=(score)
                  scoretest3=c(test1=testbm(scoreb,k=1/3)[1],test1=testbm(scoreb,k=1/2)[1],test2=testbm(scoreb,k=2/3)[1]
                  )
                  #alpha=0
                  score=as.matrix(res*ew1)
                  scoreb=score*prods(score,sn)
                  # scoreb=(score)
                  scoretest3b=c(test1=testbm(scoreb,k=1/3)[1],test1=testbm(scoreb,k=1/2)[1],test2=testbm(scoreb,k=2/3)[1]
                  )
                  #delta=0/alpha=0 and alpha=0
                  scoretest4=cbind(scoretest3,scoretest3b)
                  scoretest4=apply(scoretest4, 1,cauchypv)
                  # )
                  
                  c(A=scoretest,#delta=0
                    B=scoretest2,#delta=alpha=0
                    C=scoretest3,#delta=0/alpha=0
                    D=scoretest3b,#alpha=0
                    E=scoretest4#delta=0/alpha=0 and alpha=0
                  )
            } 
            
            lapply(seq(0,3), forsn)
            # forsn(1)
      }
      
      agetpv<-function(u,x){
            output=apply(u, 2, getpv,x)
            
            apply(output, 1, cauchypv)
      }
      
      output0 <- apply(y, 2, getpv,x)
      output0 <- lapply(seq(4), function(i) sapply(seq(output0),function(j)output0[[j]][[i]] ))
      output0 <- lapply(seq(4), function(i) as.matrix(apply(output0[[i]], 1, cauchypv)))
      c(output0,GungorLungerPesaranY(y,x,seed = seed),list(SpanningHKGL(y,x)))
}  


