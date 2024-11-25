#' Script that generates the tables for the size and plot for the power
#' 2024-11-23 - Ardia & Sessinou

##########################################################################################
##########################################################################################
## GENERATE THE TABLES FOR THE SIZE

f_sel <- function(K, N, type) {
      if (type == 1) {
            nam <- c("HK", "GL", "CCT,L=0", "CCT,L=1", "CCT,L=2")
            sel <- c("HK", "SpanningGL", "CCTdaL0", "CCTdaL1", "CCTdaL2")
            if (N == 400) {
                  nam <- c("GL", "CCT,L=0", "CCT,L=1", "CCT,L=2")
                  sel <- c("SpanningGL", "CCTdaL0", "CCTdaL1", "CCTdaL2")
            }
      }
      if (type == 3) {
            nam <- c("GRS", "F1", "BJ", "PY", "GL", "CCTa,L=0", "CCTa,L=1", "CCTa,L=2")
            sel <- c("GRS", "F1", "MK_BJ2", "PY", "mveGL", "CCTaL0", "CCTaL1", "CCTaL2")
            if (N == 400) {
                  nam <- c("GL", "PY", "CCTa,L=0", "CCTa,L=1", "CCTa,L=2")
                  sel <- c("mveGL", "Ans2", "CCTaL0", "CCTaL1", "CCTaL2")
            }
      }
      if (type == 4) {
            nam <- c("KM", "F2", "CCTd,L=0", "CCTd,L=1", "CCTd,L=2")
            sel <- c("MK_BJ1", "F2", "CCTdL0", "CCTdL1", "CCTdL2")
            if (N == 400) {
                  nam <- c("CCTd,L=0", "CCTd,L=1", "CCTd,L=2")
                  sel <- c("CCTdL0", "CCTdL1", "CCTdL2")
            }
      }
      out <- list(nam = nam, sel = sel)
      return(out)
}

# Setup in the paper
list_D <- as.list(1:12); n_D <- length(list_D)
list_K <- list(2, 10, 50, 100); n_K <- length(list_K)
list_N <- list(2, 10, 50, 100, 400); n_N <- length(list_N)
type <- 1
do_sparse <- TRUE

# Setup to ga fast and test
list_D <- as.list(1); n_D <- length(list_D)
list_K <- list(2, 10); n_K <- length(list_K)
list_N <- list(2, 10); n_N <- length(list_N)
type <- 1
do_sparse <- TRUE

tmp <- f_sel(K = 2, N = 2, type)
nam <- tmp$nam

n_T <- length(nam)
tmp_tbl <- array(NA, dim = c(n_K, n_N, n_D, length(nam)))
dimnames(tmp_tbl) <- list(paste0("K=", unlist(list_K)), 
                          paste0("N=", unlist(list_N)),
                          paste0("D=", unlist(list_D)),
                          nam)

for (iD in 1:n_D) {
      tmp1 <- c()
      for (iK in 1:n_K) {
            for (iN in 1:n_N) {
                  D <- list_D[[iD]]
                  K <- list_K[[iK]]
                  N <- list_N[[iN]]
                  if (do_sparse) {
                        nam_file <- paste("SparseWyRobStudModel", D, "K", K, "N", N, sep = "_")
                  } else {
                        nam_file <- paste("DenseWyRobStudModel", D, "K", K, "N", N, sep = "_")
                  }
                  
                  tmp <- read.csv(file = paste0("_csv/", nam_file, ".csv"), header = TRUE)
                  tmp <- tmp[,-1]
                  X <- as.matrix(tmp)
                  
                  tmp <- f_sel(K, N, type)
                  
                  tmp_tbl[iK,iN,iD,tmp$nam] <- X[9,tmp$sel]
            }
      }
}

tbl <- c()
for (iD in 1:n_D) {
      tmp <- c()
      for (iT in 1:n_T) {
            tmp <- cbind(tmp, tmp_tbl[,,iD,iT])
      }
      
      tbl <- rbind(tbl, tmp)
}

tbl
colnames(tbl) <- paste0(rep(nam, each = n_N), paste0(",N=", unlist(list_N)))
rownames(tbl) <- paste0(rep(paste0("D=",unlist(list_D)), each = n_K), paste0(",K=", unlist(list_K)))
tbl
xtable::xtable(tbl)
dim(tbl)
# remove L1
tbl <- tbl[,-c(16:20)] # setup 1
tbl <- tbl[,-c(31:35)] # setup 2
tbl <- tbl[,-c(16:20)] # setup 3
tbl

tbl2 <- as.data.frame(tbl)
for (i1 in 1:nrow(tbl)) {
      for (i2 in 1:ncol(tbl)) {
            tmp <- tbl[i1,i2]
 
            if ((tmp > 0.07) || (tmp < 0.03) || is.na(tmp)) {
                  #tbl2[i1,i2] <- format(tbl[i1,i2], digits = 3, nsmall = 3)
                  tbl2[i1,i2] <- format(100 * tbl[i1,i2], digits = 1, nsmall = 1)
            } else {
                  #tbl2[i1,i2] <- paste0("$\\grbb{", format(tbl[i1,i2], digits = 3, nsmall = 3), "}$")   
                  tbl2[i1,i2] <- paste0("$\\grbb{", format(100 * tbl[i1,i2], digits = 1, nsmall = 1), "}$")    
            }
      }
}
tbl2
xtable::xtable(tbl2)

print(xtable::xtable(tbl2), include.rownames = TRUE, 
      sanitize.text.function = function(x){x}, file = "out_tbl.txt")

stop("break here")

##########################################################################################
##########################################################################################
## GENERATE THE PLOTS FOR THE POWER

library("ggplot2")
library("gridExtra")
library("grid")
library("cowplot")

f_plot_1 <- function(X, i, K, N, type) {
      if (type == 1) {
            nam <- c("CCT,L=2", "HK", "GL")
            sel <- c("CCTdaL2", "HK", "SpanningGL")
            # N+K+1 < T (T=250) => tout passe, sinon GL, CCT, PY
            if (N == 400) {
                  nam <- c("CCT,L=2", "HK", "GL")
                  sel <- c("CCTdaL2", "NA", "SpanningGL")
            }
      }
      if (type == 2) {
            sel <- c("CCTdaL0", "CCTdaL1", "CCTdaL2", "SpanningGL")
            nam <- c("CCT,L=0", "CCT,L=1", "CCT,L=2", "GL")
      }
      if (type == 3) {
            nam <- c("CCTa,L=2", "GRS", "F1", "BJ", "PY", "GL")
            sel <- c("CCTaL2", "GRS", "F1", "MK_BJ2", "PY", "mveGL")
            if (N == 400) {
                  nam <- c("CCTa,L=2", "GRS", "F1", "BJ", "PY", "GL")
                  #nam <- c("CCTa,L=2", "PY", "GL")
                  sel <- c("CCTaL2", "NA", "NA", "NA", "Ans2", "mveGL")
            }
      }
      if (type == 4) {
            nam <- c("CCTd,L=2", "KM", "F2")
            sel <- c("CCTdL2", "MK_BJ1", "F2")
            if (N == 400) {
                  #nam <- c("CCTd,L=2", "CCTda,L=2")
                  nam <- c("CCTd,L=2", "KM", "F2")
                  sel <- c("CCTdL2", "NA", "NA")
            }
      }
      
      Xb <- cbind(X, NA)
      colnames(Xb) <- c(colnames(X), "NA")
      X <- Xb[,sel]
      colnames(X) <- nam
      rownames(X) <- seq(from = -.4, to = .4, by = .05)
      df <- reshape2::melt(X, c("alpha", "test"), value.name = "rejection")
      
      str_dgp <- c("iid-G", "iid-Std", "iid-sStd", "GARCH-G", 
                   "GARCH-Std", "GARCH-sStd", "AR-GARCH-G", 
                   "AR-GARCH-Std", "AR-GARCH-sStd", "AR-G", "AR-Std", "AR-sStd")
      
      str_title <- paste0(str_dgp[i], ";K=", K, ";N=", N)
      str_title <- paste0("K=", K, ";N=", N)
      
      p0 <- ggplot(data = df, aes(x = alpha, y = rejection, group = test)) +
            geom_line(aes(linetype = test, col = test), linewidth = 0.8) +
            geom_point(aes(shape = test, col = test), size = 2)  + 
            geom_hline(yintercept = 0.05) + 
            labs(subtitle = str_title) + 
            xlab("") + ylab("") + labs(fill = "")
      
      p1 <- p0 + theme(legend.position = "bottom", 
                       legend.title = element_blank(),
                       #legend.title = element_text(size = 15),
                       legend.text = element_text(size = 12),
                       panel.background = element_rect(fill = "white"), 
                       panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
      
      p2 <- p0 + theme(legend.title = element_blank(),
                       legend.position = "none", panel.background = element_rect(fill = "white"), 
                       panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
      
      out <- list(p0 = p0, p1 = p1, p2 = p2)
}

f_plot_2 <- function(type, do_sparse) {
      
      dgp <- 9
      
      list_KN <- list(c(K=2, N=2), c(K=2, N=50), c(K=2, N=400),
                      c(K=10, N=2), c(K=10, N=50), c(K=10, N=400),
                      c(K=100, N=2), c(K=100,N=50), c(K=100, N=400))
      
      if (do_sparse) {
            pdf(file = paste0("fig_sparse_type_", type, ".pdf"), width = 15, height = 10)      
      } else {
            pdf(file = paste0("fig_dense_type_", type, ".pdf"), width = 15, height = 10)      
      }
      p <- vector("list", 9)
      for (i in 1:9) {
            
            K <- list_KN[[i]]["K"]
            N <- list_KN[[i]]["N"]
            
            if (do_sparse) {
                  nam_file <- paste("SparseWyRobStudModel", dgp, "K", K, "N", N, sep = "_")     
            } else {
                  nam_file <- paste("DenseWyRobStudModel", dgp, "K", K, "N", N, sep = "_")
            }
            #browser()
            print(nam_file)
            tmp <- read.csv(file = paste0("_csv/", nam_file, ".csv"), header = TRUE)
            tmp <- tmp[,-1]
            X <- as.matrix(tmp)
            tmp <- f_plot_1(X, i = dgp, K = K, N = N, type = type)
            p[[i]] <- tmp
      }
      
      df_leg <- get_legend(p[[1]]$p1)
      
      grid.arrange(p[[1]]$p2, p[[2]]$p2, p[[3]]$p2, 
                   p[[4]]$p2, p[[5]]$p2, p[[6]]$p2,
                   p[[7]]$p2, p[[8]]$p2, p[[9]]$p2, 
                   arrangeGrob(nullGrob(), df_leg, nullGrob(), nrow = 1), 
                   ncol = 3, nrow = 4, heights = c(4, 4, 4, .5))
      dev.off()
}

type <- 1; do_sparse <- TRUE
f_plot_2(type, do_sparse)

type <- 3; do_sparse <- FALSE
f_plot_2(type, do_sparse)

type <- 4; do_sparse <- FALSE
f_plot_2(type, do_sparse)

stop("break here")



