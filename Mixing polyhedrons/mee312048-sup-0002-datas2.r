###Mixing Polygon Simulation (3 isotopes); v1.0, 15.2.2013 
###J.A. Smith et. al. 2013 (The University of NSW)###
###For updates and instructions, visit http://www.famer.unsw.edu.au/downloads.html

 ##Housekeeping
  rm(list=ls(all=TRUE))
  graphics.off()
  setwd("C:/R_scripts/MixPolySim")
  install.packages("geometry")  
  install.packages("ptinpoly")  
  library(geometry)
  library(ptinpoly)
  
 ##USER INPUT (a minimum of 4 sources is required)
  sources <- read.table("Sources_example_3i.csv",header=T,sep=",")  #always put in order of x(e.g. 13C), y(15N), z(iso3)
  mixture <- read.table("Mixture_example_3i.csv",header=T,sep=",")
  TEF <-  read.table("TEF_example_3i.csv", header=T,sep=",")
  its <- 1500  #specify the number of iterations ("its")
    
 ##Now RUN the simulation 
  Par_values <- array(0, c(its,(nrow(sources)*6+3)))  
  p <- array(0, c(its,(nrow(mixture))))
  for (i in 1:its) {  
    v <- array(0, c(nrow(sources),3))
    f <- array(0, c(nrow(TEF),3))
        for (j in 1:nrow(sources)) {
        v[j,1] <- rnorm(1, mean=sources[j,1], sd=sources[j,2])  
        v[j,2] <- rnorm(1, mean=sources[j,3], sd=sources[j,4])  
        v[j,3] <- rnorm(1, mean=sources[j,5], sd=sources[j,6])  
        f[j,1] <- rnorm(1, mean=TEF[j,1], sd=TEF[j,2])  
        f[j,2] <- rnorm(1, mean=TEF[j,3], sd=TEF[j,4])  
        f[j,3] <- rnorm(1, mean=TEF[j,5], sd=TEF[j,6])  
        }
    V <- v+f
    hulln <- convhulln(V, options="i") 
    hulln_v <- as.vector(hulln)
    hulln_i <- unique(hulln_v)
    hulln_i <- sort(hulln_i)
    indices <- V[hulln_i,]
    hulln_c <- convhulln(indices, options="i") #done again to match-up row indexes
    mixture_m <- as.matrix(mixture)
    P <- pip3d(indices, hulln_c, mixture_m)
    p[i,] <- P
    hulln_a <- convhulln(V, options="FA")
    poly_v <- hulln_a$vol  
    vals <- c(v[,1],v[,2],v[,3],f[,1],f[,2],f[,3],0,0,0) 
    Par_values[i,] <- vals  
    Par_values[i,ncol(Par_values)-2] <- poly_v
    Par_values[i,ncol(Par_values)-1] <- i
    Par_values[i,ncol(Par_values)] <- var(Par_values[1:i,ncol(Par_values)-2])
    if (i %% 10 ==0) cat(paste("iteration", i, "\n")) 
    }
  
 ##FIGURE 1: variance of polygon volume during simulation
   Iterations <- Par_values[,ncol(Par_values)-1]
   Variance <- Par_values[,ncol(Par_values)]
   plot(Iterations, Variance, type="n")   
   lines(Iterations, Variance, lty=1, lwd=1.5, col="blue")
 ##FIGURE 2: proportion of iterations that each consumer was inside mixing polygon
   p[p == 0] <- 1   
   p[p == -1] <- 0   
   Probabilities <- colSums(p)/its
   print(Probabilities)
   windows() 
   barplot(Probabilities, xlab="Consumer", ylab="Probability consumer is within 3-d mixing polygon", 
           ylim=c(0,1), names.arg=seq(1,nrow(mixture),by=1))
      
 ##Write data to files
   p_a <- rbind(p,Probabilities)
   write.table(p_a, file = "Consumer_Probabilities_3i.csv", sep = ",", row.names=FALSE)
   col_names <- c(rep("d13C",nrow(sources)),rep("d15N",nrow(sources)),rep("d34S",nrow(sources)),
                  rep("13C_TEF",nrow(sources)),rep("15N_TEF",nrow(sources)),rep("34S_TEF",nrow(sources)),
                  "Poly_Vol","Iteration","Variance") 
   col_nums <- c(rep(1:nrow(sources),6),0,0,0)
   col_n <- paste(col_names, col_nums)
   write.table(Par_values, file = "Parameter_Values_3i.csv", sep = ",", row.names=FALSE, col.names=col_n)
  

