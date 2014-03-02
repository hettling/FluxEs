.Rref <- function(A, tol=sqrt(.Machine$double.eps),verbose=FALSE,fractions=FALSE){
  ## A: coefficient matrix
  ## tol: tolerance for checking for 0 pivot
  ##verbose: ifTRUE, print intermediate steps
  ## fractions: try to express nonintegers as rational numbers                 
  ## Written by John Fox                       
  if (fractions) {
    mass <- require(MASS)
    if (!mass) stop("fractions=TRUE needs MASS package")
  }
  if ((!is.matrix(A)) || (!is.numeric(A)))
    stop("argument must be a numeric matrix")
  n <- nrow(A)
  m <- ncol(A)
  for (i in 1:min(c(m, n))){
    col <- A[,i]
    col[1:n < i] <- 0
    ## find maximum pivot in current column at or below current row                 
    which <- which.max(abs(col))
    pivot <- A[which, i]
    if (abs(pivot) <= tol) next     # check for 0 pivot                        
    if (which > i) A[c(i, which),] <- A[c(which, i),]  # exchange rows         
    A[i,] <- A[i,]/pivot            # pivot    
    row <- A[i,]
    A <- A - outer(A[,i], row)      # sweep    
    A[i,] <- row                    # restore current row                      
    if (verbose)
      if (fractions) print(fractions(A))
      else print(round(A,round(abs(log(tol,10)))))
  }
  for (i in 1:n)
    if (max(abs(A[i,1:m])) <= tol)
      A[c(i,n),] <- A[c(n,i),] # 0 rows to bottom                              
  if (fractions) fractions (A)
  else round(A, round(abs(log(tol,10))))
}



CalcFluxVector <- function(S, independent.fluxes=vector()){  
  ##Check for impossible combinations of independent fluxes
  row.indices <- which(rowSums(S)==0)
  flux.comb <- lapply(row.indices, function(x)names(which(S[x,]!=0)))
  for (i in seq_along(flux.comb)){
    idx <- which(flux.comb[[i]] %in% names(independent.fluxes))
    if (length(idx)>1){
      cat("Fluxes ", paste(flux.comb[[i]][idx], collapse=", "), "cannot be independent fluxes together.\n")
    }
  }
  qq <- qr(S)
  ##If no independent fluxes are given, choose automatically and initialize with random numbers
  if (length(independent.fluxes) == 0){
    independent.fluxes <- runif(ncol(S) - qq$rank, min=0, max=2000) ##caution, magic number
    names(independent.fluxes) <- setdiff(colnames(S), colnames(stats:::Thin.col(S)))
    cat("No independent fluxes given. Chosing the following independent fluxes: ", names(independent.fluxes), "\n")    
    cat("Initializing with random values, ", independent.fluxes, "\n")
  }
  else if (length(independent.fluxes) < (ncol(S) - qq$rank)){
    cat("Stoichiometry Matrix has rank ", qq$rank, " with ", ncol(S), " fluxes.\n")
    cat("You need ", ncol(S) - qq$rank, " independent fluxes to calculate the flux vector.\n")
    cat("Suggested independent fluxes : ", setdiff(colnames(S), colnames(stats:::Thin.col(S))), "\n")
    return(vector())
  }
  ##put columns of known fluxes to the end of Stoichiometry matrix
  unknown.fluxnames <- setdiff(colnames(S), names(independent.fluxes))

  S.perm <- cbind(S[,unknown.fluxnames], S[,names(independent.fluxes)])
  S.perm.echelon <- .Rref(S.perm)
  ##The echelon form for all unknown fluxes must now be equal to the identity matrix
  ##stopifnot(all(S.perm.echelon[,unknown.fluxnames] == diag(nrow=nrow(S.perm), ncol=length(unknown.fluxnames))))
  v.unknown <- as.vector(-S.perm.echelon[,names(independent.fluxes)] %*% independent.fluxes)
  names(v.unknown) <- unknown.fluxnames
  v <- c(v.unknown, independent.fluxes)[colnames(S)]
  return (v)
}



CheckSimulationOutputConsistency <- function(output, model, tolerance=1e-4){  
  ret <- TRUE
  for (s in model@species){
    if (! all(rowSums(output[,paste(s@id, 1:(2^s@carbon.count), sep="_")])-1 < tolerance)){
      cat("Simulation output of isotopomers for species ", s@id, " does not add up to one.\n")
      ret <- FALSE
    }
  }
  return(ret)
}

CalcMassIsotopomerMatrix <- function(output, model, timepoints=NA){   
  ##Calculate MS multiplets from lsoda output and store them in a matrix         
  ms_multiplet_matrix <- matrix(nrow=length(output[,1]), ncol=0)       
  for (s in model@species){  
    ##create matrix containing zeros each MS multiplet has a column, each time a row       
    m <- matrix(0, ncol=s@carbon.count+1, nrow=length(output[,1]))       
    colnames(m) <- paste(s@id, "_", 0:s@carbon.count, sep="") 
    for (c in 1:2^s@carbon.count){ ##don't consider first multiplet because it has no label  
      ##get the time series for specific isotopomer
      time_series <- output[,paste(s@id, "_", c, sep="")]   
      ##look how many carbons are labeled:         
      labeled_carbon_count <- sum(wle::binary(c-1, dim=s@carbon.count)$binary)
      m[,labeled_carbon_count+1] <- m[,labeled_carbon_count+1] + time_series     
    }      
    ms_multiplet_matrix <- cbind(ms_multiplet_matrix, m)        
  }        
  ms_multiplet_matrix <- cbind(output[,"time"], ms_multiplet_matrix)   
  colnames(ms_multiplet_matrix)[1] <- "time"
  if(length(timepoints)>0){
    indices <- apply(as.matrix(timepoints), 1, function(x) which(ms_multiplet_matrix[,1]==x))
    ms_multiplet_matrix <- ms_multiplet_matrix[indices,]
  }
  return(abs(ms_multiplet_matrix))
}

CalcExtremePathways <- function(model){
  S <- stoichiometryMatrix(model)
  ##initialize matrix
  T0 <- cbind(diag(sum(sapply(model@parameters, function(p)p@is.flux))), t(S))
  influxes <- unlist(sapply(model@parameters, function(p)if(p@external.influx)p@id))
  ##For influxes ("external fluxes constrained to be negative" in Schilling et al.), multiply
  ##  corresponding row of matrix with -1
  for (influx in influxes)
    T0[influx,] <- T0[influx,] * -1
  
  TX.before <- T0

  ##1. identify metabolites that are not associated with an exchange flux
  external.fluxnames <- unlist(sapply(model@parameters, function(p)if (p@external.influx || p@external.outflux)p@id))
  int.met <- rownames(S)[which(rowSums(S[,external.fluxnames])==0)]

  for (curr.met in int.met){
  
    ##2.Copy from TX.before all rows which have a zero in the columns for int.met
    TX.current <- matrix(nrow=0, ncol=ncol(T0))
    TX.current <- rbind(TX.current, TX.before[which(TX.before[, curr.met] == 0),])
    
    ##3. Of the remaining rows in TX.before, add all combinations that sum to 0 for the current metabolite
    all.combinations <- combn(nrow(TX.before),2)
    apply(all.combinations, 2, function(x){
      if(sum(TX.before[c(x[1], x[2]),curr.met])==0 && TX.before[x[1],curr.met]!=0){
        TX.current <<- rbind(TX.current, TX.before[x[2],] * abs(TX.before[x[1],curr.met]) + TX.before[x[1],] * abs(TX.before[x[2],curr.met]))##eq B.2
        return(0)}
      
    })  
    ##4: Check if rows can be nonnegative combination of each other, if yes, eliminate corresponding row.
    dependent.rows.indices <- vector()
    for (i in 1:nrow(TX.current)){
      a.i <- which(TX.current[i,]==0)
      for (h in 1:nrow(TX.current)){
        a.h <- which(TX.current[i,]==0)
      }
      if (identical(intersect(a.i, a.h), a.i))##Eq. B.4
        dependent.rows.indices <- append(dependent.rows.indices, i)
    }
    TX.current <- TX.current[setdiff(seq_len(nrow(TX.current)), dependent.rows.indices),] ##remove dependent rows
    
    ##5: Repeat for all mu metabolites
    TX.before <- TX.current 
  }

  TX <- TX.current##rbind(TX.current, TE)
  TE <- TX
  
  ##6: Start from the first non-zero column on the right side of the tableau (the columns which have metabolite names)
  metabolite.names <- colnames(TX)[which(colnames(TX)!="")]
  
  for (curr.met in metabolite.names){
    ##7. 
    nonzero.column.idx <- which(TX[,curr.met]!=0)
    for (idx in nonzero.column.idx){
      ##Search for corresponding row in TE
      corr.ext.row <- which(TE[,curr.met] != 0)[1]    
      ##Multiply corresponding row in TE with the value in the cell for the metabolite in the current row in TX
      TX[idx,] <- TX[idx,]  + TX[idx,curr.met] * TE[corr.ext.row,] 
    }
  }  
  extreme.pathways <- TX[,1:(no.internal.fluxes + no.external.fluxes)]
  return(extreme.pathways)
}


GetIndependentFluxes <- function(S){
  return(sort(setdiff(colnames(S), colnames(stats:::Thin.col(S)))))
}


CalcFluxVectorFromLinearBasis <- function(independent.fluxes, S=S, ...){
  independent.fluxnames <- names(independent.fluxes)
  dependent.fluxnames <- setdiff(colnames(S), independent.fluxnames)

  ##If no names for independent fluxes are given they are automatically determined
  if (length(independent.fluxnames)==0){
    independent.fluxnames <- GetIndependentFluxes(S)    
    cat("Names attribute for independent fluxes missing.\n Setting ",
        independent.fluxnames, " as independent.\n")
    if (length(independent.fluxnames)!=length(independent.fluxes))
      stop("Number of independent fluxes must be ncol(S)-rank(S).\n")
    
    names(independent.fluxes) <- independent.fluxnames
  }
  ##Permutate columns in S such that the independent fluxes are the last columns
  S.perm <- cbind(S[,dependent.fluxnames, drop=FALSE], S[,independent.fluxnames, drop=FALSE])  
  ##Calculate row-reduced echelon form of the  matrix S.perm
  S.perm.echelon <- .Rref(S.perm)
  ##Calculate dependent fluxes by matrix multiplication
  dependent.fluxes <- as.vector(-S.perm.echelon[,independent.fluxnames, drop=FALSE] %*% independent.fluxes)
  names(dependent.fluxes) <- dependent.fluxnames
  v <- c(dependent.fluxes, independent.fluxes)[colnames(S)]
  ##if (min(v)<0)
    ##cat("Warning: Linear basis yields negative flux values!\n")
  return(v)  
}

CalcFluxVectorFromLinp <- function(independent.fluxes, S=S, ...){
  dependent.fluxnames <- colnames(stats:::Thin.col(S))
  independent.fluxnames <- setdiff(colnames(S), dependent.fluxnames)
  A <- S
  colnames(A) <- colnames(S)
  A[,dependent.fluxnames] <- matrix(rep(0, length(dependent.fluxnames) * nrow(S)), ncol=length(dependent.fluxnames))
  A[,independent.fluxnames] <- diag(length(independent.fluxnames))
  b <- independent.fluxes
  E <- S
  f <- rep(0, nrow(S))
  G <- diag(nrow=ncol(S))
  h <- rep(0, ncol(S))
  sdB <- rep(1e-6, nrow(S))
  xx <- xsample(A=A, B=b, E=E, F=f, G=G, H=h, sdB=sdB, iter=2)
  return(xx$X[1,])
}


####################functions right now specific for Schilling model


simulate <- function(v){
  start.values <- unlist(CalcIsotopomerStartValues(model, labeled.input))  
  ##Get the sizes of the metabolite pools:
  pool.sizes <- sapply(model@species, function(s)s@initialConcentration)  
  output <- lsodes(start.values, times, "derivsc", parms=c(pool.sizes, v), dllname=model@id, nout=0) 
  colnames(output) <- c("time", iso.names)
  return(output)
}

chi.square <- function(v){
  if (min(v)<0 | sum(S %*% v)>1e-5){    
    return(1e6)
  }
  ##Get the start values for the simulation. External glucose is  labeled at the first carbon
  start.values <- unlist(CalcIsotopomerStartValues(model, labeled.input))  
  ##Get the sizes of the metabolite pools:
  pool.sizes <- sapply(model@species, function(s)s@initialConcentration)  
  output <- lsodes(start.values, times, "derivsc", parms=c(pool.sizes, v), dllname=model@id, nout=0) 
  colnames(output) <- c("time", iso.names)  
  ##Check if the sum of all isotopomers for each metabolite is 1 at each time point
  ##  cat("Consistent : ", CheckSimulationOutputConsistency(output, model), "\n")  
  sim.isotopomer.data <- CalcMassIsotopomerMatrix(output, model, timepoints=mock.data[,"time"])
  exp.data <- mock.data[,-which(colnames(mock.data)=="time")]
  sim.data <- sim.isotopomer.data[,colnames(exp.data)]
  if (data.sd>0)
    chisq.vec <- ((exp.data - sim.data)^2) / (data.sd^2)
  else
    chisq.vec <- ((exp.data - sim.data)^2) 
  chisquare.cost <- sum(chisq.vec)
  ##cat("Cost = ", chisquare.cost, "\n")
  return(chisquare.cost) 
}


cost <- function(independent.fluxes, calc.flux.function, ...){
  ##cat("Params  : ", independent.fluxes, "\n")
  names(independent.fluxes) <- independent.fluxnames
  v <- calc.flux.function(independent.fluxes, ...)
  return(chi.square(v))
}


plot.ensemble <- function(ens.mass.isos){
  ##Plot calculated mass isotopomers before and after the fit together with the mock data:
  par(mfrow=c(3,3))
  for (sp in model@species){
    isos <- paste(sp@id, seq_len(sp@carbon.count+1)-1, sep="_")
    plot(mock.data[,"time"], mock.data[,isos[1]], type="n", main=sp@id, ylim=c(0, 1), xlab="time (s)", ylab="mass isotopomer fraction")
    colors <- rainbow(length(isos))
    for (i in seq_along(isos)){
      iso <- isos[i]
      ##plot isos
      for (j in 1:length(ens.mass.isos)){
        lines(ens.mass.isos[[j]][,"time"], ens.mass.isos[[j]][,iso], col=colors[i])      
      }
      ##plot data
      points(mock.data[,"time"], mock.data[,iso], type="p", main=iso)
      
    }
    legend("topright", legend=isos, fill=colors)
  }
}
