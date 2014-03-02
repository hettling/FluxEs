require('rsbml')

##extended SBML Reaction class
setClass("FluxEsReaction", contains="Reaction",
         representation(reactants.atom.positions="vector",
                        products.atom.positions="vector",
                        flux.name="character"))

##extended SBML Species class
setClass("FluxEsSpecies", contains="Species",
         representation(carbon.count="numeric", external="logical", input="logical"))

##constructor for FluxEsSpecies, external defaults to FALSE
setMethod("initialize", "FluxEsSpecies", function(.Object, ..., carbon.count=numeric(), external=FALSE, input=FALSE){
  callNextMethod(.Object, carbon.count=carbon.count, external=external, input=input, ...)
})

##Method to get a binary matrix with all possible isotopomers for a FluxEsSpecies
setGeneric("isotopomerMatrix", function(object){standardGeneric("isotopomerMatrix")})
setMethod("isotopomerMatrix", "FluxEsSpecies", function(object){
  if(object@carbon.count==1)
    return(rbind(0,1))
  else
    return(t(sapply(0:(2^object@carbon.count-1), function(x)wle:::binary(x, object@carbon.count)$binary)))               
})

##extended SBML parameter class
setClass("FluxEsParameter", contains="Parameter",
         representation(is.flux="logical",
                        external.influx="logical",
                        external.outflux="logical"))

##constructor for FluxEsParameter, 'is.flux' and 'external' default to FALSE
setMethod("initialize", "FluxEsParameter", function(.Object, ..., is.flux=FALSE, external.influx=FALSE, external.outflux=FALSE){
  callNextMethod(.Object, is.flux=is.flux, external.influx=external.influx, external.outflux=external.outflux, ...)
})


##extended SBML Model class
setClass("FluxEsModel", contains="Model",
         representation(stoichiometry.matrix="matrix"))

##Method to calculate the stoichiometry matrix of a FluxEsModel
setMethod("stoichiometryMatrix", signature("FluxEsModel"), function(object, omit.external.metabolites=TRUE, external.fluxes.to.back=TRUE, ...){
  S <- callNextMethod(object, ...)
  if (omit.external.metabolites){
    S <- S[unlist(sapply(object@species, function(x)if (!x@external)x@id)), ,drop=FALSE]
  }
  if (external.fluxes.to.back){
    ext.paramnames <- unlist(sapply(object@parameters, function(x)if(x@external.influx || x@external.outflux)x@id))
    S <- S[,c(setdiff(colnames(S), ext.paramnames), ext.paramnames), drop=FALSE]
  }
  S
})

##Method to calculate the ordinary differential equantions for all isotopomer in a FluxEsModel
setGeneric("calcODEs", function(object, type="all_isos"){standardGeneric("calcODEs")})
setMethod("calcODEs", "FluxEsModel", function(object, type="all_isos"){  
  if (! type %in% c("all_isos", "MS"))
    stop("Type of ODEs must be either 'all_isos' or 'MS'.")
  if (type == "all_isos")
    return(.CalcIsotopomerODEs(object))
  if (type == "MS")
    return(.CalcMsODEs(object))
})




##HANNES: TODO: CHange Input to given vector of atoms and not to the index if the Isotopomer!!!\
##Hannes: Do it in the same way as the CalcMSstartvals
##Method to calculate the fractional starting values for simulation for each isotopomer in a FluxEs model. By default, values for the
## isotopomers are calculated from the natural abundances of 12C and 13C. The user can define zero or more input isotopomers.
## The fraction of the input isotopomer is then set to one and to zero for all other isotopomers of this metabolite, respectively.
setGeneric("CalcStartValues", function(object, labeled.input){standardGeneric("CalcStartValues")})
setMethod("CalcStartValues", "FluxEsModel",
          function(object, labeled.input=list()){

            .GetNaturalAbundances <- function(sp)
              apply(isotopomerMatrix(model@species[[sp]]), 1, function(iso){
                ##natural abundancies: 12C has probability 0.989, 13C has probability 0.011
                prod(sapply(iso, function(x)if(x==0)0.989 else 0.011))
              })
            
            ##isotopomer start values to be returned are initially set to natural abundances
            start.values <- lapply(object@species, function(s).GetNaturalAbundances(s@id))
            
            ##set names             
            labeled.input <- lapply(labeled.input, function(x){              
                names(x)=c("species", "atoms", "fraction")[1:length(x)]              
              x
            })            
            ##attach corresponding vectors to carbon atoms  for the isotopomers
            labeled.input <- lapply(labeled.input, function(li){              
              isotopomer.matrix <- isotopomerMatrix(object@species[[li$species]])
              iso <- rep(0, ncol(isotopomer.matrix))
              iso[li$atoms] <- 1
              isotopomer.idx <- unlist(sapply(seq_len(nrow(isotopomer.matrix)), function(x)if(all(isotopomer.matrix[x,]==iso))x))
              res <- rep(0, nrow(isotopomer.matrix))
              res[isotopomer.idx] <- 1
              li$isotopomer.fractions <- res
              li
            })
            ##Now the tricky part: Input isotopomers can have a certain percentage or fraction,
            ## for e.g. 1-Glc 20%, 1,2,3,4,5,6-Glc 60% and 20% natural abundance. If the fractions of
            ## one input species do not add up to one, we assume the rest is labeled according to natural
            ## abundances.
            all.input.specnames <- unique(sapply(labeled.input, function(x)x$species))
            for (sp in all.input.specnames){              
              labeled.sp.input <- lapply(labeled.input, function(x)if(x$species==sp)x) ##list with input for specific species
              labeled.sp.input <- labeled.sp.input[!sapply(labeled.sp.input, is.null)] ##remove possible NULL elements
              iso.vecs <- do.call(rbind, lapply(labeled.sp.input, function(x)x$isotopomer.fractions))              
              fractions <- unlist(sapply(labeled.sp.input, function(x)ifelse(is.null(x$fraction), NA, x$fraction)))
              ##treat the  case that no fractions are assigned, make them uniform for all input isotopomers
              if (all(is.na(fractions))) 
                fractions <- rep(1/length(labeled.sp.input), length(labeled.sp.input))
              ##treat the case that a fraction is assigned to only one of multiple input isotopomers
              if (any(is.na(fractions)) & ! all(is.na(fractions)))
                fractions[which(is.na(fractions))] <-  (1-sum(fractions, na.rm=TRUE))/length(which(is.na(fractions)))              
              iso.vecs <- iso.vecs * fractions
              ##treat the case that not all fractions add up to one
              remaining.frac <- 1.-sum(fractions)
              if (remaining.frac>0) ##assume natural abundances for the rest
                iso.vecs <- rbind(iso.vecs, .GetNaturalAbundances(sp) * remaining.frac)
              ##sum vectors to get start values for the species
              start.values[[sp]] <- colSums(iso.vecs)
            }                                                                     
            return(start.values)
  
})


##"Private" functions used by the FluxEsModel class

CalcMsStartValues <- function(model, labeled.input){
  ms.start.vals <- lapply(model@species, function(s){
    unlist(lapply(s@carbon.count, function(i) choose(i, 0:i))) * 0.989^rev(seq(0, s@carbon.count)) * 0.011^(seq(0, s@carbon.count))
  })
  all.input.specs <- sapply(labeled.input, `[[`, 1)
  for (sp in unique(all.input.specs)){
    freqs <- lapply(labeled.input, function(x)if (x[[1]]==sp)x[[3]])
    
    labeled.atoms <- lapply(labeled.input, function(x)if (x[[1]]==sp)x[[2]])
    
    ##Following (mass) isotopomers have to be set according to the frequencies:
    idx <- lapply(labeled.input, function(x)if (x[[1]]==sp)length(x[[2]])+1) 
    
    sp.startvals <- rep(0, model@species[[sp]]@carbon.count+1)
    ##set the start values according to the fractions
    for (i in 1:length(idx)){
      sp.startvals[idx[[i]]] <- sp.startvals[idx[[i]]] + freqs[[i]]
    }
    remaining.freq <- 1 - sum(unlist(freqs)) ##set the remaining fraction to natural abundance
    if (remaining.freq > 1 | remaining.freq < 0)
      error("Fractions provided for lebeled input are incorrect.")
    if (remaining.freq > 0)
      sp.startvals <- sp.startvals + ms.start.vals[[sp]] * remaining.freq 
    ms.start.vals[[sp]] <- sp.startvals
  }
  return(ms.start.vals)
}

.CalcIsotopomerODEs <- function(model){
  ins <- .GetIsotopomerInfluxTerms(model)
  outs <- .GetIsotopomerOutfluxTerms(model)  
  iso.names <- unlist(sapply(model@species,
                             function(s)sapply(1:2^s@carbon.count,
                                               function(x)paste(s@id, x, collapse="_", sep="_"))))
  odes <- lapply(iso.names, function(iso){
    ode = ""
    pool.name <- unlist(strsplit(iso, "_"))[1]
    ode <- paste(ode, ifelse(ins[[iso]]=="", "", paste("( ", ins[[iso]], " )")), sep="")
    ode <- paste(ode, ifelse(outs[[iso]]=="", "", paste(" - ( ", outs[[iso]], " )")), sep="")
    ode <- ifelse(ode=="", "", paste("(", ode, ") / ",  pool.name, sep=""))
    ode <- ifelse(ode=="", "0", ode)
    ode <- paste(iso, "=", ode)
    ode    
  })
  names(odes) <- iso.names 
  return(odes)  
}

.CalcMsODEs <- function(model){
  allspecs.out <- list()
  allspecs.in <- list()
  for (sp in model@species){
    cat("Writing ODE for species ", sp@id, "\n")
    ##cat("Species : ", sp@id, "\n")
    sp.out <- as.list(rep("", sp@carbon.count+1))
    sp.in <- as.list(rep("", sp@carbon.count+1))
    ##for (i in seq_len(sp@carbon.count + 1)){ ##all mass isotopomers, watch out, it starts from 1
    if (sp@external){
      sp.in <- as.list(rep("0", sp@carbon.count + 1))
      sp.out = as.list(rep("0", sp@carbon.count + 1))
    } else {
      
        for (r in model@reactions){
          cat("\tProcessing reaction ", r@id, "\n")
          if (sp@id %in% lapply(r@reactants, id)){
              ##recover()
              for (i in seq_len(sp@carbon.count + 1)){
                  sp.out[[i]] <- ifelse(sp.out[[i]]=="", sp.out[[i]], paste(sp.out[[i]], "+"))
                  sp.out[[i]] <- paste(sp.out[[i]], paste(sp@id, i-1, sep="_"), " * ", r@flux.name)
              }
              ##cat("Sp.out[[i]] : ", sp.out[[i]], "\n")
          }
        if (sp@id %in% lapply(r@products, id)){          
          ##if (sp@id == "squalene" && r@id == "v31")
          ##  recover()          
          reactants <- r@reactants
          ##take stoichiometry of reactants into account
          reactants.withstoich <- unlist(lapply(reactants, function(x)replicate(x@stoichiometry, x)))
          ##carbon atoms for all reactants
          possible.reactant.mass.isos <- lapply(reactants.withstoich, function(x)seq_len(model@species[[x@id]]@carbon.count+1)-1)
          cat("Making all iso combinations\n")
          all.combns <- expand.grid(possible.reactant.mass.isos)
          rowsums <- rowSums(all.combns)
          cat("\n")
          cat("All combinations calculated\n")
          
          ## Iterate over all possible mass isotopomer of the product species
          ## Select all relevant combinations of reactant mas isotopomers that can
          ##   possibly result in the given product isotopomer.
          for (i in seq_len(sp@carbon.count + 1)){  
            cat("Processing mass isotopomer ", i, "\n")            
            all.relevant.combns <- all.combns[which(rowsums==i-1), , drop=FALSE] ##all rows where the reactants can create the current mass iso
            colnames(all.relevant.combns) <- lapply(reactants.withstoich, id)
            cat("Relevant combinations calculated\n")
            cat("Dim all relevant combns : ", dim(all.relevant.combns), "\n")
            ##
            a <- data.frame(v=apply(all.relevant.combns, 1, function(x) 
                              paste(sort(paste(colnames(all.relevant.combns), x, sep="_")), collapse=" * ")))
            d <- ddply(a, .(v), summarise, new = paste(ifelse(length(v)>1, paste(length(v), "*"), ""), unique(v)))
            cat("ddply done\n")
            ##
            ##tmpmat <- as.matrix(apply(all.relevant.combns, 1, function(x)paste(colnames(all.relevant.combns), x, sep="_")))
            ##cat("tmpmat done\n")
            ##mynewin <- ifelse(sp.in[[i]]=="", sp.in[[i]], paste(sp.in[[i]], "+"))
            new.in <- paste("(", d[, "new"], ") * ", r@flux.name)
            new.in <- paste(new.in, collapse=" + ")

            sp.in[[i]] <- ifelse(sp.in[[i]]=="", sp.in[[i]], paste(sp.in[[i]], "+"))
            sp.in[[i]] <- paste(sp.in[[i]], new.in)

            ##sp.in[[i]] <- ##paste(sp.in[[i]], paste(paste("(", apply(tmpmat, 2, paste, collapse=" * "), ") * ", r@flux.name), collapse=" + "))
            cat("sp.in[[i]] : ", sp.in[[i]], "\n")                                                
            ##mynewin <- paste(mynewin, new.in)
            ##cat("mynewin : ", mynewin, "\n")
          }
        }
      }    
    }
    allspecs.out[[sp@id]] <- sp.out
    allspecs.in[[sp@id]] <- sp.in
    ##}
  } 
  ms.odes <- as.list(mapply(function(x, y){
    x <- ifelse(x=="", x, paste("(", x, ")"))
    y <- ifelse(y=="", y, paste("(", y, ")"))
    if (y=="")
      x
    else
      paste(x, y, sep=" - ")
  }, unlist(allspecs.in), unlist(allspecs.out)))
  
  pool.names <- unlist(lapply(model@species, function(x)rep(x@id, x@carbon.count+1)))
  
  ms.names <- unlist(lapply(model@species, function(x)paste(x@id, seq_len(x@carbon.count+1)-1, sep="_")))
  ms.odes <- as.list(paste(ms.names, " = ", "(", ms.odes, ") /", pool.names))
  names(ms.odes) <- ms.names
  
  input.pools <- names(which(sapply(model@species, function(s)s@input)))
  ms.odes[which(pool.names %in% input.pools)] <- paste(ms.names[which(pool.names %in% input.pools)], " = 0")
  return(ms.odes)  
}

.GetIsotopomerInfluxTerms <- function(model){
  all.influxes <- list()
  for (r.counter in seq_along(model@reactions)){
    r <- model@reactions[[r.counter]]
    ##Get the names of reactants and products, make sure to get the stoichiometries right!
    reactant.names <- unlist(lapply(r@reactants, function(rr)rep(rr@id, rr@stoichiometry)))
    product.names <- unlist(lapply(r@products, function(rr)rep(rr@id, rr@stoichiometry)))
    for (i in seq_along(product.names)){
      product.pool <- product.names[i]
      current.product.positions <- as.numeric(charToRaw(r@products.atom.positions[[i]]))-96 ##Convert letters of atom positions to integers    
      ##Get the reactants which give at least one carbon atom to the current product pool
      relevant.reactants <- vector()
      relevant.reactant.positions <- list()
      for (j in seq_along(reactant.names)){
        reactant.positions <- as.numeric(charToRaw(r@reactants.atom.positions[[j]]))-96
        if (any(reactant.positions %in% current.product.positions)){
          relevant.reactants[j] <- reactant.names[[j]]
          relevant.reactant.positions[[j]] <- reactant.positions
        }
      }   
      influx.terms <- vector()
      ##Iterate over the reactant pools and get the relevant isotopomers for each reactant
      for (j in seq_along(relevant.reactants)){
        ##if (model@species[[relevant.reactants[j]]]@input==FALSE){
          relevant.reactant.isos <- .GetRelevantIsotopomers(relevant.reactants[j], product.pool, relevant.reactant.positions[[j]], current.product.positions)
          current.influx.terms <- apply(relevant.reactant.isos, 1, function(x)paste("(", paste(relevant.reactants[j], x, collapse=" + ", sep="_"), ")"))
          ##If there are multiple reactants, the relevant isotopomer for the influx have to be multiplied.
          if(length(influx.terms)==0)
            influx.terms <- current.influx.terms
          else
            influx.terms <- paste(influx.terms, current.influx.terms, sep=" * ") 
        }
        influx.terms <- paste(influx.terms, r@flux.name, sep=" * ")
        iso.names <- paste(product.pool, seq_len(2^model@species[[product.pool]]@carbon.count), sep="_")
        for (j in seq_along(iso.names)){
          all.influxes[[iso.names[j]]] = ifelse(length(all.influxes[[iso.names[j]]])==0, influx.terms[j], paste(all.influxes[[iso.names[j]]], influx.terms[j], sep=" + "))
        }
      ##}
    }
  }
  return(all.influxes)
}

.GetRelevantIsotopomers <- function(reactant.name, product.name, reactant.positions, product.positions){
  source.matrix <- as.matrix(isotopomerMatrix(model@species[[reactant.name]])[,which(reactant.positions %in% product.positions)])
  target.matrix <- as.matrix(isotopomerMatrix(model@species[[product.name]])[,which(product.positions %in% reactant.positions)])
  result <- t(apply(target.matrix, 1, function(x){
    rets <- vector()
    for (i in seq_len(nrow(source.matrix))){
      if (all(source.matrix[i,]==x)){      
        rets <- c(rets, i)
      }
    }
    rets
  }))
  ##If there is only one relevant reactant isotopomer per product isotopomer, we must create a matrix by transposing the vector returned by apply
  if (nrow(result)==1)result <- t(result) 
  ##recover()
  return(result)
}

.GetIsotopomerOutfluxTerms <- function(model){
  ##Collect all isotopomers for the outfluxes
  all.outflux.terms <- list()  
  for (s in model@species){
    source.isos <- paste(s@id, 1:(2^s@carbon.count), sep="_")
    if (s@id %in% rownames(stoichiometryMatrix(model)) && s@input == FALSE){ ##&& p@poolType != "constant_input"){
      relevant.fluxes <- abs(stoichiometryMatrix(model)[s@id,][which(stoichiometryMatrix(model)[s@id,]<0)])    
      for (iso in source.isos){    
        if (length(relevant.fluxes)>0)
          all.outflux.terms <- c(all.outflux.terms, paste(iso, " * (", paste(paste(relevant.fluxes, names(relevant.fluxes), sep= " * "), collapse=" + "), ")"))
        else
          all.outflux.terms <- c(all.outflux.terms,"")
      }
    } else {
      for (iso in source.isos){
        all.outflux.terms <- c(all.outflux.terms,"")
      }
    }
    
  } 
  names(all.outflux.terms) <- unlist(sapply(model@species,
                                            function(s)sapply(1:2^s@carbon.count,
                                                              function(x)paste(s@id, x, collapse="_", sep="_"))))
  return(all.outflux.terms)
}

CalcStoichiometryMatrix <- function(model, omit.externals=TRUE){
  flux.names <- unique(sapply(model@reactions, function(x)x@flux.name))
  S <- matrix(0, nrow=length(model@species), ncol=length(flux.names), dimnames=list(as.vector(sapply(model@species, function(s)s@id)), flux.names))
  for (r in model@reactions){
    for (i in seq_along(r@reactants))
      S[r@reactants[[i]]@id, r@flux.name] <- S[r@reactants[[i]]@id, r@flux.name] - r@reactants[[i]]@stoichiometry
    for (i in seq_along(r@products))
      S[r@products[[i]]@id, r@flux.name] <- S[r@products[[i]]@id, r@flux.name] + r@products[[i]]@stoichiometry
  }
  if (omit.externals){
    internal.species.names <- unlist(sapply(model@species, function(x)if (!x@external)x@id))
    S <- S[internal.species.names,]
  }
  return(S)
}
