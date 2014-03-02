require('stringr')

ParseModel <- function(file.name){
  all.lines <- readLines(file.name, warn = FALSE)
  all.lines <- all.lines[-which(all.lines=="")]
  
  ##parse for key string '##' to indicate whether equations, internals or externals are described
  idx <- which(! is.na(str_match(all.lines, "##.+")))

  eq.lines <- all.lines[(idx[1]+1):(idx[2]-1)]  
  ext.lines <- all.lines[(idx[2]+1):(idx[3]-1)]
  int.lines <- all.lines[(idx[3]+1):length(all.lines)]

  ##Get fluxnames and species names and carbion transitions from equations
  all.fluxnames <- str_match(eq.lines, "(.+):")[,2]
  all.reactant.terms <- str_extract_all(unlist(str_extract_all(eq.lines, ":.+\\w+.+->")), "\\w+")
  all.product.terms <- str_extract_all(unlist(str_extract_all(eq.lines, "->.+\\w+")), "\\w+")
  
  all.reactants <- lapply(all.reactant.terms, function(x)x[seq(1, length(x), 2)])
  all.reactant.positions <- lapply(all.reactant.terms, function(x)x[seq(2, length(x), 2)])
  
  all.products <- lapply(all.product.terms, function(x)x[seq(1, length(x), 2)])
  all.product.positions <- lapply(all.product.terms, function(x)x[seq(2, length(x), 2)])
  
  ##function to create species references for reaction object from vectors of species names with taking stoichiometries into account
  MakeSpeciesRefs <- function(species){
    lapply(unique(species), function(s){
      new("SpeciesReference", id=s, species=s, stoichiometry=sum(species==s))
    })
  }
  
  ##create all reaction objects for the model
  reactions <- mapply(function(flux, reac, prod, reac.pos, prod.pos){
    new("FluxEsReaction", id=flux, flux.name=flux, reactants=reac, products=prod, reactants.atom.positions=reac.pos, products.atom.positions=prod.pos, reversible=FALSE)
  }, all.fluxnames, lapply(all.reactants, MakeSpeciesRefs), lapply(all.products, MakeSpeciesRefs), all.reactant.positions, all.product.positions)
  
  ##set reversible attribute to reversible reactions  
  ##reversible.idx <- grep("<->", eq.lines)
  ##for (ri in reversible.idx)
  ##  reactions[[ri]]@reversible=TRUE
  
  ##get species names and concentrations  
  ext.speciesconc <- as.numeric(str_match(ext.lines, "(\\w+)\\s+(.+)")[,3])
  names(ext.speciesconc) <- str_match(ext.lines, "(\\w+)\\s+(.+)")[,2]  
  int.speciesconc <- as.numeric(str_match(int.lines, "(\\w+)\\s+(.+)")[,3])
  names(int.speciesconc) <- str_match(int.lines, "(\\w+)\\s+(.+)")[,2]
  all.speciesconc <- c(int.speciesconc, ext.speciesconc)
  
  ##create all species objects for the model
  all.species.names <- unlist(c(all.reactants, all.products))
  ##species names should not begin with a number or a special character
  if (! all(sapply(sapply(all.species.names, substring, 1, 1), function(x) grepl("[[:alpha:]]", x))))
    stop("Metabolite(s) ", unique(names(which(!sapply(sapply(all.species.names, substring, 1, 1), function(x) grepl("[[:alpha:]]", x))))), " do not start with letter!\n")
  
  ##test if species names given in internals and externals agree with the species names from the equations  
  if (! all(sort(unique(all.species.names))==sort(unique(c(names(ext.speciesconc), names(int.speciesconc)))))){
    stop("Entries in #internal metabolites and #external metabolites do not match metabolites in #equations!\n")
  }
  species.pos <- unlist(c(all.reactant.positions, all.product.positions))
  
  species <- lapply(unique(all.species.names), function(s){
    new("FluxEsSpecies", id=s, name=s, carbon.count=nchar(species.pos[which(all.species.names==s)[1]]), external=s%in%names(ext.speciesconc), initialConcentration=all.speciesconc[s])
  })  
  names(species) <- unique(all.species.names)

  ##create all flux parameters for the model
  parameters <- unique(lapply(reactions, function(r){
    is.external.influx <- any(sapply(c(r@reactants), function(x)species[[x@id]]@external))    
    is.external.outflux <- any(sapply(c(r@products), function(x)species[[x@id]]@external))    
    new("FluxEsParameter", id=r@flux.name, name=r@flux.name, is.flux=TRUE, external.influx=is.external.influx, external.outflux=is.external.outflux, value=1)
  }))
  ##create model object
  return(new("FluxEsModel", species=species, reactions=reactions, parameters=parameters))
}


CreateCmodel <- function(model, odes, file.name=paste(model@id, ".c", sep="")){
  cat("#include <R.h>\n\n", file=file.name, append=FALSE)
  ## Get the parameter names and values
  ## Caution: some parameters do not have values (eg size of virtual pools...) and are not included
  parameterNames <- unlist(sapply(model@parameters, function(x)if(length(x@value)>0)x@id))
  parameterNames <- c(unlist(as.vector(sapply(model@species, function(s)s@id))), parameterNames)        
  ## Print the parameter vector        
  cat("static double p[", deparse(length(parameterNames), control = c() ), "];\n\n", file=file.name, append=TRUE)
  ## Print the parameters
  cat("// indices for parameters\n", file=file.name, append=TRUE)
  for (i in seq_along(parameterNames))
    cat("#define", parameterNames[i], i-1, "\n", file=file.name, append=TRUE)
  ## Define names for indexes for isotopomer species
  cat("// indices for isotopomer species\n", file = file.name, append = TRUE)        
  isos <- names(odes)  
  for(i in seq_along(isos))
    cat("#define", isos[i], i-1, "\n", file=file.name, append=TRUE)  
  ## Print the initializer function which has to have the same name as the dll without the extension
  cat("\nvoid ", model@id, "(void (* odeparms)(int *, double *)){\n", sep = "", file=file.name, append=TRUE)
  cat("\tint N =", length(parameterNames), ";\n", file=file.name, append=TRUE)
  cat("\todeparms(&N, p);\n}\n\n", file=file.name, append=TRUE)
  ## Print the function consisting the derivatives
  cat("void derivsc(int *neq, double *t, double *y, double *ydot, double *yout, int*ip){\n", file=file.name, append=TRUE)
  ## Substitute all names of isotopomers with y[isotopomer]
  odes <- lapply(odes, function(x)gsub("(\\w+_\\w+)", paste("y[\\1]"), x))
  ## For the left-hand side of the ode, it has to be ydot[iso] and not y[iso]
  ## We assume here that one ode only has one "=" character
  odes <- lapply(odes, function(x)gsub("y\\[(\\w+)\\]\\s+=", "ydot[\\1] =",x))
  ## Substitute all names for parameters with p[parameter]
  odes <- lapply(odes, function(x)c(paste(x, ";\n")))
  for (p.name in parameterNames){
    odes <- lapply(odes, function(ode){
      gsub(paste("\\s", p.name, "\\s", sep=""), paste("p[", p.name, "]", sep=""), ode)
    })
  }
  ## write odes to file
  counter <- 0
  for (ode in odes){
    counter <- counter + 1
    cat("\t", ode, file=file.name, append=TRUE)  
  }
  cat("}\n", file=file.name, append=TRUE)        
  system(paste("R CMD SHLIB -c --preclean ", file.name, sep = ""))
  dyn.load(paste(model@id, ".so", sep=""))    
}
