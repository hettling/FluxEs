require('deSolve')
require('plyr')

source('FluxEsFunctions.R')
source('FluxEsClasses.R')
source('FluxEsIO.R')


file.name <- "examplemodel.txt"

model <- ParseModel(file.name)
model@id <- "mymodel"
##set labeled input metabolite
model@species[["acetate"]]@input <- TRUE

labeled.input <- list(list("acetate", c(1), 1))

##Simulation of only mass isotopomers
ms.odes <- calcODEs(model, type="MS")
CreateCmodel(model, ms.odes)


ms.startvals <- unlist(CalcMsStartValues(model, labeled.input))
names(ms.startvals) <- unlist(sapply(model@species, function(x)paste(x@id, 0:x@carbon.count, sep="_")))

pool.sizes <- sapply(model@species, function(s)s@initialConcentration)

timepoints <- c(0, 1, 2, 6, 10, 29)

##vector of fluxes
v <- c(1,1,1,1,1,1)

output <- lsodes(ms.startvals, seq(0, 29, 0.1), "derivsc", parms=c(pool.sizes, v), dllname=model@id, nout=0)

##output <- ms.output[which(ms.output[,"time"] %in% timepoints),]



