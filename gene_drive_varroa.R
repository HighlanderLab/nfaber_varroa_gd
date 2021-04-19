##################
## ---- Setup ----
##################

rm(list = ls())
library('tidyverse')
library('AlphaSimR')
library('profvis')
library('reshape2')
library('viridis')
library('compiler')
library('iterators')
library('snow')
library('doSNOW')
library('foreach')
library('data.table')

###########################################
## ---- Gene drive simulation function ----
###########################################

gd <- function(input){
  
  # Make sure R doesn't screw around with scientific notation of id's
  originalOptions <- options(scipen = 999)
  on.exit(options(originalOptions))
  
  ##########################################
  ## ---- Assign all variables in input ----
  ##########################################
  
  for (i in 1:ncol(input)) {assign(names(input)[i], input[,i])}
  
  # Initialise founder population
  founderPop <- quickHaplo(nInd=(nInit+gd)*2, nChr=1, segSites=loci, genLen=0)
  
  # Set model settings
  SP <- SimParam$new(founderPop)
  SP$addTraitA(nQtlPerChr=loci)
  SP$setGender("yes_sys")
  
  ######################
  ## ---- Functions ----
  ######################
  
  # Altered AlphaSimR functions
  source(file = 'AlphaSimR_alt.R')
  
  initiateBroodcells <- function(years = 1, gaplength, gapday, broodbreak){
    broodcells <- read.delim("broodcells.csv", 
                             header=TRUE, sep=",")
    broodcells <- select(broodcells,day,newWorker,newDrone,adults)
    colnames(broodcells) <- c("day","workers","drones","adults")
    
    if (gaplength > 0){
      broodcells$workers[gapday:(gapday+gaplength-1)] <- round(broodcells$workers[gapday:(gapday+gaplength-1)]*broodbreak)
      broodcells$drones[gapday:(gapday+gaplength-1)] <- round(broodcells$drones[gapday:(gapday+gaplength-1)]*broodbreak)
    }
    
    broodcells <- broodcells[rep(1:365, years), ]
    broodcells$day <- 1:(365*years)
    return(broodcells)
  }
  
  invasionRate <- function(ratio){
    invades <- (1 + exp(-(-2.87 + 0.00385*(ratio*10000))))^-1
    return(invades)
  }
  
  invasionRateAlt <- function(brood){
    invades <- 1-exp(-(brood))
    return(invades)
  }
  
  enterBroodcells <- function(pop, nDrones, nWorkers){
    dronePops <- c()
    workerPops <- c()
    
    assignedCells <- sample(1:(nDrones+nWorkers), pop@nInd, replace=TRUE, 
                            prob = c(rep(8/9,nDrones),rep(1/9,nWorkers)))
    
    for (cell in unique(assignedCells)){
      subpop <- pop[assignedCells == cell]
      if (cell <= nDrones){
        dronePops <- c(dronePops,subpop)
      } else {
        workerPops <- c(workerPops,subpop)
      }
    }
    
    return(list(dronePops,workerPops))
  }
  
  assignMates <- function(pop){
    if ("M" %in% pop@gender & "F" %in% pop@gender){
      crossPlan <- randCrossDistr(pop, nCrosses = sum(pop@gender=="F"), simParam = SP)
      
      # log mate in second ebv column
      pop@ebv[crossPlan[,1],2] <- as.numeric(pop@id[crossPlan[,2]])
    }
    return(pop)
  }
  
  reproduce <- function(pop, brood){
    
    if ("F" %in% pop@gender & "M" %in% pop@gender){ 
      
      # Homing
      if (drive) {
        pop <- homing(pop = pop)
      }
      
      # Each female generates one haploid male offspring from an unfertilized egg
      males <- makeDH(pop[pop@gender=="F"], nDH = 1, simParam = SP)
      males@gender <- rep("M",males@nInd)
      males@mother <- pop@id[pop@gender=="F"]
      males <- prepEBVslots(males)
      
      # Determine the amount of female offspring distribution
      if (brood == "drone"){
        nProgeny <- c(0,1,2,3,4,5) # drone brood
        probs <- c(0.39,0.08,0.14,0.22,0.15,0.02) 
      } else if (brood == "worker") {
        nProgeny <- c(0,1,2,3) # worker brood
        probs <- c(0.46,0.35,0.17,0.02)
      }
      nProgeny = sample(nProgeny, sum(pop@gender=="F"), prob = probs, replace = TRUE)
      
      if (sum(nProgeny) == 0){
        return(males)
      }
      
      crossPlan <- matrix(c(pop@id[pop@gender=="F"],pop@ebv[pop@gender=="F",2]),ncol=2)
      
      females <- makeCrossDistr(pop, crossPlan = crossPlan, nProgeny = nProgeny, simParam = SP)
      females@gender <- rep("F",females@nInd)
      females <- prepEBVslots(females)
      females <- assignReproductionCycles(females)
      
      # Observe maximum amount of female offspring per cell
      if (brood == "drone"){
        if (females@nInd > (16 - males@nInd)){
          females <- mortality(females, 1 - (16 - males@nInd) / females@nInd)
        }
      } else if (brood == "worker") {
        if (females@nInd > (8 - males@nInd)){
          females <- mortality(females, 1 - (8 - males@nInd) / females@nInd)
        }
      }
      
      offspring <- mergePops(list(males,females))
      
      return(offspring)
    }
  }
  
  sortMales <- function(broodPops, males){
    broodPops <- lapply(broodPops, FUN=function(broodPop) mergePops(list(broodPop, males[males@id %in% broodPop@ebv[,2]])))
    return(broodPops)
  }
  
  mortality <- function(pop, mortality){
    survive <- (pop@ebv[,1] > 0) & (runif(pop@nInd) > mortality)
    pop <- pop[survive]
    return(pop)
  }
  
  meanIBS <- function(pop){
    genotypes <- pullQtlGeno(pop, simParam = SP)
    homozygosity <- apply(genotypes,1,FUN=function(locus) sum(locus==0 | locus == 2)/pop@nLoci)
    return(mean(homozygosity))
  }
  
  prepEBVslots <- function(pop) {
    pop@ebv <- matrix(rep(NA,3*pop@nInd), ncol = 3) # prep ebv slots
    # slot 1: reproductive cycles
    # slot 2: mate id
    # slot 3: days remaining in brood
    return(pop)
  }
  
  assignReproductionCycles <- function(pop){
    pop@ebv[pop@gender=="F",1] <- sample(1:4, sum(pop@gender=="F"), replace = TRUE)
    return(pop)
  }
  
  substractReproductionCycle <- function(pop) {
    pop@ebv[,1] = pop@ebv[,1]-1
    return(pop)
  }
  
  setBroodStatus <- function(pop, brood) {
    if (brood=="worker") {
      pop@ebv[pop@gender=="F",3] = rep(12,sum(pop@gender=="F"))
      #pop@ebv[pop@gender=="F",3] = round(rnorm(sum(pop@gender=="F"), mean = 20, sd = 1))
    } else if (brood=="drone") {
      pop@ebv[pop@gender=="F",3] = rep(14,sum(pop@gender=="F"))
      #pop@ebv[pop@gender=="F",3] = round(rnorm(sum(pop@gender=="F"), mean = 22, sd = 1))
    } else if (brood=="initial") {
      pop@ebv[pop@gender=="F",3] = rep(0,sum(pop@gender=="F"))
      #pop@ebv[pop@gender=="F",3] = sample(0,sum(pop@gender=="F"), replace = TRUE)
    }
    return(pop)
  }
  
  substractBroodStatus <- function(pop) {
    pop@ebv[,3] = pop@ebv[,3]-1
    return(pop)
  }
  
  createHaploDF <- function(pop, simParam = SP) {
    haplotype <- setDT(bind_cols(name = row.names(pullQtlHaplo(pop, simParam = simParam)),
                                 pullQtlHaplo(pop, simParam = simParam)))
    haplotype[, `:=`(id = strsplit(name, split = "_")[[1]][1],
                     haplo = strsplit(name, split = "_")[[1]][2]),
              by = name] [, name := NULL]

    haplotype[, locus1 := fifelse(test = QTL_1 == 1,
                                  yes = fifelse(test = QTL_2 == 1, yes = "GD", no = "RE"),
                                  no = fifelse(test = QTL_2 == 1, yes = "NF", no = "WT"))]
    return(haplotype[ , .(id, haplo, locus1)])
  }
  
  initPop <- function(founderPop, nInit, gd, simParam = SP) {
    
    # Initialise mated female population
    pop <- newPop(founderPop, simParam = simParam)
    pop <- prepEBVslots(pop)
    
    broodPops <- list()
    for (pair in seq(1,pop@nInd,2)){
      broodPops <- c(broodPops, pop[pair:(pair+1)])
    }
    broodPops <- lapply(broodPops, assignMates)
    pop <- mergePops(broodPops); rm(broodPops)
    
    pop <- assignReproductionCycles(pop)
    pop <- setBroodStatus(pop,"initial")
    
    if (gd > 0) { # Initialise gene drive
      pop <- editGenomeFix(pop, ind = 1:pop@nInd, chr = rep(1, pop@nLoci), segSites = 1:pop@nLoci, allele = 0, simParam = simParam)
      popGD <- pop[(nInit*2+1):pop@nInd]
      pop <- pop[1:(nInit*2)]
  
      popGDmales <- popGD[popGD@gender=="M"]
      popGDfemales <- popGD[popGD@gender=="F"]
  
      popGDmales <- editGenomeFix(popGDmales, ind = 1:popGDmales@nInd, chr = rep(1, popGDmales@nLoci), segSites = 1:popGDmales@nLoci, allele = 1, simParam = simParam)
      popGDfemales <- editHaplo(popGDfemales, ind = 1:popGDfemales@nInd, chr = rep(1, popGDfemales@nLoci), segSites = 1:popGDfemales@nLoci, allele = 1, haplotype = rep(1, popGDfemales@nLoci), simParam = simParam)
  
      pop <- mergePops(list(pop, popGDmales, popGDfemales))
    }
    
    return(pop)
  }
  
  homing <- function(pop, p_nhej = 0.02, simParam = SP) {
    haplotype <- as_tibble(createHaploDF(pop))
    
    for (row in seq(1, nrow(haplotype), 2)) {
      for (column in 3:ncol(haplotype)){
        tmpHaplo <- haplotype[row:(row+1), c(1:2, column)]
        if (!any(tmpHaplo[,3] == "GD")) {next}
        qtl <- which(3:ncol(haplotype) == column)*2
        if (any(tmpHaplo[, 3] == "GD") & any(tmpHaplo[, 3] == "WT")){
          if (runif(1) < 0.95) { #does it cut?
            if (runif(1) < 1-p_nhej) { #HDR if smaller than 0.98
              #homing
              pop <- editHaplo(pop = pop, 
                               ind = which(seq(1, nrow(haplotype), 2) == row), 
                               chr = rep(1, 2),
                               segSites = (qtl-1):qtl,
                               allele = 1, #1 because it's the gene drive
                               haplotype = which(tmpHaplo[,3] != "GD"),
                               simParam = simParam)
            } else if (runif(1) < 2/3) { #NHEJ, first check for frame shift
              #change haplotype to NF (01)
              pop <- editHaplo(pop = pop, 
                               ind = which(seq(1, nrow(haplotype), 2) == row), 
                               chr = rep(1, 2),
                               segSites = (qtl-1):qtl,
                               allele = c(0, 1), #01 for NF
                               haplotype = which(tmpHaplo[,3] != "GD"), 
                               simParam = simParam)
            }  else {
              #change haplotype to RE (10) 
              pop <- editHaplo(pop = pop, 
                               ind = which(seq(1, nrow(haplotype), 2) == row), 
                               chr = rep(1, 2),
                               segSites = (qtl-1):qtl,
                               allele = c(1, 0), #10 for NF
                               haplotype = which(tmpHaplo[,3] != "GD"), 
                               simParam = simParam)
            }
          }
        }
      }
    }
    return(pop)
  }
  
  reproductivity <- function (pop, strategy = 1){
    if (strategy == 1){ # Neutral gene drive
      return(pop)
    }
    pop_m <- pop[pop@gender == "M"]
    pop_f <- pop[pop@gender == "F"]
    if (strategy == 2){ # Male infertility
      haploDF <- createHaploDF(pop_m)
      haplo1 <- filter(haploDF, haplo == 1 & (locus1 == "GD" | locus1 == "NF")) %>%
        select(id)
      haplo2 <- filter(haploDF, haplo == 2 & (locus1 == "GD" | locus1 == "NF")) %>%
        select(id)
      if (length(haplo1$id) > 0 & length(haplo2$id) > 0){
        pop_m <- pop_m[pop_m@id != intersect(haplo1, haplo2)$id]
        pop <- mergePops(list(pop_m, pop_f))
      }
      # filter out of pop
    } else if (strategy == 3) { # Female infertility
      #horrible shit happens to the foundresses
      haploDF <- createHaploDF(pop_f)
      haplo1 <- filter(haploDF, haplo == 1 & (locus1 == "GD" | locus1 == "NF")) %>%
        select(id)
      haplo2 <- filter(haploDF, haplo == 2 & (locus1 == "GD" | locus1 == "NF")) %>%
        select(id)
      if (length(haplo1$id) > 0 & length(haplo2$id) > 0){
        pop_f <- pop_f[pop_f@id != intersect(haplo1, haplo2)$id]
        pop <- mergePops(list(pop_m, pop_f))
      }
    }
    return(pop)
  }
  
  #######################
  ## ---- Simulation ----
  #######################
  
  # Initiate broodcells backbone data
  broodcells <- initiateBroodcells(years=years, gaplength=gaplength, gapday=gapday, broodbreak=broodbreak)
  if (testing){
    broodcells <- broodcells[-(1:108),]
  }
  
  # Track population information
  var <- c('nInd','IBS','droneRep','workerRep','potentialInf',
           'nPhoretic','pPhoretic','beesPhoretic','phoreticTime',
           'varroaPerDrone','varroaPerWorker',
           'offspringPerDrone','offspringPerWorker', 
           'WT', 'GD', 'NF', 'RE', 'homoWT', 'homoGD', 'heteroGD')
  results <- array(NA,
                   dim=c(length(broodcells$day),length(var)),
                   dimnames=list(days=broodcells$day,var=var))

  pop <- initPop(founderPop, nInit, gd, SP)
  males <- pop[pop@gender=="M"]
  pop <- pop[pop@gender=="F"]
  
  for(day in 1:nrow(broodcells)){
    droneRep <- 0; workerRep <- 0; infestations <- 0
    varroaPerDrone <- 0; varroaPerWorker <- 0; 
    offspringPerDrone <- 0; offspringPerWorker <- 0
    
    nDrones <- broodcells$drones[day]
    nWorkers <- broodcells$workers[day]
    nAdults <- broodcells$adults[day]
    
    pop <- substractBroodStatus(pop)
    
    if (nDrones + nWorkers > 0){
      
      # Select females who are available for infestation
      reproductive <- pop@ebv[,3] <= 0
      
      infestations <- sum(runif(sum(reproductive)) < invasionRate((nWorkers + nDrones)/nAdults))
      #infestations <- sum(runif(reproductive) < (nWorkers + nDrones)/nAdults)
      
      if (infestations > 0 & sum(reproductive) > 0) {
        
        if (sum(reproductive) > infestations) {
          enteringBrood <- sample(which(reproductive), infestations)
          reproductive[-enteringBrood] <- FALSE
        }
        
        nonReproductive <- pop[!reproductive]
        
        # Reproductive females pick drone or worker brood
        broodPops <- enterBroodcells(pop[reproductive], nDrones, nWorkers); rm(pop)
        broodTypes <- c(rep("drone", length(broodPops[[1]])),rep("worker", length(broodPops[[2]])))
        
        droneRep <- length(broodPops[[1]])
        workerRep <- length(broodPops[[2]])
        varroaPerDrone <- mean(sapply(broodPops[[1]], FUN = function(x) x@nInd))
        varroaPerWorker <- mean(sapply(broodPops[[2]], FUN = function(x) x@nInd))
        
        broodPops <- c(broodPops[[1]], broodPops[[2]])
        
        # Sort the males into females' cells
        broodPops <- sortMales(broodPops, males)
        
        # Produce offspring
        offspring <- mapply(broodPops, FUN = reproduce, brood = broodTypes)
        broodPops <- lapply(broodPops, substractReproductionCycle)
        
        offspringPerDrone <- mean(sapply(offspring[broodTypes=="drone"], FUN = function(x) x[x@gender=="F"]@nInd))
        offspringPerWorker <- mean(sapply(offspring[broodTypes=="worker"], FUN = function(x) x[x@gender=="F"]@nInd))
        
        # Determine how long each Varroa will remain in cell + phoretic stage
        offspring <- mapply(offspring, FUN = setBroodStatus, brood = broodTypes)
        broodPops <- mapply(broodPops, FUN = setBroodStatus, brood = broodTypes)
        
        # Random mating within each cell
        offspring <- lapply(offspring, assignMates)
        
        # Emergence from brood cells
        pop <- mergePops(list(nonReproductive,
                              mergePops(broodPops),
                              mergePops(offspring))); rm(broodPops); rm(offspring)
        
        # Females emerge from brood, store males until female produces offspring
        males <- mergePops(list(males,pop[pop@gender=="M"]))
        pop <- pop[pop@gender=="F"]
      }
    } 
    
    # Some varroas die
    pop <- mortality(pop, mortality = mort)
    
    # Track population values
    nInd          <- pop@nInd
    IBS           <- meanIBS(pop)
    
    droneRep      <- droneRep
    workerRep     <- workerRep
    potentialInf  <- infestations
    varroaPerDrone <- varroaPerDrone
    varroaPerWorker <- varroaPerWorker
    
    nPhoretic     <- sum(pop@ebv[,3]<=0)
    pPhoretic     <- nPhoretic / nInd
    beesPhoretic  <- nPhoretic / nAdults
    phoreticTime  <- abs(mean(pop@ebv[pop@ebv[,3]<=0, 3]))
    
    haploDF       <- createHaploDF(pop)
    countHaplo    <- haploDF[, .N, by = locus1]
    haploDFWider  <- setDT(pivot_wider(haploDF, id_cols = id, names_from = haplo,
                                       names_prefix = "haplo", values_from = locus1))
    WT            <- ifelse(length(countHaplo[locus1 == "WT", N]/(nInd * 2)) > 0,
                            yes = countHaplo[locus1 == "WT", N]/(nInd * 2), no = 0)
    GD            <- ifelse(length(countHaplo[locus1 == "GD", N]/(nInd * 2)) > 0,
                            yes = countHaplo[locus1 == "GD", N]/(nInd * 2), no = 0)
    NF            <- ifelse(length(countHaplo[locus1 == "NF", N]/(nInd * 2)) > 0,
                            yes = countHaplo[locus1 == "NF", N]/(nInd * 2), no = 0)
    RE            <- ifelse(length(countHaplo[locus1 == "RE", N]/(nInd * 2)) > 0,
                            yes = countHaplo[locus1 == "RE", N]/(nInd * 2), no = 0)
    homoWT        <- haploDFWider[haplo1 == "WT" & haplo2=="WT", .N]
    homoGD        <- haploDFWider[haplo1 == "GD" & haplo2=="GD", .N]
    heteroGD      <- haploDFWider[xor(haplo1 == "GD", haplo2 == "GD"), .N]
    
    for (v in var) {
      results[day,v] <- get(x=v)
    }
    
    if(sum(pop@ebv[,1]) == 0 | pop@nInd > 20000){
      break
    }
  } # end of generation
  
  ###########################
  ## ---- Process output ----
  ###########################

  ## add all the parameters to the results data frame
  results <- cbind(day=broodcells$day, input, results)
  
  return(results)
}

############################################
## ---- Compile the simulation function ----
############################################

gd.comp <- cmpfun(gd)
