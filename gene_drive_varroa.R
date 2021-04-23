##################
## ---- Setup ----
##################

rm(list = ls())
library('AlphaSimR')
library('tidyverse')
library('profvis')
library('reshape2')
library('viridis')
library('compiler')
library('iterators')
library('snow')
library('doSNOW')
library('foreach')
library('data.table')
library('patchwork')

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
  
  source(file = 'AlphaSimR_alt.R') # Altered AlphaSimR functions
  
  # Initialise founder population
  founderPop <- quickHaploInbr(nInd=nInit+gd, nChr=1, segSites=loci, inbred=inbred, inbrCoef=inbrCoef)
  suppPop <- quickHaploInbr(nInd=supplement, nChr=1, segSites=loci, inbred=inbred, inbrCoef=inbrCoef)
  if (gd > 0) {
    founderPop@genMap[[1]][1:2] <- c(0,0)
    suppPop@genMap[[1]][1:2] <- c(0,0)
  }
  
  # Set model settings
  SP <- SimParam$new(founderPop)
  SP$addTraitA(nQtlPerChr=loci)
  SP$setSexes("yes_sys")
  
  ######################
  ## ---- Functions ----
  ######################
  
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
    invades <- (1 + exp(-(invasionIntercept + invasionSlope*(ratio*10000))))^-1
    return(invades)
  }
  
  enterBroodcells <- function(pop, nDrones, nWorkers){
    dronePops <- c()
    workerPops <- c()
    
    assignedCells <- sample(1:(nDrones+nWorkers), pop@nInd, replace=TRUE, 
                            prob = c(rep(DCP/(DCP+1),nDrones),rep(1/(DCP+1),nWorkers)))
    
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
    if ("M" %in% pop@sex & "F" %in% pop@sex){
      crossPlan <- randCrossDistr(pop, nCrosses = sum(pop@sex=="F"), simParam = SP)
      
      # log mate in second ebv column
      pop@ebv[crossPlan[,1],2] <- as.numeric(pop@id[crossPlan[,2]])
    }
    return(pop)
  }
  
  reproduce <- function(pop, brood, strategy){
    
    # Homing
    if (drive) {
      pop <- homing(pop = pop)
    }
    
    # Each female generates one haploid male offspring from an unfertilized egg
    infertiles <- reproductivity(pop, strategy)
    mothers <- pop[pop@sex=="F" 
                   & !(pop@id %in% infertiles) 
                   & !(pop@ebv[,2] %in% infertiles)]
    
    if (mothers@nInd > 0){
    
      males <- makeDH(mothers, nDH = 1, simParam = SP)
      males@sex <- rep("M",males@nInd)
      males@mother <- mothers@id
      males <- prepEBVslots(males)
      
      # Determine the amount of female offspring distribution
      if (brood == "drone"){
        nProgeny <- c(0,1,2,3,4,5) # drone brood
        probs <- c(0.39,0.08,0.14,0.22,0.15,0.02) 
      } else if (brood == "worker") {
        nProgeny <- c(0,1,2,3) # worker brood
        probs <- c(0.46,0.35,0.17,0.02)
      }
      nProgeny = sample(nProgeny, mothers@nInd, prob = probs, replace = TRUE)
      
      crossPlan <- matrix(c(mothers@id,mothers@ebv[,2]),ncol=2)
      
      if (sum(nProgeny) > 0){
        
        females <- makeCrossDistr(pop, crossPlan = crossPlan, nProgeny = nProgeny, simParam = SP)
        females@sex <- rep("F",females@nInd)
        females <- prepEBVslots(females)
        females <- assignReproductionCycles(females)
        
        # Observe maximum amount of female offspring per cell
        if (brood == "drone"){
          if (females@nInd > (maxDrone - males@nInd)){
            females <- mortality(females, 1 - (maxDrone - males@nInd) / females@nInd)
          }
        } else if (brood == "worker") {
          if (females@nInd > (maxWorker - males@nInd)){
            females <- mortality(females, 1 - (maxWorker - males@nInd) / females@nInd)
          }
        }
        
        females <- setBroodStatus(pop = females, brood = brood)
        offspring <- mergePops(list(males,females))
        
        return(offspring)
      }
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
    pop@ebv <- matrix(rep(NA,3*pop@nInd), ncol = 3)
    # slot 1: reproductive cycles
    # slot 2: mate id
    # slot 3: days remaining in brood
    return(pop)
  }
  
  assignReproductionCycles <- function(pop){
    pop@ebv[pop@sex=="F",1] <- sample(1:4, sum(pop@sex=="F"), replace = TRUE)
    return(pop)
  }
  
  substractReproductionCycle <- function(pop) {
    pop@ebv[,1] = pop@ebv[,1]-1
    return(pop)
  }
  
  setBroodStatus <- function(pop, brood) {
    if (brood=="worker") {
      pop@ebv[pop@sex=="F",3] = rep(12,sum(pop@sex=="F"))
    } else if (brood=="drone") {
      pop@ebv[pop@sex=="F",3] = rep(14,sum(pop@sex=="F"))
    } else if (brood=="initial") {
      pop@ebv[pop@sex=="F",3] = rep(0,sum(pop@sex=="F"))
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
    females <- newPop(founderPop, simParam = simParam)
    females@sex <- rep("F",females@nInd)
    females <- prepEBVslots(females)
    
    males <- makeDH(females, nDH = 1, simParam = simParam)
    males@sex <- rep("M",males@nInd)
    males <- prepEBVslots(males)
    
    #assign mates
    females@ebv[,2] <- as.numeric(males@id)
    
    if (gd > 0) { # Initialise gene drive
      
      femalesGD <- females[1:gd]
      malesGD <- males[1:gd]
      
      if (gdZygosity == "homozygous") {
        femalesGD <- editGenomeFix(femalesGD, ind = 1:femalesGD@nInd, chr = rep(1, 2), segSites = 1:2, allele = 1, simParam = simParam)
        malesGD <- editGenomeFix(malesGD, ind = 1:malesGD@nInd, chr = rep(1, 2), segSites = 1:2, allele = 1, simParam = simParam)
      } else if (gdZygosity == "heterozygous") {
        femalesGD <- editHaplo(femalesGD, ind = 1:femalesGD@nInd, chr = rep(1, 2), segSites = 1:2, allele = 1, haplotype = rep(1, 2), simParam = simParam)
        femalesGD <- editHaplo(femalesGD, ind = 1:femalesGD@nInd, chr = rep(1, 2), segSites = 1:2, allele = 0, haplotype = rep(2, 2), simParam = simParam)
        malesGD <- editGenomeFix(malesGD, ind = 1:malesGD@nInd, chr = rep(1, 2), segSites = 1:2, allele = 1, simParam = simParam)
      }
      
      if (nInit > 0) {
        females <- females[(gd+1):females@nInd]
        males <- males[(gd+1):males@nInd]
        
        females <- editGenomeFix(females, ind = 1:females@nInd, chr = rep(1, 2), segSites = 1:2, allele = 0, simParam = simParam)
        males <- editGenomeFix(males, ind = 1:males@nInd, chr = rep(1, 2), segSites = 1:2, allele = 0, simParam = simParam)
        
        pop <- mergePops(list(femalesGD, females, malesGD, males))
      } else {
        pop <- mergePops(list(femalesGD, malesGD))
      }
      
    } else {
      pop <- mergePops(list(females, males))
    }
    
    pop <- assignReproductionCycles(pop)
    pop <- setBroodStatus(pop,"initial")
    
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
            if (runif(1) < 1-p_nhej) { #HDR if smaller than p_nhej
              pop <- editHaplo(pop = pop, 
                               ind = which(seq(1, nrow(haplotype), 2) == row), 
                               chr = rep(1, 2),
                               segSites = (qtl-1):qtl,
                               allele = 1, 
                               haplotype = which(tmpHaplo[,3] != "GD"),
                               simParam = simParam)
            } else if (runif(1) < 2/3) { #NHEJ, first check for frame shift
              pop <- editHaplo(pop = pop, 
                               ind = which(seq(1, nrow(haplotype), 2) == row), 
                               chr = rep(1, 2),
                               segSites = (qtl-1):qtl,
                               allele = c(0, 1), #01 for NF
                               haplotype = which(tmpHaplo[,3] != "GD"), 
                               simParam = simParam)
            }  else {
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
  
  reproductivity <- function (pop, strategy){
    if (strategy == 1){ # Neutral gene drive
      return(NULL)
    }
    pop_m <- pop[pop@sex == "M"]
    pop_f <- pop[pop@sex == "F"]
    if (strategy == 2){ # Male infertility
      haploDF <- createHaploDF(pop_m)
      haplo1 <- filter(haploDF, haplo == 1 & (locus1 == "GD" | locus1 == "NF")) %>%
        select(id)
      haplo2 <- filter(haploDF, haplo == 2 & (locus1 == "GD" | locus1 == "NF")) %>%
        select(id)
      infertiles <- intersect(haplo1, haplo2)$id
      
    } else if (strategy == 3) { # Female infertility
      haploDF <- createHaploDF(pop_f)
      haplo1 <- filter(haploDF, haplo == 1 & (locus1 == "GD" | locus1 == "NF")) %>%
        select(id)
      haplo2 <- filter(haploDF, haplo == 2 & (locus1 == "GD" | locus1 == "NF")) %>%
        select(id)
      infertiles <- intersect(haplo1, haplo2)$id
    }
    return(infertiles)
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
  males <- pop[pop@sex=="M"]
  pop <- pop[pop@sex=="F"]
  
  for(day in 1:nrow(broodcells)){
    droneRep <- 0; workerRep <- 0; infestations <- 0
    varroaPerDrone <- 0; varroaPerWorker <- 0; 
    offspringPerDrone <- 0; offspringPerWorker <- 0
    
    nDrones <- broodcells$drones[day]
    nWorkers <- broodcells$workers[day]
    nAdults <- broodcells$adults[day]
    
    pop <- substractBroodStatus(pop)
    
    if (nDrones + nWorkers > 0){
      
      reproductive <- pop@ebv[,3] <= 0 # Select females who are in their dispersal phase
      
      infestations <- sum(runif(sum(reproductive)) < invasionRate((nWorkers + nDrones)/nAdults))
      
      if (infestations > 0 & sum(reproductive) > 0) {
        
        if (sum(reproductive) > infestations) {
          enteringBrood <- sample(which(reproductive), infestations)
          reproductive[-enteringBrood] <- FALSE; rm(enteringBrood)
        }
        
        nonReproductive <- pop[!reproductive]
        
        # Reproductive females pick drone or worker brood
        broodPops <- enterBroodcells(pop[reproductive], nDrones, nWorkers); rm(pop); rm(reproductive)
        broodTypes <- c(rep("drone", length(broodPops[[1]])),rep("worker", length(broodPops[[2]])))
        
        # save information about reproduction in the brood cells
        droneRep <- length(broodPops[[1]])
        workerRep <- length(broodPops[[2]])
        varroaPerDrone <- mean(sapply(broodPops[[1]], FUN = function(x) x@nInd))
        varroaPerWorker <- mean(sapply(broodPops[[2]], FUN = function(x) x@nInd))
        
        broodPops <- c(broodPops[[1]], broodPops[[2]])
        
        # Sort the males into females' cells
        broodPops <- sortMales(broodPops, males)
        
        # Produce offspring
        offspring <- mapply(broodPops, FUN = reproduce, brood = broodTypes, strategy = strategy)
        
        broodPops <- lapply(broodPops, substractReproductionCycle)
        broodPops <- mapply(broodPops, FUN = setBroodStatus, brood = broodTypes)
        broodPops <- mergePops(broodPops)
        
        if (length(compact(offspring)) > 0){
          
          offspringPerDrone <- mean(sapply(compact(offspring[broodTypes=="drone"]), FUN = function(x) x[x@sex=="F"]@nInd))
          offspringPerWorker <- mean(sapply(compact(offspring[broodTypes=="worker"]), FUN = function(x) x[x@sex=="F"]@nInd))
          offspring <- compact(offspring)
          
          # Random mating within each cell
          offspring <- lapply(offspring, assignMates)
          
          # Emergence from brood cells
          offspring <- mergePops(offspring)
          males <- mergePops(list(males,offspring[offspring@sex=="M"]))
          
          broodPops <- mergePops(list(broodPops, offspring)); rm(offspring)
        }
        
        pop <- mergePops(list(nonReproductive, broodPops)); rm(nonReproductive); rm(broodPops)
        
        # Females emerge from brood, store males until female produces offspring
        pop <- pop[pop@sex=="F"]

      }
    } 
    
    # Death due to natural causes
    pop <- mortality(pop, mortality = mort)
    
    # Death due to acaricide treatment when pop too large
    if (pop@nInd > 2405){
      pop <- mortality(pop, mortality = acaricides)
      
      if (supplement > 0){
        supplementPop <- initPop(suppPop, 0, supplement, simParam = SP)
        
        pop <- mergePops(list(pop,supplementPop[supplementPop@sex=="F"]))
        males <- mergePops(list(males,supplementPop[supplementPop@sex=="M"]))
      }
    }
    
    # Check if there's still varroas left
    if(sum(pop@ebv[,1]) == 0 | pop@nInd > 10000){
      break
    }
    
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
    
  } # end of day
  
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
