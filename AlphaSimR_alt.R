# Alternative randCross which allows Poisson distribution of offspring
randCrossDistr <- function (pop, nCrosses, balance = TRUE, parents = NULL, 
          ignoreSexes = FALSE, simParam = NULL) 
{
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  if (is.null(parents)) {
    parents = 1:pop@nInd
  }
  else {
    parents = as.integer(parents)
  }
  n = length(parents)
  if (n <= 1) {
    stop("The population must contain more than 1 individual")
  }
  if (simParam$sex == "no" | ignoreSexes) {
    crossPlan = sampHalfDialComb(n, nCrosses)
    crossPlan[, 1] = parents[crossPlan[, 1]]
    crossPlan[, 2] = parents[crossPlan[, 2]]
  }
  else {
    female = which(pop@sex == "F" & (1:pop@nInd) %in% 
                     parents)
    nFemale = length(female)
    if (nFemale == 0) {
      stop("population doesn't contain any females")
    }
    male = which(pop@sex == "M" & (1:pop@nInd) %in% parents)
    nMale = length(male)
    if (nMale == 0) {
      stop("population doesn't contain any males")
    }
    if (balance) {
      female = female[sample.int(nFemale, nFemale)]
      female = rep(female, length.out = nCrosses)
      tmp = male[sample.int(nMale, nMale, replace = TRUE)]
      n = nCrosses%/%nMale + 1
      male = NULL
      for (i in 1:n) {
        take = nMale - (i:(nMale + i - 1))%%nMale
        male = c(male, tmp[take])
      }
      male = male[1:nCrosses]
      crossPlan = cbind(female, male)
    }
    else {
      crossPlan = sampAllComb(nFemale, nMale, nCrosses)
      crossPlan[, 1] = female[crossPlan[, 1]]
      crossPlan[, 2] = male[crossPlan[, 2]]
    }
  }
  return(crossPlan)
}

# Alternative makeCross which allows Poisson distribution of offspring
makeCrossDistr <- function (pop, crossPlan, nProgeny = 1, simParam = NULL) 
{
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  if (pop@ploidy%%2L != 0L) {
    stop("You cannot cross aneuploids")
  }
  if (is.character(crossPlan)) {
    crossPlan = cbind(match(crossPlan[, 1], pop@id), match(crossPlan[, 
                                                                     2], pop@id))
    if (any(is.na(crossPlan))) {
      stop("Failed to match supplied IDs")
    }
  }
  if ((max(crossPlan) > nInd(pop)) | (min(crossPlan) < 1L)) {
    stop("Invalid crossPlan")
  }
  if (length(nProgeny) > 1) {
    crossPlan = cbind(rep(crossPlan[, 1], nProgeny), 
                      rep(crossPlan[, 2], nProgeny))
  } else if (nProgeny > 1) {
    crossPlan = cbind(rep(crossPlan[, 1], each = nProgeny), 
                      rep(crossPlan[, 2], each = nProgeny))
  }
  tmp = AlphaSimR:::cross(pop@geno, crossPlan[, 1], pop@geno, crossPlan[, 
                                                            2], simParam$femaleMap, simParam$maleMap, simParam$isTrackRec, 
              pop@ploidy, pop@ploidy, simParam$v, simParam$femaleCentromere, 
              simParam$maleCentromere, simParam$quadProb, simParam$nThreads)
  rPop = new("RawPop", nInd = nrow(crossPlan), nChr = pop@nChr, 
             ploidy = pop@ploidy, nLoci = pop@nLoci, geno = tmp$geno)
  if (simParam$isTrackRec) {
    simParam$addToRec(tmp$recHist)
  }
  return(newPop(rawPop = rPop, mother = pop@id[crossPlan[, 
                                                         1]], father = pop@id[crossPlan[, 2]], simParam = simParam))
}

makeCross2Distr <- function (females, males, crossPlan, nProgeny = 1, simParam = NULL) 
{
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  if ((females@ploidy%%2L != 0L) | (males@ploidy%%2L != 0L)) {
    stop("You can not cross aneuploids")
  }
  if (is.character(crossPlan)) {
    crossPlan = cbind(match(crossPlan[, 1], females@id), 
                      match(crossPlan[, 2], males@id))
    if (any(is.na(crossPlan))) {
      stop("Failed to match supplied IDs")
    }
  }
  if ((max(crossPlan[, 1]) > nInd(females)) | (max(crossPlan[, 
                                                             2]) > nInd(males)) | (min(crossPlan) < 1L)) {
    stop("Invalid crossPlan")
  }
  if (length(nProgeny) > 1) {
    crossPlan = cbind(rep(crossPlan[, 1], nProgeny), 
                      rep(crossPlan[, 2], nProgeny))
  } else if (nProgeny > 1) {
    crossPlan = cbind(rep(crossPlan[, 1], each = nProgeny), 
                      rep(crossPlan[, 2], each = nProgeny))
  }
  tmp = AlphaSimR:::cross(females@geno, crossPlan[, 1], males@geno, crossPlan[, 
                                                                  2], simParam$femaleMap, simParam$maleMap, simParam$isTrackRec, 
              females@ploidy, males@ploidy, simParam$v, simParam$femaleCentromere, 
              simParam$maleCentromere, simParam$quadProb, simParam$nThreads)
  rPop = new("RawPop", nInd = nrow(crossPlan), nChr = females@nChr, 
             ploidy = as.integer((females@ploidy + males@ploidy)/2), 
             nLoci = females@nLoci, geno = tmp$geno)
  if (simParam$isTrackRec) {
    simParam$addToRec(tmp$recHist)
  }
  return(newPop(rawPop = rPop, mother = females@id[crossPlan[, 
                                                             1]], father = males@id[crossPlan[, 2]], simParam = simParam))
}

editGenomeFix <- function (pop, ind, chr, segSites, allele, simParam = NULL) 
{
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  ind = unique(as.integer(ind))
  stopifnot(all(ind %in% (1:pop@nInd)))
  chr = as.integer(chr)
  segSites = as.integer(segSites)
  stopifnot(length(chr) == length(segSites))
  allele = as.integer(allele)
  stopifnot(all(allele == 0L | allele == 1L))
  allele = as.raw(allele)
  if(length(allele) == 1L){
    allele = rep(allele, length(segSites))
  }
  stopifnot(length(allele) == length(segSites))
  for (selChr in unique(chr)) {
    sel = which(chr == selChr)
    for (i in sel) {
      BYTE = (segSites[i] - 1L)%/%8L + 1L
      BIT = (segSites[i] - 1L)%%8L + 1L
      for (selInd in ind) {
        for (j in 1:pop@ploidy) {
          TMP = pop@geno[[selChr]][BYTE, j, selInd]
          TMP = rawToBits(TMP)
          TMP[BIT] = allele[i]
          TMP = packBits(TMP)
          pop@geno[[selChr]][BYTE, j, selInd] = TMP
        }
      }
    }
  }
  PHENO = pop@pheno
  EBV = pop@ebv
  pop = resetPop(pop = pop, simParam = simParam)
  pop@pheno = PHENO
  pop@ebv = EBV
  return(pop)
}

editHaplo <- function (pop, ind, chr, segSites, allele, haplotype = 1:pop@ploidy ,simParam = NULL) 
{
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  ind = unique(as.integer(ind))
  stopifnot(all(ind %in% (1:pop@nInd)))
  chr = as.integer(chr)
  segSites = as.integer(segSites)
  stopifnot(length(chr) == length(segSites))
  allele = as.integer(allele)
  stopifnot(all(allele == 0L | allele == 1L))
  allele = as.raw(allele)
  if(length(allele) == 1L){
    allele = rep(allele, length(segSites))
  }
  stopifnot(length(allele) == length(segSites))
  for (selChr in unique(chr)) {
    sel = which(chr == selChr)
    for (i in sel) {
      BYTE = (segSites[i] - 1L)%/%8L + 1L
      BIT = (segSites[i] - 1L)%%8L + 1L
      for (selInd in ind) {
        for (j in haplotype) {
          TMP = pop@geno[[selChr]][BYTE, j, selInd]
          TMP = rawToBits(TMP)
          TMP[BIT] = allele[i]
          TMP = packBits(TMP)
          pop@geno[[selChr]][BYTE, j, selInd] = TMP
        }
      }
    }
  }
  PHENO = pop@pheno
  EBV = pop@ebv
  pop = resetPop(pop = pop, simParam = simParam)
  pop@pheno = PHENO
  pop@ebv = EBV
  return(pop)
}

quickHaploInbr <- function (nInd, nChr, segSites, genLen = 1, ploidy = 2L, inbred = FALSE, inbrCoef = 0.5){
  ploidy = as.integer(ploidy)
  nInd = as.integer(nInd)
  nChr = as.integer(nChr)
  segSites = as.integer(segSites)
  if (length(segSites) == 1) 
    segSites = rep(segSites, nChr)
  if (length(genLen) == 1) 
    genLen = rep(genLen, nChr)
  nBins = segSites%/%8L + (segSites%%8L > 0L)
  centromere = genLen/2
  genMap = vector("list", nChr)
  geno = vector("list", nChr)
  for (i in 1:nChr) {
    genMap[[i]] = seq(0, genLen[i], length.out = segSites[i])
    geno[[i]] = array(sample(as.raw(0:255), nInd * ploidy * 
                               nBins[i], replace = TRUE), dim = c(nBins[i], ploidy, 
                                                                  nInd))
    if (inbred) {
      if (ploidy > 1) {
        for (j in 2:ploidy) {
          homozygous <- runif(nBins) < inbrCoef
          geno[[i]][homozygous, j, ] = geno[[i]][homozygous, 1, ]
        }
      }
    }
  }
  return(new("MapPop", nInd = nInd, nChr = nChr, ploidy = ploidy, 
             nLoci = segSites, geno = as.matrix(geno), genMap = as.matrix(genMap), 
             centromere = centromere))
}
