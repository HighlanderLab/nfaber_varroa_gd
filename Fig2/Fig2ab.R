##############################
## ---- Working directory ----
##############################

setwd(dir = '/Users/s2018147/Documents/nfaber_varroa_gd/')
source(file = 'gene_drive_varroa.R')

#############################
## ---- Parameter inputs ----
#############################

input <- list()

input$iter <- 1:10
input$years <- 3
input$nInit <-  c(1, 10, 100, 1000)
input$gd <- 0
input$gdZygosity <- "homozygous"
input$strategy <- 1
input$drive <- FALSE
input$testing <- FALSE
input$mort <- 0.005

input$invasionSlope <- 0.00385
input$invasionIntercept <- -2.87
input$DCP <- 8
input$maxDrone <- 16
input$maxWorker <- 8

input$loci <- 1000
input$inbred <- TRUE
input$inbrCoef <- 0.9

input$gaplength <- 0 
input$gapday <- 0 
input$broodbreak <- 0 
input$acaricides <- 0
input$supplement <- 0

inputs <- expand.grid(input)
# View(inputs)

#######################################
## ---- Set up parallel processing ----
#######################################

nproc <- 6
cl.tmp = makeCluster(rep('localhost',nproc), type = 'SOCK')
registerDoSNOW(cl.tmp)
getDoParWorkers()

############################
## ---- Run simulations ----
############################

# debug(gd.comp); res <- gd.comp(inputs[1,])

## run in parallel
res <- foreach(rowNum = 1:nrow(inputs), .packages = c('AlphaSimR','reshape2','tidyverse','data.table'), .verbose = TRUE) %dopar% { gd.comp(input = inputs[rowNum, ]) }
stopCluster(cl.tmp)
res <- do.call(rbind, res)

############################
## ---- Process results ----
############################

var <- c('nInd','IBS','droneRep','workerRep','potentialInf',
         'nPhoretic','pPhoretic','beesPhoretic','phoreticTime',
         'varroaPerDrone','varroaPerWorker',
         'offspringPerDrone','offspringPerWorker', 
         'WT', 'GD', 'NF', 'RE', 'homoWT', 'homoGD', 'heteroGD')
res[!is.na(res$nInd), var] = lapply(res[!is.na(res$nInd), var], FUN = function(z) { z[is.na(z)] = 0; z })
res = as_tibble(res)

res$iter <- factor(res$iter)
res$nInit <- factor(res$nInit)
res$drive <- factor(res$drive)
res$broodbreak <- factor(res$broodbreak)

res <- group_by(res, iter, nInit)

# save(res, file = "Fig2/Fig2ab.RData")
# load(file = "Fig2/Fig2ab.RData")

##################
## ---- Plots ----
##################

PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + 
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        title=element_text(size=14, hjust=0.5), 
        legend.title=element_text(size=12),
        legend.position = "bottom", 
        legend.justification = "center",
        axis.title=element_text(size=12))

p1 <- ggplot(data = res) +
  geom_line(mapping = aes(x = day, y = nInd, colour = nInit, group = interaction(nInit,iter)), size = 0.5, alpha = 0.5) +
  geom_hline(yintercept=2405, linetype = "dashed", color = "red", size = 1) +
  scale_colour_viridis(option="plasma", discrete=TRUE, name = "Initial population size", guide=guide_legend(nrow=1, title.position = "top")) +
  ylim(c(0, NA)) +
  xlab("Day") +
  ylab("Population size") +
  scale_x_continuous(breaks=seq(0,5*365,365)) +
  PaperTheme
p1

p2 <- ggplot(data = res) +
  geom_line(mapping = aes(x = day, y = IBS, colour = nInit, group = interaction(nInit,iter)), size = 0.5, alpha = 0.5) +
  scale_colour_viridis(option="plasma", discrete=TRUE, name = "Initial population size", guide=guide_legend(nrow=1, title.position = "top")) +
  ylim(c(0.5, NA)) +
  xlab("Day") +
  ylab("Mean homozygosity") +
  scale_x_continuous(breaks=seq(0,5*365,365)) +
  PaperTheme
p2

########################
## ---- Bonus plots ----
########################

plotThemVariables <- function(data, parameter, paramName, colParam = NULL, colParamName, save = FALSE){
  p <- ggplot(data = data) +
    geom_rect(mapping = aes(xmin = gapday, xmax = (gapday+gaplength), ymin = -Inf, ymax = Inf), fill="lightgrey") +
    geom_line(mapping = aes(x = day, y = get(x=parameter), colour = get(x=colParam), group = interaction(get(x=colParam),iter)), size = 0.5, alpha = 0.5) +
    scale_colour_viridis(option="plasma", discrete=TRUE, name = colParamName, guide=guide_legend(ncol=1, title.position = "top")) +
    theme_bw() +
    ylim(c(0, NA)) +
    xlab("Day") +
    ylab(paramName) +
    scale_x_continuous(breaks=seq(0,5*365,365))
  
  if (parameter=="nInd"){
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2405, ymax = Inf, alpha=0.15, fill="red")
  }
  
  if (save){
    PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + theme(legend.position = "right", strip.background = element_blank(),panel.grid = element_blank())
    ggsave(plot = p + PaperTheme, filename = paste(parameter,".png", sep=""), height = 10, width = 10 * 2, unit = "cm")
  }
  
  return(p)
}

varNames <- c('Population size','Mean homozygosity','Invaded drone cells',
              'Invaded worker cells','Potential invasions',
         'Phoretic Varroa','Phoretic Varroa (%)','Phoretic Varroa per bee',
         'Mean time Varroa spent phoretic',
         'Mean Varroa per drone cell','Mean Varroa per worker cell',
         'Mean female offspring per drone cell','Mean female offspring per worker cell',
         'Wildtype', 'Gene drive', 'Non-functional', 'Resistance', 
         'Homozygous wildtype', 'Homozygous gene drive', 'Heterozygous gene drive')

for (variable in var){
  (p <- plotThemVariables(data = res, parameter = variable, colParam = "nInit", 
                paramName = varNames[which(var==variable)], colParamName = "Initial Varroa"))
}
