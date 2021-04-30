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
input$years <- 10
input$nInit <-  100
input$gd <- 500
input$gdZygosity <- "homozygous"
input$strategy <- 1
input$drive <- TRUE
input$testing <- FALSE
input$mort <- 0.005

input$invasionSlope <- 0.00385
input$invasionIntercept <- -2.87
input$DCP <- 8
input$maxDrone <- 16
input$maxWorker <- 8

input$loci <- 2
input$inbred <- FALSE
input$inbrCoef <- 0

input$gaplength <- 0
input$gapday <- 0
input$broodbreak <- 0
input$acaricides <- c(0,0.8,0.95)
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

# save(res, file = "Fig5/Fig5S1.RData")
# load(file = "Fig5/Fig5S1.RData")

##################
## ---- Plots ----
##################

PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + 
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        title=element_text(size=14, hjust=0.5), 
        legend.title=element_text(size=12),
        legend.title.align=0.5,
        legend.position = "bottom", 
        legend.justification = "center",
        axis.title=element_text(size=12))

haploRep <- rename(res, "Acaricides" = acaricides, "GD init" = gd, "Homing" = drive, "GD/GD" = homoGD, "GD Heterozygotes" = heteroGD, "WT/WT" = homoWT) %>%
  pivot_longer(cols = `WT/WT`:`GD Heterozygotes`, names_to = "genotype", values_to = "count") 
haploRep$genotype <- factor(haploRep$genotype, levels = c("WT/WT","GD Heterozygotes","GD/GD"))

p1 <- ggplot(data = haploRep) + 
  geom_line(aes(x = day, y = count, group = interaction(genotype, iter), colour = genotype), size = 0.5, alpha = 0.5) +
  geom_hline(yintercept=2405, linetype = "dashed", color = "red", size = 1) +
  facet_grid( ~ Acaricides, labeller = labeller(.cols = label_both, .rows = label_both)) +
  scale_colour_viridis(option="plasma", end=0.9, discrete=TRUE, name = "Genotype", guide=guide_legend(ncol=3, title.position = "top")) +
  scale_x_continuous(breaks=seq(0,10*365,365), labels = 0:10) +
  ylab("Individuals") +
  xlab("Year") +
  PaperTheme
p1

haploRep <- rename(res, "Acaricides" = acaricides, "GD init" = gd, "Homing" = drive) %>% 
  pivot_longer(cols = WT:RE, names_to = "allele", values_to = "frequency")
haploRep$allele <- factor(haploRep$allele, levels = c("WT","GD","RE","NF"))

p2 <- ggplot(data = haploRep) + 
  geom_line(aes(x = day, y = frequency, group = interaction(allele, iter), colour = allele), size = 0.5, alpha = 0.5) +
  facet_grid(. ~ Acaricides, labeller = labeller(.cols = label_both, .rows = label_both)) +
  scale_colour_viridis(option="plasma", end=0.9, discrete=TRUE, name = "Allele", guide=guide_legend(ncol=4, title.position = "top")) +
  scale_x_continuous(breaks=seq(0,10*365,365), labels = 0:10) +
  ylab("Frequency") +
  xlab("Year") +
  PaperTheme
p2

p <- p1 / p2 + plot_annotation(tag_levels = 'A')
p

ggsave(plot = p + PaperTheme, path = "Fig5", filename = "Fig5S1.png", height = 20, width = 20, unit = "cm")

########################
## ---- Bonus plots ----
########################

plotThemVariables <- function(data, parameter, paramName, colParam = NULL, colParamName, save = FALSE){
  p <- ggplot(data = data) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2405, ymax = Inf, alpha=0.15, fill="red") +
    geom_rect(mapping = aes(xmin = gapday, xmax = (gapday+gaplength), ymin = -Inf, ymax = Inf), fill="lightgrey") +
    geom_line(mapping = aes(x = day, y = get(x=parameter), colour = get(x=colParam), group = interaction(get(x=colParam),iter)), size = 0.5, alpha = 0.5) +
    scale_colour_viridis(option="plasma", discrete=TRUE, name = colParamName, guide=guide_legend(ncol=1, title.position = "top")) +
    theme_bw() +
    ylim(c(0, NA)) +
    xlab("Day") +
    ylab(paramName) +
    scale_x_continuous(breaks=seq(0,5*365,365))
    
  p <- p + facet_grid(gd ~ drive, labeller = label_both)
  
  if (save){
    PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + theme(legend.position = "right", strip.background = element_blank(),panel.grid = element_blank())
    ggsave(plot = p + PaperTheme, filename = paste(parameter,".png", sep=""), height = 10, width = 10 * 2, unit = "cm")
  }
  
  return(p)
}

p <- plotThemVariables(data = res, 
                        parameter = "nInd", 
                        colParam = "inbred", 
                        paramName = "Mean homozygosity", 
                        colParamName = "Inbred", save=TRUE)
p

varNames <- c('Population size','Mean homozygosity','Invaded drone cells',
              'Invaded worker cells','Potential invasions',
         'Phoretic Varroa','Phoretic Varroa (%)','Phoretic Varroa per bee',
         'Mean time Varroa spent phoretic',
         'Mean Varroa per drone cell','Mean Varroa per worker cell',
         'Mean female offspring per drone cell','Mean female offspring per worker cell',
         'Wildtype', 'Gene drive', 'Non-funcitonal', 'Resistance', 
         'Homozygous wildtype', 'Homozygous gene drive', 'Heterozygous gene drive')

for (variable in var){
  (p <- plotThemVariables(data = res1, parameter = variable, colParam = "nInit", 
                paramName = varNames[which(var==variable)], colParamName = "Initial Varroa",
                save = TRUE))
}
