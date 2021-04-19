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
input$nInit <- 10 #c(1, 10, 100, 1000)
input$gd <- c(1,10,50) # c(10,50,100)
input$mort <- 0.005 # daily mortality
input$loci <- 2
input$drive <- c(FALSE, TRUE)
input$testing <- FALSE
input$gaplength <- 30 #c(0,7,14,30)
input$gapday <- 120 # c(100,120,150,200)
input$broodbreak <- 0.1

inputs <- expand.grid(input)
# View(inputs)

#######################################
## ---- Set up parallel processing ----
#######################################

nproc <- 7
cl.tmp = makeCluster(rep('localhost',nproc), type = 'SOCK')
registerDoSNOW(cl.tmp)
getDoParWorkers()

############################
## ---- Run simulations ----
############################

# debug(gd.comp); res <- gd.comp(inputs[21,])

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

# save(res, file = "model_1.23.RData")
# load(file = "model_1.21/model_1.21.RData")

##################
## ---- Plots ----
##################

plotThemVariables <- function(data, parameter, paramName, colParam = NULL, colParamName, save = FALSE){
  p <- ggplot(data = data) +
    geom_rect(mapping = aes(xmin = gapday, xmax = (gapday+gaplength), ymin = -Inf, ymax = Inf), fill="lightgrey") +
    geom_line(mapping = aes(x = day, y = get(x=parameter), colour = get(x=colParam), group = interaction(get(x=colParam),iter)), alpha = 0.5) +
    scale_colour_viridis(option="viridis", discrete=TRUE, name = colParamName, guide=guide_legend(ncol=1, title.position = "top")) +
    theme_bw() +
    ylim(c(0, NA)) +
    xlab("Day") +
    ylab(paramName) +
    scale_x_continuous(breaks=seq(0,5*365,365))
    
  p <- p + facet_grid(gapday ~ drive, labeller = label_both)
  
  if (save){
    PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + theme(legend.position = "right", strip.background = element_blank(),panel.grid = element_blank())
    ggsave(plot = p + PaperTheme, filename = paste(parameter,".png", sep=""), height = 10, width = 10 * 2, unit = "cm")
  }
  
  return(p)
}

p <- plotThemVariables(data = res, 
                        parameter = "GD", 
                        colParam = "broodbreak", 
                        paramName = "Population size", 
                        colParamName = "% of normal brood size")
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

res0 <- res[res$broodbreak==0,]
res1 <- res[res$broodbreak==0.1,]
res2 <- res[res$broodbreak==0.5,]

haploRep <- rename(res0, "GD init" = gd, "Homing" = drive, "GD Homozygotes" = homoGD, "GD Heterozygotes" = heteroGD, "WT Homozygotes" = homoWT) %>%
  pivot_longer(cols = `WT Homozygotes`:`GD Heterozygotes`, names_to = "genotype", values_to = "count") 
haploRep$genotype <- factor(haploRep$genotype, levels = c("WT Homozygotes","GD Heterozygotes","GD Homozygotes"))

p <- ggplot(data = haploRep) + 
  geom_rect(mapping = aes(xmin = gapday, xmax = (gapday+gaplength), ymin = -Inf, ymax = Inf), fill="lightgrey") +
  geom_line(aes(x = day, y = count, group = interaction(genotype, iter), colour = genotype)) +
  facet_grid(gapday ~ Homing, labeller = labeller(.cols = label_both, .rows = label_both)) +
  scale_colour_viridis(option="viridis", discrete=TRUE, name = "Genotype", guide=guide_legend(ncol=1, title.position = "top")) +
  scale_linetype_manual(values = c("solid","dashed")) +
  theme_bw() +
  scale_x_continuous(breaks=seq(0,5*365,365)) +
  ylab("Individuals") +
  xlab("Day") +
  ggtitle("0% brood available")
p

PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + theme(legend.position = "right", strip.background = element_blank(),panel.grid = element_blank(),plot.title = element_text(hjust = 0.5))
ggsave(plot = p + PaperTheme, filename = "genotypes0.png", height = 15, width = 10 * 2, unit = "cm")

haploRep <- rename(res0, "GD init" = gd, "Homing" = drive) %>% 
  pivot_longer(cols = WT:RE, names_to = "allele", values_to = "frequency")
haploRep$allele <- factor(haploRep$allele, levels = c("WT","GD","RE","NF"))

p <- ggplot(data = haploRep) + 
  geom_rect(mapping = aes(xmin = gapday, xmax = (gapday+gaplength), ymin = -Inf, ymax = Inf), fill="lightgrey") +
  geom_line(aes(x = day, y = frequency, group = interaction(allele, iter), colour = allele)) +
  facet_grid(gapday ~ Homing, labeller = labeller(.cols = label_both, .rows = label_both)) +
  scale_colour_viridis(option="viridis", discrete=TRUE, name = "Allele", guide=guide_legend(ncol=1, title.position = "top")) +
  theme_bw() +
  scale_x_continuous(breaks=seq(0,5*365,365)) +
  ylab("Frequency") +
  xlab("Day") +
  ggtitle("0% brood available")
p

PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + theme(legend.position = "right", strip.background = element_blank(),panel.grid = element_blank(),plot.title = element_text(hjust = 0.5))
ggsave(plot = p + PaperTheme, filename = "alleles0.png", height = 15, width = 10 * 2, unit = "cm")

