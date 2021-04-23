##############################
## ---- Working directory ----
##############################

setwd(dir = '/Users/s2018147/Documents/nfaber_varroa_gd/')
source(file = 'gene_drive_varroa.R')

#############################
## ---- Parameter inputs ----
#############################

input <- list()

input$iter <- 1
input$years <- 1
input$nInit <-  100
input$gd <- 10
input$gdZygosity <- "homozygous"
input$strategy <- 1
input$drive <- TRUE
input$testing <- FALSE
input$mort <- 0.005

input$invasionSlope <- seq(0.00185,0.00585,0.001)
input$invasionIntercept <- seq(-1.87,-3.87,-0.5)
input$DCP <- seq(0,20,4)
input$maxDrone <- seq(4,28,4)
input$maxWorker <- seq(2,14,2)

input$loci <- 2
input$inbred <- FALSE
input$inbrCoef <- 0

input$gaplength <- 0 
input$gapday <- 0 
input$broodbreak <- 0 
input$acaricides <- 0
input$supplement <- 0

inputs <- expand.grid(input)
# View(inputs)

inputs <- inputs %>%
  filter((invasionSlope==0.00385 & invasionIntercept==-2.87 & DCP==8 & maxDrone==16 & maxWorker==8) |
           (invasionSlope==0.00385 & invasionIntercept==-2.87 & DCP==8 & maxDrone==16) |
           (invasionSlope==0.00385 & invasionIntercept==-2.87 & DCP==8 & maxWorker==8) |
           (invasionSlope==0.00385 & invasionIntercept==-2.87 & maxDrone==16 & maxWorker==8) |
           (invasionSlope==0.00385 & DCP==8 & maxDrone==16 & maxWorker==8) |
           (invasionIntercept==-2.87 & DCP==8 & maxDrone==16 & maxWorker==8))

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

# save(res, file = "Fig2/Fig2S2.RData")
# load(file = "Fig2/Fig2S2.RData")

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

day365 <- res[res$day==365,]

day365$setje <- rep(1:26, each = 10)
day365$setje <- factor(day365$setje)

day365$variableRange <- case_when(
  day365$invasionSlope != 0.00385 ~ "invasionSlope",
  day365$invasionIntercept != -2.87 ~ "invasionIntercept",
  day365$DCP != 8 ~ "DCP",
  day365$maxDrone != 16 ~ "maxDrone",
  day365$maxWorker != 8 ~ "maxWorker",
  TRUE ~ "default"
)

## Duplicate the default for every variable range
defaults <- day365[day365$variableRange=="default",] %>% 
  slice(rep(1:n(), times = 5)) %>%
  dplyr::mutate(variableRange = c("invasionSlope","invasionIntercept","DCP","maxDrone","maxWorker"))
day365 <- day365[day365$variableRange!="default",]
day365 <- rbind(day365,defaults)

day365 <- group_by(day365, variableRange, setje, invasionSlope,invasionIntercept,DCP,maxDrone,maxWorker) 

## Pivot whole data frame larger
testje <- pivot_longer(day365, cols = c(invasionSlope,invasionIntercept,DCP,maxDrone,maxWorker),
                      names_to = "Variable", values_to = "Value")
testje <- pivot_longer(testje, cols = c(nInd,GD),
                      names_to = "Outcome", values_to = "Number")

## Pivot whole data frame larger
defaults <- pivot_longer(defaults, cols = c(invasionSlope,invasionIntercept,DCP,maxDrone,maxWorker),
                       names_to = "Variable", values_to = "Value")
defaults <- pivot_longer(defaults, cols = c(nInd,GD),
                       names_to = "Outcome", values_to = "Number")

stats <- testje %>%
  group_by(variableRange, setje, Variable, Value, Outcome) %>%
  summarise(mean = mean(Number),
            stDev = sd(Number), 
            count = n(), .groups = "keep") %>%
  dplyr::mutate(stErr = stDev / sqrt(count), 
                lowerCI = mean - qt(1 - (0.05 / 2), count - 1) * stErr,
                upperCI = mean + qt(1 - (0.05 / 2), count - 1) * stErr)

testje <- testje %>% 
  filter(Variable==variableRange)
stats <- stats %>% 
  filter(Variable==variableRange)
defaults <- defaults %>% 
  filter(Variable==variableRange) %>%
  filter(iter==1)

testje$Outcome <- factor(testje$Outcome, levels = c("nInd","GD"), labels = c("Population size", "GD frequency"))
testje$variableRange <- factor(testje$variableRange, 
                               levels = c("invasionSlope","invasionIntercept","DCP","maxDrone","maxWorker"), 
                               labels = c("Invasion rate slope","Invasion rate intercept","Drone cell preference","Max offspring per drone cell","Max offspring per worker cell"))
stats$Outcome <- factor(stats$Outcome, levels = c("nInd","GD"), labels = c("Population size", "GD frequency"))
stats$variableRange <- factor(stats$variableRange, 
                               levels = c("invasionSlope","invasionIntercept","DCP","maxDrone","maxWorker"), 
                               labels = c("Invasion rate slope","Invasion rate intercept","Drone cell preference","Max offspring per drone cell","Max offspring per worker cell"))
defaults$Outcome <- factor(defaults$Outcome, levels = c("nInd","GD"), labels = c("Population size", "GD frequency"))
defaults$variableRange <- factor(defaults$variableRange, 
                              levels = c("invasionSlope","invasionIntercept","DCP","maxDrone","maxWorker"), 
                              labels = c("Invasion rate slope","Invasion rate intercept","Drone cell preference","Max offspring per drone cell","Max offspring per worker cell"))
defaults$Number <- rep(c(2500,0.2), times=5)
defaults$min <- rep(c(0.00335,-3.12,6,14,7), each=2)
defaults$max <- rep(c(0.00435,-2.62,10,18,9), each=2)

p <- ggplot() +
  geom_point(data = testje, mapping = aes(x = Number, y = Value), colour = "mediumvioletred", size=3, shape=1) +
  geom_errorbarh(data = stats, mapping = aes(y = Value, xmin=lowerCI, xmax=upperCI), height=0) +
  geom_label(data = defaults, mapping = aes(label = "Default", x = Number, y = Value), size=3, label.size=0, hjust=0) +
  geom_rect(data=defaults, aes(xmin=-Inf,xmax=Inf,ymin=min,ymax=max),
            alpha=0.1) +
  facet_grid(variableRange ~ Outcome, scales = "free") + 
  xlab("Population size and GD frequency on day 365") +
  ylab("Parameter value") +
  PaperTheme
p

ggsave(plot = p, path = "Fig2", filename = "Fig2S2.png", height = 25, width = 20, unit = "cm")

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
