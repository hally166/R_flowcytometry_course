#install packages
install.packages("BiocManager")
BiocManager::install("flowCore")
BiocManager::install("flowViz")
BiocManager::install("flowAI")
library(flowCore)
library(flowViz)
library(flowAI)

#load file
myfile <- "C:/FCSfiles/8peak500v.fcs"
fcsfile <- read.FCS(myfile)
fcsfile

#explore file
exprs(fcsfile)[1:10,]
summary(fcsfile)
str(keyword(fcsfile))
summary(fcsfile[,7:24])

#clean the data
fcsfile_flowAI <- flow_auto_qc(fcsfile)
nrow(fcsfile)
nrow(fcsfile_flowAI)

#compensate data
fcsfile_comp <-compensate(fcsfile, spillover(fcsfile)[[1]])

#transform data
chnls <- colnames(fcsfile_comp[,7:24])
chnls
trans <- estimateLogicle(fcsfile_comp, chnls)
fcsfile_trans <- transform(fcsfile_comp, trans)

#basic plotting
plot(fcsfile, c("FSC-A", "SSC-A"))
plot(fcsfile, c("610/20 (561)-A", "450/50 (355)-A"))
plot(fcsfile_trans, c("610/20 (561)-A", "450/50 (355)-A"))
plot(fcsfile_trans)
plot(fcsfile_trans, "610/20 (561)-A")
plot(fcsfile_trans, "610/20 (561)-A", breaks=1024)

#groups of files - flowset
files <- list.files(path="C:/FCSfiles/", pattern=".fcs$")
fs <- read.flowSet(files, path="C:/FCSfiles/")
fs

#clean the flowset
fs_flowAI <- flow_auto_qc(fs)
head(fsApply(fs, nrow))
head(fsApply(fs_flowAI, nrow))

#compensate flowset
comp <-fsApply(fs,function(x)spillover(x)[[1]], simplify=FALSE)
fs_comp <-compensate(fs, comp)

#transform flowset
tf <- estimateLogicle(fs_comp[[1]], channels = colnames(fs[[1]][,7:24]))
fs_trans <- transform(fs_comp, tf)
fs_trans

#explore flowSet
summary(fs_trans[[1]])
head(fsApply(fs_trans, nrow))
fsApply(fs_trans, each_col, median)

#plot flowSet
plot(fs[[1]], c("610/20 (561)-A", "450/50 (355)-A"))
plot(fs_trans[[1]], c("610/20 (561)-A", "450/50 (355)-A"))

#Let's escape flowViz
BiocManager::install("ggcyto")
library(ggcyto)
autoplot(fs_trans[[1]], x="610/20 (561)-A", y="450/50 (355)-A", bins = 256)
autoplot(fs_trans, x="610/20 (561)-A", y="450/50 (355)-A", bins = 256)

#use the + symbol to change the plot paramaters.  Check the docs there are a lot of options
p <- ggcyto(fs_trans, aes(x = "610/20 (561)-A", y =  "450/50 (355)-A"))
p <- p + geom_hex(bins = 128)
p

#Gating - called filters
BiocManager::install("flowWorkspace")
library(flowWorkspace)

gs <- GatingSet(fs_trans)
rg1 <- rectangleGate("FSC-H"=c(100000, Inf), filterId="NonDebris")
add(gs, rg1, parent = "root")
getNodes(gs)
recompute(gs)
autoplot(gs,x = 'FSC-H', y = 'SSC-H', "NonDebris", bins = 256)

rg2 <- rectangleGate("FSC-H"=c(100000, 150000),"FSC-W"=c(50000, 75000))
add(gs, rg2, parent = "NonDebris", name = "singlets")
getNodes(gs)
recompute(gs)
autoplot(gs,x = 'FSC-H', y = 'FSC-W', "singlets", bins = 256)

p <- ggcyto(fs_trans, aes(x = "FSC-H", y = 'FSC-W'))+ geom_hex(bins = 256)
g <- getGate(gs, "singlets")
p <- p + geom_gate(g)
p

#if you need to remove a filter -- don't use this now unless you need to
Rm('singlets', gs)
recompute(gs)

#look at the gated data
plot(gs)
getStats(gs)
getStats(gs, "singlets", "percent")

autoplot(gs[[1]])
fs_singlets <- getData(gs, "/NonDebris/singlets")
fsApply(fs_singlets, each_col, median)

#automatic gating -- clear the R environment by pressing the bruch button to the right
library(flowCore)
library(flowWorkspace)
library(openCyto)

files <- list.files(path="C:/FCSfiles/", pattern=".fcs$")
fs <- read.flowSet(files, path="C:/FCSfiles/")
tf <- estimateLogicle(fs[[1]], channels = colnames(fs[[1]][,7:24]))
fs_trans <- transform(fs, tf)

gs <- GatingSet(fs_trans)

#gate the main population of events
thisData <- getData(gs)
nonDebris_gate <- fsApply(thisData, function(fr) openCyto:::.flowClust.2d(fr, channels = c("FSC-H","SSC-H")))
add(gs, nonDebris_gate, parent = "root", name = "nonDebris")
recompute(gs)
autoplot(gs,x = 'FSC-H', y = 'SSC-H', "nonDebris", bins = 256)

#gate the singlets
thisData <- getData(gs, "nonDebris") #get parent data
singlet_gate <- fsApply(thisData, function(fr) openCyto:::.singletGate(fr, channels =c("FSC-H", "FSC-W")))
add(gs, singlet_gate, parent = "nonDebris", name = "singlets")
recompute(gs)
autoplot(gs,x = 'FSC-H', y = 'FSC-W', "singlets", bins = 256) 

#don't use autoplot to control the plot layout
p1 <- ggcyto(gs, aes(x = 'FSC-H', y = 'SSC-H')) + geom_hex(bins = 128) + geom_gate("nonDebris") + geom_stats(nudge_y=50000)
p1
p2 <- ggcyto(gs, aes(x = 'FSC-H', y = 'FSC-W')) + geom_hex(bins = 128) + geom_gate("singlets") + geom_stats()
p2

#export data
library(gridExtra)
grid.arrange(as.ggplot(p1), as.ggplot(p2), nrow = 2)
g <- arrangeGrob(as.ggplot(p1), as.ggplot(p2), nrow = 2)
ggsave(file="plots.png", g)
write.csv(getPopStats(gs), "stats.csv")
write.csv(exprs(thisData[[1]]), "stats2.csv")

#tSNE clustering
install.packages("Rtsne")
library(Rtsne)
set.seed(42) # Sets seed for reproducibility
thisData <- getData(gs, "singlets") #get parent data
tsne_out <- Rtsne(as.matrix(exprs(thisData[[2]][,7:24])),perplexity=50) # Run TSNE
plot(tsne_out$Y, col=exprs(thisData[[2]][,7])) # Plot the result
