install.packages("BiocManager")
BiocManager::install("flowCore")
BiocManager::install("flowViz")
BiocManager::install("flowWorkspace")
BiocManager::install("openCyto")
BiocManager::install("ggcyto")
library(flowCore)
library(flowViz)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(gridExtra)

setwd("C:/yourDATA/")

#load files
files <- list.files(path="C:/yourDATA/", pattern=".fcs$")
fs <- read.flowSet(files, path="C:/yourDATA/")

#compensate
comp <-fsApply(fs,function(x)spillover(x)[[1]], simplify=FALSE)
fs_comp <-compensate(fs, comp)

#transform
tf <- estimateLogicle(fs_comp[[1]], channels = colnames(fs[[1]][,8:13]))
fs_trans <- transform(fs_comp, tf)

#create gating set
gs <- GatingSet(fs_trans)

#gate the main population of events
thisData <- getData(gs)
nonDebris_gate <- fsApply(thisData, function(fr) openCyto:::.mindensity(fr, channels = "FSC-A"))
add(gs, nonDebris_gate, parent = "root", name = "nonDebris")
recompute(gs)

#gate the singlets
thisData <- getData(gs, "nonDebris") #get parent data
singlet_gate <- fsApply(thisData, function(fr) openCyto:::.singletGate(fr, channels =c("FSC-A", "FSC-W")))
add(gs, singlet_gate, parent = "nonDebris", name = "singlets")
recompute(gs)

#gate the CD45+
thisData <- getData(gs, "singlets") #get parent data
CD45_gate <- fsApply(thisData, function(fr) openCyto:::.mindensity(fr, channels ="APC-A"))
add(gs, CD45_gate, parent = "singlets", name = "CD45")
recompute(gs)

#gate the CD3 then the CD4/8
thisData <- getData(gs, "CD45") #get parent data
CD3_gate <- fsApply(thisData, function(fr) openCyto:::.mindensity(fr, channels ="FITC-A"))
add(gs, CD3_gate, parent = "CD45", name = "CD3")
recompute(gs)

thisData <- getData(gs, "CD3") #get parent data
Tcell_quad <- fsApply(thisData, function(fr) openCyto:::quadGate.tmix(fr, channels = c("PerCP-Cy5-5-A", "APC-Cy7-A"), K = 3))
add(gs, Tcell_quad, parent = "CD3", name = c("CD4+","CD8+","--","++"))
recompute(gs)

#export the data
f2<- autoplot(gs,x = 'FSC-A', y = 'SSC-A', "nonDebris", bins = 256)
f3<- autoplot(gs,x = 'FSC-A', y = 'FSC-W', "singlets", bins = 256)
f4<- autoplot(gs, x = 'APC-A', y = 'FSC-A', "CD45", bins = 256)
f5<- autoplot(gs, x = 'FITC-A', y = 'FSC-A', "CD3", bins = 256)
f6<- autoplot(gs, x = "PerCP-Cy5-5-A", y = "APC-Cy7-A", c("CD4+","CD8+","--","++"), bins = 256)

g <- grid.arrange(as.ggplot(f2), as.ggplot(f3),as.ggplot(f4), as.ggplot(f5), as.ggplot(f6),nrow = 2)
ggsave(file="plots.png", g)
write.csv(getPopStats(gs), "stats.csv")
