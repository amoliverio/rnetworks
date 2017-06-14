# Soil Network for top taxa

library(ggraph)
library(igraph)
library(dplyr)
library(reshape2)

setwd("/Users/Angela/desktop/manu_500soil/full_dataset/65_cutoff_v2/random_network_generation/test_graph_generation/")

#------------------------------------------------------------------------------

# import with node phyla names
edges <- read.csv("../ggraph_edges.csv")
nodes<- read.csv('../ggraph_nodes_wMetadata_v2.csv')

# generate random graphs with gen_rand_edges.R function (erdos.renyi model)
source("FUN_gen_rand_edges.R")
rand_edges <- gen_rand_edges(x=1000, n=270, pm=3646) # with obs. # nodes + edges

# build node dataframe with ENV factors randomly assigned by sampling
random_nodes <- data.frame(matrix(NA, ncol = 2, nrow = 270))
random_nodes$id <- seq(1,270, by = 1)
random_nodes$ENV <- random_nodes$X2
random_nodes$X1 <- NULL
random_nodes$X2 <- NULL
ENV <- as.vector(nodes$cluster_info_reduced)

# loop to assign nodes ENV factors in propotion from RF
rand_nodes <- list()

for(i in 1:1000){
    random_nodes$ENV <- sample(ENV)
    rand_nodes[[i]] <- random_nodes
}

# merge ENV info with edges list
full_graphs <- list()
pairwise_freq <- list()
within_group_edges <- vector("numeric", length=1000)

for(i in 1:1000){
    
    full <- as.data.frame(rand_edges[[i]])
    
    # for source edges
    full$id <- full$V1
    full <- merge(full, rand_nodes[[i]], "id")
    full$source_ENV <- full$ENV
    full$ENV <- NULL
    
    # for target edges
    full$id <- full$V2
    full <- merge(full, rand_nodes[[i]], "id")
    full$target_ENV <- full$ENV
    full$ENV <- NULL
    
    # remove ids
    full$id1 <- NULL
    full$id2 <- NULL
    
    # put obs into the list
    full_graphs[[i]] <- full
    
    # and pairwise tables
    table <- as.data.frame(as.matrix(table(full$source_ENV, full$target_ENV)))
    pairwise_freq[[i]] <- table
    
    # and total within group obs (subset of the 3646)
    full_matched <- full[full$source_ENV==full$target_ENV,]
    within_group_edges[i]<- length(full_matched$V1)
}

rm(full_matched)

# no. of edges that match - across groups
within <- as.data.frame(within_group_edges)
ggplot(within, aes(within_group_edges)) + geom_density() +
    theme_bw() + geom_vline(xintercept=1766, color="red")

# How robust is each ecological cluster?
# extract info from list and prep for visualization
pairwise_df <- data.frame(matrix(vector(mode = 'numeric',length = 1), nrow = 49, ncol = 1000))

for(i in 1:1000){
    p1 <- as.data.frame(pairwise_freq[[i]])
    p2 <- p1[with(p1, order(Var1, Var2)), ]
    pairwise_df[,i] <- p2$Freq
}

avg_edge_connections <- as.data.frame(rowMeans(pairwise_df))
avg_edge_connections$Var1 <- p2$Var1
avg_edge_connections$Var2 <- p2$Var2

write.csv(avg_edge_connections, 'avg_edge_connections_v2.csv')

pairwise_df$Var1 <- p2$Var1
pairwise_df$Var2 <- p2$Var2

pairwise_within_filt <- filter(pairwise_df, Var2 == Var1)
rownames(pairwise_within_filt) <- pairwise_within_filt$Var1
pairwise_within_filt$Var1 <- NULL
pairwise_within_filt$Var2 <- NULL
pairwise_within_filt_T <- as.data.frame(t(pairwise_within_filt))
pairwise_df <- melt(pairwise_within_filt_T)
pairwise_df$Var1 <- NULL

rm(rand_edges, random_nodes, pairwise_freq, full_graphs)

#------------------------------------------------------------------------------
# visualize with histograms for each and then make a multiplot
#------------------------------------------------------------------------------

# Make vectors for each
Dry_forests <- pairwise_within_filt_T$`Dry-forests`
Dry_forests <- as.data.frame(Dry_forests)

Drylands <- pairwise_within_filt_T$Drylands
Drylands <- as.data.frame(Drylands)

High_pH <- pairwise_within_filt_T$`High_pH`
High_pH <- as.data.frame(High_pH)

Low_pH <- pairwise_within_filt_T$`Low_pH`
Low_pH <- as.data.frame(Low_pH)

Low_productivity <- pairwise_within_filt_T$`Low_productivity`
Low_productivity <- as.data.frame(Low_productivity)

Multi_factor <- pairwise_within_filt_T$Multi_factor
Multi_factor <- as.data.frame(Multi_factor)

Undefined <- pairwise_within_filt_T$Undefined
Undefined <- as.data.frame(Undefined)

# Observed within cluster values
df_obs = 19
dl_obs = 125 
hph_obs = 1353
lph_obs = 78
lp_obs = 30
mf_obs = 139
uf_obs = 22

# Calculate p-values for each cluster (pvalues for permutation tests)
# add 1 to numerator and denominator to account for misestimation of p-value
# as per Phipson and Smyth 2010, "Permutation P-values should never be zero"

# the proporiton of permuatations with larger value
DF_pval <- (sum((Dry_forests$Dry_forests) >= df_obs)+1) / (length(Dry_forests$Dry_forests) + 1)
DF_label <- as.vector(paste("p <", round(DF_pval, digits=5)))

DL_pval <- (sum((Drylands$Drylands) >= dl_obs)+1) / (length(Drylands$Drylands) + 1)
DL_label <- as.vector(paste("p <", round(DL_pval, digits=5)))

HPH_pval <- (sum((High_pH$High_pH) >= hph_obs)+1) / (length(High_pH$High_pH) + 1)
HPH_label <- as.vector(paste("p <", round(HPH_pval, digits=5)))

LPH_pval <- (sum((Low_pH$Low_pH) >= lph_obs)+1) / (length(Low_pH$Low_pH) + 1)
LPH_label <- as.vector(paste("p <", round(LPH_pval, digits=5)))

LP_pval <- (sum((Low_productivity$Low_productivity) >= lp_obs)+1) / (length(Low_productivity$Low_productivity) + 1)
LP_label <- as.vector(paste("p <", round(LP_pval, digits=5)))

MF_pval <- (sum((Multi_factor$Multi_factor) >= mf_obs)+1) / (length(Multi_factor$Multi_factor) + 1)
MF_label <- as.vector(paste("p <", round(MF_pval, digits=5)))

Uf_pval <- (sum((Undefined$Undefined) >= uf_obs)+1) / (length(Undefined$Undefined) + 1)
Uf_label <- as.vector(paste("p =", round(Uf_pval, digits=5)))

# Make plots for each
DF <- ggplot(data= Dry_forests, aes(Dry_forests)) + 
    geom_histogram(bins=20,aes(y = ..density..), fill = "gray") +
    stat_function(fun=dnorm, args=list(mean=mean(Dry_forests$Dry_forests),
                                       sd=sd(Dry_forests$Dry_forests))) +
    geom_vline(xintercept = df_obs, color="red", linetype=2) +
    annotate("text", x=3, y=.20, label= DF_label) +
    theme_classic()

DL <- ggplot(data= Drylands, aes(Drylands)) + 
    geom_histogram(bins=50,aes(y = ..density..), fill = "gray") +
    stat_function(fun=dnorm, args=list(mean=mean(Drylands$Drylands),
                                       sd=sd(Drylands$Drylands))) +
    geom_vline(xintercept = dl_obs, color="red", linetype=2) +
    annotate("text", x=50, y=.08, label= DL_label) +
    theme_classic()

HPH <- ggplot(data= High_pH, aes(High_pH)) + 
    geom_histogram(bins=50,aes(y = ..density..), fill = "gray") +
    stat_function(fun=dnorm, args=list(mean=mean(High_pH$High_pH),
                                       sd=sd(High_pH$High_pH))) +
    geom_vline(xintercept = hph_obs, color="red", linetype=2) +
    annotate("text", x=600, y=.02, label= HPH_label) +
    theme_classic()

LPH <- ggplot(data= Low_pH, aes(Low_pH)) + 
    geom_histogram(bins=50,aes(y = ..density..), fill = "gray") +
    stat_function(fun=dnorm, args=list(mean=mean(Low_pH$Low_pH),
                                       sd=sd(Low_pH$Low_pH))) +
    geom_vline(xintercept = lph_obs, color="red", linetype=2) +
    annotate("text", x=25, y=.09, label= LPH_label) +
    theme_classic()

LP <- ggplot(data= Low_productivity, aes(Low_productivity)) + 
    geom_histogram(bins=20,aes(y = ..density..), fill = "gray") +
    stat_function(fun=dnorm, args=list(mean=mean(Low_productivity$Low_productivity),
                                       sd=sd(Low_productivity$Low_productivity))) +
    geom_vline(xintercept = lp_obs, color="red", linetype=2) +
    annotate("text", x=2, y=.15, label= LP_label) +
    theme_classic()

MF <- ggplot(data= Multi_factor, aes(Multi_factor)) + 
    geom_histogram(bins=50,aes(y = ..density..), fill = "gray") +
    stat_function(fun=dnorm, args=list(mean=mean(Multi_factor$Multi_factor),
                                       sd=sd(Multi_factor$Multi_factor))) +
    geom_vline(xintercept = mf_obs, color="red", linetype=2) +
    annotate("text", x=60, y=.05, label= MF_label) +
    theme_classic()

Uf <- ggplot(data= Undefined, aes(Undefined)) + 
    geom_histogram(bins=50,aes(y = ..density..), fill = "gray") +
    stat_function(fun=dnorm, args=list(mean=mean(Undefined$Undefined),
                                       sd=sd(Undefined$Undefined))) +
    geom_vline(xintercept = uf_obs, color="red", linetype=2) +
    annotate("text", x=50, y=.03, label= Uf_label) +
    theme_classic()

source("/Users/Angela/Desktop/research/R_functions/multiplot.R")
multiplot(HPH, LP, DL, DF, LPH, cols=2)
multiplot(HPH, LP, DL, DF, LPH, MF, Uf, cols=3)
