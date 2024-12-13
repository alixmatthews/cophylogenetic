# R CODE FOR PARAFIT ANALYSIS
# updated 10 December 2024 by Alix Matthews





# ************************************* #
#### SET UP ####
# ************************************* #


##### 1) SET WORKING DIRECTORY #####
setwd("/Users/alix22792/Desktop/parafit-testing")


##### 2) LOAD LIBRARIES #####
library(ape) # version 5.8
library(phangorn) # version 2.12.2
library(parallel) # version 4.3.1
library(paco) # version 0.4.2
library(phytools) # version 2.3-0, can also be downloaded from GitHub: https://github.com/liamrevell/phytools


##### 3) SET ALPHA LEVEL #####
p_level <- 0.0549





# ************************************* #
#### LOAD AND CLEAN UP TREES ####
# ************************************* #


##### 1) SYMBIONT ("PARASITE") TREE #####
# this is a nexus format. Otherwise, use "read.tree"
original.p.tree <- read.nexus("/Users/alix22792/Desktop/parafit-testing/fake-parasite.nexus")

# view tree (if desired)
plot(original.p.tree)


#### REMOVE SYMBIONT OUTGROUPS OR DUPLICATE TAXA
# read in tip labels (to avoid any typos on drop.tip() step)
original.p.tree$tip.label


#### DROP UNWANTED TIPS
trimmed.p.tree <- drop.tip(original.p.tree, c('P10'))

# view trimmed tree, if desired
plot(trimmed.p.tree)

# Change tip labels on tree, if desired
# If you do this, these new labels should match what is in the matrix file
# first, read the remaining tips in to see the correct order needed for relabeling
trimmed.p.tree$tip.label

# relabel them in correct order, if desired
new_tiplabels.p.tree <- c("Parasite_sp1", "Parasite_sp2", "Parasite_sp3", "Parasite_sp4", "Parasite_sp5", "Parasite_sp6", "Parasite_sp7", "Parasite_sp8", "Parasite_sp9")

# assign new labels to tree, if desired
trimmed.p.tree$tip.label <- new_tiplabels.p.tree

# view relabeled tree, if desired
plot(trimmed.p.tree)



##### 2) HOST TREE #####
# this is a nexus format. Otherwise, use "read.tree"
original.h.tree <- read.nexus("/Users/alix22792/Desktop/parafit-testing/fake-host.nexus")

# view tree (if desired)
plot(original.h.tree)


#### REMOVE HOST OUTGROUPS OR DUPLICATE TAXA
# read in tip labels (to avoid any typos on drop.tip() step)
original.h.tree$tip.label


#### DROP UNWANTED TIPS
trimmed.h.tree <- drop.tip(original.h.tree, c("H37","H36"))

# view trimmed tree, if desired
plot(trimmed.h.tree)

# Change tip labels on tree, if desired
# If you do this, these new labels should match what is in the matrix file
# first, read the remaining tips in to see the correct order needed for relabeling
trimmed.h.tree$tip.label

# relabel them in correct order, if desired
new_tiplabels.trimmed.h.tree <- c("Host_sp32", "Host_sp31", "Host_sp30", "Host_sp33", "Host_sp35", "Host_sp34", "Host_sp2", "Host_sp4", "Host_sp5", "Host_sp3", "Host_sp6", "Host_sp7", "Host_sp8", "Host_sp9", "Host_sp10", "Host_sp14", "Host_sp16", "Host_sp12", "Host_sp13", "Host_sp15", "Host_sp25", "Host_sp1", "Host_sp24", "Host_sp21", "Host_sp22", "Host_sp18", "Host_sp17", "Host_sp23", "Host_sp20", "Host_sp19", "Host_sp29", "Host_sp26", "Host_sp27", "Host_sp28", "Host_sp11")

# assign new labels to tree, if desired
trimmed.h.tree$tip.label <- new_tiplabels.trimmed.h.tree

# view relabeled tree, if desired
plot(trimmed.h.tree)





# ************************************* #
#### LOAD IN HOST-SYMBIONT MATRIX ####
# ************************************* #


##### 1) READ IN MATRIX #####
##### IMPORTANT: make sure labels match with every tree
hpmatrix <- as.matrix(read.table("/Users/alix22792/Desktop/parafit-testing/fake-matrix.txt", header=T))

# view to check it, if desired
# View(hpmatrix)


##### 2) CALCULATE PATRISTIC DISTANCES AND ROTATE HEADINGS TO FIT MATRIX #####
# symbionts
symbiont.patristic.dist <- as.matrix(cophenetic(trimmed.p.tree))
symbiont.patristic.dist <- symbiont.patristic.dist[colnames(hpmatrix), colnames(hpmatrix)]

# hosts
host.patristic.dist <- as.matrix(cophenetic(trimmed.h.tree))
host.patristic.dist <- host.patristic.dist[rownames(hpmatrix), rownames(hpmatrix)]





# ************************************* #
#### RUN PARAFIT ####
# ************************************* #


##### 1) LOOP MULTIPLE PARAFIT RUNS AND PRINT RESULTS #####
# set up the loop
runParaFit <- function(seed) {
    set.seed(seed)  # set a unique seed for each parafit run to ensure randomness

    # run parafit
    result <- parafit(
        host.patristic.dist,
        symbiont.patristic.dist,
        hpmatrix,
        nperm = 999,
        test.links = TRUE,
        correction = 'cailliez' # 'lingoes' is another correction option, the default is 'none'
    )

    # extract key metrics
    global_value <- result$ParaFitGlobal
    p_global_value <- result$p.global
    p_F1 <- result$link.table[, 4]  # column 4 corresponds to p.F1
    p_F2 <- result$link.table[, 6]  # column 6 corresponds to p.F2

    # return a list containing all extracted key metrics values
    return(list(global = global_value, p.global = p_global_value, p.F1 = p_F1, p.F2 = p_F2, link.table = result$link.table))
}


# set the number of runs and prep the seeds (usually 100, but can use fewer for testing)
num_runs <- 100
seeds <- sample(1:10000, num_runs)


# run the loop (the parallelized ParaFit with unique seeds)
# this step is the the most CPU-intensive
# change mc.cores value to something that fits your machine
res <- mclapply(seeds, runParaFit, mc.preschedule = TRUE, mc.cores = 10)


# convert results into a structured format
global_values <- sapply(res, function(x) x$global)
p_global_values <- sapply(res, function(x) x$p.global)


# calculate the mean ParaFitGlobal and p.global values
mean_ParaFitGlobal <- mean(global_values)
mean_p_global <- mean(p_global_values)


# extract link table from the first result (this assumes link structure is consistent across runs, which it should be)
link_names <- cbind(
    rownames(hpmatrix)[res[[1]]$link.table[, 1]],  # row names of links (hosts)
    colnames(hpmatrix)[res[[1]]$link.table[, 2]]  # column names of links (parasites)
)


# !VERY IMPORTANT!
# check if R correctly ordered host-symbiont links (this links to p-values)
# the output of this should be your expected host-symbiont associations
link_names


# create a matrix for p.F1 and p.F2 across runs
p_F1_matrix <- do.call(rbind, lapply(res, function(x) x$p.F1))
p_F2_matrix <- do.call(rbind, lapply(res, function(x) x$p.F2))


# adjust p-values using BH correction
# NOTE: if p_F*_matrix values are relatively homogeneous and/or large across runs, it is likely that they will all be ranked closely and BH adjustment could lead to the same value for all runs. If the p_F*_matrix values are quite different and there is lots of variation across runs, you will see different adjusted p-values for each run
p_F1_adj <- apply(p_F1_matrix, 2, p.adjust, method = "BH")
p_F2_adj <- apply(p_F2_matrix, 2, p.adjust, method = "BH")


# calculate mean adjusted p-values for each individual link
p_F1_adj_means <- colMeans(p_F1_adj)
p_F2_adj_means <- colMeans(p_F2_adj)


# combine results into a table for readability
p_means_table <- data.frame(
    Host = link_names[, 1],
    Symbiont = link_names[, 2],
    Mean_p_F1_adj = p_F1_adj_means,
    Mean_p_F2_adj = p_F2_adj_means
)


# output the mean ParaFitGlobal and p.global values
cat("Mean ParaFitGlobal:", mean_ParaFitGlobal, "| Mean p.global:", mean_p_global, "\n")


# output the table of adjusted p-values for each host-symbiont link
print(p_means_table)


##### 2) OPTIONAL DETAILED REPORT #####
# you can ignore warnings... it's about adding the columns names to the file for each run.

# save detailed ParaFit results to a log file, if you want to cross-check
log_file <- "parafit_results_log.txt"
writeLines("ParaFit Results Log\n", log_file)  # start by writing a header to the file


# initialize a flag to control column names
write_column_names <- TRUE

for (i in seq_along(res)) {
    result <- res[[i]]
    link_table <- result$link.table

    # get host and symbiont names from the hpmatrix based on indices in link_table
    host_names <- rownames(hpmatrix)[link_table[, 1]]
    symbiont_names <- colnames(hpmatrix)[link_table[, 2]]

    # append global test results
    cat(
        sprintf(
            "\nRun %d\nGlobal Test:\nParaFitGlobal = %f, p.global = %f\n",
            i, result$global, result$p.global
        ),
        file = log_file,
        append = TRUE
    )

    # append link table results with species names instead of numeric indices
    cat("\nHost-Parasite Links:\n", file = log_file, append = TRUE)
    write.table(
        cbind(Host = host_names, Symbiont = symbiont_names, link_table),
        file = log_file,
        append = TRUE,
        quote = FALSE,
        sep = "\t",
        col.names = TRUE,
        row.names = FALSE
    )
}





# ************************************* #
#### TANGLEGRAM PLOTTING ####
# ************************************* #


##### 1) INITIALIZE COLORS/LINE TYPES #####
# initialize vectors to store colors and line types for h-p links
# this uses the adjusted p-values, which might not vary across some runs after correction (see line 206)
# you can also replace p_F1_adj with p_F1_matrix (non-adjusted p-values) if you'd rather display uncorrected p-values
p_colors <- rep(NA, ncol(p_F1_adj))
p_type <- rep(NA, ncol(p_F1_adj))


# Loop over each column of p_F1_adj
for (i in 1:ncol(p_F1_adj)) {
    # Get the p-values for this column
    p_vals <- p_F1_adj[, i]

    # Set default values for color and line type
    color <- "red"  # Default color for mixed significant/non-significant values
    lty <- "dotdash"   # Default line type for mixed values

    # Check if all p-values in this column are above the threshold (not significant)
    if (all(p_vals > p_level)) {
        color <- "gray"
        lty <- "longdash"
    }
    # Check if all p-values in this column are below the threshold (significant)
    else if (all(p_vals <= p_level)) {
        color <- "blue"
        lty <- "solid"
    }

    # Store the color and line type for this column
    p_colors[i] <- color
    p_type[i] <- lty
}


##### 2) PLOT THE TANGLEGRAM #####
# this minimizes "link crossing" so the tanglegram is less tangled

ph1 <- trimmed.h.tree
ph2 <- trimmed.p.tree
assoc_links <- link_names
cophy <- cophylo(ph1, ph2, assoc=assoc_links, rotate=T)

# set as PDF and size of printing area, adjust as needed
pdf("fake-tanglegram.pdf", height = 30, width = 45)
cophyloplot_min<-plot(cophy, link.lty=p_type, link.col=p_colors, link.lwd=8, fsize=3.5, lwd = 6, pts=F)
dev.off()



# NOTE: If label names are too long.... do this....

# EMMANUEL PARADIS AT http://grokbase.com/t/r/r-sig-phylo/12c3djkx3h/modification-of-figures-produced-using-cophyloplot
# FIRST MAKE A COPY OF THE INTERNAL FUNCTION
# titi<-plotCophylo2
# then do fix(titi) and modify what you want (here look for calls to
# text(... cex = ...), save and close. Do fix(cophyloplot) and change
# "plotCophylo2(..." by "titi(...", save and close.
# fix(titi) #
# if (show.tip.label) {
#   text(a[1:N.tip.x, ], cex = 0.5, font = font, pos = 4, labels = x$tip.label)
#   text(b2[1:N.tip.y, ], cex = 0.5, font = font, pos = 2,
#        labels = y$tip.label)
# }
#
# titi
# fix(cophyloplot) # Do fix(cophyloplot) and change "plotCophylo2(..." by "titi(...", save and close.


