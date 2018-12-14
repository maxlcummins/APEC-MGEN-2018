library(dplyr)
library(ggtree)
library(phytools)
library(readr)
library(magrittr)

#read genotype data
data <- read_csv("~/Dropbox/Doctorate/Manuscripts/BigChook/Figures/BigChook3_simple.csv")

#read APEC metadata
APEC_Origins_BigChook <- read_csv("~/Dropbox/Doctorate/Manuscripts/BigChook/Figures/ST117_Origins.csv")

#clean column names
colnames(APEC_Origins_BigChook) <- c('name','date_submittted', 'origin')

#read tree file
tree <- read.tree(file = "~/Dropbox/Doctorate/Manuscripts/BigChook/Figures/CC117_clean_noout.tree")

#Here a column called 'phyloname' is generated through joining strain names, sequence type and serotype
data$phyloname <- paste(data$name, "_", data$OH_type, sep = "")

#remove non-CC117 strains
data <- filter(data, ST %in% c(117,4045))

#Here columns are reordered
data <- data %>% select(name, phyloname, ST, phylogroup, O_type, H_type, OH_type, everything())

#join the genotype data with the metadata
genotype_metadata <- left_join(data, APEC_Origins_BigChook)

#Trim everything after an underscore in tip names in the tree to remove _R1.fastq.gz.
#NOTE - if you modify these script ensure your sample names don't have an underscore in them!
tree$tip.label <- gsub(pattern = "_.*", "", tree$tip.label)

#root midpoint
tree <- midpoint.root(tree)

#midpoint root the tree...
#tree <- root(tree, outgroup = "MRSN346647")

#Removes the outgroup from the tree - it complicates the viewing of our branches
#tree <- drop.tip(tree, "MRSN346647")

#change ST column to have 'ST' preceding the sequence type
gsub(pattern = "(.*)", "ST\\1", x = data$ST) -> data$ST

#trim phyloname in tree tip labels to just sample names to allow matching later on
gsub(pattern = "_ST.*", "", x = tree$tip.label) -> tree$tip.label

#pull out the tip labels from the tree so we can filter in the next step
as.data.frame(tree$tip.label) -> tips

#change the column name of tips to 'name'
colnames(tips) <- 'name'

#join a list of phylonames from the ARIBAlord output to the list of tip labels in the tree
newnames <- left_join(tips,data[,c('name','phyloname')])

#Rename the reference tip
#NOTE: This name change is index based, be sure that you are changing the right name!
newnames[7,2] <- "2009C-3133_O119:H4"

#replace tip labels in the tree with phylonames
newnames$phyloname -> tree$tip.label

#creates a list containing only the samples that appear in the tree
#we will use this  in a few lines in combination with semi-join to remove ARIBAlord rows that are missing from the tree
newnames$phyloname -> tips$name

#rename this column to phyloname
colnames(tips) <- 'phyloname'

#Reduce the dataframe with phylogenetic and genotypic data to just phylonames and genetic data
genotype <- data %>% select(name, phyloname, matches("^(r|i|v|r|p)_"))

#remove ariba data from samples not included in the tree
data <- semi_join(data, tips)
genotype <- semi_join(genotype, tips)

#Order columns by their names
data <- data[order(data$name),]
genotype <- genotype[order(genotype$name),]

#Here columns are reordered
data <- data %>% select(phyloname, name, everything())

#Provides the names and colors that will be used to annotate the dates.
state <- c("QLD", "NSW", "VIC", "WA", NA)
state_cols <- c("#754000",  #brown, QLD
                "#0015a7",  #blue, NSW
                "#00ad01",  #green, VIC
                "#b6004e",  #red, WA
                "black")    #black, Unknown

#Define colors for gene-type dependent coloring of gene hits
colorgenotype <- c("N" = "white", "I" = "#8dd3c7", "R" = "#bebada", "V" = "#fb8072", "P" = "#80b1d3", "REF" = "black")

#reorder the genotype and genotype_metadata tables so that phyloname is in the first column
#This is needed by ggtree to bind the metadata to the tree and color the tiplabels
data <- data %>% select(phyloname, name, everything())
genotype_metadata <- genotype_metadata %>% select(phyloname, name, everything())

#generate a tree that colors the tip labels based on the phylogroups they are in:
ptree <- ggtree(tree) %<+%
  genotype_metadata +
  geom_tiplab(size = 2.5, hjust = -.125, aes(color = origin), align = TRUE) +
  geom_text2(aes(subset = !isTip, label=label), size = 1, hjust=-.125) +
  scale_color_manual(breaks = c(state), values = c(state_cols), na.value = "black")


#Save rownames of phyloname so we can reapply them after they are lost in subsequent steps
genotype$phyloname -> namesave

#splits the gene hit types to separate dfs to be processed to allow different colours for different gene types
i <- genotype %>% select(starts_with("i_"))
r <- genotype %>% select(starts_with("r_"))
v <- genotype %>% select(starts_with("v_"))
p <- genotype %>% select(starts_with("p_"))

#remove any remaining columns with no hits
i <- i[,colSums(i) > 0]
r <- r[,colSums(r) > 0]
v <- v[,colSums(v) > 0]
p <- p[,colSums(p) > 0]

#change gene hits to a different number based on gene type
i[i > 0] <- "I"
r[r > 0] <- "R"
v[v > 0] <- "V"
p[p > 0] <- "P"

#combines processed sub-dfs
genotype <- cbind(i,r,v,p)

#orders the columns based on their names
genotype <- genotype[ , order(names(genotype))]

# #trim off the r_ etc from start of each gene hit
colnames(genotype) <- gsub(pattern = "^[r,i,v,p]_","",colnames(genotype), perl = TRUE)

#assign rownames
rownames(genotype) <- namesave

#duplicate genotype to a new dataframe called genotype_binary
genotype_binary <- genotype

#convert letters indicating presence of a gene of a certain type back to 1
#we will use 'genotype_binary' to allow the generation of sums for each gene
genotype_binary[genotype_binary == "I"] <- 1
genotype_binary[genotype_binary == "R"] <- 1
genotype_binary[genotype_binary == "V"] <- 1
genotype_binary[genotype_binary == "P"] <- 1

#Define colors for gene-type dependent coloring of gene hits
colorgenotype <- c("N" = "white", "I" = "#8dd3c7", "R" = "#bebada", "V" = "#fb8072", "P" = "#80b1d3", "REF" = "black")

#trim the gene name prefixes indicating their gene type
colnames(genotype_binary) <- gsub(pattern = "^[r,i,v,p]_","",colnames(genotype_binary), perl = TRUE)

#creates a df to be used in generation of new colnames with gene sums
old_new_colnames <- rbind(colnames(genotype),colnames(genotype))

#Sum each hit for each gene
genesums <- as.data.frame(colSums(data.matrix(genotype_binary)))

#rename genesum column to 'sum'
colnames(genesums) <- 'sum'

#paste together the new colnames and assign to our df with old and new names
genesums$rowsumcat <- paste(rownames(genesums), " (", genesums$sum , "/", nrow(genotype),")", sep="")

#duplicates our genotype table for a second figure containing gene sums
genotype_wsums <- genotype

#assigns new colnames to hit_table4
colnames(genotype_wsums) <- genesums$rowsumcat

#We need to create a row for the reference genome so that we can grey out the genotype box
#This row must have the correct rowname
REF_rep <- rep(NA, ncol(genotype_wsums))
genotype_wsums <- rbind(genotype_wsums, REF_rep)
rowname_replace <- c(rownames(genotype_wsums))
rowname_replace <- c(rowname_replace[1:23],"2009C−3133_ST117_O119:H4")
rownames(genotype_wsums) <- rowname_replace

#change cases where a gene was not detected to the representation of such with an N isntead of a zero
genotype_wsums[genotype_wsums == 0] <- "N"

#Generates the tree
a <- gheatmap(p = ptree, data = genotype_wsums,
              colnames_offset_y = -0.4,
              font.size = 1.5,
              hjust = 0,
              colnames_position = "top",
              colnames_angle = 90,
              width = 3,
              color = "black",
              offset = 0.25) + 
  scale_fill_manual(values = colorgenotype, na.value = 'grey') +
  theme(legend.position = "none")

#Plots the tree
a

