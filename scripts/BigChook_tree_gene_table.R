library(dplyr)
library(ggtree)
library(phytools)
library(readr)
library(magrittr)

#read genotype data
data <- read_csv("~/Dropbox/Doctorate/Manuscripts/BigChook/Figures/BigChook3_simple.csv")

#read tree file
tree <- read.tree(file = "~/Dropbox/Doctorate/Manuscripts/BigChook/Figures/BigChook3.tree")

#Here tip names are changed to reflect their new strain names for publishing
tree$tip.label <- gsub(pattern = "A", replacement = "AVC", tree$tip.label)

#Here a column called 'phyloname' is generated through joining strain names, sequence type and serotype
data$phyloname <- paste(data$name, "_ST", data$ST, "_", data$OH_type, sep = "")

#Here columns are reordered
data <- data %>% select(name, phyloname, ST, phylogroup, O_type, H_type, OH_type, everything())

#Trim everything after an underscore in tip names in the tree to remove _R1.fastq.gz.
#NOTE - if you modify these script ensure your sample names don't have an underscore in them!
tree$tip.label <- gsub(pattern = "_.*", "", tree$tip.label)

#midpoint root the tree...
tree <- midpoint.root(tree = tree)

#change STs to have 'ST'
gsub(pattern = "(.*)", "ST\\1", x = data$ST) -> data$ST

#trim phyloname in tree tip labels to just sample names to allow matching later on
gsub(pattern = "_ST.*", "", x = tree$tip.label) -> tree$tip.label

#generate and sort a list of STs
STs <- sort(unique(data$ST))

#generate a list of sample colors equal to length of ST list
ST_cols <- rainbow(n=length(STs))

#define colours associated with phylogroups
c("A","B1","B2","D") -> phylogroups
c("#cc545e",
  "#64a860",
  "#9970c1",
  "#b98d3e") -> phylocols

#pull out the tip labels from the tree so we can filter in the next step
as.data.frame(tree$tip.label) -> tips

#change the column name of tips to 'name'
colnames(tips) <- 'name'

#join a list of phylonames from the ARIBAlord output to the list of tip labels in the tree
newnames <- left_join(tips,data[,c('name','phyloname')])

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

#generate a tree that colors the tip labels based on the phylogroups they are in:
#red - A
#green - B1
#purple - B2
#gold/mustard - D
BC_tree <- ggtree(tree, branch.length = "none") %<+%
  data +
  geom_tiplab(offset = 0, size = 1.6, aes(color = phylogroup)) +
  scale_color_manual(breaks = phylogroups, values = phylocols)

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
colorgenotype <- c("N" = "white", "I" = "#8dd3c7", "R" = "#bebada", "V" = "#fb8072", "P" = "#80b1d3")

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

#plot heatmap alongside BC_tree
#color columns based on the gene-types with which they are associated:
#teal - instertion/integrase
#blue - plasmid
#purple - resistance
#red - virulence
# Note: This figure will show the sum of genes in the column names.
#       You can replace genotype_wsums to genotype and it will remove the sums from the column names
c <- gheatmap(p = BC_tree, data = genotype_wsums,
              #colnames_offset_y = .0000000001,
              font.size = 1,
              #colnames_offset_x = 0.1,
              colnames_position = "top",
              hjust = 0,
              colnames_angle = 90,
              width = 12,
              color = "black",
              offset = 28) +
  scale_fill_manual(values = colorgenotype) +
  theme(legend.position = "none")
