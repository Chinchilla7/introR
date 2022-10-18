#PART 1
#1.1
getwd()

#R graphics engine version 15 is not supported by this version of RStudio. 
#The Plots tab will be disabled until a newer version of RStudio is installed. 
#MacOS version is 10.13.6, which could only install rstudio version 1.3.1093

3+5

x <- 3
y <- 5
y
x + y
number <- x + y

#exercise
x <- 5
y <- 10
number <- x + y
#
#1.2
# Create a numeric vector and store the vector as a variable called 'glengths'
glengths <- c(4.6, 3000, 50000)
glengths

# Create a character vector and store the vector as a variable called 'species'
species <- c("ecoli", "human", "corn")
species

# Forget to put quotes around corn
species <- c("ecoli", "human", corn)

#exercise 
combined <- c(glengths, species)
#
# Create a character vector and store the vector as a variable called 'expression'
expression <- c("low", "high", "medium", "high", "low", "medium", "high")

# Turn 'expression' vector into a factor
expression <- factor(expression)

#exercise 
samplegroup <- c("CTL", "CTL", "CTL", "KO", "KO", "KO", "OE", "OE", "OE")
samplegroup <- factor(samplegroup)
#
# Create a data frame and store it as a variable called 'df'
df <- data.frame(species, glengths)
df

#exercise 
titles <- c("Catch-22", "Pride and Prejudice", "Nineteen Eighty Four")
pages <- c(453, 432, 328)
favorite_books <- data.frame(titles, pages)
#
list1 <- list(species, df, number)

#exercise 
list2 <- list(species, glengths, number)
#
#1.3
function_name(input)
getwd()

glengths <- c(glengths, 90) # adding at the end	
glengths <- c(30, glengths) # adding at the beginning

sqrt(81)

sqrt(glengths)

round(3.14159)

#exercise 
mean(glengths)
test <- c(1, NA, 2, 3, NA, 4)
mean(test, na.rm=TRUE)
sort(glengths, decreasing = TRUE)
#

?round
args(round)


round(3.14159, digits=2)
round(3.14159, 2)

name_of_function <- function(argument1, argument2) {
    statements or code that does something
    return(something)
}

square_it <- function(x) {
    square <- x * x
    return(square)
}
square_it(5)

#exercise 
multiply_it <- function(x,y) {
  product <- x * y
  return(product)
}
#
#1.4
?read.csv

metadata <- read.csv(file="data/mouse_exp_design.txt")

#exercise 
proj_summary <- read.table(file = "data/project-summary.txt", header = TRUE, row.names = 1)
proj_summary
#
metadata
head(metadata)

#exercise 
class(glengths)
#"numeric"
class(metadata)
#"data.frame"

summary(proj_summary)
#Median :0.005345 
length(samplegroup)
#9
dim(proj_summary)
#9 8
str(rownames(metadata))
#chr [1:12] "sample1" "sample2" "sample3" "sample4" "sample5" "sample6" "sample7" "sample8" "sample9" ...
length(colnames(proj_summary))
#8
#
#Practice
temp_conv <- function(temp_f) {
  temp_c = (temp_f - 32) * 5 / 9
  temp_k = temp_c + 273.15
  return (temp_k)
}
round(temp_conv(70), digits = 1)
#
#PART 2
#2.1
age <- c(15, 22, 45, 52, 73, 81)
age[5]
age[-5]
age[c(3,5,6)]   ## nested
# OR
## create a vector first then select
idx <- c(3,5,6) # create vector of the elements of interest
age[idx]

age[1:4]

#exercise 1
alphabets <- c("C", "D", "X", "L", "F")
alphabets[c(1,2,5)]
alphabets[-3]
alphabets[5:1]
#

age
age > 50

age > 50 | age < 18
age
age[age > 50 | age < 18]  ## nested
# OR
## create a vector first then select
idx <- age > 50 | age < 18
age[idx]

which(age > 50 | age < 18)
age[which(age > 50 | age < 18)]  ## nested
# OR
## create a vector first then select
idx_num <- which(age > 50 | age < 18)
age[idx_num]

expression[expression == "high"]    ## This will only return those elements in the factor equal to "high"

#exercise
samplegroup != "KO"
samplegroup[samplegroup != "KO"]
#

expression
str(expression)

expression <- factor(expression, levels=c("low", "medium", "high"))     # you can re-factor a factor 
str(expression)

#exercise
samplegroup
str(samplegroup)
samplegroup <- factor(samplegroup, levels=c("KO","CTL","OE"))
str(samplegroup)
#

#2.2
sessionInfo() #Print version information about R, the OS and attached or loaded packages
# OR
search() #Gives a list of attached packages

install.packages("ggplot2")
library(ggplot2)
search()

#exercise
install.packages("tidyverse")

#2.3
# Extract value 'Wt'
metadata[1, 1]
# Extract value '1'
metadata[1, 3] 
# Extract third row
metadata[3, ] 
# Extract third column
metadata[ , 3]  
# Extract third column as a data frame
metadata[ , 3, drop = FALSE] 
# Dataframe containing first two columns
metadata[ , 1:2] 
# Data frame containing first, third and sixth rows
metadata[c(1,3,6), ] 
# Extract the celltype column for the first three samples
metadata[c("sample1", "sample2", "sample3") , "celltype"] 

# Check column names of metadata data frame
colnames(metadata)
# Check row names of metadata data frame
rownames(metadata)

# Extract the genotype column
metadata$genotype 
# Extract the first five values/elements of the genotype column
metadata$genotype[1:5]

#exercise
metadata[c("sample2", "sample8"), c("genotype", "replicate")]
metadata$replicate[c(4,9)]
metadata[, "replicate", drop = FALSE]
#

metadata$celltype == "typeA"
logical_idx <- metadata$celltype == "typeA"
metadata[logical_idx, ]
which(metadata$celltype == "typeA")

idx <- which(metadata$celltype == "typeA")
metadata[idx, ]
idx <- which(metadata$replicate > 1)

which(metadata$replicate > 1)
metadata[idx, ]

metadata[which(metadata$replicate > 1), ]
sub_meta <- metadata[which(metadata$replicate > 1), ]

#exercise
idx <- which(metadata$genotype=="KO")
metadata[idx, ]
#

list1[[2]]

comp2 <- list1[[2]]
class(comp2)

list1[[1]]
list1[[1]][1]

#exercise
random <- list(metadata, age, list1, samplegroup, number)
random[4]
#

names(list1) 
# Name components of the list
names(list1) <- c("species", "df", "number")
names(list1)

# Extract 'df' component
list1$df

#exercise
names(random) <- c("metadata", "age", "list1", "samplegroup", "number")
random$age

#2.4
rpkm_data <- read.csv("data/counts.rpkm.csv")
head(rpkm_data)
ncol(rpkm_data)
nrow(metadata)

vector1 %in% vector2

A <- c(1,3,5,7,9,11)   # odd numbers
B <- c(2,4,6,8,10,12)  # even numbers

# test to see if each of the elements of A is in B	
A %in% B

A <- c(1,3,5,7,9,11)   # odd numbers
B <- c(2,4,6,8,1,5)  # add some odd numbers in 
A %in% B
A[intersection]
any(A %in% B)
all(A %in% B)

#exercise
B %in% A
B[B %in% A]
#

A <- c(10,20,30,40,50)
B <- c(50,40,30,20,10)  # same numbers but backwards 

# test to see if each element of A is in B
A %in% B

# test to see if each element of A is in the same position in B
A == B

# use all() to check if they are a perfect match
all(A == B)

x <- rownames(metadata)
y <- colnames(rpkm_data)
all(x %in% y)
all(rownames(metadata) %in% colnames(rpkm_data))
x == y
all(x == y)

#exercise
important_genes<- c("ENSMUSG00000083700", "ENSMUSG00000080990", "ENSMUSG00000065619", "ENSMUSG00000047945", "ENSMUSG00000081010", "ENSMUSG00000030970")
important_genes %in% rownames(rpkm_data)
ris <- rownames(rpkm_data) %in% important_genes
ext <- rpkm_data[ris, ]
ris2 <- which(rownames(rpkm_data) %in% important_genes)
ext2 <- rpkm_data[ris2, ]
ext3 <- rpkm_data[important_genes, ]

#2.5
teaching_team <- c("Jihe", "Mary", "Meeta", "Radhika")
# Extracting values from a vector
teaching_team[c(2, 4)] 
teaching_team
teaching_team[c(4, 2)] 
# Extracting all values and reordering them
teaching_team[c(4, 2, 1, 3)]
# Saving the results to a variable
reorder_teach <- teaching_team[c(4, 2, 1, 3)] 

#exercise
first <- c("A","B","C","D","E")
second <- c("B","D","E","A","C")  # same letters but different order
second[c(4,1,5,2,3)]
#

match(first,second)
reorder_idx <- match(first,second) 
second[reorder_idx]
second_reordered <- second[reorder_idx]  

first <- c("A","B","C","D","E")
second <- c("D","B","A")  # remove values
match(first,second)
second[match(first, second)]

# Check row names of the metadata
rownames(metadata)

# Check the column names of the counts data
colnames(rpkm_data)

genomic_idx <- match(rownames(metadata), colnames(rpkm_data))
genomic_idx
# Reorder the counts data frame to have the sample names in the same order as the metadata data frame
rpkm_ordered  <- rpkm_data[ , genomic_idx]
# View the reordered counts
View(rpkm_ordered)
all(rownames(metadata) == colnames(rpkm_ordered))

#exercise
subset_rpkm <- rpkm_ordered[ , -c(2,9)]

mc <- match(colnames(subset_rpkm), rownames(metadata))
metadata[mc, ]
#

#2.6
mean(rpkm_ordered$sample1)
library(purrr)  # Load the purrr

samplemeans <- map_dbl(rpkm_ordered, mean) 

# Named vectors have a name assigned to each element instead of just referring to them as indices ([1], [2] and so on)
samplemeans

# Check length of the vector before adding it to the data frame
length(samplemeans)

age_in_days <- c(40, 32, 38, 35, 41, 32, 34, 26, 28, 28, 30, 32)    	

# Add the new vector as the last column to the new_metadata dataframe
new_metadata <- data.frame(metadata, samplemeans, age_in_days) 

# Take a look at the new_metadata object
View(new_metadata)

#practice
animals <- read.csv("data/animals.txt")
class(animals)
nrow(animals)
ncol(animals)

animals[1,1]
animals[c(2,5),]
animals[which(animals$speed > 50), "color", drop =F]
animals$color[which(animals$color == "Grey")] <- "Gray"
animals_list <- list(animals$speed, animals$color)
names(animals_list) <- colnames(animals)

ctrl_samples <- data.frame(row.names = c("sample3", "sample10", "sample8", "sample4", "sample15"), date = c("01/13/2018", "03/15/2018", "01/13/2018", "09/20/2018","03/15/2018"))
length(which(rownames(ctrl_samples) %in% rownames(proj_summary)))
proj_summary_ctrl <- proj_summary[which(rownames(proj_summary) %in% rownames(ctrl_samples)),]
m <- match(rownames(proj_summary_ctrl), rownames(ctrl_samples))
proj_summary_ctrl <- cbind(proj_summary_ctrl, batch=ctrl_samples[m,])

proj_summary_noctl <- proj_summary$treatment[which(proj_summary$treatment != "control")]
keep <- map_lgl(proj_summary_noctl, is.numeric)
proj_summary_noctl <- proj_summary_noctl[,keep]

#EXERCISE 3
#3.1
## load the new_metadata data frame into your environment from a .RData object
load("data/new_metadata.RData")

# this data frame should have 12 rows and 5 columns
View(new_metadata)

library(ggplot2)

ggplot(new_metadata)
ggplot(new_metadata) +
  geom_point()
#Error in `check_required_aesthetics()`:
#! geom_point requires the following missing aesthetics: x and y
#Run `rlang::last_error()` to see where the error occurred.

ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans))

ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype)) 
  
ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,
                 shape=celltype)) 

ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,
                 shape=celltype), size=2.25) 
ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,
                 shape=celltype), size=3.0) +
  theme_bw() 

ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,
                 shape=celltype), size=2.25) +
  theme_bw() +
  theme(axis.title = element_text(size=rel(1.5)))

#exercise
ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,
                 shape=celltype), size=2.25) +
  theme_bw() +
  theme(axis.title = element_text(size=rel(1.5))) +
  xlab("Age (days)") +
  ylab("Mean expression") + 
  ggtitle("Exercise plot") +
  theme(plot.title=element_text(hjust=0.5)) #title is centered, ~20 theme layers can be added?

#3.2
name_of_function <- function(arguments) {
    statements or code that does something
}

theme_bw() +
theme(axis.title=element_text(size=rel(1.5))) +
theme(plot.title=element_text(size=rel(1.5), hjust=0.5))
#OR
theme_bw() +
theme(axis.title=element_text(size=rel(1.5)), plot.title=element_text(size=rel(1.5), hjust=0.5))

personal_theme <- function(){
  theme_bw() +
  theme(axis.title=element_text(size=rel(1.5))) +
  theme(plot.title=element_text(size=rel(1.5), hjust=0.5))
}

ggplot(new_metadata) +
  geom_point(aes(x=age_in_days, y=samplemeans, color=genotype, shape=celltype), size=rel(3.0)) +
  xlab("Age (days)") +
  ylab("Mean expression") +
  ggtitle("Expression with Age") +
  personal_theme()

  #3.3
  ggplot(new_metadata) +
  geom_boxplot(aes(x=genotype, y=samplemeans, fill=celltype)) +
  ggtitle("Genotype differences in average gene expression") +
  xlab("Genotype") +
  ylab("Mean expression")  +
    theme_bw() +
    theme(axis.title = element_text(size = rel(1.25))) +
    theme(plot.title=element_text(hjust = 0.5, size = rel(1.5)))

new_metadata$genotype <- factor(new_metadata$genotype, levels = c("Wt", "KO"))
ggplot(new_metadata) +
  geom_boxplot(aes(x=new_metadata$genotype, y=samplemeans, fill=celltype)) +
  ggtitle("Genotype differences in average gene expression") +
  xlab("Genotype") +
  ylab("Mean expression")  +
    theme_bw() +
    theme(axis.title = element_text(size = rel(1.25))) +
    theme(plot.title=element_text(hjust = 0.5, size = rel(1.5)))

ggplot(new_metadata) +
  geom_boxplot(aes(x=genotype, y=samplemeans, fill=celltype)) +
  ggtitle("Genotype differences in average gene expression") +
  xlab("Genotype") +
  ylab("Mean expression")  +
  scale_color_manual(values=c("purple","orange")) + #no change observed
    theme_bw() +
    theme(axis.title = element_text(size = rel(1.25))) +
    theme(plot.title=element_text(hjust = 0.5, size = rel(1.5)))

ggplot(new_metadata) +
  geom_boxplot(aes(x=genotype, y=samplemeans, fill=celltype)) +
  ggtitle("Genotype differences in average gene expression") +
  xlab("Genotype") +
  ylab("Mean expression")  +
  scale_fill_manual(values=c("purple","orange")) + #change observed, different from scale_color_manual
    theme_bw() +
    theme(axis.title = element_text(size = rel(1.25))) +
    theme(plot.title=element_text(hjust = 0.5, size = rel(1.5)))
#difference is scale_color_manual() works with scatter plot, scale_fill_#manual() works with box plot 

#3.4
# Save a data frame to file
write.csv(sub_meta, file="data/subset_meta.csv")
?write.csv

# Save a vector to file
write(glengths, file="data/genome_lengths.txt")
?write
# Save a vector to file as a single column
write(glengths, file="data/genome_lengths.txt", ncolumns = 1)

## Open device for writing
pdf("figures/scatterplot.pdf")
## Make a plot which will be written to the open device, in this case the temp file created by pdf()/png()
ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,
                 shape=celltype), size=rel(3.0)) 
## Closing the device is essential to save the temporary file created by pdf()/png()
dev.off()

#3.5
 sessionInfo()  #This time it is not interchangeable with search()
#exercise
 # Create vector of work days
work_days <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")
# Create a function to round the output of the sum function
round_the_sum <- function(x){
  return(round(sum(x)))
}
# Create a function to add together three numbers
add_numbers <- function(x,y,z){
  sum(x,y,z)
}

add_numbers(5,9,6)

 #Error: package or namespace load failed for 'Seurat' in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]): there is no package called 'multtest'
BiocManager::install('multtest')
install.packages('Seurat')

#You would like to ask for help on an online forum. To do this you want the users of the forum to reproduce your problem, so you want to provide them as much relevant information and data as possible.
#You want to provide them with the list of packages that you currently have loaded, the version of R, your OS and package versions. Use the appropriate function(s) to obtain this information.
sessionInfo()
#You want to also provide a small data frame that reproduces the error (if working with a large data frame, you’ll need to subset it down to something small). For this exercse use the data frame df, and save it as an RData object called df.RData.
save(df, file = "data/df.RData")
#What code should the people looking at your help request should use to read in df.RData?
load(file="data/df.RData")

#3.6
library(tidyverse)
#exercise
random_numbers <- c(81, 90, 65, 43, 71, 29)
random_numbers %>% mean()
random_numbers %>% 
  mean() %>% 
  round(digits = 3)
#

# Read in the functional analysis results
functional_GO_results <- read_delim(file = "data/gprofiler_results_Mov10oe.txt", delim = "\t" )

# Take a look at the results
functional_GO_results

# Return only GO biological processes
bp_oe <- functional_GO_results %>%
  filter(domain == "BP")
View(bp_oe)

#exercise
bp_oe <- bp_oe %>% 
  filter(relative.depth > 4)
#

# Selecting columns to keep
bp_oe <- bp_oe %>%
  select(term.id, term.name, p.value, query.size, term.size, overlap.size, intersection)

View(bp_oe)

# Order by adjusted p-value ascending
bp_oe <- bp_oe %>%
  arrange(p.value)

# Provide better names for columns
bp_oe <- bp_oe %>% 
  dplyr::rename(GO_id = term.id, 
                GO_term = term.name)
#exercise
bp_oe <- bp_oe %>% 
  dplyr::rename(genes = intersection)
#
# Create gene ratio column based on other columns in dataset
bp_oe <- bp_oe %>%
  mutate(gene_ratio = overlap.size / query.size)

#exercise
bp_oe <- bp_oe %>% 
  mutate(term_percent = overlap.size / term.size)
#

#Practice
animals_tb <- animals %>%
        rownames_to_column(var = "animal_names") %>%
        as_tibble()

ggplot(animals_tb) +
        geom_point(aes(x = animal_names, y = speed), color = "purple") +
        theme_bw() +
        ggtitle("Speed Comparisons Between Animals") + 
        ylab("Speed (km/h)") +
        xlab("Animal") +
        theme(plot.title=element_text(hjust=0.5))

names_ordered_by_speed <- animals_tb %>% arrange(speed) %>% pull(animal_names)

animals_tb$animal_names <- factor(animals_tb$animal_names, 
                                  levels = names_ordered_by_speed)

ggplot(animals_tb) +
        geom_point(aes(x = animal_names, y = speed), color = "purple") +
        theme_bw() +
        ggtitle("Speed Comparisons Between Animals") + 
        ylab("Speed (km/h)") +
        xlab("Animal") +
        theme(plot.title=element_text(hjust=0.5))

pdf("results/animals_by_speed_scatterplot.pdf")
ggplot(animals_tb) +
        geom_point(aes(x = animal_names, y = speed), color = "purple") +
        theme_bw() +
        ggtitle("Speed Comparisons Between Animals") + 
        ylab("Speed (km/h)") +
        xlab("Animal") +
        theme(plot.title=element_text(hjust=0.5))

dev.off()

animals_gray_tan <- animals_tb %>% 
  filter(color == "Gray" | color == "Tan") %>%
  arrange(speed)
write.csv(animals_gray_tan,
          file = "results/animals_tb_ordered.csv",
          quote = FALSE)

#Practice 
sex <- c("M", "F","M", "F","M", "F","M", "F","M", "F","M", "F")
stage <- c("I", "II", "II","I", "II", "II","I", "II", "II","I", "II", "II")
treatment <- c("A","A","A","A", "B",  "B", "B", "B","P","P","P","P")
myc <- c("2343","457","4593","9035","3450","3524","958","1053","8674","3424","463","5105")

meta <- data.frame(sex, stage, treatment, myc)
rownames(meta) <- paste("sample", 1:12, sep="")

meta[ , c(1,3)]
meta[c(5,7,9,10), 3]
filter(meta, treatment == "P")
filter(meta, myc > 5000) %>% select(stage, treatment)
meta[, -3]
meta[-7:-9, ]
meta [1:6, ]
pre_treatment <- c(T, F, F, F, T, T, F, T, F, F, T, T)

cbind(pre_treatment, meta)
colnames(meta) <- c("A", "B", "C", "D")

list_hw <- list(glengths, df, number)
list_hw[[2]]
list_hw[[1]][3]
list_hw[[2]][2, 1]
names(list_hw) <- c("genome_lengths","genomes","record")

meta <- read.delim("Mov10_full_meta.txt", sep="\t", row.names=1)
data <- read.delim("normalized_counts.txt", sep="\t", row.names=1)

expression <- data["MOV10", ]
class(expression)
expression <- as.numeric(expression)
class(expression)

df <- data.frame(meta, expression)

ggplot(df) +
  geom_jitter(aes(x= sampletype, y= expression, color = sampletype)) +
  theme_bw() +
  ggtitle("Expression of MOV10") +
  xlab(NULL) +
  ylab("Normalized counts") +
  theme(legend.position = "none",
        plot.title=element_text(hjust=0.5, size=rel(1.5)),
        axis.text=element_text(size=rel(1.25)),
        axis.title=element_text(size=rel(1.5)),
        axis.text.x=element_text(angle=45, hjust=1))















