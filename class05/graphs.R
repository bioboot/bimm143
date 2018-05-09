#' ---
#' title: "Bioinformatics Class 5"
#' author: "Barry Grant"
#' output:
#'   html_document:
#'     code_folding: hide
#'     keep_md: true
#' ---


#' # Class 5 graphs
#' We can enter text and comments using the hash/pound symbol (#) at the start of a line

# boxplot
boxplot( rnorm(1000,0) )

hist( rnorm(1000,0) )

summary( rnorm(1000,0) )

# Flip my boxplot
boxplot( rnorm(1000,0), horizontal = TRUE )


# Read first data file.
baby <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)

# Point and line plot
plot(baby, type="b", pch=15, cex=1.5, main="Some main title" )


# 1B ----
# Read another tab seperated file
feat <- read.table("bimm143_05_rstats/feature_counts.txt", 
                   sep="\t", header = TRUE)

# Make a barplot
par(mar=c(5,11,4,2))
barplot( feat[,2], horiz=TRUE,  
         names.arg=feat[,1], las=1 )


## Section 2. ----
file <- "bimm143_05_rstats/male_female_counts.txt"
mf_counts <- read.table(file, sep="\t", header=TRUE)

# Make a barplot
barplot(mf_counts$Count, col=rainbow( nrow(mf_counts)) )

## color by male and female
mycols=c("blue2","red2")
barplot(mf_counts$Count, col=mycols)

## Section2B. ----
file_2b <- "bimm143_05_rstats/up_down_expression.txt"
up.down <- read.delim(file_2b)

plot(up.down$Condition1, up.down$Condition2, col=up.down$State)


# How many genes are up and down regulated?
table(up.down$State)


## Plot of mtcars
par(mar=c(2,11,2,2))
ord <- order(mtcars$mpg)
barplot(mtcars$mpg[ord], names.arg = rownames(mtcars)[ord], 
        horiz = TRUE, las=2, col=rainbow(nrow(mtcars)))

## Dot plots are a reasonable substitute for bar plots.
dotchart(mtcars$mpg[ord], labels = rownames(mtcars)[ord], xlab="MPG")


#' ### Add a session info
#' It is a good practice to add a session info at the end of your document. It will increase reproducibility and costs only one line of code
 
sessionInfo()




## Section 2C.

map.colors <- function (value,high.low,palette) {
  proportion <- ((value-high.low[1])/(high.low[2]-high.low[1]))
  index <- round ((length(palette)-1)*proportion)+1
  return (palette[index])
}

map.colors2 <- function (value,high.low,palette) {
  proportion <- ((value-high.low[1])/(high.low[2]-high.low[1]))
  index <- round ((length(palette)-1)*proportion)+1
  return (palette[index])
}
