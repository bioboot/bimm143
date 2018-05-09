# Class 6 Bioinformatics

# Functions

add <- function(x, y=1) {
  # Sum x and y inputs
  x + y
}

# My second function
rescale <- function(x) {
  rng <-range(x, na.rm=TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

# Test on a small example where you know the answer
rescale(1:10)

# How would you get your function to work here…
rescale( c(1,2,NA,3,10) )

# What should your function do here?
#rescale( c(1,10,"string") )

## Next generation rescaling - WOW!
rescale2 <- function(x, na.rm=TRUE, plot=FALSE) {
  if(na.rm) {
    rng <-range(x, na.rm=na.rm)
  } else {
    rng <-range(x)
  }
  print("Hello")
  
  answer <- (x - rng[1]) / (rng[2] - rng[1])

  print("is it me you are looking for?")
  
  if(plot) { 
    plot(answer, typ="b", lwd=4) 
  }
  print("I can see it in ...")
}

rescale3 <- function(x, na.rm=TRUE, plot=FALSE) {
  if(na.rm) {
    rng <-range(x, na.rm=na.rm)
  } else {
    rng <-range(x)
  }
  print("Hello")
  
  answer <- (x - rng[1]) / (rng[2] - rng[1])
  
  print("is it me you are looking for?")
  
  if(plot) { 
    plot(answer, typ="b", lwd=4) 
  }
  print("I can see it in ...")

  return(answer)
  
}


# Lets test rescale2
rescale2( c(1,2,NA,3,10) )

## Hands on Lab sheet for class 6! ----

# Section 2B.
## One time only package install
#install.packages("bio3d")

# Use the bio3d package!
library(bio3d)
s1 <- read.pdb("4AKE")  # kinase with drug
s1

s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")

## Lets start to make this into a function
x <- "4AKE"
s <- read.pdb(x)
