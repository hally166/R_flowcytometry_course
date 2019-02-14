# Wellcome Sanger Institute Cytometry Core Facility: 
# Introduction to R for Flow Cytometry
## Christopher Hall

### Contents
1. Working Directory
2. Console Window
3. Objects
4. Data types
5. Vectors
6. Matrices
7. Lists
8. Data Frames
9. Scripts
10. Stats
11. Plots
12. Programming
13. For Loops
14. While Loops
15. If
16. Apply
17. Packages
18. Flow Cytometry
19. Plotting
20. Group of files
21. Gating
22. Autogating
23. Clustering
24. Exporting your data
25. Help
26. Analyse your data
27. Code
28. Results

R is case sensitive and the red indented text in this handout is code.
## Working Directory
This is where R looks for data and stores the output.  To see where it is type into the console:
```R
getwd()
```
To change the working directory type:
```R
setwd(“c:/myfiles/”)
```
Note that forward slashes are used instead of backslashes. Backslashes mean something in R, they are ‘escape characters’. For example, \n = new line, \t = tab
## Console Window
Think of it as a calculator.  Type commands and receive answers.  Finish each command by pressing enter.  Try:
```R
1+1
6*8
5+6*8
4/(4-2)
```
Pay attention to the order of calculations. You can use the up arrow to bring back previous inputs.
Brackets – orders (power) – division – multiplication – addition - subtraction
## Objects
Objects are things that mean other things ;) Type:
```R
myname <- “Chris”
myname
```
Strings, i.e. text needs to be enclosed in double quotes. Now try:
```R
a <- 1
b <- 2
a+b
a+b+myname
```
myname is not a number!  To remove an object type:
```R
rm(myname)
```
With objects you can be clever and make your life easier, for example:
```R
radius <- 10
height <- 15
vol_cyl <- pi*radius^2*height
vol_cyl
```
Or even better you could make a function to do the calculation by just typing vol_cyl:
```R
vol_cyl <- function(radius, height){
  volume <- pi*(radius^2)*height
  return(volume)
```
Then typing either of these commands will give you the cylinder volume: 
```R
vol_cyl(radius, height)
vol_cyl(10, 15)
```
You now have the volume of a cylinder. You can change radius and height as you wish and still use the vol_cyl object to calculate the volume.  pi is a built in constant. 
## Data types
**Character**: "words", “writing”
**Numeric**: 2, 15.5
**Integer**: 2L (the L tells R to store this as an integer, i.e. a whole number)
**Logical**: TRUE, FALSE
**Complex*: 1+4i (complex number)
```R
class(radius)
```
### Vectors
A collection of elements of a similar data type.
```R
MyVector <- c(10,11,12,13,14,15)
MyVector
```
The c() allows you to input multiple values at ones.  Another example:
```R
EasyVector <- c(10:14)
EasyVector
```
You can examine and change a vector using:
```R
MyVector[3]
length[MyVector]
EasyVector <- c(EasyVector, 15)
EasyVector 
```
Try this and think about the result:
```R
MyVector*EasyVector
```
### Matrices
A matrix is a vector with dimensions. 
```R
MyMatrix <- matrix(1:50, nrow=10, ncol=5) 
MyMatrix
```
### Lists
An ordered collection of elements that can be different data types.
```R
MyList <- c(name=“Chris”, number=a, matrix=MyMatrix)
MyList[[2]]
MyList[[“number”]]
```
### Data Frames
Are the basically tables.  Each column can have a different data type.
```R
x <- c(10,20,33,44)
y <- c("blue", "green", "purple", NA)
z <- c(FALSE,TRUE,TRUE,FALSE)
mydata <- data.frame(x,y,z)
mydata
```
**head()** - first 6 rows
**tail()** - last 6 rows
**dim()** - dimensions of data frame
**nrow()** - number of rows
**ncol()** - number of columns
**str()** - structure of data frame - name, type and preview of data in each column
**names()** or colnames() - both show the names attribute for a data frame
There are other data types, Google is your friend.
## Scripts
A script is a reusable list of commands, think text file.  Click File>New>R Script.
This will make an empty file where you can type in your commands.  To run them either click Ctrl+Enter on a line, or highlight the section and press Run.
Stats
R contains a number of example data sets.  Type data() to see them.
To load a dataset use data(Indometh), to see the data type Indometh, and to see a description use help(Indometh).
The dollar sign ($) allows you to select columns in data frames such as this. Try:
mean(Indometh$conc)
You can replace mean with mode, sd, max, min etc.  For help try help(max).
Help() is your friend.
## Plots
There are many ways to plot in R.  This is just one way.  Better ways include using ggplot2 and, for flow cytometry, flowViz or ggcyto.  More on this later. 
```R
hist(Indometh$conc)
hist(Indometh$conc, col=”red”, breaks=10)
boxplot(Indometh$conc)
plot(Indometh$time, Indometh$conc)
```
To format your table try adding:
•	xlim, ylim: range of the x-axis and the y-axis
•	xlab, ylab: labels for the x-axis and the y-axis
•	main: main title of the graph
•	pch: symbol type. See ?points.
•	cex: symbol size
•	col: symbol colour by either colour name or index. See colors().
```R
plot(Indometh$time, Indometh$conc, xlab="Time", ylab="Concentration", main="Pharmacokinetics of Indomethacin", col="red", pch=13)
```
Export your plot using the export button in the plot window. Or by typing this and checking your working directory:
```R
dev.copy(png,"myplot.png",width=8,height=6,units="in",res=100)
dev.off()
```
Always finish with dev.off() when exporting a plot by command line.
## Programming
### For Loops
‘For loops’ allow you to repeat sections of codes for a predetermined set of times.  The code to br run is enclosed in curley brackets{}, for example:
```R
for(i in 1:10){print(i)}
```
### While Loops
‘While loops’ repeat a section of code as long as a condition is met.  For example:
```R
countdown <- 10
while(countdown>0){
print(countdown)
countdown <- countdown-1}
```
### If
If statements respond to certain values.
```R
countdown <- 10
while(countdown>0){
  print(countdown)
  countdown <- countdown-1
  if (countdown == 0){ print("Takeoff")}
}
```
Pay attention to the curly brackets!
## Apply
The apply family of functions enable you to iterate through and to manipulate matrices, lists, arrays, and data frames. We care about the apply functions because they are used to handle groups of flow cytometry files in R.  There are a few different ones: lapply(), sapply(), maaply() and more.  Google is your friend, again.  We will use fsApply later.
## Packages
R has many packages which are pre written sets of R functions that can be used on your data.  For example ggplot2 has sophisticated plotting tools and flowCore allows you to analyse flow cytometry data.  They are written by the R community and can be found in two places.  The general purpose R packages are on CRAN and biological R packages are on BioConductor.  To load them:
```R
install.packages("ggplot2")
```
Or for Bioconductor:
```R
install.packages("BiocManager")
BiocManager::install(“flowCore”)
```
You may be asked which packages to update, you normally want to type ‘a’ and enter.
You must then apply the downloaded packages to your R installation using:
```R
library(ggplot2)
```
