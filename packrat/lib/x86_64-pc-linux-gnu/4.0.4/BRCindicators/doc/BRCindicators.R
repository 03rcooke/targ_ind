## ----knitr_setup, echo=FALSE--------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      cache = TRUE, 
                      tidy.opts = list(width.cutoff = 60),
                      tidy = FALSE)

## ----libraries----------------------------------------------------------------
require(snowfall)


## ----set_up, eval=FALSE-------------------------------------------------------
## library(devtools)
## install_github('biologicalrecordscentre/sparta')


## ----create_some_data---------------------------------------------------------
# Create data
n <- 8000 # size of dataset
nyr <- 50 # number of years in data
nSamples <- 200 # set number of dates
nSites <- 100 # set number of sites
set.seed(125) # set a random seed

# Create somes dates
first <- as.Date(strptime("1950/01/01", "%Y/%m/%d")) 
last <- as.Date(strptime(paste(1950+(nyr-1),"/12/31", sep=''), "%Y/%m/%d")) 
dt <- last-first 
rDates <- first + (runif(nSamples)*dt)

# taxa are set semi-randomly
taxa_probabilities <- seq(from = 0.1, to = 0.7, length.out = 26)
taxa <- sample(letters, size = n, TRUE, prob = taxa_probabilities)

# sites are visited semi-randomly
site_probabilities <- seq(from = 0.1, to = 0.7, length.out = nSites)
site <- sample(paste('A', 1:nSites, sep=''), size = n, TRUE, prob = site_probabilities)

# the date of visit is selected semi-randomly from those created earlier
time_probabilities <- seq(from = 0.1, to = 0.7, length.out = nSamples)
time_period <- sample(rDates, size = n, TRUE, prob = time_probabilities)

myData <- data.frame(taxa, site, time_period)


## ----library_sparta-----------------------------------------------------------
# Load the sparta package
library(sparta)


## ----format_data--------------------------------------------------------------
# Preview of my data
head(myData)

# First format our data
formattedOccData <- formatOccData(taxa = myData$taxa,
                                  site = myData$site,
                                  survey = myData$time_period)

# Here we are going to use the package snowfall to parallelise
library(snowfall)

# I have 4 cpus on my PC so I set cpus to 4
# when I initialise the cluster
sfInit(parallel = TRUE, cpus = 4)

# Export my data to the cluster
sfExport('formattedOccData')

# I create a function that takes a species name and runs my model
occ_mod_function <- function(taxa_name){
  
  library(sparta)
  
  # Note that this will write you results to your computer
  # the location is set to your user folder
  occ_out <- occDetFunc(taxa_name = as.character(taxa_name),
                        n_iterations = 200,
                        burnin = 15, 
                        occDetdata = formattedOccData$occDetdata,
                        spp_vis = formattedOccData$spp_vis,
                        write_results = TRUE,
                        output_dir = '~/Testing_indicator_pipe',
                        seed = 123)  
} 

# I then run this in parallel
system.time({
para_out <- sfClusterApplyLB(unique(myData$taxa), occ_mod_function)
})

# Stop the cluster
sfStop()

# We can see all the files this has created
list.files('~/Testing_indicator_pipe')


## ----eval = FALSE-------------------------------------------------------------
## library(devtools)
## install_github(repo = 'biologicalrecordscentre/BRCindicators')


## ----set_up_run---------------------------------------------------------------
library(BRCindicators)

# All we have to supply is the directory where out data is saved
# You will note this is the 'output_dir' passed to sparta above.
trends_summary <- summarise_occDet(input_dir = '~/Testing_indicator_pipe')

# Lets see the summary
head(trends_summary[,1:5])


## ----trends_summary-----------------------------------------------------------
trends_summary[1:3, 'a'] <- NA
trends_summary[1:5, 'b'] <- NA
trends_summary[2:4, 'c'] <- 1000
trends_summary[45:50, 'd'] <- NA

# Let's have a look at these changes
head(trends_summary[,1:5])
tail(trends_summary[,1:5])


## ----rescaled_indicator-------------------------------------------------------
# Let's run this data through our scaling function (all defaults used)
rescaled_trends <- rescale_species(Data = trends_summary)

# Here's the result
head(rescaled_trends[,c('year', 'indicator', 'a', 'b', 'c', 'd')])
tail(rescaled_trends[,c('year', 'indicator', 'a', 'b', 'c', 'd')])


## ----create_confidence_intervals----------------------------------------------
# This function takes just the species columns
scaled_species <- rescaled_trends[,!colnames(rescaled_trends) %in% c('year', 'indicator')]
indicator_CIs <- bootstrap_indicator(Data = scaled_species)

# Returned are the CIs for our indicator
head(indicator_CIs)


## ----smoothing_indicator------------------------------------------------------
# The smoothing function takes the indicator values
smoothed_indicator <- GAM_smoothing(rescaled_trends[,'indicator'])

# In this example there is little support for a non-linear trend and 
# so the line almost linear
plot(x = rescaled_trends[,'year'], y = rescaled_trends[,'indicator'])
lines(x = rescaled_trends[,'year'], y = smoothed_indicator, col = 'red')

# But if our indicator did support a non-linear trend it might look 
# like this
eg_indicator <- jitter(sort(rnorm(50)), amount = 0.5)
eg_smoothed <- GAM_smoothing(eg_indicator)
plot(x = 1:50, y = eg_indicator)
lines(x = 1:50, y = eg_smoothed, col = 'red')


## ----plot_indicator-----------------------------------------------------------
# Plot our indicator.
plot_indicator(indicator = rescaled_trends[,'indicator'],
               smoothed_line = smoothed_indicator,
               CIs = indicator_CIs)


## ----BMAdata------------------------------------------------------------------
# Here is an example dataset for the BMA method
data <- data.frame(species = rep(letters, each = 50),
                   year = rep(1:50, length(letters)), 
                   index = runif(n = 50 * length(letters), min = 0, max = 1), 
                   se = runif(n = 50 * length(letters), min = 0.01, max = .1))
head(data)


## ----runBMA-------------------------------------------------------------------
bma_indicator <- bma(data)


## ----runBMAparameters---------------------------------------------------------
bma_indicator2 <- bma(data,
                     parallel = TRUE,
                     n.iter = 500,
                     m.scale = 'log10')


## ----BMAresults---------------------------------------------------------------
head(bma_indicator)


## ----BMAplot------------------------------------------------------------------
plot_indicator(indicator = bma_indicator[,'Index'],
               CIs = bma_indicator[,c(3,4)])


## ----msi1, fig.height=4-------------------------------------------------------
# Create some example data in the format required
nyr = 20
species = rep(letters, each = nyr)
year = rev(rep(1:nyr, length(letters)))

# Create an index value that increases with time
index = rep(seq(50, 100, length.out = nyr), length(letters))
# Add randomness to species
index = index * runif(n = length(index), 0.7, 1.3)
# Add correlated randomness across species, to years
index = index * rep(runif(0.8, 1.2, n = nyr), length(letters))

se = runif(n = nyr * length(letters), min = 10, max = 20)

data <- data.frame(species, year, index, se)

# Our species are decreasing
plot(data$year, data$index)

# Species index values need to be 100 in the base year. Here I use
# the first year as my base year and rescale to 100. The standard error
# in the base year should be 0.
min_year <- min(data$year)

for(sp in unique(data$species)){

  subset_data <- data[data$species == sp, ]
  multi_factor <- 100 / subset_data$index[subset_data$year == min_year]
  data$index[data$species == sp] <- data$index[data$species == sp] * multi_factor
  data$se[data$species == sp] <- data$se[data$species == sp] * multi_factor
  data$se[data$species == sp][1] <- 0

}

# Our first year is now indexed at 100
plot(data$year, data$index)

# Alternativly I could read in data from a csv
# data <- read.csv('path/to/my/data.csv')

# Run the MSI function
msi_out <- msi(data)

head(msi_out$CV)

# I can capture the output figures too
# pdf('test.pdf')
#   msi_out <- msi(data)
# dev.off()


## ----msi2, fig.height=4-------------------------------------------------------
msi_out <- msi(data,
               nsim = 500, # The number of Mote Carlo simulations
               SEbaseyear = 10, # The year to index on
               plotbaseyear = 15, # The year to set as 100 in plots
               index_smoot = 'INDEX', # plotbaseyear uses MSI not trend
               span = 0.7, # 'wigglyness' of line, between 0 and 1
               lastyears = 5, # last X years of time series for short-term trends
               maxCV = 10, # maximum allowed Coefficient of Variation 
               changepoint = 10, # compare trends before and after this year
               truncfac = 8, # max year-to-year index ratio
               TRUNC = 5, #set all indices below TRUNC to this
               plot = TRUE # should the plots be returned?)
               )


## ----msi3---------------------------------------------------------------------
# The returned object has 2 elements
head(msi_out$results)


## ----msi4---------------------------------------------------------------------
# The returned object has 2 elements
msi_out$trends

# I could write this as a csv too
# write.csv(msi_out$trends, file = 'path/to/my/output.csv')


## ----msi_span, fig.height=4---------------------------------------------------
for(i in c(0.3, 0.5, 0.7)){ # use a range of values for span
  
  msi_out <- msi(data, span = i, # span is set to i
                 nsim = 200, plot = FALSE)
  print( # print makes the plot visible in the for loop
    plot(msi_out, title = paste('MSI - span =', i)) # plot
  )
  
}


## ----lambda_1-----------------------------------------------------------------
# number of species
nsp = 50

# number of years
nyr = 40

#number of iterations
iter = 500

# Build a random set of data
myArray <- array(data = rnorm(n = nsp*nyr*iter,
                               mean = 0.5,
                               sd = 0.1),
                  dim = c(nsp, nyr, iter),
                  dimnames = list(paste0('SP',1:nsp),
                                  1:nyr,
                                  1:iter))

# Ensure values are bounded by 0 and 1
myArray[myArray > 1] <- 1
myArray[myArray < 0] <- 0

str(myArray)


## ----lambda_2-----------------------------------------------------------------
# Run the lambda_interpolation method on this data
myIndicator <- lambda_indicator(myArray)

# Plot the indicator
plot_indicator(myIndicator$summary[,'indicator'],
               myIndicator$summary[,c('lower' ,'upper')])



## ----lambda_3-----------------------------------------------------------------
myIndicator <- lambda_indicator(myArray,
                                index = 1, # Set the index value to 1 not 100
                                year_range = c(30,40), # Set year range
                                threshold_yrs = 5) # set a threshold
plot_indicator(myIndicator$summary[,'indicator'],
               myIndicator$summary[,c('lower' ,'upper')])


## ----cache=TRUE---------------------------------------------------------------
# I call my function 'run_pipeline' and the only arguement it
# takes is the directory of sparta's output
run_pipeline <- function(input_dir){

  require(sparta)
  require(BRCindicators)
  
  # Create the trends summary
  trends_summary <- summarise_occDet(input_dir = input_dir)

  # Rescale the values and get the indicator values
  # Here I set the index to 1 and change the value limits
  rescaled_trends <- rescale_species(Data = trends_summary,
                                     index = 1,
                                     max = 100,
                                     min = 0.001)
  
  # Bootstrap the indicator to get CIs
  scaled_species <- rescaled_trends[,!colnames(rescaled_trends) %in% c('year', 'indicator')]
  # This time I set the iterations to twice the default and 
  # use custom confidence intervals
  indicator_CIs <- bootstrap_indicator(Data = scaled_species,
                                       CI_limits = c(0.25, 0.75),
                                       iterations = 20000)
  
  # Get the smoothed indicator line
  smoothed_indicator <- GAM_smoothing(rescaled_trends[,'indicator'])
  
  # This time I specify the years and index value
  plot_indicator(indicator = rescaled_trends[,'indicator'],
                 year = rescaled_trends[,'year'],
                 index = 1,
                 CIs = indicator_CIs,
                 smoothed_line = smoothed_indicator)
  
  ## I'll return all my data  
  return(cbind(smoothed_indicator, indicator_CIs, as.data.frame(trends_summary)))
 }


## ---- in_one------------------------------------------------------------------
# Now we can run the pipeline in one line, like a boss
indicator_data <- run_pipeline(input_dir = '~/Testing_indicator_pipe')

head(indicator_data)

