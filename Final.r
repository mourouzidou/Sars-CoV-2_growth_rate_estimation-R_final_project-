# Read the files
observed_data <- readLines("ms_obs_final.out")
sim_data <- readLines("ms_sim_final.out")
pars_used <- readLines('pars_final.txt')

#   Convert the observed data into a matrix
observed_matrix <- matrix(unlist(observed_data), 
                          nrow = length(observed_data), byrow = TRUE)

# _____________________________________________________
#                SIMULATION DATA PREPARATION
#
# Initialize a list to store the simulated datasets as matrices
sim_matrices <- list()

# Loop through each dataset in the datasets list
for (dataset in datasets) {
  # Retrieve the dataset data (second element)
  dataset_data <- dataset$data
  dataset_matrix <- matrix(unlist(dataset_data), 
                           nrow = length(dataset_data), byrow = TRUE)
    sim_matrices[[length(sim_matrices) + 1]] <- dataset_matrix
}

#________________________________________________________





#   Initialize an empty list to store the datasets
#   and the growth parameter used for each
datasets <- list()

#   Initialize variables
current_dataset <- NULL
line_count <- 0
dataset_index <- 1

#   loop that iterates over the lines and 
#   seperates the datasets when an empty line is met
for (line in sim_data) {
  if (line == "") {
    if (line_count > 0) {
      # Store the current dataset in the list
      datasets[[dataset_index]] <- list(name = pars_used[dataset_index], 
                                        data = current_dataset)
      current_dataset <- NULL
      line_count <- 0
      dataset_index <- dataset_index + 1 
    }
  } else {
    # Append the line to the current dataset
    current_dataset <- c(current_dataset, line)
    line_count <- line_count + 1
  }
}

# store the last dataset
datasets[[dataset_index]] <- list(name = pars_used[dataset_index], 
                                  data = current_dataset)

# calulate number of sequences function (rows)
n_calc <- function(matrix_data){
  return(nrow(matrix_data))
}

# calculate number of positions function (characters of a sequence)
# all the sequences have the same length, 
# thus we simply use nchar for the 1st row 
S_calc <- function(matrix_data) {
  return(nchar(matrix_data[1,]))
}

#...........................................................................
#                    Calculate K statistic function
#...........................................................................

calculate_K <- function(matrix_data) {
  d <- 0
  n <- n_calc(matrix_data)
  combs <- (n * (n - 1)) / 2
  
  for (i in 1:(n - 1)) {
    seq1 <- matrix_data[i, ]
    for (j in (i + 1):n) {
      seq2 <- matrix_data[j, ]
      diff_count <- sum(utf8ToInt(seq1) != utf8ToInt(seq2))
      d <- d + diff_count
    }
  }
  
  k <- d / combs
  return(k)
}



#...........................................................................
#                    Calculate W statistic function
#...........................................................................


#       calculate a1 function
a1_calc <- function(matrix_data) {
  a1 <- 0
  n <- n_calc(matrix_data)
  for (s in 1:(n - 1)) {
    a <- 1/s
    a1 <- a1 + a
  }
  return(a1)
}



#      calculate W statistic using a1 and S 
W_calc <- function(matrix_data) {
  return(S_calc(matrix_data)/a1_calc(matrix_data))
}



#...........................................................................
#                    Calculate Tajima's D statistic function
#...........................................................................

#       function for a2 calculation 
a2_calc <- function(matrix_data) {        
  a2 <- 0
  for (s in 1:(n_calc(matrix_data) - 1)) {
    a <- 1/s**2
    a2 <- a2 + a
  }
  return(a2)
}

#      function for b1 calculation
b1_calc <- function(matrix_data) {
  n = n_calc(matrix_data)
  return((n + 1)/(3*(n-1)))
}

#      function for b2 calculation
b2_calc <- function(matrix_data){
  n = n_calc(matrix_data)
  return((2*(n**2 + n + 3))/(9*n*(n-1)))
}

#      function for c1 calculation
c1_calc <- function(matrix_data){
  b1 <- b1_calc(matrix_data)
  a1 <- a1_calc(matrix_data)
  return(b1 - 1/a1)
}

#       function for c2 calculation
c2_calc <- function(matrix_data){
  b2 <- b2_calc(matrix_data)
  a1 <- a1_calc(matrix_data)
  a2 <- a2_calc(matrix_data)
  n <- n_calc(matrix_data)
  return(b2 - ((n+2)/(a1*n)) + (a2/(a1^2)))
}

#       function for e1 calculation
e1_calc <- function(matrix_data){
  c1 <- c1_calc(matrix_data)
  a1 <- a1_calc(matrix_data)
  return(c1/a1)
}

#       function for e2 calculation
e2_calc <- function(matrix_data) {
  c2 <- c2_calc(matrix_data)
  a1 <- a1_calc(matrix_data)
  a2 <- a2_calc(matrix_data)
  return(c2/(a1^2 + a2))
  
}

#   Finally calculate Tajima's D

D_calc <- function(matrix_data){
  K <- calculate_K(matrix_data)
  W <- W_calc(matrix_data)
  e1 <- e1_calc(matrix_data)
  e2 <- e2_calc(matrix_data)
  S <- S_calc(matrix_data)
  return((K-W)/sqrt((e1*S)+ (e2*S)*(S-1)))
  
}

#__________________________________________________________________________
# calculate the K, W, Tajima's D statistic for observed dataset (K0, W0, D0)
#___________________________________________________________________________
  
W0 <- W_calc(observed_matrix)    # 29.46951
D0 <- D_calc(observed_matrix)    # -1.845542
K0 <- calculate_K(observed_matrix) # 14.28245

#__________________________________________________________________________
#
#         Calculate K, W, D statistics for every simulated dataset
#__________________________________________________________________________

# Initialize variables
k_sim_values <- numeric(length(sim_matrices))
D_sim_values <- numeric(length(sim_matrices))
W_sim_values <- numeric(length(sim_matrices))

# Loop through each matrix in the sim_matrices list
for (i in 1:length(sim_matrices)) {
  matrix_data <- sim_matrices[[i]]
  
  
  k_sim_values[i] <- calculate_K(matrix_data)
  D_sim_values[i] <- D_calc(matrix_data)
  W_sim_values[i] <- W_calc(matrix_data)
}

#__________________________________________________________________________
#
#                 Normalize the statistic vectors
#__________________________________________________________________________

#   calculate the mean values of the statistic matrices
k_mean <- mean(k_sim_values)
d_mean <- mean(D_sim_values)
w_mean <- mean(W_sim_values)

#  calculate the std values for the statistic matrices
k_std <- sd(k_sim_values)
d_std <- sd(D_sim_values)
w_std <- sd(W_sim_values)

# create new normalized vectors for all the statistics
k_normalized <- (k_sim_values - k_mean)/k_std
d_normalized <- (D_sim_values - d_mean)/d_std
w_normalized <- (W_sim_values - w_mean)/w_std

# normalize also the observed statistics using the same parameters (mean and sd)
k0_norm <- (K0 - k_mean)/k_std
d0_norm <- (D0 - d_mean)/d_std
w0_norm <- (W0 - w_mean)/w_std


#___________________________________________________________________
#                 Statistic comparison using euclidean distance 
#___________________________________________________________________
# calculate the euclidean distances between the observed 
#normalized and each of the normalized simulated datasets
euc_distances = sqrt((k_normalized - k0_norm)^2 +
                       (d_normalized - d0_norm)^2 +
                       (w_normalized - w0_norm)^2)

# Find the 500 smallest distances and keep their indexes
smallest_indexes <- order(euc_distances, decreasing = FALSE)[1:500]

# get the 500 smallest distances by their indexes
smallest_distances <- euc_distances[smallest_indexes]

# get the growth parameters that generated the 500 datasets 
#with the smallest distance from the observed dataset
best_growth <- pars_used[smallest_indexes]

best_growth <- as.numeric(best_growth) 
mean(best_growth)    #110.1053
median(best_growth)  #106.7937






#_______________________________________________
#_______________PLOTS___________________________

#________________________
# Construct the histogram
#________________________
hist(best_growth, 
     main = "Distribution of Growth Parameters",
     xlab = "Growth Parameter",
     ylab = "Frequency",
     col = "skyblue",
     border = "white",
     breaks = 20)

# Add a density line
lines(density(best_growth), col = "blue", lwd = 2)

# Add labels for mean and median
abline(v = mean(best_growth), col = "red", lwd = 2, lty = 2)
text(mean(best_growth), 30, "Mean", col = "red", pos = 3)

abline(v = median(best_growth), col = "green", lwd = 2, lty = 2)
text(median(best_growth), 40, "Median", col = "green", pos = 3)


#____________________________
# Construct the density plot
#____________________________
plot(density(best_growth),
     main = "Density Plot of Growth Parameters",
     xlab = "Growth Parameter",
     ylab = "Density",
     col = "steelblue",
     lwd = 2)

# Add a rug plot
rug(best_growth, col = "gray")

# Add labels for mean and median
abline(v = mean(best_growth), col = "red", lwd = 2, lty = 2)
text(mean(best_growth), 0.015, "Mean", col = "red", pos = 1)

abline(v = median(best_growth), col = "green", lwd = 2, lty = 2)
text(median(best_growth), 0.017, "Median", col = "green", pos = 1)



