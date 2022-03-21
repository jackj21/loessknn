library(R.methodsS3) 
library(ggplot2)
  # Your function will have the following inputs.
  # 
  # * x - a numeric input vector
  # * y - a numeric response
  #
  # Note span and degree are shown with their default values. (Read about this in the description)
  # * degree should be 1 or 2 only
  # * span can be any value in interval (0, 1) non-inclusive.
  #
  # If show.plot = TRUE then you must show a plot of either the final fit

myloess <- function(x, y, span = 0.5, degree = 1, show.plot = TRUE){
  
  # Return error messages when the input values are invalid.
  if (!is.numeric(x) || !is.vector(x)) stop("x should be a numeric input vector")
  if (!is.numeric(y) || !is.vector(y)) stop("y should be a numeric response")
  if (degree != 1 && degree != 2) stop("Degree should be 1 or 2 only")
  if (span >= 1 || span <= 0) stop("Span can be any value in interval (0,1) non-inclusive")
  if (!is.logical(show.plot)) stop("show.plot should be either TRUE or FALSE")
  
  # N_total: total number of points in the data set
  N_total <- length(x)
  # n_points: number of points in each window in a vector
  n_points <- round(N_total*span)
  # Win_total: total number of windows
  Win_total <- 0
  # SSE: Error Sum of Squares (Tells us how good of a fit we had).
  SSE <- 0
  
  y_fitted <- c()
  closest_x <- c()
  closest_y <- c()
  #
  for (N in 1:N_total){
    
    # Find the closest n_points
    x_temp <- x
    y_temp <- y
    prev_win <- closest_x
    
    ## Find the closest one point, delete the point, repeat n_points times.
    for (n in 1:n_points){
      index <- which.min(abs(x_temp - x[N]))
      closest_x[n] <- x_temp[index]
      closest_y[n] <- y_temp[index]
      x_temp <- x_temp[-index]
      y_temp <- y_temp[-index]
    }
    
    # A new window is created if previous window is not identical to this one
    if(!setequal(prev_win,closest_x)){
      Win_total <- Win_total + 1
    }
    
    # Calculate the distance
    distance <- abs(closest_x - x[N])
    
    # Calculate the scaled distance(distance divided by the farthest point)
    scaled_distance <- distance/distance[length(distance)]
    
    # Calculate the weight of the points
    wts <- (1-abs(scaled_distance)^3)^3
    
    # If degree is 1 use (y~x)
    # If degree is 2 use (y ~ x + x^2)
    if (degree == 1) {
      wls_model <- lm(closest_y ~ closest_x, weights=wts)
      y_fitted[N] <- wls_model$coefficient[1] + x[N]*wls_model$coefficient[2]
    }
    
    if (degree == 2) {
      wls_model <- lm(closest_y ~ poly(closest_x, 2, raw=TRUE), weights=wts)
      y_fitted[N]<- wls_model$coefficient[1] + x[N]*wls_model$coefficient[2] + x[N]^2*wls_model$coefficients[3]
    }
    
    SSE <- SSE + (mean(y)-y_fitted[N])^2
  }
  
  # Create the loessplot
  df <- data.frame(x,y_fitted)
  loessplot <- ggplot(df,aes(x,y_fitted)) + geom_line()
  loessplot <- loessplot + geom_point(aes(x,y)) + ggtitle(span)
  
  # If show.plot = TRUE, show the plot
  if (show.plot == TRUE) {
    loessplot
  }
  
  # Create the returned list
  final_list <- list(span, degree, N_total, Win_total, n_points, SSE, loessplot)
  names(final_list) <- c("span", "degree", "N_total", "Win_total", "n_points","SSE","loessplot")
  return(final_list)
}


# Your function should return a named list containing the following:
# span: proportion of data used in each window (controls the bandwidth)
# degree: degree of polynomial
# N_total: total number of points in the data set
# Win_total: total number of windows
# n_points: number of points in each window in a vector
# SSE: Error Sum of Squares (Tells us how good of a fit we had).
# loessplot: An object containing the ggplot so that we can see the plot later. 
#  We want this even if show.plot = FALSE
#  Note: you are NOT allowed to simply use stat_smooth() or geom_smooth() to have it automatically do LOESS.
#  You should use geom_line() or similar to plot your final the LOESS curve.

# Make sure you can access the objects properly using the $ notation.

# Your function will have the following inputs similar to what you would find with the
#  knn() function
#
# * train - matrix or data frame of training set cases
# * test - matrix or data frame of test set cases.  
#     (A vector will be interpreted as a row vector for a single case.)
# * y_train - Either a numeric vector, or factor vector for the responses in the training set
# * y_test - Either a numeric vector, or factor vector for the responses in the testing set
# * k - number of neighbors considered, the default value is 3
#
# If weighted = TRUE, then your function must used the distance weighted kNN as described above,
#  otherwise it should do the default knn method.


mykNN <- function(train, test, y_train, y_test, k = 3, weighted = TRUE){
  
  # Function to calculate mode to get most occurrences
  mode <- function(x) {
    u <- unique(x)
    u[which.max(tabulate(match(x,u)))]
  }
  
  # Function to calculate accuracy of predictions from confusion matrix
  accuracy_cm <- function(x) { sum(diag(x) / sum(rowSums(x))) * 100}
  
  # Create new vectors for distance, weights, yhat, accuracy
  dist <- c()
  yhat <- c()
  levels <- levels(y_test)
  
  # classification problem
  if (is.factor(y_test)) {
    
    # default kNN 
    if (weighted == FALSE) {
    
      # Iterate through points in test and train and calculate distance
      for (i in 1:nrow(test)) {
        for (j in 1:nrow(train)) {
          l2 <- sqrt(sum((train[j,] - test[i,])^2))    # Calculate Euclidean distance
          dist[j] <- l2                              # Assign distance to dist vec
          
        } 

        neighbors <- data.frame(dist, y_train)
        n <- neighbors[order(dist),]
        nn <- n[1:k,]
        
        yhat[i] <- mode(nn$y_train)
        
      }
      
      yhat <- factor(levels[yhat], levels)
      
      confusion <- table(yhat, y_test)  # Confusion matrix
      accuracy <- accuracy_cm(confusion)
      error <- 100 - accuracy                        # Error rate of predictions
      
      classification <- list(yhat, accuracy, error, confusion, k)
      names(classification) <- c("yhat", "accuracy", "error", "confusion", "k")
      
      return(classification)
    }
  
    
  # distance weighted kNN 
  else {
    # Iterate through points in test and train and calculate distance
    for (i in 1:nrow(test)) {
      for (j in 1:nrow(train)) {
        l2 <- sqrt(sum((train[j,] - test[i,])^2))    # Calculate Euclidean distance
        dist[j] <- l2                              # Assign distance to dist vec
        
      } 
      weights <- 1 / dist 
      neighbors <- data.frame(dist, weights, y_train)
      n <- neighbors[order(dist),]
      nn <- n[1:k,]
      df <- aggregate(nn$weights, list(nn$y_train), sum)
      order(df[,2], decreasing=TRUE)
      
      yhat[i] <- df[1,1]

    }
    
    yhat <- factor(levels[yhat], levels)

    confusion <- table(yhat, y_test)             # Confusion matrix
    accuracy <- accuracy_cm(confusion)              # Accuracy rate of predictions
    error <- 100 - accuracy                        # Error rate of predictions
    classification <- list(yhat, accuracy, error, confusion, k)
    names(classification) <- c("yhat", "accuracy", "error", "confusion", "k")
    
    return(classification)
    
  }
    
  } # end classification 
    
  # regression problem - same thing as classification just average values
  else {
    # default kNN 
    if (weighted == FALSE) {
      
      # Iterate through points in test and train and calculate distance
      for (i in 1:nrow(test)) {
        for (j in 1:nrow(train)) {
          l2 <- sqrt(sum((train[j,] - test[i,])^2))    # Calculate Euclidean distance
          dist[j] <- l2                              # Assign distance to dist vec
          
        } 
        
        neighbors <- data.frame(dist, y_train)
        n <- neighbors[order(dist),]
        nn <- n[1:k,]
        
        yhat[i] <- mean(nn$y_train)                  # Assign mean as predicted point to yhat
        
      }
      
      residual <- y_test - yhat
      sse <- sum(residual^2)
      
      regression <- list(yhat, residual, sse, k)
      names(regression) <- c("yhat", "residual", "SSE", "k")
      
      return(regression)
    }
    
    
    # distance weighted kNN 
    else {
      # Iterate through points in test and train and calculate distance
      for (i in 1:nrow(test)) {
        for (j in 1:nrow(train)) {
          l2 <- sqrt(sum((train[j,] - test[i,])^2))    # Calculate Euclidean distance
          dist[j] <- l2                              # Assign distance to dist vec
          
        } 
        weights <- 1 / dist 
        neighbors <- data.frame(dist, weights, y_train)
        n <- neighbors[order(dist),]
        nn <- n[1:k,]
        
        yhat[i] <- mean(nn$y_train)                 # Assign mean as predicted point to yhat
        
      }
      residual <- y_test - yhat
      sse <- sum(residual^2)
      regression <- list(yhat, residual, sse, k)
      names(regression) <- c("yhat", "residual", "SSE", "k")
      
      return(regression)
      
    }
   

  }
}
# If you are doing classification, then your function must return:
#  * A factor vector (yhat) for the predicted categories for the testing data
#  * The accuracy of the classification
#  * The error rate = 1 - accuracy
#  * A confusion matrix
#  * The value of k used

# If you are doing regression, then your function must return:
#  * A numeric vector (yhat) for the predicted responses for the testing data
#  * The residual vector
#  * The SSE
#  * The value of k used
