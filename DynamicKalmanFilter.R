## This is the main project file, it is where the code needs to be added



## preparing data frames

## the following data frame will save your Fingerprint estimates 
## please notice the column names

## the following data frame will save your kalaman estimates 
## please notice the column names

kalmanEstimates <- data.frame(
  measurement_x=double(),
  measurement_y=double(),
  estimate_x=double(),
  estimate_y=double(),
  deltaT = character(),
  TimeStamp = character()
  )


## put the Measurement point table with measured locations of the points here

 iterations <- nrow(MeasurementTable)

 #Normalization Matrix
 H = matrix(  c(1,0,0,0,
                0,1,0,0), # the data elements 
              nrow=2,              # number of rows 
              ncol=4,              # number of columns 
              byrow = TRUE)        # fill matrix by rows 
 
 #this is the initial state
 initial_xy = matrix(  c(0,0), # the data elements 
                       nrow=2,              # number of rows 
                       ncol=1,              # number of columns 
                       byrow = TRUE)        # fill matrix by rows
 
  #this is the covariance matrix, covariance taken as 1
 P = matrix( c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1), # the data elements 
             nrow=4,              # number of rows 
             ncol=4,              # number of columns 
             byrow = TRUE)        # fill matrix by rows 
 
 #this is the external motion matrix, we dont have any so its 0, it should be set accordingly
 u = matrix(  c(0,0,0,0), # the data elements 
              nrow=4,              # number of rows 
              ncol=1,              # number of columns 
              byrow = TRUE)        # fill matrix by rows 
 
 x = matrix(  c(initial_xy[1,1],initial_xy[2,1],0,0), # the data elements 
              nrow=4,              # number of rows 
              ncol=1,              # number of columns 
              byrow = TRUE)        # fill matrix by rows 
 
 
 
 
 for(i in 1:iterations){
   ## measurement matrix
   z = matrix(c (MeasurementTable$x[i],
              MeasurementTable$y[i]), # the data elements
          nrow=2,              # number of rows
         ncol=1,              # number of columns
        byrow = TRUE)        # fill matrix by rows
   
   ## change in time between each reading
   ## Add package 'lubridate'
   DeltaT <- as.numeric(ymd_hms(MeasurementTable$TimeStamp[i+1]) - ymd_hms(MeasurementTable$TimeStamp[i]))
   
   ## State transition matrix
   F = matrix( c(1,0,DeltaT,0,0,1,0,DeltaT,0,0,1,0,0,0,0,1), # the data elements 
             nrow=4,              # number of rows 
             ncol=4,              # number of columns 
             byrow = TRUE)        # fill matrix by rows 
   
   ##the accuracy depends on the problem and is ususally provided by the sensor with each reading.
   ## measurement noise matrix
   R = matrix( c(MeasurementTable[i]$accuracy_x,0,0,MeasurementTable[i]$accuracy_y), # the data elements
           nrow=2,              # number of rows
          ncol=2,              # number of columns
         byrow = TRUE)        # fill matrix by rows
   
   #Calculating Kalman Gain
   S = H %*% P %*% t(H)+ R 
   K = P%*%t(H)%*%ginv(S)
   
   #updating the measurements
   x = x + (K %*% (z - H %*% x))
   Y = z - H %*% x
   
   # computing the estimate and covraiance
   #time update
   x = (F %*% x) + u
   P = F %*% P %*% t(F)
   
   ## setting data for the kalamn estimates
 kalman_values <- data.frame(
   measurement_x = z[1,1],
   measurement_y= z[2,1],
   estimate_x=Y[1,1],
   estimate_y=Y[2,1],
   deltaT =DeltaT,
   TimeStamp = MeasurementTable[i]$TimeStamp
 )
 
 ## adding data to the kalman estimates
 kalmanEstimates <- rbind(kalmanEstimates,  kalman_values)
 
 }
 
 