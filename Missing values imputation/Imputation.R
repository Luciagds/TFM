# This code loads NOAA data using the worldmet R package and imputes 
# missing values using a user-defined function that computes the IDW
# method, followed by either ARIMA or Moving Average imputation (R's
# package imputeTS). A statistical analysis is performed to select 
# the optimal method.

# Lucía García-Duarte Sáenz

memory.limit(6400000)
options(warn=-1)

library(dplyr)
library(tidyverse)
library(plyr)
library(naniar)
library(xts)
library(fpp2)
library(cowplot)
library(reshape2)
library(data.table)
library(imputeTS)

####################################################################
##################         GET DATA           ######################
####################################################################

# Load station names
stations = read.csv('NOAA_summary_opt1.csv', header = TRUE, colClasses= 'character')
stations = stations$code

# Create sequence of data
start <- as.POSIXct("1996-01-01 00:00:00", format="%Y-%m-%d %H:%M:%S")
interval <- 60
end <- as.POSIXct("2021-03-28 21:00:00", format="%Y-%m-%d %H:%M:%S")
date=seq(from=start, by=interval*60, to=end)
days_seq=tibble(date)

# Load only a station
library(worldmet) # NOAA
data = importNOAA(
    code = stations[1],
    year = 1996:2021,
    hourly = TRUE,
    n.cores = 4,
    quiet = FALSE,
    path = NA
  )
data = data[!duplicated(data$date,fromLast = TRUE),] # remove duplicates
mydata <- data[c('code','date','longitude','latitude','ws','wd','air_temp','visibility','dew_point','RH')] # select variables of interest

# Join with my data and convert to xts 
mydata = days_seq %>% left_join(mydata)
mydata[which(is.na(mydata$code)),]$code = stations[1]
mydata[which(is.na(mydata$latitude)),]$latitude = as.numeric(latitude[1])
mydata[which(is.na(mydata$longitude)),]$longitude = as.numeric(longitude[1])

# Store the data
alldata <- rbind.fill(mydata)

# Start loop
for (i in 2:37){
  #print(i)
  
  library(worldmet) # NOAA
  data = importNOAA(
      code = stations[i],
      year = 1996:2021,
      hourly = TRUE,
      n.cores = 1,
      quiet = FALSE,
      path = NA
    )
  
  data = data[!duplicated(data$date,fromLast = TRUE),]
  mydata <- data[c('code','date','longitude','latitude','ws','wd','air_temp','visibility','dew_point','RH')]

  # join with my data and convert to xts
  mydata = days_seq %>% left_join(mydata)
  mydata[which(is.na(mydata$code)),]$code = stations[i]
  mydata[which(is.na(mydata$latitude)),]$latitude = as.numeric(latitude[i])
  mydata[which(is.na(mydata$longitude)),]$longitude = as.numeric(longitude[i])

  alldata <- rbind.fill(alldata,mydata)
}

# Save file
write.csv(alldata,'imputed_all_vars.csv')


####################################################################
##################     SPATIAL IMPUTATION     ######################
####################################################################

### INVERSE DISTANCE WEIGHTING METHOD (IDWM)

# Define distance and imputation functions
get_nn <- function(mystat, stations, limit){
  # function to compute the neighbours and their distance 
  #      mystat: the station to process
  #      stations: the list of stations to compare with
  #      limit: maximum distance to consider a neighbor (in km) 
  
  # Compute distances between mystat and the other stations
  r = 6371 # earth radius (km)
  dist = matrix(0, length(stations$code))
  s = which(stations$code==mystat)
  
  # Loop over the stations
  for (i in 1:length(stations$code)){
    if (mystat != stations$code[i]){
      cc = sin(pi/180*(stations$lat[s] - stations$lat[i])/2)^2 + cos(pi/180*stations$lat[s])*cos(pi/180*stations$lat[i])*sin(pi/180*(stations$lon[s] - stations$lon[i])/2)^2
      dist[i] = 2*r*asin(sqrt(cc)) # haversine formula
    }
  }
  
  # Sort the stations by increasing distance and select neighbours
  dist.indx = order(dist)
  dist = dist[dist.indx]
  pos = which(dist < limit)[-1] # remove the station itself (mystat dist = 0)
  distances = as.data.frame(cbind(dist[pos],dist.indx[pos])) # select neigbours by maximum ratio
  colnames(distances) = c('dist','dist.indx')
  return (distances)
}
invdistImpute <- function(data, stations,limit) {
  # function to impute using inverse distance method 
  # data must be of the form: 
  #      each column is a station with the values of a specific variable, and the date column
  
  data_imp = copy(data)
  
  # Loop over the stations (each data column)
  for (i in 1:nrow(stations)){ 
    #print(i)
    myNA = which(is.na(data[i]))
    mystat = stations$code[i]
    dist = get_nn(mystat,stations,limit=limit)
    
    # Loop over the missing values of the mystat
    for (j in myNA){
      
      # Check that there are some stations within 100km radius
      if (length(dist$dist) != 0){
        values = data[j,dist$dist.indx]
        indeces = which(!is.na(values))
        
        # Check that the stations do not have missing values
        if (length(indeces)!= 0){
          values = values[indeces]
          val = as.double(values)%*%(dist$dist[indeces])^(-2)/sum((dist$dist[indeces])^(-2)) 
          
          # QUALITY CONTROL
          # Select range (1 day, data point in the middle)
          if (j<=12){
            indx = seq(1,25)
          } else if (j>=length(days_seq$date)-12){
            indx = seq(length(days_seq$date)-24,length(days_seq$date))
          } else {
            indx = seq(j-12,j+12)
          }
          # Check max/min values
          if (val>max(data[indx, i], na.rm = T)&(max(data[indx, i], na.rm = T)!=-Inf)){
            val = max(data[indx, i], na.rm = T)
          }
          if (val<min(data[indx, i], na.rm = T)&(min(data[indx, i], na.rm = T)!=Inf)){
            val = min(data[indx, i], na.rm = T)
          }
          
          # Impute
          data_imp[j,i] = val
        }
      }
    }
  }
  return(data_imp)
}

# Load data
alldata = read.csv('imputed_all_vars.csv', header = TRUE)

# Make stations and data_idw have the same order of the stations
stations = read.csv('NOAA_summary_opt1.csv', header = TRUE)
stations = stations %>%
  arrange(code)

limit = 100

######################### AIR_TEMP ######################### 
# Perform imputation
data_idw = copy(alldata)
data_idw = as.data.frame(alldata[, names(alldata) %in% c("air_temp","code",'date')])
statsNA(data_idw$air_temp) # check %missings
data_idw = spread(data_idw, code, air_temp)
data_idw = data_idw[,c('date',stations$code)]
data_idw$date = NULL

my_var = copy(data_idw)
statistics = as.data.frame(colSums(my_var, na.rm = T)/nrow(my_var))
colnames(statistics) = c('mean_air_temp')
statistics$sd_air_temp = apply(my_var,2, function (x) sd(x, na.rm = TRUE))
statistics$median_air_temp = apply(my_var,2, function (x) median(x, na.rm = TRUE))
statistics$min_air_temp = apply(my_var,2, function (x) min(x, na.rm = TRUE))
statistics$max_air_temp = apply(my_var,2, function (x) max(x, na.rm = TRUE))

data_imputed_air_temp = invdistImpute(data_idw,stations,limit=limit)
mv_air_temp=c()
for (i in 1:37){
  mv_air_temp[i] = which(is.na(data_imputed_air_temp[i])==FALSE)[1] 
}
data_imputed_air_temp$date = days_seq$date
data_imputed_air_temp=data_imputed_air_temp[,c(38,1:37)]


######################### WS ######################### 
# Perform imputation
data_idw = copy(alldata)
data_idw = as.data.frame(alldata[, names(alldata) %in% c("ws","code",'date')])
statsNA(data_idw$ws) # check %missings
data_idw = spread(data_idw, code, ws)
data_idw = data_idw[,c('date',stations$code)]
data_idw$date = NULL

my_var = copy(data_idw)
statistics$mean_ws = colSums(my_var, na.rm = T)/nrow(my_var)
statistics$sd_ws = apply(my_var,2, function (x) sd(x, na.rm = TRUE))
statistics$median_ws = apply(my_var,2, function (x) median(x, na.rm = TRUE))
statistics$min_ws = apply(my_var,2, function (x) min(x, na.rm = TRUE))
statistics$max_ws = apply(my_var,2, function (x) max(x, na.rm = TRUE))

data_imputed_ws = invdistImpute(data_idw,stations,limit=limit)
mv_ws=c()
for (i in 1:37){
  mv_ws[i] = which(is.na(data_imputed_ws[i])==FALSE)[1] 
}
data_imputed_ws$date = days_seq$date
data_imputed_ws=data_imputed_ws[,c(38,1:37)]

######################### WD ######################### 
# Perform imputation
data_idw = copy(alldata)
data_idw = as.data.frame(alldata[, names(alldata) %in% c("wd","code",'date')])
statsNA(data_idw$wd) # check %missings
data_idw = spread(data_idw, code, wd)
data_idw = data_idw[,c('date',stations$code)]
data_idw$date = NULL

my_var = copy(data_idw)
statistics$mean_wd = colSums(my_var, na.rm = T)/nrow(my_var)
statistics$sd_wd = apply(my_var,2, function (x) sd(x, na.rm = TRUE))
statistics$median_wd = apply(my_var,2, function (x) median(x, na.rm = TRUE))
statistics$min_wd = apply(my_var,2, function (x) min(x, na.rm = TRUE))
statistics$max_wd = apply(my_var,2, function (x) max(x, na.rm = TRUE))

data_imputed_wd = invdistImpute(data_idw,stations,limit=limit)
mv_wd=c()
for (i in 1:37){
  mv_wd[i] = which(is.na(data_imputed_wd[i])==FALSE)[1] 
}
data_imputed_wd$date = days_seq$date
data_imputed_wd=data_imputed_wd[,c(38,1:37)]
#data_imputed_wd[1:10945,21] <- NA

######################### VISIBILITY ######################### 
# Perform imputation
data_idw = copy(alldata)
data_idw = as.data.frame(alldata[, names(alldata) %in% c("visibility","code",'date')])
statsNA(data_idw$visibility) # check %missings
data_idw = spread(data_idw, code, visibility)
data_idw = data_idw[,c('date',stations$code)]
data_idw$date = NULL

my_var = copy(data_idw)
statistics$mean_visibility = colSums(my_var, na.rm = T)/nrow(my_var)
statistics$sd_visibility = apply(my_var,2, function (x) sd(x, na.rm = TRUE))
statistics$median_visibility = apply(my_var,2, function (x) median(x, na.rm = TRUE))
statistics$min_visibility = apply(my_var,2, function (x) min(x, na.rm = TRUE))
statistics$max_visibility = apply(my_var,2, function (x) max(x, na.rm = TRUE))

data_imputed_visibility = invdistImpute(data_idw,stations,limit=limit)
mv_visibility=c()
for (i in 1:37){
  mv_visibility[i] = which(is.na(data_imputed_visibility[i])==FALSE)[1] 
}
data_imputed_visibility$date = days_seq$date
data_imputed_visibility=data_imputed_visibility[,c(38,1:37)]

######################### RH ######################### 
# Perform imputation
data_idw = copy(alldata)
data_idw = as.data.frame(alldata[, names(alldata) %in% c("RH","code",'date')])
statsNA(data_idw$RH) # check %missings
data_idw = spread(data_idw, code, RH)
data_idw = data_idw[,c('date',stations$code)]
data_idw$date = NULL

my_var = copy(data_idw)
statistics$mean_RH = colSums(my_var, na.rm = T)/nrow(my_var)
statistics$sd_RH = apply(my_var,2, function (x) sd(x, na.rm = TRUE))
statistics$median_RH = apply(my_var,2, function (x) median(x, na.rm = TRUE))
statistics$min_RH = apply(my_var,2, function (x) min(x, na.rm = TRUE))
statistics$max_RH = apply(my_var,2, function (x) max(x, na.rm = TRUE))

data_imputed_RH = invdistImpute(data_idw,stations,limit=limit)
mv_RH=c()
for (i in 1:37){
  mv_RH[i] = which(is.na(data_imputed_RH[i])==FALSE)[1] 
}
data_imputed_RH$date = days_seq$date
data_imputed_RH=data_imputed_RH[,c(38,1:37)]

######################### DEW_POINT ######################### 
# Perform imputation
data_idw = copy(alldata)
data_idw = as.data.frame(alldata[, names(alldata) %in% c("dew_point","code",'date')])
statsNA(data_idw$dew_point) # check %missings
data_idw = spread(data_idw, code, dew_point)
data_idw = data_idw[,c('date',stations$code)]
data_idw$date = NULL

my_var = copy(data_idw)
statistics$mean_dew_point = colSums(my_var, na.rm = T)/nrow(my_var)
statistics$sd_dew_point = apply(my_var,2, function (x) sd(x, na.rm = TRUE))
statistics$median_dew_point = apply(my_var,2, function (x) median(x, na.rm = TRUE))
statistics$min_dew_point = apply(my_var,2, function (x) min(x, na.rm = TRUE))
statistics$max_dew_point = apply(my_var,2, function (x) max(x, na.rm = TRUE))

data_imputed_dew_point = invdistImpute(data_idw,stations,limit=limit)
mv_dew_point=c()
for (i in 1:37){
  mv_dew_point[i] = which(is.na(data_imputed_dew_point[i])==FALSE)[1] 
}
data_imputed_dew_point$date = days_seq$date
data_imputed_dew_point=data_imputed_dew_point[,c(38,1:37)]

write.csv(statistics,'statistics_raw_mean_sd.csv')


# Check %missings
data_air_temp_melted <- melt(data_imputed_air_temp, id=c("date"))
colnames(data_air_temp_melted) <- c('date','code','air_temp')
statsNA(data_air_temp_melted$air_temp)

data_ws_melted <- melt(data_imputed_ws, id=c("date"))
colnames(data_ws_melted) <- c('date','code','ws')
statsNA(data_ws_melted$ws)

data_wd_melted <- melt(data_imputed_wd, id=c("date"))
colnames(data_wd_melted) <- c('date','code','wd')
statsNA(data_wd_melted$wd)

data_visibility_melted <- melt(data_imputed_visibility, id=c("date"))
colnames(data_visibility_melted) <- c('date','code','visibility')
statsNA(data_visibility_melted$visibility)

data_RH_melted <- melt(data_imputed_RH, id=c("date"))
colnames(data_RH_melted) <- c('date','code','RH')
statsNA(data_RH_melted$RH)

data_dew_point_melted <- melt(data_imputed_dew_point, id=c("date"))
colnames(data_dew_point_melted) <- c('date','code','dew_point')
statsNA(data_dew_point_melted$dew_point)


####################################################################
##################     TEMPORAL IMPUTATION    ######################
####################################################################

myImp <- function(data,dates){
  for (i in 2:38){
    print(i)
    mydata <- data[i] 
    mydata = cbind(dates, mydata)
    colnames(mydata) = c('date','variable')
    mydata <-  as.xts(mydata$variable, order.by = mydata$date) 
    data[i] = as.data.frame(na_kalman(mydata, model = "auto.arima", smooth=TRUE))$V1 # ARIMA 
    #data[i] = as.data.frame(na.ma(mydata, k=4))$V1                                  # MOVING AVERAGE

  }
  return (data)
}

######################### AIR_TEMP ######################### 
data_imputed_air_temp = myImp(data_imputed_air_temp,days_seq)
# QC (quality control)
data_idw_var = copy(alldata)
data_idw_var = as.data.frame(alldata[, names(alldata) %in% c("air_temp","code",'date')])
data_idw_var = spread(data_idw_var, code, air_temp)
data_idw_var = data_idw_var[,c('date',stations$code)]
data_idw_var$date = NULL
mymax = apply(data_idw_var,2, function (x) max(x, na.rm = TRUE))
names(mymax) = NULL
mymin = apply(data_idw_var,2, function (x) min(x, na.rm = TRUE))
names(mymin) = NULL
# considers date is included
for (i in 2:38){
  data_imputed_air_temp[1:mv_air_temp[i-1]-1,i] = data_imputed_air_temp[mv_air_temp[i-1],i]
  indx1 = which(data_imputed_air_temp[i] > mymax[i])
  data_imputed_air_temp[indx1,i] = mymax[i]
  indx2 = which(data_imputed_air_temp[i] < mymin[i])
  data_imputed_air_temp[indx2,i] = mymin[i]
}

######################### WS ######################### 
data_imputed_ws = myImp(data_imputed_ws,days_seq)
# QC (quality control)
data_idw_var = copy(alldata)
data_idw_var = as.data.frame(alldata[, names(alldata) %in% c("ws","code",'date')])
data_idw_var = spread(data_idw_var, code, ws)
data_idw_var = data_idw_var[,c('date',stations$code)]
data_idw_var$date = NULL
mymax = apply(data_idw_var,2, function (x) max(x, na.rm = TRUE))
names(mymax) = NULL
mymin = apply(data_idw_var,2, function (x) min(x, na.rm = TRUE))
names(mymin) = NULL
for (i in 2:38){
  data_imputed_ws[1:mv_ws[i-1]-1,i] = data_imputed_ws[mv_ws[i-1],i]
  indx1 = which(data_imputed_ws[i] > mymax[i])
  data_imputed_ws[indx1,i] = mymax[i]
  indx2 = which(data_imputed_ws[i] < mymin[i])
  data_imputed_ws[indx2,i] = mymin[i]
}
######################### WD ######################### 
data_imputed_wd = myImp(data_imputed_wd,days_seq)
# QC (quality control)
data_idw_var = copy(alldata)
data_idw_var = as.data.frame(alldata[, names(alldata) %in% c("wd","code",'date')])
data_idw_var = spread(data_idw_var, code, wd)
data_idw_var = data_idw_var[,c('date',stations$code)]
data_idw_var$date = NULL
mymax = apply(data_idw_var,2, function (x) max(x, na.rm = TRUE))
names(mymax) = NULL
mymin = apply(data_idw_var,2, function (x) min(x, na.rm = TRUE))
names(mymin) = NULL
# considers date is included
for (i in 2:38){
  data_imputed_wd[1:mv_wd[i-1]-1,i] = data_imputed_wd[mv_wd[i-1],i]
  indx1 = which(data_imputed_wd[i] > mymax[i])
  data_imputed_wd[indx1,i] = mymax[i]
  indx2 = which(data_imputed_wd[i] < mymin[i])
  data_imputed_wd[indx2,i] = mymin[i]
}

######################### VISIBILITY ######################### 
data_imputed_visibility = myImp(data_imputed_visibility,days_seq)
# QC (quality control)
data_idw_var = copy(alldata)
data_idw_var = as.data.frame(alldata[, names(alldata) %in% c("visibility","code",'date')])
data_idw_var = spread(data_idw_var, code, visibility)
data_idw_var = data_idw_var[,c('date',stations$code)]
data_idw_var$date = NULL
mymax = apply(data_idw_var,2, function (x) max(x, na.rm = TRUE))
names(mymax) = NULL
mymin = apply(data_idw_var,2, function (x) min(x, na.rm = TRUE))
names(mymin) = NULL
# considers date is included
for (i in 2:38){
  data_imputed_visibility[1:mv_visibility[i-1]-1,i] = data_imputed_visibility[mv_visibility[i-1],i]
  indx1 = which(data_imputed_visibility[i] > mymax[i])
  data_imputed_visibility[indx1,i] = mymax[i]
  indx2 = which(data_imputed_visibility[i] < mymin[i])
  data_imputed_visibility[indx2,i] = mymin[i]
}

######################### RH ######################### 
data_imputed_RH = myImp(data_imputed_RH,days_seq)
# QC (quality control)
data_idw_var = copy(alldata)
data_idw_var = as.data.frame(alldata[, names(alldata) %in% c("RH","code",'date')])
data_idw_var = spread(data_idw_var, code, RH)
data_idw_var = data_idw_var[,c('date',stations$code)]
data_idw_var$date = NULL
mymax = apply(data_idw_var,2, function (x) max(x, na.rm = TRUE))
names(mymax) = NULL
mymin = apply(data_idw_var,2, function (x) min(x, na.rm = TRUE))
names(mymin) = NULL
# considers date is included
for (i in 2:38){
  data_imputed_RH[1:mv_RH[i-1]-1,i] = data_imputed_RH[mv_RH[i-1],i]
  indx1 = which(data_imputed_RH[i] > mymax[i])
  data_imputed_RH[indx1,i] = mymax[i]
  indx2 = which(data_imputed_RH[i] < mymin[i])
  data_imputed_RH[indx2,i] = mymin[i]
}

######################### DEW_POINT ######################### 
data_imputed_dew_point = myImp(data_imputed_dew_point,days_seq)
# QC (quality control)
data_idw_var = copy(alldata)
data_idw_var = as.data.frame(alldata[, names(alldata) %in% c("dew_point","code",'date')])
data_idw_var = spread(data_idw_var, code, dew_point)
data_idw_var = data_idw_var[,c('date',stations$code)]
data_idw_var$date = NULL
mymax = apply(data_idw_var,2, function (x) max(x, na.rm = TRUE))
names(mymax) = NULL
mymin = apply(data_idw_var,2, function (x) min(x, na.rm = TRUE))
names(mymin) = NULL
# considers date is included
for (i in 2:38){
  data_imputed_dew_point[1:mv_dew_point[i-1]-1,i] = data_imputed_dew_point[mv_dew_point[i-1],i]
  indx1 = which(data_imputed_dew_point[i] > mymax[i])
  data_imputed_dew_point[indx1,i] = mymax[i]
  indx2 = which(data_imputed_dew_point[i] < mymin[i])
  data_imputed_dew_point[indx2,i] = mymin[i]
}


# Check %missings
data_air_temp_melted <- melt(data_imputed_air_temp, id=c("date"))
colnames(data_air_temp_melted) <- c('date','code','air_temp')
statsNA(data_air_temp_melted$air_temp)

data_ws_melted <- melt(data_imputed_ws, id=c("date"))
colnames(data_ws_melted) <- c('date','code','ws')
statsNA(data_ws_melted$ws)

data_wd_melted <- melt(data_imputed_wd, id=c("date"))
colnames(data_wd_melted) <- c('date','code','wd')
statsNA(data_wd_melted$wd)

data_visibility_melted <- melt(data_imputed_visibility, id=c("date"))
colnames(data_visibility_melted) <- c('date','code','visibility')
statsNA(data_visibility_melted$visibility)

data_RH_melted <- melt(data_imputed_RH, id=c("date"))
colnames(data_RH_melted) <- c('date','code','RH')
statsNA(data_RH_melted$RH)

data_dew_point_melted <- melt(data_imputed_dew_point, id=c("date"))
colnames(data_dew_point_melted) <- c('date','code','dew_point')
statsNA(data_dew_point_melted$dew_point)

# Save files
write.csv(data_imputed_air_temp,'data_imputed_air_temp.csv')
write.csv(data_imputed_ws,'data_imputed_ws.csv')
write.csv(data_imputed_wd,'data_imputed_wd.csv')
write.csv(data_imputed_visibility,'data_imputed_visibility.csv')
write.csv(data_imputed_RH,'data_imputed_RH.csv')
write.csv(data_imputed_dew_point,'data_imputed_dew_point.csv')

write.csv(data_air_temp_melted,'data_air_temp_melted.csv')
write.csv(data_ws_melted,'data_ws_melted.csv')
write.csv(data_wd_melted,'data_wd_melted.csv')
write.csv(data_visibility_melted,'data_visibility_melted.csv')
write.csv(data_RH_melted,'data_RH_melted.csv')
write.csv(data_dew_point_melted,'data_dew_point_melted.csv')


####################################################################
##################     STATISTICAL ANALYSIS    #####################
####################################################################

# Load files
data_imputed_air_temp=read.csv('temporal imputed data/MOVING AVERAGE/data_imputed_air_temp.csv')
data_imputed_ws=read.csv('temporal imputed data/MOVING AVERAGE/data_imputed_ws.csv')
data_imputed_wd=read.csv('temporal imputed data/MOVING AVERAGE/data_imputed_wd.csv')
data_imputed_visibility=read.csv('temporal imputed data/MOVING AVERAGE/data_imputed_visibility.csv')
data_imputed_RH=read.csv('temporal imputed data/MOVING AVERAGE/data_imputed_RH.csv')
data_imputed_dew_point=read.csv('temporal imputed data/MOVING AVERAGE/data_imputed_dew_point.csv')

data_air_temp_melted=read.csv('data_air_temp_melted.csv')
data_ws_melted=read.csv('data_ws_melted.csv')
data_wd_melted=read.csv('data_wd_melted.csv')
data_visibility_melted=read.csv('data_visibility_melted.csv')
data_RH_melted=read.csv('data_RH_melted.csv')
data_dew_point_melted=read.csv('data_dew_point_melted.csv')


####### STATISTICS SUMMARY

my_var = copy(data_imputed_air_temp)
my_var$date = NULL
my_var$X = NULL
statistics = as.data.frame(colSums(my_var, na.rm = T)/nrow(my_var))
colnames(statistics) = c('mean_air_temp')
statistics$sd_air_temp = apply(my_var,2, function (x) sd(x, na.rm = TRUE))
statistics$median_air_temp = apply(my_var,2,function (x) median(x, na.rm = T))
statistics$min_air_temp = apply(my_var,2,function (x) min(x, na.rm = T))
statistics$max_air_temp = apply(my_var,2,function (x) max(x, na.rm = T))

my_var = copy(data_imputed_ws)
my_var$date = NULL
my_var$X = NULL
statistics$mean_ws = colSums(my_var, na.rm = T)/nrow(my_var)
statistics$sd_ws = apply(my_var,2, function (x) sd(x, na.rm = TRUE))
statistics$median_ws = apply(my_var,2,function (x) median(x, na.rm = T))
statistics$min_ws = apply(my_var,2,function (x) min(x, na.rm = T))
statistics$max_ws = apply(my_var,2,function (x) max(x, na.rm = T))

my_var = copy(data_imputed_wd)
my_var$date = NULL
my_var$X = NULL
statistics$mean_wd = colSums(my_var, na.rm = T)/nrow(my_var)
statistics$sd_wd = apply(my_var,2, function (x) sd(x, na.rm = TRUE))
statistics$median_wd = apply(my_var,2,function (x) median(x, na.rm = T))
statistics$min_wd = apply(my_var,2,function (x) min(x, na.rm = T))
statistics$max_wd = apply(my_var,2,function (x) max(x, na.rm = T))

my_var = copy(data_imputed_visibility)
my_var$date = NULL
my_var$X = NULL
statistics$mean_visibility = colSums(my_var, na.rm = T)/nrow(my_var)
statistics$sd_visibility = apply(my_var,2, function (x) sd(x, na.rm = TRUE))
statistics$median_visibility = apply(my_var,2,function (x) median(x, na.rm = T))
statistics$min_visibility = apply(my_var,2,function (x) min(x, na.rm = T))
statistics$max_visibility = apply(my_var,2,function (x) max(x, na.rm = T))

my_var = copy(data_imputed_RH)
my_var$date = NULL
my_var$X = NULL
statistics$mean_RH = colSums(my_var, na.rm = T)/nrow(my_var)
statistics$sd_RH = apply(my_var,2, function (x) sd(x, na.rm = TRUE))
statistics$median_RH = apply(my_var,2,function (x) median(x, na.rm = T))
statistics$min_RH = apply(my_var,2,function (x) min(x, na.rm = T))
statistics$max_RH = apply(my_var,2,function (x) max(x, na.rm = T))

my_var = copy(data_imputed_dew_point)
my_var$date = NULL
my_var$X = NULL
statistics$mean_dew_point = colSums(my_var, na.rm = T)/nrow(my_var)
statistics$sd_dew_point = apply(my_var,2, function (x) sd(x, na.rm = TRUE))
statistics$median_dew_point = apply(my_var,2,function (x) median(x, na.rm = T))
statistics$min_dew_point = apply(my_var,2,function (x) min(x, na.rm = T))
statistics$max_dew_point = apply(my_var,2,function (x) max(x, na.rm = T))

write.csv(statistics,'statistics_comparison_ALL.csv')
kk = read.csv('statistics_raw.csv')



############# PLOTS
# Change data by the MOVING AVERAGE data to obtain its results
data_imputed_air_temp=read.csv('temporal imputed data/ARIMA + QUALITY CONTROL/data_imputed_air_temp.csv')
data_imputed_ws=read.csv('temporal imputed data/ARIMA + QUALITY CONTROL/data_imputed_ws.csv')
data_imputed_wd=read.csv('temporal imputed data/ARIMA + QUALITY CONTROL/data_imputed_wd.csv')
data_imputed_visibility=read.csv('temporal imputed data/ARIMA + QUALITY CONTROL/data_imputed_visibility.csv')
data_imputed_RH=read.csv('temporal imputed data/ARIMA + QUALITY CONTROL/data_imputed_RH.csv')
data_imputed_dew_point=read.csv('temporal imputed data/ARIMA + QUALITY CONTROL/data_imputed_dew_point.csv')

data_air_temp_melted <- melt(data_imputed_air_temp, id=c("date"))
data_air_temp_melted$X = NULL
colnames(data_air_temp_melted) <- c('date','code','air_temp_imp')
data_air_temp_melted = data_air_temp_melted %>% filter(code != "X")
data_air_temp_melted$code = substring(data_air_temp_melted$code, 2)

data_ws_melted <- melt(data_imputed_ws, id=c("date"))
data_ws_melted$X = NULL
colnames(data_ws_melted) <- c('date','code','ws')
data_ws_melted = data_ws_melted %>% filter(code != "X")
data_ws_melted$code = substring(data_ws_melted$code, 2)

data_wd_melted <- melt(data_imputed_wd, id=c("date"))
data_wd_melted$X = NULL
colnames(data_wd_melted) <- c('date','code','wd')
data_wd_melted = data_wd_melted %>% filter(code != "X")
data_wd_melted$code = substring(data_wd_melted$code, 2)

data_visibility_melted <- melt(data_imputed_visibility, id=c("date"))
data_visibility_melted$X = NULL
colnames(data_visibility_melted) <- c('date','code','visibility')
data_visibility_melted = data_visibility_melted %>% filter(code != "X")
data_visibility_melted$code = substring(data_visibility_melted$code, 2)

data_RH_melted <- melt(data_imputed_RH, id=c("date"))
data_RH_melted$X = NULL
colnames(data_RH_melted) <- c('date','code','RH')
data_RH_melted = data_RH_melted %>% filter(code != "X")
data_RH_melted$code = substring(data_RH_melted$code, 2)

data_dew_point_melted <- melt(data_imputed_dew_point, id=c("date"))
data_dew_point_melted$X = NULL
colnames(data_dew_point_melted) <- c('date','code','dew_point')
data_dew_point_melted = data_dew_point_melted %>% filter(code != "X")
data_dew_point_melted$code = substring(data_dew_point_melted$code, 2)



data_idw_var = copy(alldata)
data_idw_var = as.data.frame(alldata[, names(alldata) %in% c("air_temp","code",'date')])
data_idw_var = spread(data_idw_var, code, air_temp)
#data_idw_var = data_idw_var[,c('date',stations$code)]
data_idw_var <- melt(data_idw_var, id=c("date"))
colnames(data_idw_var) <- c('date','code','air_temp')
var_data = cbind(data_air_temp_melted, data_idw_var$air_temp)
colnames(var_data) <- c('date','code','air_temp_imp','air_temp')
windows()
ggplot(var_data) + 
  geom_density(aes(x=air_temp_imp, fill = 'air_temp_imp', color = 'air_temp_imp'), alpha=0.1, adjust = 2) +  
  geom_density(aes(x=air_temp, fill = 'air_temp', color = 'air_temp'), alpha=0.1, adjust = 2) + 
  facet_wrap( ~ code, ncol=6, scales = "free") + labs(x = 'air temperature (ºC)', fill = "Data:") + theme(legend.position = c(0.9, 0.05),legend.direction = "horizontal") + guides(color=FALSE)

data_idw_var = copy(alldata)
data_idw_var = as.data.frame(alldata[, names(alldata) %in% c("ws","code",'date')])
data_idw_var = spread(data_idw_var, code, ws)
data_idw_var <- melt(data_idw_var, id=c("date"))
colnames(data_idw_var) <- c('date','code','ws')
var_data = cbind(data_ws_melted, data_idw_var$ws)
colnames(var_data) <- c('date','code','ws_imp','ws')
windows()
ggplot(var_data) + 
  geom_density(aes(x=ws_imp, fill = 'ws_imp', color = 'ws_imp'), alpha=0.1, adjust = 2) +  
  geom_density(aes(x=ws, fill = 'ws', color = 'ws'), alpha=0.1, adjust = 2) + 
  facet_wrap( ~ code, ncol=6, scales = "free") + labs(x = 'wind speed (m/s)', fill = "Data:") + theme(legend.position = c(0.9, 0.05),legend.direction = "horizontal") + guides(color=FALSE)

data_idw_var = copy(alldata)
data_idw_var = as.data.frame(alldata[, names(alldata) %in% c("wd","code",'date')])
data_idw_var = spread(data_idw_var, code, wd)
data_idw_var <- melt(data_idw_var, id=c("date"))
colnames(data_idw_var) <- c('date','code','wd')
var_data = cbind(data_wd_melted, data_idw_var$wd)
colnames(var_data) <- c('date','code','wd_imp','wd')
windows()
ggplot(var_data) + 
  geom_density(aes(x=wd_imp, fill = 'wd_imp', color = 'wd_imp'), alpha=0.1, adjust = 2) +  
  geom_density(aes(x=wd, fill = 'wd', color = 'wd'), alpha=0.1, adjust = 2) + 
  facet_wrap( ~ code, ncol=6, scales = "free") + labs(x = 'wind direction (ºN)', fill = "Data:") + theme(legend.position = c(0.9, 0.05),legend.direction = "horizontal") + guides(color=FALSE)

data_idw_var = copy(alldata)
data_idw_var = as.data.frame(alldata[, names(alldata) %in% c("visibility","code",'date')])
data_idw_var = spread(data_idw_var, code, visibility)
data_idw_var <- melt(data_idw_var, id=c("date"))
colnames(data_idw_var) <- c('date','code','visibility')
var_data = cbind(data_visibility_melted, data_idw_var$visibility)
colnames(var_data) <- c('date','code','visibility_imp','visibility')
windows()
ggplot(var_data) + 
  geom_density(aes(x=visibility_imp, fill = 'visibility_imp', color = 'visibility_imp'), alpha=0.1, adjust = 2) +  
  geom_density(aes(x=visibility, fill = 'visibility', color = 'visibility'), alpha=0.1, adjust = 2) + 
  facet_wrap( ~ code, ncol=6, scales = "free") + labs(x = 'visibility (m)', fill = "Data:") + theme(legend.position = c(0.9, 0.05),legend.direction = "horizontal") + guides(color=FALSE)

data_idw_var = copy(alldata)
data_idw_var = as.data.frame(alldata[, names(alldata) %in% c("RH","code",'date')])
data_idw_var = spread(data_idw_var, code, RH)
data_idw_var <- melt(data_idw_var, id=c("date"))
colnames(data_idw_var) <- c('date','code','RH')
var_data = cbind(data_RH_melted, data_idw_var$RH)
colnames(var_data) <- c('date','code','RH_imp','RH')
windows()
ggplot(var_data) + 
  geom_density(aes(x=RH_imp, fill = 'RH_imp', color = 'RH_imp'), alpha=0.1, adjust = 2) +  
  geom_density(aes(x=RH, fill = 'RH', color = 'RH'), alpha=0.1, adjust = 2) + 
  facet_wrap( ~ code, ncol=6, scales = "free") + labs(x = 'relative humidity (%)', fill = "Data:") + theme(legend.position = c(0.9, 0.05),legend.direction = "horizontal") + guides(color=FALSE)

data_idw_var = copy(alldata)
data_idw_var = as.data.frame(alldata[, names(alldata) %in% c("dew_point","code",'date')])
data_idw_var = spread(data_idw_var, code, dew_point)
data_idw_var <- melt(data_idw_var, id=c("date"))
colnames(data_idw_var) <- c('date','code','dew_point')
var_data = cbind(data_dew_point_melted, data_idw_var$dew_point)
colnames(var_data) <- c('date','code','dew_point_imp','dew_point')
windows()
ggplot(var_data) + 
  geom_density(aes(x=dew_point_imp, fill = 'dew_point_imp', color = 'RH_imp'), alpha=0.1, adjust = 2) +  
  geom_density(aes(x=dew_point, fill = 'dew_point', color = 'RH'), alpha=0.1, adjust = 2) + 
  facet_wrap( ~ code, ncol=6, scales = "free") + labs(x = 'dew point temperature (ºC)', fill = "Data:") + theme(legend.position = c(0.9, 0.05),legend.direction = "horizontal") + guides(color=FALSE)


# Statistics for all data together
library(moments) # skewness

my_var = copy(alldata)
my_var$code = NULL
my_var$X = NULL
statistics = as.data.frame(colSums(my_var, na.rm = T)/nrow(my_var))
colnames(statistics) = c('mean_raw')
statistics$sd_raw = apply(my_var,2, function (x) sd(x, na.rm = TRUE))
statistics$skewness_raw = apply(my_var,2, function (x) skewness(x, na.rm = TRUE))
statistics$kurtosis_raw = apply(my_var,2, function (x) kurtosis(x, na.rm = TRUE))

melted_data = as.data.frame(data_ws_melted$ws)
melted_data$wd = data_wd_melted$wd
melted_data$air_temp = data_air_temp_melted$air_temp
melted_data$visibility = data_visibility_melted$visibility
melted_data$dew_point = data_dew_point_melted$dew_point
melted_data$RH = data_RH_melted$RH

statistics$mean_MA = colSums(melted_data, na.rm = T)/nrow(melted_data)
statistics$sd_MA = apply(melted_data,2, function (x) sd(x, na.rm = TRUE))
statistics$skewness_MA = apply(melted_data,2, function (x) skewness(x, na.rm = TRUE))
statistics$kurtosis_MA = apply(melted_data,2, function (x) kurtosis(x, na.rm = TRUE))





