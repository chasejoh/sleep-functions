# read.ag = MOCA function used to read and clean Actigraph data (alternative to read_ag in MOCAfunctions)
read.ag = function(filepath, ENMO_calibrate = F, calibration_file = F){
  
  check_data = fread(filepath,header = T,skip = 10, nrows = 1)[1,]
  if(is.numeric(check_data[1,])){
    file_data = fread(filepath,header = F, skip = 10, stringsAsFactors = F)
  } else {
    file_data = fread(filepath,header = T, skip = 10, stringsAsFactors = F)
  }
  
  ag_header = read.csv(filepath,header = F,stringsAsFactors = F, nrows = 10)
  device_serial = str_split(ag_header[2,],'Number: ')[[1]][2]
  start = str_split(ag_header[3,],'Time ')[[1]][2]
  date = str_split(ag_header[4,],'Date ')[[1]][2]
  frequency = as.numeric(str_split(str_split(ag_header[1,],'at ')[[1]][3],' Hz')[[1]][1])
  
  epoch_full = str_split(ag_header[5,],' 00:')[[1]][2]
  epoch_temp = str_split(epoch_full,':')
  epoch = (as.numeric(epoch_temp[[1]][1])*60) + (as.numeric(epoch_temp[[1]][2]))
  
  file_length = nrow(file_data)
  date_time_start = mdy_hms(paste(date,start, sep = ' '))
  
  if(epoch < 1){
    # For Raw data
    file_data = file_data %>% mutate(`Accelerometer X` = as.numeric(`Accelerometer X`),
                                     `Accelerometer Y` = as.numeric(`Accelerometer Y`),
                                     `Accelerometer Z` = as.numeric(`Accelerometer Z`))
    
    date_time_end = date_time_start + (file_length/frequency)
    Timestamp = seq(from = date_time_start,to = date_time_end, by = 1/frequency)
    Timestamp = Timestamp[1:length(Timestamp)-1]
    
    file_data = mutate(file_data, 
                       VM = sqrt(`Accelerometer X`^2 + `Accelerometer Y`^2 + `Accelerometer Z`^2), 
                       VMcorrG = abs(sqrt(`Accelerometer X`^2 + `Accelerometer Y`^2 + `Accelerometer Z`^2)-1))
    
    if(ENMO_calibrate == T){
      C = g.calibrate(filepath,use.temp = F, printsummary=F)
      
      if(C$offset[1] == 0 & C$scale[1] == 1){
        
        device_serial_index = str_which(calibration_file$Serial,device_serial)
        
        if(length(device_serial_index) == 0){
          
        } else {
          C$offset[1] = calibration_file$Offset_X[device_serial_index]
          C$offset[2] = calibration_file$Offset_Y[device_serial_index]
          C$offset[3] = calibration_file$Offset_Z[device_serial_index]
          
          C$scale[1] = calibration_file$Scale_X[device_serial_index]
          C$scale[2] = calibration_file$Scale_Y[device_serial_index]
          C$scale[3] = calibration_file$Scale_Z[device_serial_index]
          
        }
      }
      
      file_data = mutate(file_data, 
                         calX = `Accelerometer X`*C$scale[1] + C$offset[1], 
                         calY = `Accelerometer Y`*C$scale[2] + C$offset[2], 
                         calZ = `Accelerometer Z`*C$scale[3] + C$offset[3],
                         ENMO = sqrt(calX^2 + calY^2 + calZ^2)-1)
      
      file_data = mutate(file_data, ENMO = ifelse(ENMO < 0,0,ENMO))
      
    }
    
  } else {
    # For Count data
    date_time_end = date_time_start + (file_length*epoch)
    Timestamp = seq(from = date_time_start, to = date_time_end, by = epoch)
    Timestamp = Timestamp[1:length(Timestamp)-1]
    # If Count data does not have column names when being read in 
    
    file_data = mutate(file_data, VM = sqrt(Axis1^2 + Axis2^2 + Axis3^2))
  }
  
  Dates = as.character(date(Timestamp))
  Hour = as.character(hour(Timestamp)) %>% str_pad(width = 2, side = 'left',pad = '0')
  Minute = as.character(minute(Timestamp)) %>% str_pad(width = 2, side = 'left',pad = '0')
  Second = as.character(second(Timestamp)) %>% str_pad(width = 2, side = 'left',pad = '0')
  
  ag_data = cbind(Dates,paste(Hour,Minute,Second,sep = ':'),file_data, stringsAsFactors = F)
  
  if(epoch<1){
    if(ENMO_calibrate == T){
      colnames(ag_data) = c('Date','Time','AxisX','AxisY','AxisZ', 'VM', 'VMcorrG', 'CalibratedX','CalibratedY','CalibratedZ','ENMO')
    } else {
      colnames(ag_data) =c('Date','Time','Axis1','Axis2','Axis3','Steps', 'Lux', 'Inclinometer_Off', 'Inclinometer_Stand', 'Inclinometer_Sit', 'Inclinometer_Lying', 'VM')
    }
  }else{
    colnames(ag_data) = c('Date','Time','Axis1','Axis2','Axis3','Steps', 'Lux', 'Inclinometer_Off', 'Inclinometer_Stand', 'Inclinometer_Sit', 'Inclinometer_Lying', 'VM')
  }
  
  return(ag_data)
}

# ag_epochr = MOCA function used to change epoch length
ag_epochr = function(ag_data_1sec,epoch = 60){
  rows = nrow(ag_data_1sec)
  
  new_rows = floor(rows/epoch)
  
  new_rows_index = (seq(1:new_rows)*epoch)+1
  new_rows_index = c(1,new_rows_index)
  new_rows_index = new_rows_index[-length(new_rows_index)]
  
  Date = ag_data_1sec$Date[new_rows_index]
  Time = ag_data_1sec$Time[new_rows_index]
  count_data = data.frame(Axis1 = rep(0,new_rows), Axis2 = rep(0,new_rows), Axis3 = rep(0,new_rows))
  
  for (i in 1:new_rows){
    
    if(i+1 <= new_rows){
      count_data$Axis1[i] = sum(ag_data_1sec$Axis1[new_rows_index[i]:(new_rows_index[i+1]-1)]) 
      count_data$Axis2[i] = sum(ag_data_1sec$Axis2[new_rows_index[i]:(new_rows_index[i+1]-1)]) 
      count_data$Axis3[i] = sum(ag_data_1sec$Axis3[new_rows_index[i]:(new_rows_index[i+1]-1)]) 
      
    } else {
      count_data$Axis1[i] = sum(ag_data_1sec$Axis1[new_rows_index[i]:(new_rows_index[i]+epoch-1)]) 
      count_data$Axis2[i] = sum(ag_data_1sec$Axis2[new_rows_index[i]:(new_rows_index[i]+epoch-1)]) 
      count_data$Axis3[i] = sum(ag_data_1sec$Axis3[new_rows_index[i]:(new_rows_index[i]+epoch-1)])
      
      # Handles excess if not exact 1 hour length
      # if((new_rows_index[i]+epoch-1) < rows){
      #   excess$Axis1 = sum(ag_data_1sec$Axis1[(new_rows_index[i]+epoch-1):rows])
      #   excess$Axis2 = sum(ag_data_1sec$Axis2[(new_rows_index[i]+epoch-1):rows])
      #   excess$Axis3 = sum(ag_data_1sec$Axis3[(new_rows_index[i]+epoch-1):rows])
      # }
    }
    
  }
  
  epoch_data = data.frame(Date = Date, Time = Time, count_data)
  epoch_data$VM = sqrt(epoch_data$Axis1^2 + epoch_data$Axis2^2 + epoch_data$Axis3^2)
  
  return(epoch_data)
}
