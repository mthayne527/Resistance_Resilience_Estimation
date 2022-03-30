# Script to estimate resistance and resilience from disturbances
while (!is.null(dev.list())){dev.off()}
# rm(list=ls())
lake = "lake_name"
#suggested packages
library(data.table)
library(ggplot2)
library(lubridate)
library(dismo)
library(gbm)
library(plyr)
library(dplyr)
library(stringr)


### Set your working directory to the location, with the sample.csv file and the resistance_resilience_support_script.R file
setwd("working_document_destination")

### Load the script with the support functions
source("SupportScript_File.R")

### Download data frame
allFiles = list.files(pattern = "*.csv")
allFiles <- str_extract(allFiles, '.*(?=\\.csv)')
allFiles = as.vector(allFiles)

### Create data frame with date, response variables, and predictors for BRT modelling
full_df = full_df[,c("variable_1","variable_2"),]


### Define the variables to determine resistance and resilience for (Important: these should also appear in full_df)
stability_variables = c("varibale_1", "variable_2")

### Do you want to start from the beginning, or pick of where you left before

# If startFromBeginning is FALSE, where do you want to start? (Don't make mistakes in the names here)
start_file = "first_disturbance_name"#for example wind disturbances

# Do you want to start from the beginning or start with the named start file above
startFromBeginning = F

### Make an empty data frame for filling resistance and resilience values
df_stability = NULL

### Settings for estimation of resistance resilience
min_width = 72 # Minimum width of the rolling recovery window (number of measurements)

### For-loop 

stop=F
rownumber=0

# Test:
#i = "df_event_2" ; j = "o2_adj"
####################################
#Stability Function#################
####################################
for(i in allFiles){
  
  # If you don't want to start at the beginning of allFiles, this lets you skip loops until you reach the event you want to start with
  if(!startFromBeginning){
    if(i == start_file){
      startFromBeginning = T
    }else{
      next
    }
  }
  # Load data
  df = get(i)
  
  for(j in stability_variables){
    
    rownumber = rownumber + 1
    
    # If the stability variable is not in the dataframe, skip
    if(!(j %in% colnames(df))){
      next
    }
    ### Step 1: Define disturbance charachteristics and pre-event conditions
    df = data.frame(df)
    
    time_series_start = df$date[1]
    time_wind_peak = df$date[which.max(df$wind)]
    time_first_half = full_df[full_df$date >= time_series_start & full_df$date <= time_wind_peak,]
    time_wind_start = time_first_half$date[which.min(time_first_half$wind)]
    time_series_start = time_wind_start - days(3)
    time_series_end = last(df$date)
    time_second_half = full_df[full_df$date >= time_wind_peak & full_df$date <= time_series_end,]
    time_wind_end = time_second_half$date[which.min(time_second_half$wind)]
    time_series_end = time_wind_end + days(3)
    df = data.frame(full_df[full_df$date >= time_series_start & full_df$date <= time_series_end,])
    
    df = df[, colSums(is.na(df)) != nrow(df)]#Remove columns with only NA's
    
    #df$variable[df$variable < value] <- NA #Option for removing outliers
    
    window_pre_event = c(time_series_start, time_wind_start)
    
    pre_event_level = mean(df[df$date >= window_pre_event[1] & df$date <= window_pre_event[2], j], na.rm = T)
    
    sd_pre_event_level = sd(df[df$date >= window_pre_event[1] & df$date <= window_pre_event[2], j], na.rm = T)
    
    
    ### Step 2: Determine peak effect
    
    # Find peak in the first number of hours after the peak in shear stress (peak in shear stress occurs 3 days after start dataframe)
    window_peak_event = c(df$date[1] + days(3), df$date[nrow(df)])
    
    max_variable = max(df[df$date >= window_peak_event[1] & df$date <= window_peak_event[2], j], na.rm = T)
    min_variable = min(df[df$date >= window_peak_event[1] & df$date <= window_peak_event[2], j], na.rm = T)
    # If the absolute difference between pre_event and max_variable is largest, max_variable is peak_effect
    # If the absolute difference with min_variable  is largest, min_variable is peak_effect
    if(abs(max_variable-pre_event_level) >=  abs(min_variable-pre_event_level)){
      peak_effect = max_variable
    }else{
      peak_effect = min_variable
    }
    
    # Identify the date/time at which the peak occurs. Second line is to resolve any draws.
    timeOfPeak = df$date[which(df[,j] == peak_effect)]
    timeOfPeak = timeOfPeak[length(timeOfPeak)]
    
    ### Step 3: Determine return level
    
    range_post_event = c(timeOfPeak, df$date[nrow(df)])
    
    # Define window of post-event conditions
    df=data.table(df)
    
    window_post_event= defineWindow(df,j , range_post_event[1], range_post_event[2], minimum_width = min_width)
    
    df=data.frame(df)
    
    post_event_level = mean(df[df$date >= window_post_event[1] & df$date <= window_post_event[2], j], na.rm = T)
    
    ### Step 4: Calculate resistance and resilience
    resistance = calc_resistance(pre_event_level, peak_effect)
    resilience = calc_resilience(pre_event_level, peak_effect, post_event_level)
    
    ### Step 5: Plot the result
    make_plot(df, i, j, pre_event_level, sd_pre_event_level, window_pre_event,
              timeOfPeak, peak_effect, post_event_level, window_post_event, resistance, resilience)
    
    
    ### Step 6: Ask user for approval. Register resistance, resilience and pre-event conditions in data frame
    
    outputOK = F
    
    while(outputOK==F){
      userInput = readline("Is the analysis for this perturbation correct? (y/n). Type 'stop' to end the analysis. ")
      
      if(tolower(userInput) %in% c("y","yes")){
        outputOK=T
        ggsave(paste0("file_destination", i, " - ", j, ".png"))
        
        if(is.null(df_stability)){
          df_stability = fill_df_stability(df, window_pre_event, window_post_event,
                                           i, j, resistance, resilience)
        }else{
          df_stability = rbindlist(list(df_stability, fill_df_stability(df, window_pre_event, window_post_event,
                                                                        i, j, resistance, resilience)),
                                   use.names = T, fill=T, idcol = NULL)
        }
        
      }else if(tolower(userInput) %in% c("n","no")){
        
        # Manually set pre-, peak, and post- conditions.
        
        ### Pre-event
        userInput = readline("Are you happy about the pre-event window? (y/n) ")
        if(tolower(userInput) %in% c("n","no")){
          userInput = readline("What time is the start of the pre-event condition? Format yyyy-mm-dd HH:MM:SS ")
          
          window_pre_event[1] = as.POSIXct(userInput)
          
          userInput = readline("What time is the end of the pre-event condition? Format yyyy-mm-dd HH:MM:SS ")
          
          window_pre_event[2] = as.POSIXct(userInput)
          
        }
        
        ### Peak effect
        userInput = readline("Are you happy about location of the peak effect? (y/n) ")
        if(tolower(userInput) %in% c("n","no")){
          
          userInput = readline("Do you want to see the results of a boosted regression tree plot? (y/n) ")
          
          if(tolower(userInput) %in% c("y", "yes")){
            call_BRT(df, j)
            
            userInput = readline("Is the main response positive or negative? (positive/negative/x) ")
            
            if(tolower(userInput) == "positive"){
              peak_effect = max_variable
              
              timeOfPeak = df$date[which(df[,j] == peak_effect)]
              timeOfPeak = timeOfPeak[length(timeOfPeak)]
              
            }else if(tolower(userInput) == "negative"){
              peak_effect = min_variable
              
              timeOfPeak = df$date[which(df[,j] == peak_effect)]
              timeOfPeak = timeOfPeak[length(timeOfPeak)]
              
            }else{
              # Nothing
            }
          }
          
          make_plot(df, i, j, pre_event_level, sd_pre_event_level, window_pre_event,
                    timeOfPeak, peak_effect, post_event_level, window_post_event, resistance, resilience)
          
          userInput = readline("Do you want to specify window on where to find the peak? (y/n) ")
          
          if(tolower(userInput) %in% c("y", "yes")){
            
            userInput = readline("When does the window start? yyyy-mm-dd HH:MM:SS ")
            
            window_peak_event[1] = as.POSIXct(userInput)
            
            userInput = readline("When does the window end? yyyy-mm-dd HH:MM:SS ")
            
            window_peak_event[2] = as.POSIXct(userInput)
            
            max_variable = max(df[df$date >= window_peak_event[1] & df$date <= window_peak_event[2], j], na.rm = T)
            min_variable = min(df[df$date >= window_peak_event[1] & df$date <= window_peak_event[2], j], na.rm = T)
            
            # If the absolute difference between pre_event and max_variable is largest, max_variable is peak_effect
            # If the absolute difference with min_variable  is largest, min_variable is peak_effect
            if(abs(max_variable-pre_event_level) >=  abs(min_variable-pre_event_level)){
              peak_effect = max_variable
            }else{
              peak_effect = min_variable
            }
            
            # Identify the date/time at which the peak occurs. Second line is to resolve any draws.
            timeOfPeak = df$date[which(df[,j] == peak_effect)]
            timeOfPeak = timeOfPeak[length(timeOfPeak)]
            
            make_plot(df, i, j, pre_event_level, sd_pre_event_level, window_pre_event,
                      timeOfPeak, peak_effect, post_event_level, window_post_event, resistance, resilience)
            
          }
          
        }
        
        ### Post-event
        userInput = readline("Are you happy about the post-event window? (y/n) ")
        if(tolower(userInput) %in% c("n","no")){
          
          userInput = readline("Do you want to extend the post-event window? (y/n) ")
          
          if(tolower(userInput) %in% c("y", "yes")){
            userInput = readline("With how many hours do you want to extend the window? ")
            
            # If the user entered a negative number: cut part of the data frame
            if(as.numeric(userInput) <=0){
              cutOff = df$date[nrow(df)] - hours(abs(as.numeric(userInput)))
              
              df = df[df$date <= cutOff,]
            }else{
              # If the user entered a positive number: add part of the full data frame to df
              
              rangeToAdd = c(df$date[nrow(df)], df$date[nrow(df)] + hours(as.numeric(userInput)))
              
              full_df_toAdd = full_df[full_df$date > rangeToAdd[1] & full_df$date <= rangeToAdd[2],]
              
              temp_df=rbind.fill(df, full_df_toAdd)
              df = temp_df[,colnames(temp_df) %in% colnames(df)]
              rm(temp_df)
              
            }
            
            range_post_event = c(timeOfPeak, df$date[nrow(df)])
            
            # Define window of post-event conditions
            df=data.table(df)
            
            window_post_event= defineWindow(df,j , range_post_event[1], range_post_event[2], minimum_width = min_width)
            
            df=data.frame(df)
            
          }else{
            userInput = readline("What time is the start of the post-event condition? Format yyyy-mm-dd HH:MM:SS ")
            
            window_post_event[1] = as.POSIXct(userInput)
            
            userInput = readline("What time is the end of the post-event condition? Format yyyy-mm-dd HH:MM:SS ")
            
            window_post_event[2] = as.POSIXct(userInput)
          }
          
        }
        
        # Pre-event level
        pre_event_level = mean(df[df$date >= window_pre_event[1] & df$date <= window_pre_event[2], j], na.rm = T)
        
        sd_pre_event_level = sd(df[df$date >= window_pre_event[1] & df$date <= window_pre_event[2], j], na.rm = T)
        
        #  Post-event level
        post_event_level = mean(df[df$date >= window_post_event[1] & df$date <= window_post_event[2],j], na.rm = T)
        
        # Resistance and resilience
        resistance = calc_resistance(pre_event_level, peak_effect)
        resilience = calc_resilience(pre_event_level, peak_effect, post_event_level)
        
        make_plot(df, i, j, pre_event_level, sd_pre_event_level, window_pre_event,
                  timeOfPeak, peak_effect, post_event_level, window_post_event, resistance, resilience)
        
        userInput = readline("Is it good now? (y/n) ")
        
        if(tolower(userInput) %in% c("y","yes")){
          outputOK=T
          
          ggsave(paste0("file_destination", i, " - ", j, ".png"))
          
          if(is.null(df_stability)){
            df_stability = fill_df_stability(df, window_pre_event, window_post_event,
                                             i, j, resistance, resilience)
          }else{
            df_stability = rbindlist(list(df_stability, fill_df_stability(df, window_pre_event, window_post_event,
                                                                          i, j, resistance, resilience)),
                                     use.names = T, fill=T, idcol = NULL)
          }
          
        }else{
          print("Let's try again...")
        }
        
        
        
      }else if(tolower(userInput) == "stop"){
        outputOK=T
        stop=TRUE
        print(paste("Exit for-loop, stop analysis. Stopped at:",i,j))
      }else{
        print("Answer not y/n/stop. Type again. ")
      }
      
    }
    if(stop){break}
  }
  if(stop){break}
  ggsave(paste0("file_destination", i, " - ", j, ".png"))
  
}