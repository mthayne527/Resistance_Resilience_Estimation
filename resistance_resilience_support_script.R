# "Support script" for resistance_resilience_estimation.R
# Defines functions to: 
# - calculate resistance
# - calculate resilience
# - define the width of the rolling window to determine the pre-disturbance condition


# Calculate resistance
calc_resistance <- function(pre_disturbance, peak_effect){
  return(1 - (2 * abs(pre_disturbance - peak_effect)) / (pre_disturbance + abs(pre_disturbance - peak_effect)))
}

# Calculate resilience
calc_resilience <- function(pre_disturbance, peak_effect, return_level){
  
  return((2 * abs(pre_disturbance - peak_effect)) / (abs(pre_disturbance - peak_effect) + abs(pre_disturbance - return_level)) - 1)
}

# # Testing
# data_frame = df_replicate; stability_variable = "Blue"; first_day = range_predisturbance[1]; last_day = range_predisturbance[2]; minimum_width = 3

# Define the range over which the pre- or post_disturbance conditions should be calculated
defineWindow <- function(data_frame, stability_variable, first_day, last_day, minimum_width){
  data_frame = copy(data_frame)
  
  data_frame = data_frame[date >= first_day & date <= last_day, c("date", stability_variable), with=FALSE]
  
  # Try to minimize SE
  SE_rolling_window=data.table(Start=as.numeric(),End=as.numeric(),SE=as.numeric())
  
  # Start at the full width of the range
  current_range = c(1,nrow(data_frame))
  
  # Do a rolling window starting at full width and moving towards minimum_width, calculating standard error (SE) for each step
  # If the full width is shorter than the min_width, return the same date as the peak:
  if(current_range[2]-current_range[1]<=(minimum_width-1)){
    return(c(first_day,first_day))
  }else{
    
    while(current_range[2]-current_range[1]>=(minimum_width-1)){
      # Calculate SE
      SE = var(data_frame[current_range[1]:current_range[2]][[stability_variable]],na.rm = T)/nrow(data_frame[current_range[1]:current_range[2]])
      
      # Enter it in the data frame
      SE_rolling_window=rbindlist(list(SE_rolling_window,
                                       data.table(Start=current_range[1],End=current_range[2],SE=SE)))
      
      
      if(current_range[2]==nrow(data_frame)){
        # Set back to the start and make window smaller
        current_range=c(1,(current_range[2]-current_range[1]))
      }else{
        # Move the window one step forward
        current_range=current_range+1
      }
    }
    
    # Take the row with the minimum SE
    minimumCV=which(SE_rolling_window$SE == min(SE_rolling_window$SE,na.rm = T))
    
    SE_rolling_window = SE_rolling_window[minimumCV]
    # If there are multiple rows, pick the widest window. 
    if(nrow(SE_rolling_window)>1){
      SE_rolling_window$width=SE_rolling_window$End-SE_rolling_window$Start+1
      SE_rolling_window=SE_rolling_window[which(SE_rolling_window$width == max(SE_rolling_window$width,na.rm = T))]
      
      # If there are then still multiple rows with the same width, pick the one with the lowest average stability_variable
      if(nrow(SE_rolling_window)>1){
        SE_rolling_window$average=-999
        for(i in 1:nrow(SE_rolling_window)){
          SE_rolling_window$average[i]=mean(data_frame[SE_rolling_window$Start[i]:SE_rolling_window$End[i]][[stability_variable]],na.rm = T)
        }
        
        SE_rolling_window=SE_rolling_window[which(SE_rolling_window$average == min(SE_rolling_window$average,na.rm = T))]
        
        # If there are then still multiple rows, it won't matter what range you pick. Select the first one.
      }
    }
    
    return(c(data_frame$date[SE_rolling_window$Start[1]],data_frame$date[SE_rolling_window$End[1]]))
  }
  
}

fill_df_stability <- function(df, window_pre_event, window_post_event,
                              name_df, stability_variable, resistance, resilience){
  #extract disturbance characteristics
  time_series_start = df$date[1]
  time_wind_peak = df$date[which.max(df$wind)]
  time_first_half = full_df[full_df$date >= time_series_start & full_df$date <= time_wind_peak,]
  time_wind_start = time_first_half$date[which.min(time_first_half$wind)]
  time_series_start = time_wind_start - days(3)
  time_series_end = last(df$date)
  time_second_half = full_df[full_df$date >= time_wind_peak & full_df$date <= time_series_end,]
  time_wind_end = time_second_half$date[which.min(time_second_half$wind)]
  time_series_end = time_wind_end + days(3)
  disturb_chars = full_df[full_df$date >= time_wind_start & full_df$date <= time_wind_end,]#dataframe to extract disturbance characteristics
  
  #option to calculate difference in SD in response variables between pre and post disturbance conditions
  df_pre = df[df$date >= window_pre_event[1] & df$date <= window_pre_event[2],2:ncol(df)]
  stdev_pre = sd(df_pre[,colnames(df_pre)==stability_variable], na.rm = T)
  
  df_post = df[df$date >= window_post_event[1] & df$date <= window_post_event[2],2:ncol(df)]
  stdev_post = sd(df_post[,colnames(df_post)==stability_variable], na.rm = T)
  
  #create dataframe for filling in resistance and resilience values
  stability_df = t(data.frame(colMeans(df_pre)))
  colnames(stability_df) = colnames(df)[2:ncol(df)]
  rownames(stability_df) = NULL
  
  stability_df = data.frame(stability_df)
  
  stability_df$Event = name_df
  stability_df$Stability_variable = stability_variable
  stability_df$Resistance = resistance
  stability_df$Resilience = resilience
  
  stability_df = stability_df[, c((ncol(stability_df)-3):ncol(stability_df),1:(ncol(stability_df)-4))]
  
  stability_df$stdev_pre = stdev_pre
  stability_df$stdev_post = stdev_post
  
  ###option to extract disturbance characteristics
  # stability_df$duration = nrow(disturb_chars)
  # #stability_df$acc_rain = sum(disturb_chars$rain)
  # #stability_df$avg_rain = mean(disturb_chars$rain)
  # stability_df$avg_wind = mean(storm_chars$wind)
  
  
  return(stability_df)
}

###BRT modelling for helping identify the peak in the response variable
call_BRT <- function(df, stability_variable, printLR = FALSE){
  #Learning rate for boosted regression
  LearningRate<-c(0.8192,0.4096,0.2048,0.1024,0.0512,0.0256,0.0128,0.0064,0.0032,0.0016,0.0008,0.0004,0.0002,0.0001,0.000005,0.0000025,0.00000175,0.000000875,0.0000004375,0.0000001,0.00000001,0.000000001,0.0000000001)
  
  df = df[, colSums(is.na(df)) != nrow(df)]#Remove columns with only NA's
  df = df[complete.cases(df), ]#cannot have NA in response column
  allColsWithVariables = 2:ncol(df)
  colNrStabilityVariable = which(colnames(df) == stability_variable) # Select column with the stability variable
  otherVariables = allColsWithVariables[allColsWithVariables != colNrStabilityVariable] # These are all other columns (and excluding the first column)
  
  #For loop for running boosted regression to determine response in the stability variable for df
  for(lr in LearningRate){
    tryCatch({
      if(printLR){print(lr)} 
      
      #fit the brt which explains resistance variability
      set.seed(123)
      BRT_result <- gbm.step(
        data=df,
        gbm.x = otherVariables,
        gbm.y = colNrStabilityVariable,
        family = "gaussian", 
        tree.complexity = 2,
        learning.rate = lr,
        max.trees=10000,
        bag.fraction = c(0.5,0.6,0.7),
        silent=T,
        plot.main=F,
        plot.folds=F)  
      
      if(is.null(BRT_result)==FALSE){
        if(BRT_result$n.trees>1000) break}
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  # Show plot
  nplots = length(otherVariables)
  gbm.plot(BRT_result, 
           n.plots = nplots,smooth = T,
           plot.layout = c(3,3), 
           common.scale = F)
  
}

###Make resistance and resilience quantification plots
make_plot <- function(df, name_df, stability_variable, 
                      pre_event_level, sd_pre_event_level, window_pre_event,
                      timeOfPeak, peak_effect, 
                      post_event_level, window_post_event,
                      resistance, resilience){
  
  p = ggplot(data = df)+
    geom_line(aes_string(x="date", y=stability_variable))+
    geom_hline(aes(yintercept = pre_event_level + 2*sd_pre_event_level), linetype = "dotted", size = 1)+
    geom_hline(aes(yintercept = pre_event_level - 2*sd_pre_event_level), linetype = "dotted", size = 1)+
    geom_segment(aes(x = window_pre_event[1], y = pre_event_level, xend = window_pre_event[2], yend = pre_event_level,colour="Pre-event"), size=3)+
    geom_segment(aes(x = timeOfPeak-hours(4), y = peak_effect, xend = timeOfPeak+hours(4), yend = peak_effect,colour="Peak Effect"), size=3)+
    geom_segment(aes(x = window_post_event[1], y = post_event_level, xend = window_post_event[2], yend = post_event_level, colour = "Post-event"), size=3)+
    scale_colour_manual(name="Legend",
                        values = c(`Pre-event` = "Red", `Peak Effect` = "Blue", `Post-event` = "Green"))+
    labs(title = name_df,
         subtitle = paste0("Resistance: ",round(resistance,3),";    Resilience: ",round(resilience,3)))+
    theme_light()
  
  print(p)
  
}