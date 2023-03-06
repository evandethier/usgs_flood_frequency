#### IMPORT PACKAGES ####
library(data.table)
library(dataRetrieval)
library(ggplot2)

#### IMPORT K FACTOR TABLE, SITE SELECTION TABLES ####    
# Import K factor table for lookup value related to skewness, recurrence interval
# From Haan, 1977
Kfactor_lookup <- fread(paste0('usgs_flood_frequency_imports/', 'Recurrence_interval_Kfactor_lookup.csv'))

# Import HCDN station table
hcdn_simple <- fread(paste0('usgs_flood_frequency_imports/','hcdn_simple.csv'))
hcdn_simple <- hcdn_simple[,':='(site_no = ifelse(site_no < 1e7, paste0('0',site_no), paste0(site_no)))]
<<<<<<< HEAD



# # USGS stations with POC data
#   setwd("/Users/evandethier/Documents/Dartmouth/IreneProject/IreneInsideWork/Irene_landsat/earthengine/earthengineDams/usa-hydrology/")
#   station_df_poc <- read_csv('annual_p00689_q_avg_station_info.csv')

# Of imported options, select table to find recurrence intervals for
station_df_1 <- setDT(station_df_hcdn)
station_info <- readNWISsite(station_df_1$site_no)

setwd("/Users/evandethier/Documents/Dartmouth/IreneProject/IreneInsideWork/Irene_landsat/earthengine/earthengineDams/usa-hydrology/")

# Select recurrence interval columns
recur_cols <- colnames(Kfactor_lookup)[-1]
# Add new columns for discharge corresponding to given recurrence interval
station_df_2 <- setDT(station_info)

# Loop through sites with concavity data (not all sites have discharge data)
for(i in 1:nrow(station_df_2)){
  site_sel <- station_df_2$site_no[i] # select site for analysis
  Qpeak_list_sel <- readNWISpeak(site_sel) # download peak flow data for given site
  # Download average flow data for given site
  Qavg_list_sel <- readNWISstat(siteNumbers = site_sel, parameterCd = '00060', statReportType="annual") 
  Qavg_annual <- NA
  if(!is.null(Qavg_list_sel)){
    Qavg_annual <- mean(Qavg_list_sel$mean_va, na.rm = T) * 0.0283168}
  
  if(nrow(Qpeak_list_sel)>10){ # only proceed with analysis if > 10-year peak flow record exists
    # some sites just have gage height (no peak flow) or have peaks of 0. find those.
    Qpeak_list_sel <- setDT(Qpeak_list_sel)[!is.na(peak_va) & peak_va > 0]
    if(nrow(Qpeak_list_sel)>10){ # reject sites with only gage height
      Qpeak_list_sel[,Q_cms:= Qpeak_list_sel$peak_va * 0.0283168] # convert Q to cms
      
      # Log Pearson III calculation columns
      Qpeak_list_sel <- Qpeak_list_sel[,
                                       log_Q_cms := log(Qpeak_list_sel$Q_cms)][,
                                                                               rank := rank(Q_cms)][, # rank among peak flows
                                                                                                    return_period := (nrow(Qpeak_list_sel)+1)/rank][,
                                                                                                                                                    exceed_prob := 1/return_period]
      Qpeak_var <- var(Qpeak_list_sel$log_Q_cms) # variance
      Qpeak_skew <- round(skewness(Qpeak_list_sel$log_Q_cms), 1) # skew
      Qpeak_mean <- mean(Qpeak_list_sel$log_Q_cms) # log mean
      
      if(Qpeak_skew > 3){ # some skews are outside K table range. Identify and replace with max table value.
        Qpeak_skew <- 3
      } else if(Qpeak_skew < -3){
        Qpeak_skew <- -3
      }
      # Look up K value from K factor table 
      K_factor_sel <- Kfactor_lookup[which(Qpeak_skew == Kfactor_lookup$Cw), ..recur_cols]
      
      # Calculate Q for recurrence intervals (1, 2, 5, 10, 25, 50, 100, 200 yrs)
      Recurrence_Qs <- exp(Qpeak_mean + K_factor_sel * sqrt(Qpeak_var))
      # Add Q for selected recurrence intervals to row of dataframe with that data
      station_df_2[i,c(recur_cols) := Recurrence_Qs][i, # add recurrence data
                                                     Qavg := Qavg_annual # add annual average Q (cms)
      ]
      # [Optional] print calculated line for checking
      print(cbind(station_df_2[i, c('site_no', 'station_nm')], Recurrence_Qs, Qavg_annual))
    }}}

# Convert drainage area, mi2 to km2
station_df_2[,drainage_area_km2 := drain_area_va * 2.58999]

# Write file to disk
fwrite(station_df_2, file = 'hcdn_rivers_recurrence.csv')

# Select column for plotting
Q_recur_plots <- list()
for(i in 1:4){
  # Select plot recurrence interval
  col_sel_plot <- c('2 yr', '5 yr', '10 yr', '25 yr')[i]
  col_sel_save <- c('2yr', '5yr', '10yr', '25yr')[i]
  Q_recur_runoff_plot <- ggplot(na.omit(station_df_2, cols = col_sel_plot),
                                # %>% subset(substr(huc_cd, start = 1, stop = 2) == '01'), 
                                aes(x = dec_long_va, y = dec_lat_va, 
                                    color = !!sym(col_sel_plot)/drainage_area_km2 * 86.4)) + 
    geom_map(data = us_ca, map = us_ca, 
             aes(map_id = id, x = long, y = lat),
             color = "grey30", fill = 'white', lwd = 0.25) +
    geom_point(size = 1.5) + 
    theme_evan +
    scale_x_continuous(limits = c(-160,-66)) + 
    scale_y_continuous(limits = c(25,65)) + 
    scale_color_gradientn(limits = c(0,200), colors = c('#D1B660', '#049CBF','#0468BF'), oob = squish) +
    theme(legend.position = 'right') +
    labs(
      x = 'Longitude',
      y = 'Latitude', 
      # color = expression(paste('Runoff (m'^'3'*'/s / km'^'2'*')')),
      color = expression(paste('Runoff (mm d'^'-1'*')')),
      title = paste0(col_sel_plot, ' recurrence runoff')
    ) 
  
  # Save plot
  ggsave(Q_recur_runoff_plot, filename = paste0('HCDN_runoff_', col_sel_save,'_recurrence.pdf'), 
         width = 7, height = 5, useDingbats = F)
  Q_recur_plots[[i]] <- Q_recur_runoff_plot
}

# Save aggregate plot of multiple recurrences
ggsave(ggarrange(plotlist = Q_recur_plots, align = 'hv', labels = c('A','B','C','D'),
                 common.legend = T), 
       width = 8, height = 6,
       filename = 'HCDN_runoff_recurrence_combined.pdf', useDingbats = F, onefile = F)


=======
>>>>>>> fb711bc724153e27e79568e81588e2dcc91d4410
