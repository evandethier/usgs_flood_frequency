#### IMPORT PACKAGES ####
library(data.table)
library(dataRetrieval)
library(ggplot2)
library(lubridate)
library(e1071)
library(smwrBase)

#### IMPORT K FACTOR TABLE, SITE SELECTION TABLES ####    
# Import K factor table for look-up value related to skewness, recurrence interval
# From Haan, 1977
Kfactor_lookup <- fread(paste0('usgs_flood_frequency_imports/', 'Recurrence_interval_Kfactor_lookup.csv'))
recur_cols <- colnames(Kfactor_lookup)[-1]
# Import HCDN station table
hcdn_simple <- fread(paste0('usgs_flood_frequency_imports/','hcdn_simple.csv'))
hcdn_simple <- hcdn_simple[,':='(site_no = ifelse(site_no < 1e7, paste0('0',site_no), paste0(site_no)))]


#### DOWNLOAD PEAK FLOW DATA ####
# Select site
# Androscoggin: 
site_sel <- hcdn_simple[grepl(casefold('WHITE RIVER AT WEST HARTFORD, VT'), casefold(station_nm))]
site_no_sel <- site_sel[,site_no]

# Download peak flow data from USGS for selected site
Q_peak_sel_download <- data.table(readNWISpeak(site_no_sel)) # download peak flow data for given site

# Rename columns
Q_peak_sel <- Q_peak_sel_download[,':='(Q_cms = peak_va * 0.0283168,
                             gage_height_m = gage_ht/3.28,
                             month = month(peak_dt),
                             water_year = calcWaterYear(peak_dt))]

# Select columns
Q_peak_sel <- Q_peak_sel[,.(agency_cd, site_no, peak_dt, water_year, month, Q_cms, gage_height_m)]
# Plot data
ggplot(Q_peak_sel, aes(x = water_year, y = Q_cms)) + 
  geom_point() +
  labs(
    x = 'Year',
    y = 'Annual peak discharge (m3/s)'
  )

# Plot distribution
ggplot(Q_peak_sel, aes(x = Q_cms)) + 
  geom_histogram(bins = 20) +
  labs(
    x = 'Year',
    y = 'Annual peak discharge (m3/s)'
  ) +
  scale_x_log10()

#### FLOOD FREQUENCY ANALYSIS ####
# Remove NA values
Q_peak_sel <- Q_peak_sel[!is.na(Q_cms)]

#### SIMPLE RECURRENCE INTERVAL CALCULATION ####
# Log-transformed discharge, rank among peak flows
Q_peak_sel <- Q_peak_sel[,':='(log_Q_cms = log(Q_cms),
                               rank = rank(-Q_cms))]

# Count number of years in the record
years_of_record <- nrow(Q_peak_sel)

# Plot flow by rank
ggplot(Q_peak_sel, aes(x = rank, y = Q_cms)) + 
  geom_point() +
  labs(
    x = 'Rank',
    y = 'Annual peak discharge (m3/s)'
  )

# Calculate return period (also called recurrence interval)
# and exceedence probability: probability that discharge of given size will be exceeded in a year
Q_peak_sel <- Q_peak_sel[,':='(return_period = (years_of_record + 1)/rank)]
Q_peak_sel <- Q_peak_sel[,':='(exceed_prob = 1/return_period)]

# Plot exceedence probability for each discharge
ggplot(Q_peak_sel, aes(x = return_period, y = Q_cms)) + 
  geom_point() +
  labs(
    x = 'Return period (years)',
    y = 'Annual peak discharge (m3/s)'
  ) +
  scale_x_log10() +
  scale_y_log10()


#### LOG PEARSON III RECURRENCE INTERVAL CALCULATION ####
Q_peak_sel <- Q_peak_sel[,':='(return_period = (nrow(peakQ)+1)/rank)]

# Make vector of log-transformed discharge
peakQ_ln <- Q_peak_sel[,log_Q_cms]
# Take mean, std. deviation, and skewness of log-transformed discharge
peakQ_ln_avg <- mean(peakQ_ln, na.rm = T)
peakQ_ln_sd <- sd(peakQ_ln, na.rm = T)
peakQ_ln_skew <- skewness(peakQ_ln, na.rm = T)

# Make a vector of recurrences to calculate
recurrences_to_calculate <- seq(0.01,0.999,0.001)

# Log Pearson III -- recurrence calculation
logPearsonIII <- qlpearsonIII(recurrences_to_calculate, peakQ_ln_avg, peakQ_ln_sd, peakQ_ln_skew)
logPearsonIII_qpeaks <- data.table(probabilities = recurrences_to_calculate,
                                   return_period = 1/(1-recurrences_to_calculate), 
                                   Q_cms = logPearsonIII)

# Standard deviations away from mean
distribution_min_4SDs <- peakQ_ln_avg - 4*peakQ_ln_sd
distribution_max_4SDs <- peakQ_ln_avg + 4*peakQ_ln_sd

SD8_range <- seq(distribution_min_4SDs,distribution_max_4SDs,0.01)

logPearsonIII_density_sel <- data.table(
  log_Q_cms = SD8_range,
  probabilities = dpearsonIII(x = SD8_range, mean = peakQ_ln_avg, sd = peakQ_ln_sd, skew = peakQ_ln_skew)
)


density_sel <- data.table(log_Q_cms = SD8_range,
                      probabilities = dnorm(SD8_range, peakQ_ln_avg, peakQ_ln_sd)
)

ggplot() +
  geom_histogram(data = Q_peak_sel, aes(x = log_Q_cms, y = ..density..)) + 
  geom_line(data = density_sel, aes(x = SD8_range, y = probabilities, color = 'Log-normal')) +
  geom_line(data = logPearsonIII_density_sel, aes(x = log_Q_cms, y = probabilities, color = 'Log Pearson III')) +
  geom_vline(xintercept = peakQ_ln_avg, lty = 'dashed') +
  scale_color_manual(values = c('Log-normal' = 'blue', 'Log Pearson III' = 'orange')) +
  # stat_function(fun = dnorm, n = 101, args = list(mean = peakQ_ln_avg, sd = peakQ_ln_sd)) +
  labs(
    x = 'ln(Q)',
    y = 'Density'
  )

# Calculate recurrence intervals based on log-normal distribution
log_normal <- qlnorm(recurrences_to_calculate, peakQ_ln_avg, peakQ_ln_sd)
log_normal_qpeaks <- data.table(return_period = 1/(1-recurrences_to_calculate), Q_cms = log_normal)

# Plot recurrence intervals based on different approaches
recurrence_plot <- ggplot(Q_peak_sel, aes(x = return_period, y = Q_cms)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F, aes(color = 'Linear regression'), lwd = 0.5, fullrange = T) +
  geom_line(data = logPearsonIII_qpeaks, aes(color = 'Log Pearson III'), lwd = 0.5) +
  geom_line(data = log_normal_qpeaks, aes(color = 'Normal'), lwd = 0.5) +
  scale_color_manual(values = c('Normal' = 'blue', 'Log Pearson III' = 'orange', 'Linear regression' = 'red')) +
  theme_classic() +
  scale_x_log10() +
  labs(
    x = 'Recurrence',
    y = 'Discharge (1,000s cfs)',
    color = 'Model'
  )


