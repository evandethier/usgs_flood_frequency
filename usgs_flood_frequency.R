#### IMPORT PACKAGES ####
library(data.table)
library(dataRetrieval)
library(ggplot2)
library(lubridate)
library(patchwork)

# New packages
library(e1071)
library(smwrBase)
library(ggpubr)


#### THEMES ####
season_facet <- theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linetype = 'dashed',color = 'grey70'),
    panel.grid.major.x = element_blank(),
    # panel.grid = element_blank(),
    legend.position = 'none',
    panel.border = element_rect(size = 0.5),
    text = element_text(size=8),
    axis.text = element_text(size = 8), 
    plot.title = element_text(size = 9),
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0, margin = margin(0,0,0,0, unit = 'pt'))
  )


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
# site_sel <- hcdn_simple[grepl(casefold('CRYSTAL RIVER AB AVALANCHE C'), casefold(station_nm))]
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

# Compute maximum, minimum, and average peak flow
Q_peak_min <- min(Q_peak_sel[,Q_cms], na.rm = T)
Q_peak_max <- max(Q_peak_sel[,Q_cms], na.rm = T)
Q_peak_mean <- mean(Q_peak_sel[,Q_cms], na.rm = T)

# Average annual flow
Q_annual_avg <- data.table(readNWISstat(site_no_sel,
                                parameterCd = '00060',
                                 statReportType = 'annual',
                                 statType = 'mean'))

Q_cms_total_annual_avg_m3 <- mean(Q_annual_avg[,mean_va]) * 0.0283168 * 60 * 60 * 24 * 365

# What percent of the total water that flows through the river comes in the peak flow on an average year?
# What about for the peak year?
# Assume: peak flow is sustained for 1 day (unlikely, but a decent approximation)
percent_total_water_avg <- Q_peak_mean * 60 * 60 * 24/Q_cms_total_annual_avg_m3 * 100
percent_total_water_max <- Q_peak_max * 60 * 60 * 24/Q_cms_total_annual_avg_m3 * 100

# Analyze peak flow timing
peak_flow_timing <- Q_peak_sel[,.(N_events = .N,
                                  fraction_of_events = .N/nrow(peak_flow_timing)),
                               by = .(month)]

# Plot peak flow timing: number of events per month
peak_flow_timing_plot <- ggplot(peak_flow_timing, aes(x = month, y = N_events)) + 
  geom_bar(stat = 'identity') +
  season_facet +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c('Jan.','Apr.','Jul.', 'Oct.')) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
  labs(
    x = 'Month',
    y = 'N Floods'
  )

# Another way of doing this doesn't require a separate table
peak_flow_timing_plot <- ggplot(Q_peak_sel, aes(x = month)) + 
  geom_bar() +
  season_facet +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c('Jan.','Apr.','Jul.', 'Oct.')) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
  labs(
    x = 'Month',
    y = 'N Floods'
  )



# Plot data
peak_flow_timeseries <- ggplot(Q_peak_sel, aes(x = water_year, y = Q_cms)) + 
  geom_point() +
  scale_y_continuous(limits = c(0, Q_peak_max)) +
  # scale_y_log10(limits = c(Q_peak_min, Q_peak_max)) +
  season_facet + 
  labs(
    x = 'Year',
    y = 'Annual peak discharge (m3/s)'
  )

# Plot distribution
peak_flow_distribution <- ggplot(Q_peak_sel, aes(x = Q_cms)) + 
  # geom_histogram(aes(y = ..density..), bins = 15) +
  geom_histogram(bins = 15) +
  scale_x_continuous(limits = c(0, Q_peak_max), sec.axis = dup_axis()) +
  # scale_x_log10(limits = c(Q_peak_min, Q_peak_max)) +
  season_facet + 
  labs(
    x = 'Annual peak discharge (m3/s)',
    # y = 'Density'
    y = 'Count'
  ) +
  # scale_x_log10() +
  rotate()

# Combined plot
combined_peak_flow_plot <- 
  peak_flow_timeseries + 
  peak_flow_distribution + 
    theme(axis.title.y.left = element_blank(),
          axis.text.y.left = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  plot_layout(widths = c(1, 0.3)) +
  plot_annotation(tag_levels = 'a')

print(combined_peak_flow_plot)

combined_peak_flow_plot_with_timing <- 
  peak_flow_timing_plot /
  combined_peak_flow_plot +
  plot_layout(heights = c(0.7,1)) +
  plot_annotation(tag_levels = 'a')
print(combined_peak_flow_plot_with_timing)

#### FLOOD FREQUENCY ANALYSIS ####
#### Prepare data for analysis
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
probabilities_to_calculate <- seq(0.01,0.999,0.001)

# Log Pearson III -- recurrence calculation
logPearsonIII <- qlpearsonIII(probabilities_to_calculate, peakQ_ln_avg, peakQ_ln_sd, peakQ_ln_skew)
logPearsonIII_qpeaks <- data.table(probabilities = probabilities_to_calculate,
                                   return_period = 1/(1-probabilities_to_calculate), 
                                   Q_cms = logPearsonIII)

# Standard deviations away from mean
distribution_min_4SDs <- peakQ_ln_avg - 4*peakQ_ln_sd
distribution_max_4SDs <- peakQ_ln_avg + 4*peakQ_ln_sd

# Set up a vector of discharge values spanning +/- 4 standard deviations from the log-transformed mean
# We will use this vector to model our flood probabilities
SD8_range <- seq(distribution_min_4SDs,distribution_max_4SDs,0.01)



# Compute probability density function for log-transformed discharges using normal distribution
density_sel <- data.table(log_Q_cms = SD8_range,
                      probabilities = dnorm(SD8_range, peakQ_ln_avg, peakQ_ln_sd)
)

# Compute probability density function for log-transformed discharges using Log Pearson III distribution
logPearsonIII_density_sel <- data.table(
  log_Q_cms = SD8_range,
  probabilities = dpearsonIII(x = SD8_range, mean = peakQ_ln_avg, sd = peakQ_ln_sd, skew = peakQ_ln_skew)
)

# Log-normal and Log Pearson III distributions on flood histogram
observed_vs_modeled_distributions <- ggplot() +
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

print(observed_vs_modeled_distributions)

# Calculate recurrence intervals based on log-normal distribution
# This takes in a vector of recurrences, an average, and a standard 
log_normal <- qlnorm(probabilities_to_calculate, peakQ_ln_avg, peakQ_ln_sd)
log_normal_qpeaks <- data.table(return_period = 1/(1-probabilities_to_calculate), Q_cms = log_normal)

# Plot recurrence intervals based on different approaches
recurrence_plot <- ggplot(Q_peak_sel, aes(x = return_period, y = Q_cms)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F, aes(color = 'Linear regression'), lwd = 0.5, fullrange = T) +
  geom_line(data = logPearsonIII_qpeaks, aes(color = 'Log Pearson III'), lwd = 0.5) +
  geom_line(data = log_normal_qpeaks, aes(color = 'Normal'), lwd = 0.5) +
  scale_color_manual(values = c('Normal' = 'blue', 'Log Pearson III' = 'orange', 'Linear regression' = 'red')) +
  theme_classic() +
  scale_x_log10() +
  # scale_y_log10() +
  labs(
    x = 'Recurrence',
    y = 'Discharge (1,000s cfs)',
    color = 'Model'
  )
print(recurrence_plot)


#### PLOT FINAL X-YEAR FLOOD ESTIMATES ####
# Get discharge for the 10, 50, 100, 500, and 1000 year floods. 
# Recall that the probability of an event occurring is 1/Recurrence Interval
# So for the 10-year event we want the probability of 1/10
X_year_floods <- logPearsonIII_qpeaks[probabilities %in% (1-1/c(10,50,100,500, 1000))]

# Add to the peak flow timeseries plot
peak_flow_timeseries_annotated <- peak_flow_timeseries + 
  geom_hline(data = X_year_floods, aes(yintercept = Q_cms, color = factor(return_period)), lwd = 1) +
  scale_color_manual(values = c('blue','green','yellow','orange','red')) +
  scale_y_continuous(limits = c(0, max(X_year_floods$Q_cms))) +
  geom_text(data = X_year_floods, 
            aes(x = -Inf, y = Q_cms, label = paste0(round(return_period), '-yr')),
            hjust = -0.25, vjust = -0.5)


# Combined plot
combined_peak_flow_plot_annotated <- 
  peak_flow_timeseries_annotated + 
  peak_flow_distribution + 
  scale_x_continuous(limits = c(0, max(X_year_floods$Q_cms)), sec.axis = dup_axis()) +
  theme(axis.title.y.left = element_blank(),
        axis.text.y.left = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  plot_layout(widths = c(1, 0.3)) +
  plot_annotation(tag_levels = 'a')

print(combined_peak_flow_plot_annotated)

combined_peak_flow_plot_with_timing_annotated <- 
  peak_flow_timing_plot /
  combined_peak_flow_plot_annotated +
  plot_layout(heights = c(0.7,1)) +
  plot_annotation(tag_levels = 'a')
print(combined_peak_flow_plot_with_timing_annotated)

