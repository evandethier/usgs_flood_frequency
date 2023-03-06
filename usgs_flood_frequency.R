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
