library(arrow)
library(dplyr)      
library(ggplot2)     
library(lubridate) 
library(tidyr)  
library(patchwork)
library(corrplot)
library(sp)          # spatial objects
library(spacetime)   # space-time data classes
library(gstat)       # variograms, kriging
library(sf)
library(terra)
library(ggrepel)
library(automap)
library(viridis)
library(viridisLite)

### PREPROCESSING

# --- 1. Read and clean data ---
df_anom <- read_parquet("df_illgraben_anomaly_scores.parquet")
coords <- read.delim("coords.txt", sep = "|", comment.char = "#", 
                     header = FALSE, stringsAsFactors = FALSE)

col_names <- c("Network","Station","Location","Channel","Latitude","Longitude","Elevation",
               "Depth","Azimuth","Dip","SensorDescription","Scale","ScaleFreq",
               "ScaleUnits","SampleRate","StartTime","EndTime")
coords <- coords[-1, ]
names(coords) <- col_names

coords <- coords %>%
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude),
         Elevation = as.numeric(Elevation)) %>%
  arrange(Station, desc(StartTime)) %>%
  group_by(Station) %>% slice(1) %>% ungroup() %>%
  select(Station, Latitude, Longitude, Elevation)

# Convert to MN95/LV95
stations_sf <- st_as_sf(coords, coords = c("Longitude","Latitude"), crs = 4326) %>%
  st_transform(2056) %>%
  mutate(X_mn95 = st_coordinates(.)[,1],
         Y_mn95 = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  select(Station, X_mn95, Y_mn95, Elevation)

# Merge coordinates with anomaly data
df_anom_coords <- df_anom %>%
  left_join(stations_sf, by = c("station" = "Station")) %>%
  filter(!is.na(X_mn95) & !is.na(Y_mn95))

# Bin time to 1-min intervals
df_anom_binned <- df_anom_coords %>%
  mutate(time_bin = floor_date(time, "1 minute")) %>%
  group_by(time_bin, station) %>%
  summarise(
    anomaly_score = mean(anomaly_score, na.rm = TRUE),
    X_mn95 = first(X_mn95),
    Y_mn95 = first(Y_mn95),
    .groups = "drop"
  )

cat("Number of binned observations:", nrow(df_anom_binned), "\n")
# Loop through increasing date ranges and time variogram computation

# Starting date
start_date <- as.POSIXct("2022-07-01")
end_date <- as.POSIXct("2022-07-15")

# Calculate number of days
n_days <- as.numeric(difftime(end_date, start_date, units = "days"))

# Storage for results
timing_results <- data.frame(
  n_days = integer(),
  end_date = character(),
  n_observations = integer(),
  n_stations = integer(),
  n_timepoints = integer(),
  computation_time_sec = numeric(),
  stringsAsFactors = FALSE
)

# Define temporal lags outside loop (constant)
tlags <- seq(0, 60, by = 5)

cat("\n=== Starting Variogram Timing Tests ===\n\n")

# Loop through each day increment
for (i in 1:n_days) {
  current_end <- start_date + days(i)
  
  cat(sprintf("Testing: %s to %s (%d day%s)\n", 
              format(start_date, "%Y-%m-%d"), 
              format(current_end, "%Y-%m-%d"),
              i,
              ifelse(i > 1, "s", "")))
  
  # Filter data for current date range
  df_focused <- df_anom_binned %>%
    filter(time_bin >= start_date & time_bin < current_end)
  
  # Keep only stations with complete data
  complete_stations <- df_focused %>%
    group_by(station) %>%
    summarise(n_obs = n()) %>%
    filter(n_obs == length(unique(df_focused$time_bin))) %>%
    pull(station)
  
  df_complete <- df_focused %>%
    filter(station %in% complete_stations)
  
  # Check if we have enough data
  if (nrow(df_complete) == 0 || length(complete_stations) < 2) {
    cat("  -> Insufficient data, skipping...\n\n")
    next
  }
  
  # Create spatial points
  sp_unique <- df_complete %>%
    distinct(station, X_mn95, Y_mn95) %>%
    arrange(station)
  
  sp_unique_coords <- SpatialPoints(
    coords = sp_unique[, c("X_mn95", "Y_mn95")],
    proj4string = CRS("+proj=somerc +lat_0=46.9524055555556 +lon_0=7.43958333333333 +k_0=1 +x_0=2600000 +y_0=1200000 +ellps=bessel +units=m +no_defs")
  )
  
  # Prepare data for STFDF
  stidf_data <- df_complete %>%
    arrange(time_bin, station) %>%
    select(anomaly_score)
  
  # Create STFDF object
  stidf_obj <- STFDF(
    sp = sp_unique_coords,
    time = sort(unique(df_complete$time_bin)),
    data = stidf_data
  )
  
  cat(sprintf("  Stations: %d | Time points: %d | Total obs: %d\n",
              length(complete_stations),
              length(unique(df_complete$time_bin)),
              nrow(df_complete)))
  
  # TIME THE VARIOGRAM COMPUTATION
  start_time <- Sys.time()
  
  vv <- suppressWarnings(variogram(
    object = anomaly_score ~ 1,
    data = stidf_obj,
    width = 2838,
    cutoff = 6000,
    tlags = tlags,
    tunit = "mins"
  ))
  
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat(sprintf("  ⏱️  Computation time: %.2f seconds (%.2f minutes)\n\n", 
              elapsed, elapsed/60))
  
  # Store results
  timing_results <- rbind(timing_results, data.frame(
    n_days = i,
    end_date = format(current_end, "%Y-%m-%d"),
    n_observations = nrow(df_complete),
    n_stations = length(complete_stations),
    n_timepoints = length(unique(df_complete$time_bin)),
    computation_time_sec = elapsed
  ))
}

cat("\n=== Summary of Results ===\n")
print(timing_results)

# Plot the timing results
ggplot(timing_results, aes(x = n_days, y = computation_time_sec)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 3) +
  labs(
    title = "Variogram Computation Time vs. Data Duration",
    x = "Number of Days",
    y = "Computation Time (seconds)",
    subtitle = paste0("Start date: ", format(start_date, "%Y-%m-%d"))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40")
  )

# Also plot observations vs time
ggplot(timing_results, aes(x = n_observations, y = computation_time_sec)) +
  geom_line(color = "coral", linewidth = 1) +
  geom_point(color = "coral", size = 3) +
  labs(
    title = "Variogram Computation Time vs. Number of Observations",
    x = "Total Observations",
    y = "Computation Time (seconds)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14)
  )