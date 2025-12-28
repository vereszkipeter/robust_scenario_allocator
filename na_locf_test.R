# Test na.locf with merged xts objects of different lengths
library(xts)
library(zoo)

# Series 1: Ends Nov 30
s1_dates <- seq(as.Date("2025-09-30"), as.Date("2025-11-30"), by = "month")
s1_values <- 1:3
series1 <- xts(s1_values, order.by = s1_dates)
colnames(series1) <- "S1"

# Series 2: Ends Dec 31
s2_dates <- seq(as.Date("2025-09-30"), as.Date("2025-12-31"), by = "month")
s2_values <- 10:13
series2 <- xts(s2_values, order.by = s2_dates)
colnames(series2) <- "S2"

# Merge
merged_series <- merge(series1, series2, all = TRUE)
print("Merged Series:")
print(merged_series)

# Apply na.locf
filled_series <- na.locf(merged_series)
print("Filled Series (na.locf):")
print(filled_series)
