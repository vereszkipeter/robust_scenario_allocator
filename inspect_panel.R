
library(targets)
panel <- tar_read(panel_monthly_full)
print(summary(panel))
print(tail(panel))
