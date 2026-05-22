library(dplyr)
library(tidyr)

save_dir <- normalizePath("~/R/ppd,bidifacplus.simul")

files <- list.files(
  save_dir,
  pattern = "^rep_.*\\.rds$",
  full.names = TRUE
)

raw_check <- bind_rows(lapply(files, readRDS))

table(raw_check$method, is.na(raw_check$error))

## rank summary
rank_summary <- raw_check %>%
  filter(is.na(error)) %>%
  group_by(method, d, ang, ang_deg) %>%
  summarise(
    nrep = n_distinct(rep),
    mean_rJ = mean(rJ, na.rm = TRUE),
    mean_rI1 = mean(rI1, na.rm = TRUE),
    mean_rI2 = mean(rI2, na.rm = TRUE),
    prop_correct_J = mean(rJ == 3, na.rm = TRUE),
    prop_correct_I1 = mean(rI1 == 3, na.rm = TRUE),
    prop_correct_I2 = mean(rI2 == 6, na.rm = TRUE),
    mean_time = mean(time, na.rm = TRUE),
    .groups = "drop"
  )

## error mean/sd
error_summary <- raw_check %>%
  filter(is.na(error)) %>%
  group_by(method, d, ang, ang_deg) %>%
  summarise(
    nrep = n_distinct(rep),
    eJ1 = sprintf("%.4f (%.4f)", mean(eJ1, na.rm = TRUE), sd(eJ1, na.rm = TRUE)),
    eJ2 = sprintf("%.4f (%.4f)", mean(eJ2, na.rm = TRUE), sd(eJ2, na.rm = TRUE)),
    eI1 = sprintf("%.4f (%.4f)", mean(eI1, na.rm = TRUE), sd(eI1, na.rm = TRUE)),
    eI2 = sprintf("%.4f (%.4f)", mean(eI2, na.rm = TRUE), sd(eI2, na.rm = TRUE)),
    .groups = "drop"
  )

print(error_summary, n = Inf)

## 각 setting별 rJ table
rank_table_rJ <- raw_check %>%
  filter(is.na(error)) %>%
  count(method, d, ang_deg, rJ) %>%
  pivot_wider(
    names_from = rJ,
    values_from = n,
    values_fill = 0
  ) %>%
  arrange(method, d, ang_deg)

print(rank_table_rJ, n = Inf)