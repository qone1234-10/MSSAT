library(dplyr)
 
save_dir <- normalizePath("~/R/ppd,bidifacplus.3block.simul")
ang_label <- "pi8"
files <- list.files(
  save_dir,
  pattern = paste0("^rep_.*_", ang_label, "\\.rds$"),
  full.names = TRUE
)
cat("saved reps:", length(files), "\n")

raw_check <- bind_rows(lapply(files, readRDS))

print(table(raw_check$method, is.na(raw_check$error)))

raw_check %>%
  filter(is.na(error)) %>%
  group_by(method) %>%
  summarise(
    nrep = n_distinct(rep),
    mean_rG = mean(rG),
    mean_r12 = mean(r12),
    mean_r13 = mean(r13),
    mean_r23 = mean(r23),
    mean_rI1 = mean(rI1),
    mean_rI2 = mean(rI2),
    mean_rI3 = mean(rI3),
    mean_time = mean(time),
    total_time = sum(time),
    .groups = "drop"
  )

rank_long <- raw_check %>%
  dplyr::filter(is.na(error)) %>%
  dplyr::select(method, rG, r12, r13, r23, rI1, rI2, rI3) %>%
  tidyr::pivot_longer(
    cols = c(rG, r12, r13, r23, rI1, rI2, rI3),
    names_to = "structure",
    values_to = "rank_hat"
  )

raw_check %>%
  dplyr::filter(is.na(error)) %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    eG  = sprintf("%.4f (%.4f)", mean(eG,  na.rm = TRUE), sd(eG,  na.rm = TRUE)),
    e12 = sprintf("%.4f (%.4f)", mean(e12, na.rm = TRUE), sd(e12, na.rm = TRUE)),
    e13 = sprintf("%.4f (%.4f)", mean(e13, na.rm = TRUE), sd(e13, na.rm = TRUE)),
    e23 = sprintf("%.4f (%.4f)", mean(e23, na.rm = TRUE), sd(e23, na.rm = TRUE)),
    eI1 = sprintf("%.4f (%.4f)", mean(eI1, na.rm = TRUE), sd(eI1, na.rm = TRUE)),
    eI2 = sprintf("%.4f (%.4f)", mean(eI2, na.rm = TRUE), sd(eI2, na.rm = TRUE)),
    eI3 = sprintf("%.4f (%.4f)", mean(eI3, na.rm = TRUE), sd(eI3, na.rm = TRUE)),
    time = sprintf("%.4f (%.4f)", mean(time, na.rm = TRUE), sd(time, na.rm = TRUE)),
    .groups = "drop"
  )

rank_long %>%
  dplyr::group_by(method, structure) %>%
  dplyr::summarise(
    rank_mean = mean(rank_hat, na.rm = TRUE),
    rank_sd   = sd(rank_hat, na.rm = TRUE),
    rank_summary = sprintf("%.2f (%.2f)", rank_mean, rank_sd),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = structure,
    values_from = rank_summary
  )

