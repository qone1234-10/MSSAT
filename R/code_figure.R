## This script creates the figure
## Starts from the merged raw CSV containing MATLAB + R methods.
## Run simulation_combine.R
library(dplyr)
library(readr)

source("simulation_combine.R")

preprocess_combined <- function(raw) {
  raw %>%
    mutate(
      method = as.character(method),
      
      ## method name
      method = recode(
        method,
        "ajive" = "AJIVE",
        "divas" = "DIVAS",
        "mssat" = "MSSAT",
        "MSSAT" = "MSSAT",
        "PPD" = "PPD",
        "BIDIFAC+" = "BIDIFAC+",
        .default = method
      ),
      
      ## angle
      ang_deg = case_when(
        !is.na(suppressWarnings(as.numeric(ang_deg))) ~ suppressWarnings(as.numeric(ang_deg)),
        !is.na(suppressWarnings(as.numeric(ang)))     ~ round(suppressWarnings(as.numeric(ang)) * 180 / pi, 1),
        TRUE ~ NA_real_
      ),
      
      ang_label = case_when(
        !is.na(ang_deg) ~ paste0("theta = ", ang_deg, " deg"),
        TRUE ~ NA_character_
      ),
      
      ## time 
      time = case_when(
        source == "R" & !is.na(suppressWarnings(as.numeric(time))) ~ suppressWarnings(as.numeric(time)),
        source == "MATLAB" & !is.na(suppressWarnings(as.numeric(runtime))) ~ suppressWarnings(as.numeric(runtime)),
        !is.na(suppressWarnings(as.numeric(time))) ~ suppressWarnings(as.numeric(time)),
        !is.na(suppressWarnings(as.numeric(runtime))) ~ suppressWarnings(as.numeric(runtime)),
        TRUE ~ NA_real_
      ),
      
      ## error
      error = case_when(
        is.logical(error) ~ error,
        is.na(error) ~ FALSE,
        TRUE ~ as.logical(error)
      ),
      
      ## true rank
      true_rJ  = suppressWarnings(as.numeric(true_rJ)),
      true_rI1 = suppressWarnings(as.numeric(true_rI1)),
      true_rI2 = suppressWarnings(as.numeric(true_rI2)),
      
      correct_rJ  = ifelse(is.na(true_rJ),  NA, rJ  == true_rJ),
      correct_rI1 = ifelse(is.na(true_rI1), NA, rI1 == true_rI1),
      correct_rI2 = ifelse(is.na(true_rI2), NA, rI2 == true_rI2),
      
      eJ_mean = rowMeans(cbind(eJ1, eJ2), na.rm = TRUE),
      eI_mean = rowMeans(cbind(eI1, eI2), na.rm = TRUE),
      eTotal_mean = rowMeans(cbind(eTotal1, eTotal2), na.rm = TRUE),
      eTotal_mean = ifelse(is.nan(eTotal_mean), NA_real_, eTotal_mean)
    ) %>%
    mutate(
      method = factor(
        method,
        levels = c("AJIVE", "DIVAS", "MSSAT", "PPD", "BIDIFAC+")
      ),
      ang_label = factor(
        ang_label,
        levels = paste0("theta = ", sort(unique(ang_deg[!is.na(ang_deg)])), " deg")
      )
    )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(cowplot)
})

base_dir <- "~/R"
in_csv <- file.path(base_dir, "combined_2block_all_methods_raw.csv")
out_pdf <- file.path(base_dir, "Rplot_combined_2block_all_methods.pdf")
out_png <- file.path(base_dir, "Rplot_combined_2block_all_methods.png")

if (!file.exists(path.expand(in_csv))) {
  stop("Input CSV not found: ", in_csv)
}

raw <- read_csv(
  "~/R/combined_2block_all_methods_raw.csv",
  show_col_types = FALSE
)

raw <- preprocess_combined(raw)

raw <- raw %>%
  mutate(d = 2 * as.numeric(d))

table(raw$source, raw$method, useNA = "ifany")

raw %>%
  count(source, method)

raw %>%
  group_by(source, method) %>%
  summarise(
    n = n(),
    time_na = sum(is.na(time)),
    time_mean = mean(time, na.rm = TRUE),
    .groups = "drop"
  )

raw <- raw %>%
  mutate(
    method = as.character(method),
    method = recode(
      method,
      "ajive" = "AJIVE",
      "divas" = "DIVAS",
      "mssat" = "MSSAT",
      "MSSAT" = "MSSAT",
      .default = method
    ),
    
    time_plot = case_when(
      source == "R"      & !is.na(time)    ~ as.numeric(time),
      source == "MATLAB" & !is.na(runtime) ~ as.numeric(runtime),
      !is.na(time)                         ~ as.numeric(time),
      !is.na(runtime)                      ~ as.numeric(runtime),
      TRUE                                 ~ NA_real_
    )
  )

raw[c(1, 10000), ]
table(raw$source, raw$method, useNA = "ifany")

raw[c(1,10000),]

cat("Input CSV:", normalizePath(path.expand(in_csv), mustWork = FALSE), "\n")
cat("Rows:", nrow(raw), "\n")
cat("Columns:\n")
print(names(raw))

## -----------------------------
## Helpers
## -----------------------------
first_existing_col <- function(dat, candidates, required = TRUE) {
  hit <- intersect(candidates, names(dat))
  if (length(hit) == 0) {
    if (required) stop("None of these columns exist: ", paste(candidates, collapse = ", "))
    return(NULL)
  }
  hit[1]
}

as_num_col <- function(dat, col) {
  if (is.null(col)) return(rep(NA_real_, nrow(dat)))
  as.numeric(dat[[col]])
}

## angle column: merged file usually has radian ang, but old plotting file used ang_deg.
if (!"ang_deg" %in% names(raw)) {
  if ("ang" %in% names(raw)) {
    raw <- raw %>% mutate(ang_deg = as.numeric(ang) * 180 / pi)
  } else {
    stop("No angle column found. Expected either ang_deg or ang.")
  }
}

## method labels/order
raw <- raw %>%
  mutate(
    method = recode(
      as.character(method),
      "mssat" = "MSSAT",
      "MSSAT" = "MSSAT",
      "ajive" = "AJIVE",
      "AJIVE" = "AJIVE",
      "divas" = "DIVAS",
      "DIVAS" = "DIVAS",
      "bidifac" = "BIDIFAC",
      "BIDIFAC" = "BIDIFAC",
      "PPD" = "PPD",
      "BIDIFAC+" = "BIDIFAC+",
      .default = as.character(method)
    ),
    d = as.numeric(d),
    ang_deg = as.numeric(ang_deg),
    ang_label = factor(
      paste0("theta = ", round(ang_deg, 1), " deg"),
      levels = paste0("theta = ", round(sort(unique(ang_deg)), 1), " deg")
    )
  ) %>%
  filter(method != "BIDIFAC") %>%
  mutate(method = factor(method, levels = c("AJIVE", "DIVAS", "PPD", "BIDIFAC+", "MSSAT")))

## -----------------------------
## Build long raw metric table
## -----------------------------
## Rank metrics
rank_long <- raw %>%
  transmute(
    rep, d, ang_deg, ang_label, method,
    joint = as.numeric(rJ),
    ind1 = as.numeric(rI1),
    ind2 = as.numeric(rI2)
  ) %>%
  pivot_longer(
    cols = c(joint, ind1, ind2),
    names_to = "structure",
    values_to = "value"
  ) %>%
  mutate(metric_type = "rank", block = NA_character_)

## Error column names from the merged file.
eJ1_col <- first_existing_col(raw, c("eJ1", "err_joint1", "err.joint1"), required = FALSE)
eJ2_col <- first_existing_col(raw, c("eJ2", "err_joint2", "err.joint2"), required = FALSE)
eI1_col <- first_existing_col(raw, c("eI1", "err_ind1", "err.ind1"), required = FALSE)
eI2_col <- first_existing_col(raw, c("eI2", "err_ind2", "err.ind2"), required = FALSE)
eT1_col <- first_existing_col(raw, c("eTotal1", "err_overall1", "err.overall1"), required = FALSE)
eT2_col <- first_existing_col(raw, c("eTotal2", "err_overall2", "err.overall2"), required = FALSE)

## Stop only if all error blocks are absent.
if (all(vapply(list(eJ1_col, eJ2_col, eI1_col, eI2_col, eT1_col, eT2_col), is.null, logical(1)))) {
  stop("No error columns found. Expected eJ1/eJ2/eI1/eI2/eTotal1/eTotal2 or matching alternatives.")
}

err_wide <- raw %>%
  mutate(
    joint_block1      = as_num_col(., eJ1_col),
    joint_block2      = as_num_col(., eJ2_col),
    individual_block1 = as_num_col(., eI1_col),
    individual_block2 = as_num_col(., eI2_col),
    overall_block1    = as_num_col(., eT1_col),
    overall_block2    = as_num_col(., eT2_col)
  ) %>%
  transmute(
    rep, d, ang_deg, ang_label, method,
    joint_block1, joint_block2,
    individual_block1, individual_block2,
    overall_block1, overall_block2
  )

err_long <- err_wide %>%
  pivot_longer(
    cols = c(joint_block1, joint_block2, individual_block1, individual_block2, overall_block1, overall_block2),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(
    metric_type = "error",
    structure = case_when(
      str_detect(metric, "^joint_") ~ "joint",
      str_detect(metric, "^individual_") ~ "individual",
      str_detect(metric, "^overall_") ~ "overall",
      TRUE ~ NA_character_
    ),
    block = case_when(
      str_detect(metric, "block1$") ~ "block1",
      str_detect(metric, "block2$") ~ "block2",
      TRUE ~ NA_character_
    )
  ) %>%
  select(rep, d, ang_deg, ang_label, method, metric_type, structure, block, value)

## Add average-by-block error rows, matching old code's block == "avg" convention.
err_avg <- err_long %>%
  filter(metric_type == "error", structure %in% c("joint", "individual", "overall"), block %in% c("block1", "block2")) %>%
  group_by(rep, d, ang_deg, ang_label, method, metric_type, structure) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(block = "avg")

## Time metric: use row-wise unified time_plot.
## Important: do NOT choose one global time column such as runtime,
## because MATLAB rows use runtime while PPD/BIDIFAC+ rows use time.
## The global first_existing_col() approach drops R-method times.
if (!"time_plot" %in% names(raw)) {
  raw <- raw %>%
    mutate(
      time_plot = case_when(
        source == "R" & "time" %in% names(.) & !is.na(suppressWarnings(as.numeric(time))) ~ suppressWarnings(as.numeric(time)),
        source == "MATLAB" & "runtime" %in% names(.) & !is.na(suppressWarnings(as.numeric(runtime))) ~ suppressWarnings(as.numeric(runtime)),
        "time" %in% names(.) & !is.na(suppressWarnings(as.numeric(time))) ~ suppressWarnings(as.numeric(time)),
        "runtime" %in% names(.) & !is.na(suppressWarnings(as.numeric(runtime))) ~ suppressWarnings(as.numeric(runtime)),
        "elapsed" %in% names(.) & !is.na(suppressWarnings(as.numeric(elapsed))) ~ suppressWarnings(as.numeric(elapsed)),
        "time_sec_file" %in% names(.) & !is.na(suppressWarnings(as.numeric(time_sec_file))) ~ suppressWarnings(as.numeric(time_sec_file)),
        TRUE ~ NA_real_
      )
    )
}

cat("\nTime check by source/method before plotting:\n")
print(as.data.frame(raw %>%
                      group_by(source, method) %>%
                      summarise(
                        n = n(),
                        time_non_na = sum(!is.na(time_plot)),
                        time_mean = mean(time_plot, na.rm = TRUE),
                        .groups = "drop"
                      )))

time_long <- raw %>%
  transmute(
    rep, d, ang_deg, ang_label, method,
    value = suppressWarnings(as.numeric(time_plot)),
    metric_type = "time_plot",
    structure = "time_plot",
    block = NA_character_
  ) %>%
  filter(!is.na(value), value > 0)

metric_raw <- bind_rows(
  rank_long %>% select(rep, d, ang_deg, ang_label, method, metric_type, structure, block, value),
  err_long,
  err_avg,
  time_long
)

## Summarise across Monte Carlo replicates.
df <- metric_raw %>%
  group_by(method, d, ang_deg, ang_label, metric_type, structure, block) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = sum(!is.na(value)),
    .groups = "drop"
  ) %>%
  mutate(
    method = factor(as.character(method), levels = c("AJIVE", "DIVAS", "PPD", "BIDIFAC+", "MSSAT"))
  )

cat("\nSummary rows by metric_type/structure/block:\n")
print(as.data.frame(df %>% count(metric_type, structure, block, name = "n_rows")))

## -----------------------------
## Plot settings from code_figure_2.R
## -----------------------------
method_cols <- c(
  "AJIVE"    = "#d95f02",
  "DIVAS"    = "#7570b3",
  "PPD"      = "#66a61e",
  "BIDIFAC+" = "#a6761d",
  "MSSAT"    = "#e7298a"
)

theme_paper <- function() {
  theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      strip.background = element_rect(fill = "grey95", colour = "grey80"),
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
}

plot_band <- function(dat, ylab, title_txt, log_y = FALSE) {
  dat <- dat %>% arrange(method)
  p <- ggplot(dat, aes(x = d, y = mean, color = method, fill = method, group = method)) +
    geom_ribbon(
      aes(ymin = pmax(mean - sd, 0), ymax = mean + sd),
      alpha = 0.18,
      colour = NA,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    facet_wrap(~ ang_label, nrow = 1) +
    scale_color_manual(
      values = method_cols,
      breaks = names(method_cols),
      drop = FALSE,
      guide = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(linewidth = 1.2, shape = 16, size = 3, alpha = 1)
      )
    ) +
    scale_fill_manual(values = method_cols, breaks = names(method_cols), drop = FALSE, guide = "none") +
    labs(x = "Signal strength", y = ylab, title = title_txt, color = NULL) +
    theme_paper()
  
  if (log_y) p <- p + scale_y_log10()
  p
}

plot_band_with_baseline <- function(dat, ylab, title_txt, baseline_df = NULL, baseline_col = NULL, log_y = FALSE) {
  dat <- dat %>% arrange(method)
  p <- ggplot(dat, aes(x = d, y = mean, color = method, fill = method, group = method)) +
    geom_ribbon(
      aes(ymin = pmax(mean - sd, 0), ymax = mean + sd),
      alpha = 0.18,
      colour = NA,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8)
  
  if (!is.null(baseline_df) && !is.null(baseline_col)) {
    p <- p + geom_hline(
      data = baseline_df,
      aes(yintercept = .data[[baseline_col]]),
      linetype = "dashed",
      linewidth = 0.6,
      colour = "black"
    )
  }
  
  p <- p +
    facet_wrap(~ ang_label, nrow = 1) +
    scale_color_manual(
      values = method_cols,
      breaks = names(method_cols),
      drop = FALSE,
      guide = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(linewidth = 1.2, shape = 16, size = 3, alpha = 1)
      )
    ) +
    scale_fill_manual(values = method_cols, breaks = names(method_cols), drop = FALSE, guide = "none") +
    labs(x = "Signal strength", y = ylab, title = title_txt, color = NULL) +
    theme_paper()
  
  if (log_y) p <- p + scale_y_log10()
  p
}

## True ranks. If the merged file has true_rJ/true_rI1/true_rI2, use them.
## Otherwise use the same rule as code_figure_2.R.
if (all(c("true_rJ", "true_rI1", "true_rI2") %in% names(raw)) && any(!is.na(raw$true_rJ))) {
  baseline_rank <- raw %>%
    group_by(ang_deg, ang_label) %>%
    summarise(
      true_joint = dplyr::first(na.omit(true_rJ)),
      true_ind1 = dplyr::first(na.omit(true_rI1)),
      true_ind2 = dplyr::first(na.omit(true_rI2)),
      .groups = "drop"
    )
} else {
  baseline_rank <- df %>%
    distinct(ang_deg, ang_label) %>%
    mutate(
      true_joint = ifelse(abs(ang_deg - 0) < 1e-8, 6, 3),
      true_ind1  = ifelse(abs(ang_deg - 0) < 1e-8, 0, 3),
      true_ind2  = ifelse(abs(ang_deg - 0) < 1e-8, 3, 6)
    )
}

## ===================== A. Joint rank =====================
df_joint_rank <- df %>% filter(metric_type == "rank", structure == "joint")

p1 <- plot_band_with_baseline(
  df_joint_rank,
  ylab = "Estimated joint rank",
  title_txt = "A. Joint rank",
  baseline_df = baseline_rank,
  baseline_col = "true_joint"
)

## ===================== B. Individual rank =====================
df_ind_rank <- df %>%
  filter(metric_type == "rank", structure %in% c("ind1", "ind2")) %>%
  mutate(block = recode(structure, "ind1" = "Block 1", "ind2" = "Block 2"))

baseline_rank_long <- baseline_rank %>%
  select(ang_label, true_ind1, true_ind2) %>%
  pivot_longer(cols = c(true_ind1, true_ind2), names_to = "tmp", values_to = "yint") %>%
  mutate(block = recode(tmp, "true_ind1" = "Block 1", "true_ind2" = "Block 2"))

p2 <- ggplot(df_ind_rank, aes(x = d, y = mean, color = method, fill = method, group = method)) +
  geom_ribbon(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), alpha = 0.18, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  geom_hline(data = baseline_rank_long, aes(yintercept = yint), linetype = "dashed", linewidth = 0.6, colour = "black") +
  facet_grid(block ~ ang_label, scales = "free_y") +
  scale_color_manual(
    values = method_cols,
    breaks = names(method_cols),
    drop = FALSE,
    guide = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = 1.2, shape = 16, size = 3, alpha = 1))
  ) +
  scale_fill_manual(values = method_cols, breaks = names(method_cols), drop = FALSE, guide = "none") +
  labs(x = "Signal strength", y = "Estimated individual rank", title = "B. Individual rank", color = NULL) +
  theme_paper()

## ===================== C. Joint reconstruction error =====================
df_err_joint_block <- df %>%
  filter(metric_type == "error", structure == "joint", block %in% c("block1", "block2")) %>%
  mutate(block = recode(block, "block1" = "Block 1", "block2" = "Block 2")) %>%
  arrange(method)

p3 <- ggplot(df_err_joint_block, aes(x = d, y = mean, color = method, fill = method, group = method)) +
  geom_ribbon(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), alpha = 0.18, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_grid(block ~ ang_label, scales = "free_y") +
  scale_color_manual(
    values = method_cols,
    breaks = names(method_cols),
    drop = FALSE,
    guide = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = 1.2, shape = 16, size = 3, alpha = 1))
  ) +
  scale_fill_manual(values = method_cols, breaks = names(method_cols), drop = FALSE, guide = "none") +
  labs(x = "Signal strength", y = "Relative Frobenius error", title = "C. Joint reconstruction error", color = NULL) +
  theme_paper()

## ===================== D. Individual reconstruction error =====================
df_err_ind_block <- df %>%
  filter(metric_type == "error", structure == "individual", block %in% c("block1", "block2")) %>%
  mutate(
    block = recode(block, "block1" = "Block 1", "block2" = "Block 2"),
    ## Same default caps as code_figure_2.R. Adjust if needed.
    y_cap = ifelse(block == "Block 1", 4, 2.5),
    mean_raw = mean,
    ymin = pmax(mean - sd, 0),
    ymax = mean + sd,
    mean_plot = pmin(mean, y_cap),
    ymin_plot = pmin(ymin, y_cap),
    ymax_plot = pmin(ymax, y_cap),
    capped = mean_raw > y_cap
  ) %>%
  arrange(method)

p4 <- ggplot(df_err_ind_block, aes(x = d, y = mean_plot, color = method, fill = method, group = method)) +
  geom_ribbon(
    aes(ymin = ymin_plot, ymax = ymax_plot),
    alpha = 0.18,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(data = df_err_ind_block %>% filter(!capped), size = 1.8) +
  geom_point(data = df_err_ind_block %>% filter(capped), aes(y = mean_plot), shape = 24, size = 2.4, stroke = 0.8) +
  facet_grid(block ~ ang_label, scales = "free_y") +
  scale_color_manual(
    values = method_cols,
    breaks = names(method_cols),
    drop = FALSE,
    guide = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = 1.2, shape = 16, size = 3, alpha = 1))
  ) +
  scale_fill_manual(values = method_cols, breaks = names(method_cols), drop = FALSE, guide = "none") +
  labs(x = "Signal strength", y = "Relative Frobenius error", title = "D. Individual reconstruction error", color = NULL) +
  theme_paper()

## ===================== E. Log-time =====================
df_time <- df %>% filter(metric_type == "time_plot")
if (nrow(df_time) > 0) {
  p5 <- plot_band(df_time, ylab = "Time (sec, log scale)", title_txt = "E. Computation time", log_y = TRUE)
} else {
  p5 <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No computation time column found") +
    theme_void() +
    labs(title = "E. Computation time")
}

## Remove legends from panels; manual compact legend below.
p1 <- p1 + labs(x = "Signal strength") + theme(legend.position = "none")
p2 <- p2 + labs(x = "Signal strength") + theme(legend.position = "none")
p3 <- p3 + labs(x = "Signal strength") + theme(legend.position = "none")
p4 <- p4 + labs(x = "Signal strength") + theme(legend.position = "none")
p5 <- p5 + labs(x = "Signal strength") + theme(legend.position = "none")

## Main layout. This defines the missing main_plot object in the uploaded file.
main_plot <- cowplot::plot_grid(
  cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(0.85, 1.15), align = "v"),
  cowplot::plot_grid(p3, p4, p5, ncol = 1, rel_heights = c(1.15, 1.15, 0.85), align = "v"),
  ncol = 2,
  rel_widths = c(1, 1),
  align = "hv"
)

legend_order <- c("MSSAT", "AJIVE", "BIDIFAC+", "DIVAS", "PPD")
legend_manual <- tibble(
  method = factor(legend_order, levels = legend_order),
  label = legend_order
) %>%
  mutate(
    label_width = nchar(label) * 0.105,
    item_width  = 0.30 + label_width,
    gap = 0.48,
    x0 = c(0, cumsum(head(item_width + gap, -1))),
    y = 1
  )
xmax <- max(legend_manual$x0 + legend_manual$item_width)

legend_plot <- ggplot(legend_manual, aes(y = y)) +
  geom_segment(aes(x = x0, xend = x0 + 0.30, yend = y, color = method), linewidth = 0.9) +
  geom_point(aes(x = x0 + 0.15, color = method), size = 2.0) +
  geom_text(aes(x = x0 + 0.40, label = label), hjust = 0, size = 3.2) +
  scale_color_manual(values = method_cols[legend_order], guide = "none") +
  coord_cartesian(xlim = c(-0.03, xmax + 0.03), ylim = c(0.92, 1.08), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

legend_compact <- cowplot::plot_grid(NULL, legend_plot, NULL, nrow = 1, rel_widths = c(0.39, 0.22, 0.39))

final_plot <- cowplot::plot_grid(main_plot, legend_compact, ncol = 1, rel_heights = c(1, 0.028))

## Save
if (capabilities("cairo")) {
  ggsave(filename = out_pdf, plot = final_plot, width = 15, height = 12, units = "in", device = cairo_pdf)
} else {
  ggsave(filename = out_pdf, plot = final_plot, width = 15, height = 12, units = "in")
}
ggsave(filename = out_png, plot = final_plot, width = 15, height = 12, units = "in", dpi = 400)

cat("\nSaved PDF:", normalizePath(path.expand(out_pdf), mustWork = FALSE), "\n")
cat("Saved PNG:", normalizePath(path.expand(out_png), mustWork = FALSE), "\n")
cat("Done.\n")