## =========================================================
## Combine MATLAB (.mat) AJIVE/DIVAS/MSSAT results and R PPD/BIDIFAC+ results
##
## Default paths:
##   MATLAB dir: ~/R/sim_dang_raw_full_parallel_24
##   R csv:      ~/R/ppd,bidifacplus.simul/raw_results_grid_ppd_bidifacplus.csv
## Output:
##   ~/R/combined_2block_all_methods_raw.csv
##
## Run:
##   source("make_combined_sim_csv_fixed.R")
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(R.matlab)
})

## -------------------------
## Paths
## -------------------------
mat_dir <- path.expand("~/R/sim_dang_raw_full_parallel_24")

r_csv_candidates <- c(
  path.expand("~/R/ppd,bidifacplus.simul/raw_results_grid_ppd_bidifacplus.csv"),
  path.expand("~/R/raw_results_grid_ppd_bidifacplus.csv"),
  path.expand("~/R/combined_R_raw_results_grid_ppd_bidifacplus.csv")
)

out_csv <- path.expand("~/R/combined_2block_all_methods_raw.csv")

## -------------------------
## MATLAB struct helpers
## R.matlab reads MATLAB struct fields through attr(x, "dimnames")[[1]],
## not through names(x).
## -------------------------
field_names <- function(s) {
  dn <- attr(s, "dimnames")
  if (is.null(dn) || length(dn) < 1 || is.null(dn[[1]])) character(0) else dn[[1]]
}

get_field <- function(s, field) {
  dn <- field_names(s)
  idx <- match(field, dn)
  if (is.na(idx)) {
    stop("No field: ", field, "; available fields: ", paste(dn, collapse = ", "))
  }
  s[[idx]]
}

get_field_or_null <- function(s, field) {
  dn <- field_names(s)
  idx <- match(field, dn)
  if (is.na(idx)) return(NULL)
  s[[idx]]
}

as_vec <- function(x) as.vector(x)

## -------------------------
## Read one MATLAB .mat result file
## -------------------------
mat_to_df_one <- function(file) {
  x <- R.matlab::readMat(file)
  
  if (!"rep.method.result" %in% names(x)) {
    stop(
      "No top-level variable rep.method.result in ", basename(file),
      "; top-level names: ", paste(names(x), collapse = ", ")
    )
  }
  
  s <- x[["rep.method.result"]]
  
  rep_idx  <- as.integer(get_field(s, "rep.idx")[1, 1])
  method   <- as.character(get_field(s, "method")[1, 1])
  d_vals   <- as.numeric(get_field(s, "d.vals"))
  ang_vals <- as.numeric(get_field(s, "ang.vals"))
  result   <- get_field(s, "result")
  
  ## top-level elapsed is one scalar per file
  elapsed_raw <- get_field_or_null(s, "elapsed")
  elapsed_file <- if (is.null(elapsed_raw)) NA_real_ else as.numeric(elapsed_raw[1, 1])
  
  ## result fields from your debug output:
  ## rJ, rI1, rI2, init.rank1, init.rank2, time, err.overall, err.joint, err.ind, error.msg
  rJ         <- get_field(result, "rJ")
  rI1        <- get_field(result, "rI1")
  rI2        <- get_field(result, "rI2")
  init_rank1 <- get_field(result, "init.rank1")
  init_rank2 <- get_field(result, "init.rank2")
  runtime    <- get_field(result, "time")
  eTotal     <- get_field(result, "err.overall")
  eJ         <- get_field(result, "err.joint")
  eI         <- get_field(result, "err.ind")
  
  grid <- expand.grid(
    d_idx = seq_along(d_vals),
    ang_idx = seq_along(ang_vals)
  )
  
  expected_n <- length(d_vals) * length(ang_vals)
  
  ## Dimensions should be 5 x 4 for ranks/time and 5 x 4 x 2 for errors.
  if (length(as_vec(rJ)) != expected_n) {
    stop("Unexpected rJ length in ", basename(file), ": got ", length(as_vec(rJ)),
         ", expected ", expected_n)
  }
  
  out <- grid %>%
    mutate(
      source = "MATLAB",
      rep = rep_idx,
      method = method,
      d = d_vals[d_idx],
      ang = ang_vals[ang_idx],
      
      rJ = as_vec(rJ),
      rI1 = as_vec(rI1),
      rI2 = as_vec(rI2),
      init_rank1 = as_vec(init_rank1),
      init_rank2 = as_vec(init_rank2),
      
      eJ1 = as_vec(eJ[, , 1]),
      eJ2 = as_vec(eJ[, , 2]),
      eI1 = as_vec(eI[, , 1]),
      eI2 = as_vec(eI[, , 2]),
      eTotal1 = as_vec(eTotal[, , 1]),
      eTotal2 = as_vec(eTotal[, , 2]),
      
      runtime = as_vec(runtime),
      time_sec_file = elapsed_file
    ) %>%
    select(
      source, rep, method, d, ang,
      rJ, rI1, rI2,
      init_rank1, init_rank2,
      eJ1, eJ2, eI1, eI2, eTotal1, eTotal2,
      runtime, time_sec_file
    )
  
  out
}

## -------------------------
## Read all MATLAB files
## -------------------------
read_all_matlab <- function(mat_dir) {
  if (!dir.exists(mat_dir)) {
    stop("MATLAB directory does not exist: ", mat_dir)
  }
  
  all_mat <- list.files(mat_dir, pattern = "\\.mat$", full.names = TRUE)
  files <- list.files(
    mat_dir,
    pattern = "^rep_[0-9]{3}_(ajive|divas|jive2)_result\\.mat$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  cat("MATLAB directory:", normalizePath(mat_dir, mustWork = FALSE), "\n")
  cat("Total .mat files:", length(all_mat), "\n")
  cat("Matched MATLAB result files:", length(files), "\n")
  
  if (length(files) == 0) {
    cat("First files in directory:\n")
    print(head(list.files(mat_dir), 30))
    stop("No MATLAB result files matched expected pattern.")
  }
  
  dfs <- vector("list", length(files))
  failed <- character(0)
  
  for (i in seq_along(files)) {
    dfs[[i]] <- tryCatch(
      mat_to_df_one(files[[i]]),
      error = function(e) {
        failed <<- c(failed, paste0(basename(files[[i]]), ": ", conditionMessage(e)))
        NULL
      }
    )
  }
  
  if (length(failed) > 0) {
    cat("\nFailed MATLAB files:", length(failed), "\n")
    cat(paste(head(failed, 20), collapse = "\n"), "\n")
    if (length(failed) > 20) cat("...\n")
  }
  
  out <- bind_rows(dfs)
  
  if (nrow(out) == 0) {
    stop("All MATLAB files failed. See failure messages above.")
  }
  
  out
}

## -------------------------
## Read R result CSV robustly
## -------------------------
find_r_csv <- function(candidates) {
  candidates <- unique(path.expand(candidates))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    cat("R CSV candidates checked:\n")
    print(candidates)
    stop("No R result CSV found. Edit r_csv_candidates near the top of this file.")
  }
  hit[[1]]
}

rename_if_exists <- function(df, old, new) {
  if (old %in% names(df) && !(new %in% names(df))) {
    names(df)[names(df) == old] <- new
  }
  df
}

standardize_r_results <- function(df) {
  if (!"source" %in% names(df)) df$source <- "R"
  
  df <- df %>%
    rename_if_exists("rep_idx", "rep") %>%
    rename_if_exists("replicate", "rep") %>%
    rename_if_exists("d_val", "d") %>%
    rename_if_exists("ang_val", "ang") %>%
    rename_if_exists("angle", "ang") %>%
    rename_if_exists("r_J", "rJ") %>%
    rename_if_exists("r_I1", "rI1") %>%
    rename_if_exists("r_I2", "rI2") %>%
    rename_if_exists("e_J1", "eJ1") %>%
    rename_if_exists("e_J2", "eJ2") %>%
    rename_if_exists("e_I1", "eI1") %>%
    rename_if_exists("e_I2", "eI2") %>%
    rename_if_exists("e_Total1", "eTotal1") %>%
    rename_if_exists("e_Total2", "eTotal2") %>%
    rename_if_exists("err_joint1", "eJ1") %>%
    rename_if_exists("err_joint2", "eJ2") %>%
    rename_if_exists("err_ind1", "eI1") %>%
    rename_if_exists("err_ind2", "eI2") %>%
    rename_if_exists("err_overall1", "eTotal1") %>%
    rename_if_exists("err_overall2", "eTotal2")
  
  needed <- c(
    "source", "rep", "method", "d", "ang",
    "rJ", "rI1", "rI2",
    "init_rank1", "init_rank2",
    "eJ1", "eJ2", "eI1", "eI2", "eTotal1", "eTotal2",
    "runtime", "time_sec_file"
  )
  
  for (v in needed) {
    if (!(v %in% names(df))) df[[v]] <- NA
  }
  
  df %>%
    mutate(source = "R") %>%
    select(all_of(needed), everything())
}

read_r_results <- function(r_csv_candidates) {
  file <- find_r_csv(r_csv_candidates)
  cat("R CSV:", normalizePath(file, mustWork = FALSE), "\n")
  df <- readr::read_csv(file, show_col_types = FALSE)
  standardize_r_results(df)
}

## -------------------------
## Main
## -------------------------
cat("============================================================\n")
cat("Combining MATLAB and R simulation results\n")
cat("============================================================\n")

mat_df <- read_all_matlab(mat_dir)
r_df <- read_r_results(r_csv_candidates)

cat("\nRows read from MATLAB:", nrow(mat_df), "\n")
cat("Rows read from R:", nrow(r_df), "\n")

combined <- bind_rows(mat_df, r_df) %>%
  mutate(
    method = recode(
      as.character(method),
      "ajive" = "AJIVE",
      "divas" = "DIVAS",
      "mssat" = "MSSAT",
      "PPD" = "PPD",
      "BIDIFAC+" = "BIDIFAC+",
      .default = as.character(method)
    ),
    true_rJ = if ("true_rJ" %in% names(.)) true_rJ else NA,
    true_rI1 = if ("true_rI1" %in% names(.)) true_rI1 else NA,
    true_rI2 = if ("true_rI2" %in% names(.)) true_rI2 else NA,
    correct_rJ = ifelse(is.na(rJ) | is.na(true_rJ), NA, rJ == true_rJ),
    correct_rI1 = ifelse(is.na(rI1) | is.na(true_rI1), NA, rI1 == true_rI1),
    correct_rI2 = ifelse(is.na(rI2) | is.na(true_rI2), NA, rI2 == true_rI2),
    eJ_mean = rowMeans(cbind(eJ1, eJ2), na.rm = TRUE),
    eI_mean = rowMeans(cbind(eI1, eI2), na.rm = TRUE),
    eTotal_mean = rowMeans(cbind(eTotal1, eTotal2), na.rm = TRUE)
  ) %>%
  arrange(source, method, rep, d, ang)

readr::write_csv(combined, out_csv)

combined %>%
  dplyr::count(source, method, name = "n") %>%
  as.data.frame()

cat("\nQuick check table(source):\n")
print(table(combined$source))

cat("\nDone.\n")