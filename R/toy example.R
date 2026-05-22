## ==========================================================
## 3-block toy example
## ==========================================================

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(Matrix)
library(cowplot)
library(grid)

source("mssat.R")

## ---------------------------
## Helper functions
## ---------------------------

## Random orthonormal columns (p x r)
rand_orth <- function(p, r) {
  Z <- matrix(rnorm(p * r), p, r)
  qr.Q(qr(Z))
}

## Orthonormal basis orthogonal to given Q
rand_orth_orthogonal <- function(n, Q, r_perp) {
  if (r_perp == 0) return(matrix(0, n, 0))
  P_perp <- diag(n) - Q %*% t(Q)
  Z <- matrix(rnorm(n * r_perp), n, r_perp)
  qr.Q(qr(P_perp %*% Z))
}

## Simple rank-1 joint score (contrast vector)
make_joint_score <- function(n) {
  v <- c(rep(1, n/2), rep(-1, n - n/2))
  v / sqrt(sum(v^2))
}

## Convert matrix to long-form data frame
mat_to_df <- function(M, title) {
  as.data.frame(M) |>
    mutate(row = row_number()) |>
    pivot_longer(-row, names_to = "col", values_to = "val") |>
    mutate(col = as.integer(gsub("V","", col)),
           panel = title)
}

## Grouped pattern generator (for individual visibility)
group_vec <- function(n, g = 3) {
  labs <- rep(1:g, length.out = n)
  z <- scale(t(model.matrix(~ factor(labs) - 1)))  # g x n
  t(rand_orth(n, g) %*% z)[,1]
}

## ---------------------------
## Design (sizes, ranks, strengths)
## ---------------------------

set.seed(1)

n  <- 200
p1 <- 100; p2 <- 200; p3 <- 1200
r1 <- 2; r2 <- 2; r3 <- 3

d1 <- sqrt(n) * c(0.8, 0.4)
d2 <- sqrt(n) * c(200, 100)
d3 <- sqrt(n) * c(16, 12, 6)

sig1 <- 0.4; sig2 <- 120; sig3 <- 6

## ---------------------------
## Joint & individual structures
## ---------------------------

vJ <- c(rep(-1,n/2),rep(1,n/2))
vJ <- vJ/norm(vJ,type = "2")

VI1 <- c(rep(-1,n/5),rep(1,n/5),rep(0,n/5),rep(1,n/5),rep(-1,n/5))
VI1 <- VI1/norm(VI1,type = "2")

VI2 <- c(rep(-1,n/4),rep(1,n/4),rep(-1,n/4),rep(1,n/4))
VI2 <- VI2/norm(VI2,type = "2")

VI3 <- matrix(nrow = n, ncol = 2)
VI3[,1] <- c(rep(-1,n/5),rep(1,n/5),rep(0,n/5),rep(-1,n/5),rep(1,n/5))
VI3[,2] <- c(rep(-1,n/4),rep(1,n/4),rep(1,n/4),rep(-1,n/4))

VI3[,1] <- VI3[,1]/norm(VI3[,1],type = "2")
VI3[,2] <- VI3[,2]/norm(VI3[,2],type = "2")

U1 <- rand_orth(p1, r1)
U1 <- cbind(c(rep(0,p1*2/4),rep(1,p1*2/4)),
            c(rep(0,p1*2/4),rep(1,p1*1/4),rep(0,p1*1/4)))

U2 <- rand_orth(p2, r2)
U2 <- cbind(c(rep(0,p2*1/5),rep(1,p2*3/5),rep(0,p2*1/5)),
            c(rep(0,p2*1/5),rep(0,p2*2/5),rep(1,p2*2/5)))

U3 <- rand_orth(p3, r3)
U3 <- cbind(c(rep(1,p3*2/6),rep(0,p3*2/6),rep(0,p3*2/6)),
            c(rep(0,p3*2/6),rep(0,p3*2/6),rep(1,p3*2/6)),
            c(rep(0,p3*1/6),rep(1,p3*1/6),rep(0,p3*4/6)))

V1 <- cbind(vJ, VI1)
V2 <- cbind(vJ, VI2)
V3 <- cbind(vJ, VI3)

D1 <- diag(d1, r1)
D2 <- diag(d2, r2)
D3 <- diag(d3, r3)

J1 <- U1[,1,drop=FALSE] %*% d1[1] %*% t(V1[,1,drop=FALSE])
I1 <- if (r1 > 1) U1[,2:r1] %*% diag(d1[2:r1], nrow = r1-1) %*% t(V1[,2:r1]) else matrix(0, p1, n)

J2 <- U2[,1,drop=FALSE] %*% d2[1] %*% t(V2[,1,drop=FALSE])
I2 <- if (r2 > 1) U2[,2:r2] %*% diag(d2[2:r2], nrow = r2-1) %*% t(V2[,2:r2]) else matrix(0, p2, n)

J3 <- U3[,1,drop=FALSE] %*% d3[1] %*% t(V3[,1,drop=FALSE])
I3 <- if (r3 > 1) U3[,2:r3] %*% diag(d3[2:r3], nrow = r3-1) %*% t(V3[,2:r3]) else matrix(0, p3, n)

E1 <- matrix(rnorm(p1*n, 0, sig1), p1, n)
E2 <- matrix(rnorm(p2*n, 0, sig2), p2, n)
E3 <- matrix(rnorm(p3*n, 0, sig3), p3, n)

X1 <- J1 + I1 + E1
X2 <- J2 + I2 + E2
X3 <- J3 + I3 + E3

## ---------------------------
## Color range utility
## ---------------------------

rng_block <- function(...) {
  mats <- list(...)
  m <- max(sapply(mats, function(M) max(abs(M))))
  c(-m, m)
}

rX1 <- rng_block(X1); rJ1 <- rng_block(J1); rI1 <- rng_block(I1); rE1 <- rng_block(E1)
rX2 <- rng_block(X2); rJ2 <- rng_block(J2); rI2 <- rng_block(I2); rE2 <- rng_block(E2)
rX3 <- rng_block(X3); rJ3 <- rng_block(J3); rI3 <- rng_block(I3); rE3 <- rng_block(E3)

scale_for <- function(block, type) {
  if (block == "Block 1") {
    lim <- switch(type, "Observed"=rX1,"Joint"=rJ1,"Individual"=rI1,"Noise"=rE1)
  } else if (block == "Block 2") {
    lim <- switch(type, "Observed"=rX2,"Joint"=rJ2,"Individual"=rI2,"Noise"=rE2)
  } else {
    lim <- switch(type, "Observed"=rX3,"Joint"=rJ3,"Individual"=rI3,"Noise"=rE3)
  }
  scale_fill_gradient2(low="#313695", mid="white", high="#A50026",
                       midpoint=0, limits=lim, oob=scales::squish)
}

## ---------------------------
## Data frame for ggplot
## ---------------------------

df1 <- bind_rows(
  mat_to_df(X1, "Observed (X₁)"),
  mat_to_df(J1, "Joint (J₁)"),
  mat_to_df(I1, "Individual (I₁)"),
  mat_to_df(E1, "Noise (E₁)")
) |> mutate(block = "Block 1")

df2 <- bind_rows(
  mat_to_df(X2, "Observed (X₂)"),
  mat_to_df(J2, "Joint (J₂)"),
  mat_to_df(I2, "Individual (I₂)"),
  mat_to_df(E2, "Noise (E₂)")
) |> mutate(block = "Block 2")

df3 <- bind_rows(
  mat_to_df(X3, "Observed (X₃)"),
  mat_to_df(J3, "Joint (J₃)"),
  mat_to_df(I3, "Individual (I₃)"),
  mat_to_df(E3, "Noise (E₃)")
) |> mutate(block = "Block 3")

df <- bind_rows(df1, df2, df3) |>
  mutate(panel_type = case_when(
    grepl("^Observed", panel) ~ "Observed",
    grepl("^Joint", panel) ~ "Joint",
    grepl("^Individual", panel) ~ "Individual",
    TRUE ~ "Noise"
  ))

## ---------------------------
## Normalize coordinates for uniform panel sizes
## ---------------------------

df_norm <- df %>%
  group_by(block, panel_type) %>%
  mutate(
    row_norm = (row - min(row)) / (max(row) - min(row)),
    col_norm = (col - min(col)) / (max(col) - min(col))
  ) %>%
  ungroup()


library(ggplot2)
library(patchwork)
library(cowplot)
library(grid)

## ---------------------------
common_theme <- theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.8, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.margin = margin(t = -6, unit = "pt"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    plot.margin = margin(0, 0, 0, 0)
  )

## ---------------------------
scale_for_block <- function(block) {
  if (block == "Block 1") lim <- rng_block(X1, J1, I1, E1)
  else if (block == "Block 2") lim <- rng_block(X2, J2, I2, E2)
  else lim <- rng_block(X3, J3, I3, E3)
  
  min_val <- as.integer(round(lim[1]))
  max_val <- as.integer(round(lim[2]))
  
  scale_fill_gradient2(
    low = "#313695", mid = "white", high = "#A50026",
    midpoint = 0, limits = lim, oob = scales::squish,
    guide = guide_colorbar(
      title = NULL,
      barwidth = 3.2,
      barheight = 0.3,
      ticks.colour = "black",
      label.position = "bottom",
      label.theme = element_text(size = 8),
      breaks = c(min_val, max_val),
      labels = c(min_val, max_val),
      nbin = 10
    )
  )
}

## ---------------------------
draw_block_row <- function(dfi, block_label) {
  block_name <- unique(dfi$block)[1]
  types <- c("Observed", "Joint", "Individual", "Noise")
  
  plist <- lapply(types, function(tp) {
    ggplot(filter(dfi, panel_type == tp),
           aes(x = col_norm, y = row_norm, fill = val)) +
      geom_raster(interpolate = FALSE) +
      scale_for_block(block_name) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_fixed(ratio = 1.8) +
      common_theme
  })
  
  plots_row <- wrap_plots(plist, nrow = 1)
  
  label_grob <- ggdraw() +
    draw_label(block_label, angle = 90, fontface = "bold",
               size = 11, x = 0.99, hjust = 0.2, vjust = 0.5)
  
  cowplot::plot_grid(
    label_grob, plots_row,
    ncol = 2, rel_widths = c(0.03, 1),
    align = "v", axis = "tb"
  )
}

## ---------------------------
row1 <- draw_block_row(filter(df_norm, block == "Block 1"), "Block 1")
row2 <- draw_block_row(filter(df_norm, block == "Block 2"), "Block 2")
row3 <- draw_block_row(filter(df_norm, block == "Block 3"), "Block 3")

grid_blocks <- (row1 / row2 / row3) + plot_layout(heights = c(1, 1, 1))

## ---------------------------
col_labels <- ggdraw() +
  draw_label("Observed",   x = 0.19, y = 0.3, fontface = "bold", size = 11) +
  draw_label("Joint",      x = 0.405, y = 0.3, fontface = "bold", size = 11) +
  draw_label("Individual", x = 0.62, y = 0.3, fontface = "bold", size = 11) +
  draw_label("Noise",      x = 0.84, y = 0.3, fontface = "bold", size = 11)

## ---------------------------
final_plot <- plot_grid(
  col_labels,   # column labels
  grid_blocks,  # main plots
  ncol = 1,
  rel_heights = c(0.035, 1)
)

final_plot

## ---------------------------
## Run MSSAT
## ---------------------------
mssat_result <- mssat(list(X1,X2,X3), center=FALSE, alpha=0.05,
                      rank_method = "given", r_init = c(2, 2, 3),
                      sigma_method = "given", sigma_init = (c(2,15,4)^2), verbose=FALSE)


## ---------------------------
## Reconstruction
## ---------------------------
reconstruct_block <- function(Xk, blk) {
  P_VAk <- blk$P_VAk
  P_VAk_given_VJ <- blk$P_VAk_given_VJ
  
  Jk <- Xk %*% P_VAk_given_VJ
  Ik <- Xk %*% (P_VAk - P_VAk_given_VJ)
  Ek <- Xk - Xk %*% P_VAk
  
  list(
    X = Xk,
    J = Jk,
    I = Ik,
    E = Ek
  )
}

rec1 <- reconstruct_block(X1, mssat_result$blocks[[1]])
rec2 <- reconstruct_block(X2, mssat_result$blocks[[2]])
rec3 <- reconstruct_block(X3, mssat_result$blocks[[3]])

## ---------------------------
rng_block <- function(...) range(c(...), na.rm = TRUE)

## ---------------------------
df_build <- function(rec, block_name) {
  mats <- list(Observed = rec$X, Joint = rec$J, Individual = rec$I, Residual = rec$E)
  do.call(rbind, lapply(names(mats), function(nm) {
    df <- as.data.frame(as.table(mats[[nm]]))
    names(df) <- c("row", "col", "val")
    df$panel_type <- nm
    df$block <- block_name
    df
  }))
}

df_norm <- rbind(
  df_build(rec1, "Block 1"),
  df_build(rec2, "Block 2"),
  df_build(rec3, "Block 3")
)
df_norm$row_norm <- as.numeric(df_norm$row)
df_norm$col_norm <- as.numeric(df_norm$col)

## ---------------------------
common_theme <- theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.8, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.margin = margin(t = -6, unit = "pt"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    plot.margin = margin(0, 0, 0, 0)
  )

## ---------------------------
scale_for_block <- function(block) {
  lim <- switch(block,
                "Block 1" = rng_block(rec1$X, rec1$J, rec1$I, rec1$E),
                "Block 2" = rng_block(rec2$X, rec2$J, rec2$I, rec2$E),
                "Block 3" = rng_block(rec3$X, rec3$J, rec3$I, rec3$E))
  
  min_val <- as.integer(round(lim[1]))
  max_val <- as.integer(round(lim[2]))
  
  scale_fill_gradient2(
    low = "#313695", mid = "white", high = "#A50026",
    midpoint = 0, limits = lim, oob = scales::squish,
    guide = guide_colorbar(
      title = NULL,
      barwidth = 3.2,
      barheight = 0.3,
      ticks.colour = "black",
      label.position = "bottom",
      label.theme = element_text(size = 8),
      breaks = c(min_val, max_val),
      labels = c(min_val, max_val),
      nbin = 10
    )
  )
}

## ---------------------------
draw_block_row <- function(dfi, block_label) {
  block_name <- unique(dfi$block)[1]
  types <- c("Observed", "Joint", "Individual", "Residual")
  
  ## block-specific aspect ratio
  p_k <- switch(block_name,
                "Block 1" = nrow(X1),
                "Block 2" = nrow(X2),
                "Block 3" = nrow(X3))
  n_k <- ncol(X1) 
  ratio_k <- 1.8 / (p_k / n_k)  
  
  plist <- lapply(types, function(tp) {
    ggplot(subset(dfi, panel_type == tp),
           aes(x = col_norm, y = row_norm, fill = val)) +
      geom_raster(interpolate = FALSE) +
      scale_for_block(block_name) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_fixed(ratio = ratio_k) + 
      common_theme
  })
  
  plots_row <- wrap_plots(plist, nrow = 1)
  
  label_grob <- ggdraw() +
    draw_label(block_label, angle = 90, fontface = "bold",
               size = 11, x = 0.99, hjust = 0.2, vjust = 0.5)
  
  cowplot::plot_grid(
    label_grob, plots_row,
    ncol = 2, rel_widths = c(0.03, 1),
    align = "v", axis = "tb"
  )
}

## ---------------------------
row1 <- draw_block_row(subset(df_norm, block == "Block 1"), "Block 1")
row2 <- draw_block_row(subset(df_norm, block == "Block 2"), "Block 2")
row3 <- draw_block_row(subset(df_norm, block == "Block 3"), "Block 3")

grid_blocks <- (row1 / row2 / row3) + plot_layout(heights = c(1, 1, 1))

## ---------------------------
col_labels <- ggdraw() +
  draw_label("Observed",   x = 0.19, y = 0.3, fontface = "bold", size = 11) +
  draw_label("Joint",      x = 0.405, y = 0.3, fontface = "bold", size = 11) +
  draw_label("Individual", x = 0.62,  y = 0.3, fontface = "bold", size = 11) +
  draw_label("Residual",   x = 0.84,  y = 0.3, fontface = "bold", size = 11)

## ---------------------------
final_plot <- plot_grid(
  col_labels,   # column labels
  grid_blocks,  # main plots
  ncol = 1,
  rel_heights = c(0.035, 1)
)

final_plot