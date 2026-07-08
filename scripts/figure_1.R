## =============================================================================
## Manhattan Plot Grid (4 models x 3 phenotypes) for TWAS Results
##
## Generates a 4x3 grid of Manhattan plots comparing TWAS association results
## across four expression models (PrediXcan, Epithelial, Fibroblast, Adipocyte)
## and three mammographic density phenotypes (Dense Area, Non-Dense Area,
## Percent Density). Points are colored by whether the association was
## detected using cell-type-specific (CTS) or non-specific (NS) models.
##
## Input: one CSV per (model, phenotype) combination, each containing at least
##   file path: MiXcan2/Data/figure_1
##   the columns:
##     CHR        - chromosome number
##     BP         - base-pair position (used to order variants within a chromosome)
##     SNP        - variant ID (used to label significant points)
##     P_nonlog   - raw (non-log) p-value; plotted as -log10(P_nonlog) on the y-axis
##     p_adjusted - multiple-testing-adjusted p-value; determines which points are
##                  called significant (p_adjusted <= 0.05), which get labeled, and
##                  where the dashed significance-threshold line is drawn
##     CTS        - flag (>= 0.5) marking a cell-type-specific association
##     NS         - flag (>= 0.5) marking a non-specific (non-cell-type-specific) association
##   
## Output: a single combined PNG figure (MiXcan2/Data/figure_1/fig.png)
## =============================================================================

library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(tidyr)

## ---------------------------------------------------------------------------
## 0) Config -- EDIT THESE for your local setup
## ---------------------------------------------------------------------------

data_dir   <- "data"                                  # folder holding input CSVs
output_png <- "figure_1.png"

# path_map[[model]][[phenotype]] = "path/to/file.csv"
path_map <- list(
  "PrediXcan" = list(
    "DA"  = file.path(data_dir, "predixcan_DA.csv"),
    "NDA" = file.path(data_dir, "predixcan_NDA.csv"),
    "PD"  = file.path(data_dir, "predixcan_PD.csv")
  ),
  "Epithelial" = list(
    "DA"  = file.path(data_dir, "epithelial_DA.csv"),
    "NDA" = file.path(data_dir, "epithelial_NDA.csv"),
    "PD"  = file.path(data_dir, "epithelial_PD.csv")
  ),
  "Fibroblast" = list(
    "DA"  = file.path(data_dir, "fibroblast_DA.csv"),
    "NDA" = file.path(data_dir, "fibroblast_NDA.csv"),
    "PD"  = file.path(data_dir, "fibroblast_PD.csv")
  ),
  "Adipose" = list(
    "DA"  = file.path(data_dir, "adipose_DA.csv"),
    "NDA" = file.path(data_dir, "adipose_NDA.csv"),
    "PD"  = file.path(data_dir, "adipose_PD.csv")
  )
)

models      <- c("PrediXcan", "Epithelial", "Fibroblast", "Adipose")
row_labels  <- dplyr::recode(models, Adipose = "Adipocyte")   # display label per row
outcomes    <- c("DA", "NDA", "PD")
outcome_map <- c("Dense Area" = "DA", "Non-Dense Area" = "NDA", "Percent Density" = "PD")

sig_alpha  <- 0.05                    # p_adjusted threshold for labeling/highlighting SNPs
y_breaks   <- seq(2.5, 12.5, by = 2.5)

point_colors <- c("CTS" = "#ff7f0e", "NS" = "#1e90ff", "Other" = "grey70")

## ---------------------------------------------------------------------------
## 1) Load & prepare a single result file for plotting
##    (adds cumulative genomic position for a standard Manhattan x-axis)
## ---------------------------------------------------------------------------

load_and_prepare <- function(file) {
  df <- readr::read_csv(file, show_col_types = FALSE) %>%
    mutate(
      CHR          = as.integer(CHR),
      BP           = as.numeric(BP),
      minus_log10p = -log10(P_nonlog),
      source = case_when(
        CTS >= 0.5 ~ "CTS",
        NS  >= 0.5 ~ "NS",
        TRUE       ~ "Other"
      )
    ) %>%
    filter(!is.na(CHR)) %>%
    arrange(CHR, BP) %>%
    group_by(CHR) %>%
    mutate(BP_index = row_number()) %>%
    ungroup()
  
  chrom_offsets <- df %>%
    group_by(CHR) %>%
    summarise(max_bp = max(BP_index, na.rm = TRUE), .groups = "drop") %>%
    arrange(CHR) %>%
    mutate(offset_val = lag(cumsum(max_bp + ifelse(CHR >= 20, 1000, 200)), default = 0))
  
  df <- df %>%
    left_join(chrom_offsets, by = "CHR") %>%
    mutate(cum_pos = BP_index + offset_val)
  
  axis_df <- df %>%
    group_by(CHR) %>%
    summarise(center = mean(cum_pos, na.rm = TRUE), .groups = "drop")
  
  list(data = df, axis_df = axis_df)
}

## ---------------------------------------------------------------------------
## 2) Single-panel Manhattan plot
## ---------------------------------------------------------------------------

manhattan_plot <- function(prep, title = "", label_alpha = sig_alpha,
                           x_breaks, x_labels, x_lim, y_lim,
                           show_y = TRUE, show_title = TRUE, show_legend = FALSE) {
  d       <- prep$data
  sig     <- d %>% filter(!is.na(p_adjusted), p_adjusted <= label_alpha)
  nonsig  <- d %>% filter(is.na(p_adjusted) | p_adjusted > label_alpha)
  
  # Draw a dashed significance-threshold line midway between the lowest
  # significant point and the highest non-significant point, if both exist.
  threshold_y <- NA
  if (nrow(sig) > 0 && nrow(nonsig) > 0) {
    lowest_sig  <- min(sig$minus_log10p, na.rm = TRUE)
    next_nonsig <- nonsig %>%
      filter(minus_log10p < lowest_sig) %>%
      summarise(val = max(minus_log10p, na.rm = TRUE)) %>%
      pull(val)
    if (is.finite(next_nonsig)) threshold_y <- (lowest_sig + next_nonsig) / 2
  }
  
  p <- ggplot(d, aes(x = cum_pos, y = minus_log10p, color = source)) +
    geom_point(size = 1.5, alpha = 0.4) +
    scale_color_manual(name = "Cell Type", values = point_colors,
                       labels = names(point_colors)) +
    ggrepel::geom_text_repel(
      data = sig, aes(label = SNP),
      size = 5.5, fontface = "bold", nudge_y = 0.6, max.overlaps = 50,
      box.padding = 0.6, point.padding = 0.4, segment.angle = 10,
      show.legend = FALSE
    ) +
    scale_x_continuous(labels = x_labels, breaks = x_breaks, limits = x_lim,
                       expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = y_lim, expand = c(0, 0), breaks = y_breaks) +
    labs(
      title = if (show_title) title else NULL,
      x = "Chromosome",
      y = expression(-log[10](italic(p)))
    ) +
    theme(
      legend.position  = if (show_legend) "bottom" else "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(color = "black", linewidth = 0.5),
      axis.ticks       = element_line(color = "black", linewidth = 0.5),
      axis.text        = element_text(color = "black", size = 14),
      axis.title       = element_text(color = "black", size = 16),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 22),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
  
  if (!show_y) {
    p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
  }
  if (!is.na(threshold_y)) {
    p <- p + geom_hline(yintercept = threshold_y, linetype = "dashed", color = "red")
  }
  
  p
}

## ---------------------------------------------------------------------------
## 3) Load every (model, phenotype) file
## ---------------------------------------------------------------------------

preps <- list()
for (out in outcomes) {
  for (mod in models) {
    file <- path_map[[mod]][[out]]
    key  <- paste(mod, out, sep = "||")
    
    if (is.null(file) || !is.character(file)) {
      message("Skipping ", key, ": no file path configured")
      next
    }
    if (!file.exists(file)) {
      message("Skipping ", key, ": file not found at ", file)
      next
    }
    preps[[key]] <- load_and_prepare(file)
  }
}

if (length(preps) == 0) stop("No input files were found. Check `data_dir` and `path_map`.")

## ---------------------------------------------------------------------------
## 4) Recompute chromosome offsets globally so all panels share one x-axis
## ---------------------------------------------------------------------------

chr_widths <- bind_rows(lapply(preps, function(pr) pr$data)) %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP_index, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(
    padding    = ifelse(CHR >= 20, 400, 200),
    offset_val = lag(cumsum(max_bp + padding), default = 0),
    center     = offset_val + max_bp / 2
  )

for (nm in names(preps)) {
  preps[[nm]]$data <- preps[[nm]]$data %>%
    dplyr::select(-cum_pos, -offset_val) %>%
    left_join(chr_widths %>% dplyr::select(CHR, offset_val), by = "CHR") %>%
    mutate(cum_pos = BP_index + offset_val)
  preps[[nm]]$axis_df <- chr_widths %>% dplyr::select(CHR, center)
}

global_x_lim <- c(0, max(chr_widths$offset_val + chr_widths$max_bp))
x_breaks     <- chr_widths$center
x_labels     <- chr_widths$CHR

global_y_lim <- c(
  0,
  ceiling(max(sapply(preps, function(pr) max(pr$data$minus_log10p, na.rm = TRUE)))) + 1
)

## ---------------------------------------------------------------------------
## 5) Build one panel per (model, phenotype) cell of the grid
## ---------------------------------------------------------------------------

nrows <- length(models)
ncols <- length(outcomes)

plot_list <- list()
for (row_idx in seq_along(models)) {
  for (col_idx in seq_along(outcomes)) {
    key <- paste(models[row_idx], outcomes[col_idx], sep = "||")
    pr  <- preps[[key]]
    if (is.null(pr)) next
    
    p <- manhattan_plot(
      pr,
      title       = names(outcome_map)[col_idx],
      label_alpha = sig_alpha,
      x_breaks    = x_breaks,
      x_labels    = x_labels,
      x_lim       = global_x_lim,
      y_lim       = global_y_lim,
      show_title  = (row_idx == 1)
    )
    plot_list[[length(plot_list) + 1]] <- p
  }
}

# Tighten margins; only keep y-axis title on the leftmost column
plot_list_tight <- lapply(seq_along(plot_list), function(i) {
  p <- plot_list[[i]] + theme(plot.margin = margin(1, 1, 1, 1))
  col <- ((i - 1) %% ncols) + 1
  if (col != 1) p <- p + theme(axis.title.y = element_blank())
  p
})

## ---------------------------------------------------------------------------
## 6) Assemble the full figure: title, row labels, grid, shared legend
## ---------------------------------------------------------------------------

# Shared legend, pulled from one panel and styled up for the combined figure
legend_source <- plot_list[[ncols + 1]]
legend <- cowplot::get_legend(
  legend_source +
    theme(
      legend.position    = "right",
      legend.text        = element_text(size = 16),
      legend.title       = element_text(size = 18, face = "bold"),
      legend.key.size    = unit(2, "lines"),
      legend.box.margin  = margin(0, 0, 0, 0),
      legend.background  = element_rect(fill = "white", color = NA)
    )
)

grid_4x3 <- cowplot::plot_grid(
  plotlist    = plot_list_tight,
  nrow        = nrows,
  ncol        = ncols,
  align       = "hv",
  rel_widths  = rep(4, ncols),
  rel_heights = rep(1, nrows),
  greedy      = TRUE
)

# Row labels (model names), rotated vertically along the left edge
label_x  <- 0.8
left_col <- cowplot::ggdraw() + theme_void()
for (i in seq_len(nrows)) {
  left_col <- left_col +
    cowplot::draw_label(
      row_labels[i],
      x = label_x, y = 1 - (i - 0.5) / nrows + 0.03,
      hjust = 1, vjust = 0.5, fontface = "bold", size = 22, angle = 90
    )
}

grid_with_left <- cowplot::plot_grid(
  left_col, grid_4x3,
  ncol = 2, rel_widths = c(0.10, 0.90), align = "h"
)

title <- cowplot::ggdraw() +
  draw_label(
    "Manhattan Plots of TWAS Performed across DA, NDA, PD Phenotypes using Epithelial, Fibroblast, Adipocyte, and PrediXcan Models",
    fontface = "bold", size = 24, x = 0.5, hjust = 0.5
  )

plots_with_legend <- cowplot::plot_grid(
  grid_with_left, legend, ncol = 2, rel_widths = c(0.9, 0.1)
)

final_grid <- cowplot::plot_grid(
  title, plots_with_legend, ncol = 1, rel_heights = c(0.08, 1)
)

final_grid <- ggdraw(final_grid) +
  theme(plot.background = element_rect(fill = "white", color = NA))

## ---------------------------------------------------------------------------
## 7) Save
## ---------------------------------------------------------------------------

ggsave(output_png, final_grid, width = 26, height = 18, dpi = 600)
=
