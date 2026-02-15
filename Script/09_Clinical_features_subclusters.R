#=========================Script Description=================================
# This script is used for clinical features visualization
#===========================Loading Packages=============================
suppressMessages(library("optparse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("gridExtra"))
suppressMessages(library("ggthemes"))
suppressMessages(library("scales"))
suppressMessages(library("circlize"))
#===========================Function Definition=============================
ha_colors <- list(
  sex = c(Female = "#E78AC3", Male = "#66C2A5"),
  condition = c(bulbar = "#FEE090", spinal = "#FDAE61"),
  progression_group = c(SP = "#1B9E77", IP = "#D95F02", FP = "#7570B3"),
  cluster = c(alpha = "#FC8D62", beta = "#8DA0CB", theta = "#A6D854")
)

make_df_for_plot <- function(se_obj, cluster_csv) {
  clus <- read.csv(cluster_csv, stringsAsFactors = FALSE)
  colnames(clus)[2] <- "k2"
  df <- as.data.frame(SummarizedExperiment::colData(se_obj))
  df$k2 <- clus$k2[match(df$label, clus$patid)]
  df
}

plot_cat_by_k2 <- function(
    df, k2_col = "k2", var_col = "condition",
    ha_colors = NULL, output_dir = ".", file_name = NULL,
    linesize = 0.5, table_base_size = 12, return_plot_only = FALSE) {
  
  stopifnot(all(c(k2_col, var_col) %in% names(df)))
  d <- subset(df, !is.na(df[[k2_col]]) & !is.na(df[[var_col]]))
  d[[k2_col]]  <- factor(d[[k2_col]])
  d[[var_col]] <- factor(d[[var_col]])
  
  tab <- table(k2 = d[[k2_col]], var = d[[var_col]])
  chi_try <- suppressWarnings(chisq.test(tab, correct = FALSE))
  if (any(chi_try$expected < 5)) {
    test_res  <- fisher.test(tab)
    test_type <- "Fisher's exact test"
  } else {
    test_res  <- chi_try
    test_type <- "Chi-squared test"
  }
  
  df_plot <- as.data.frame(tab)
  colnames(df_plot) <- c("cluster", "variable", "Freq")
  
  coeff <- data.frame(
    Method   = test_type,
    `p-value` = sprintf("%.3f", test_res$p.value),
    check.names = FALSE
  )
  tbl_grob <- tableGrob(coeff, rows = NULL,
                        theme = ttheme_minimal(base_size = table_base_size))
  
  if (!is.null(ha_colors) && var_col %in% names(ha_colors)) {
    pal <- ha_colors[[var_col]]
    pal <- pal[levels(d[[var_col]])]
  } else {
    nlev <- nlevels(d[[var_col]])
    if (nlev > 8) stop("Too many levels for default palette; please supply ha_colors.")
    pal <- RColorBrewer::brewer.pal(8, "Set2")[seq_len(nlev)]
  }
  
  #ttl <- sprintf("%s distribution by %s", var_col, k2_col)
  ylab <- sprintf("proportions of %s", var_col)
  
  p <- ggplot(df_plot, aes(fill = variable, y = Freq, x = cluster)) +
    geom_bar(position = "fill", stat = "identity", color = "white", width = 0.7) +
    annotation_custom(tbl_grob, ymin = 0, ymax = 0.2) +
    scale_fill_manual(values = pal, name = var_col) +
    scale_y_continuous(labels = percent_format(accuracy = 10),
                       expand = expansion(mult = c(0, 0.08))) +
    #labs(title = ttl, x = k2_col, y = ylab) +
    labs( x = k2_col, y = ylab) +
    theme_classic(base_size = 14) +
    theme(panel.background = element_rect(size = linesize),
          plot.title = element_text(hjust = 0),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
  
  if (!return_plot_only) {
    if (is.null(file_name)) {
      file_name <- sprintf("k2_vs_%s.svg", var_col)
    }
    save_path <- file.path(output_dir, file_name)
    ggsave(save_path, plot = p, width = 4.5, height = 4.2, dpi = 300)
  } else {
    save_path <- NULL
  }
  
  result <- list(plot = p, table = tab, proportions = prop.table(tab, margin = 1),
                test_type = test_type, test = test_res)
  if (!return_plot_only) result$save_path <- save_path
  result
}


.choose_test_two_groups <- function(y, g, mode = c("auto","parametric","nonparametric")) {
  mode <- match.arg(mode)
  g <- droplevels(factor(g))
  if (nlevels(g) != 2) stop("choose_test_two_groups: g must have exactly 2 groups.")
  grp <- split(y, g)
  n1 <- sum(!is.na(grp[[1]])); n2 <- sum(!is.na(grp[[2]]))
  
  if (mode == "parametric") {
    # force parametric path
    lev <- tryCatch(car::leveneTest(y ~ g), error = function(e) NULL)
    if (!is.null(lev) && lev$`Pr(>F)`[1] > 0.05) {
      return(list(name = "Student's t-test", test = t.test(y ~ g, var.equal = TRUE)))
    } else {
      return(list(name = "Welch's t-test",   test = t.test(y ~ g)))
    }
  }
  if (mode == "nonparametric") {
    return(list(name = "Wilcoxon rank-sum",  test = wilcox.test(y ~ g)))
  }
  
  # auto: need both groups n>=3 to trust Shapiro; otherwise Wilcoxon
  if (n1 < 3 || n2 < 3) {
    return(list(name = "Wilcoxon rank-sum",  test = wilcox.test(y ~ g)))
  }
  
  sw1 <- tryCatch(shapiro.test(grp[[1]]), error = function(e) NULL)
  sw2 <- tryCatch(shapiro.test(grp[[2]]), error = function(e) NULL)
  normal1 <- !is.null(sw1) && sw1$p.value > 0.05
  normal2 <- !is.null(sw2) && sw2$p.value > 0.05
  
  if (normal1 && normal2) {
    lev <- tryCatch(car::leveneTest(y ~ g), error = function(e) NULL)
    if (!is.null(lev) && lev$`Pr(>F)`[1] > 0.05) {
      list(name = "Student's t-test", test = t.test(y ~ g, var.equal = TRUE))
    } else {
      list(name = "Welch's t-test",   test = t.test(y ~ g))
    }
  } else {
    list(name = "Wilcoxon rank-sum",  test = wilcox.test(y ~ g))
  }
}

.format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) return(formatC(p, format="e", digits=2))
  sprintf("%.3f", p)
}

plot_cont_by_k2 <- function(
    df,
    k2_col    = "k2",
    var_col   = "age",
    ha_colors = NULL,          # for cluster fill
    output_dir = ".",
    file_name  = NULL,
    linesize = 0.5,
    table_base_size = 12,
    test_pref = c("auto","parametric","nonparametric"),
    return_plot_only = FALSE
){
  test_pref <- match.arg(test_pref)
  stopifnot(all(c(k2_col, var_col) %in% names(df)))
  d <- subset(df, !is.na(df[[k2_col]]) & !is.na(df[[var_col]]))
  d[[k2_col]]  <- droplevels(factor(d[[k2_col]]))
  if (nlevels(d[[k2_col]]) != 2) stop("plot_cont_by_k2: k2 must have exactly 2 groups for this auto test.")
  d[[var_col]] <- suppressWarnings(as.numeric(d[[var_col]]))
  if (all(is.na(d[[var_col]]))) stop(sprintf("Variable '%s' is not numeric.", var_col))
  
  # auto choose test
  tst <- .choose_test_two_groups(d[[var_col]], d[[k2_col]], mode = test_pref)
  test_name <- tst$name
  pval <- tst$test$p.value
  
  # annotation table
  coeff <- data.frame(
    Method   = test_name,
    `p-value` = .format_p(pval),
    check.names = FALSE
  )
  tbl_grob <- gridExtra::tableGrob(coeff, rows = NULL,
                                   theme = gridExtra::ttheme_minimal(base_size = table_base_size))
  
  # palette for cluster fill
  if (!is.null(ha_colors) && "cluster" %in% names(ha_colors)) {
    pal <- ha_colors[["cluster"]]
    pal <- pal[levels(d[[k2_col]])]
  } else {
    nlev <- nlevels(d[[k2_col]])
    pal <- RColorBrewer::brewer.pal(8, "Set2")[seq_len(nlev)]
  }
  
  #ttl  <- sprintf("%s by %s", var_col, k2_col)
  ylab <- var_col
  
  p <- ggplot2::ggplot(d, ggplot2::aes(y = .data[[var_col]], x = .data[[k2_col]], fill = .data[[k2_col]])) +
    ggplot2::geom_violin(trim = FALSE, color = "black", fill = "white") +
    ggplot2::geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.8,
                          color = "black", ggplot2::aes(fill = .data[[k2_col]])) +
    ggplot2::scale_fill_manual(values = pal, name = k2_col) +
    #ggplot2::labs(title = ttl, x = k2_col, y = ylab) +
    ggplot2::labs( x = k2_col, y = ylab) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(panel.background = ggplot2::element_rect(size = linesize),
                   plot.title = ggplot2::element_text(hjust = 0),
                   axis.text = ggplot2::element_text(size = 12),
                   axis.title = ggplot2::element_text(size = 14),
                   legend.text = ggplot2::element_text(size = 12),
                   legend.title = ggplot2::element_text(size = 14))
  
  # place the annotation table near the bottom
  y_min <- min(d[[var_col]], na.rm = TRUE)
  y_span <- diff(range(d[[var_col]], na.rm = TRUE))
  p <- p + ggplot2::annotation_custom(tbl_grob, ymin = y_min, ymax = y_min + 0.25*y_span)
  
  if (!return_plot_only) {
    if (is.null(file_name)) file_name <- sprintf("k2_vs_%s.svg", var_col)
    save_path <- file.path(output_dir, file_name)
    ggplot2::ggsave(save_path, plot = p, width = 4.8, height = 4.2, dpi = 300)
  } else {
    save_path <- NULL
  }
  
  result <- list(plot = p, test_name = test_name, p_value = pval)
  if (!return_plot_only) result$save_path <- save_path
  result
}
# ===========================Command Parameters Setting=============================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "02_Missing_Inspection_subclusters/norm_imp_MinProb.rds",
              help = "SummarizedExperiment object after normalize and imputation."),
  make_option(c("--cluster", "-c"),
              type = "character",default = "08_Clustering_als/cluster_assignments_2.csv",
              help = "Cluster assignment CSV file (include: patid, k2)."),
  make_option(c("--output", "-o"),
              type = "character", default = "09_Clinical_features_subclusters",
              help = "Output directory path."),
  make_option(c("--var", "-v"),
              type = "character", default = NULL,
              help = "Single variable column name to compare against k2 (deprecated, use --vars)."),
  make_option(c("--vars", "-V"),
              type = "character", default = NULL,
              help = "Comma-separated variable column names to compare against k2 (e.g., 'age,age_at_onset,Nfl,pNFh')."),
  make_option(c("--merge", "-m"),
              type = "character", default = "FALSE",
              help = "Whether to merge multiple plots into one page (TRUE/FALSE)."),
  make_option(c("--row", "-r"),
              type = "integer", default = NULL,
              help = "Number of rows for merged plot grid (optional, auto-calculated if not specified)."),
  make_option(c("--col", "-C"),
              type = "integer", default = NULL,
              help = "Number of columns for merged plot grid (optional, auto-calculated if not specified)."),
  make_option(c("--seed", "-e"),
              type = "integer", default = 9,
              help = "Random seed.")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
#============================================================================
if (is.null(opt$output)) {
  print("NO OUTPUT PATH SUPPLIED,current directory will be used!")
  output_dir <- getwd()
} else {
  output_dir <- opt$output
  if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
  }
}

if (is.null(opt$input)) {
  stop("Please provide the cleaned SummarizedExperiment object after normalize and imputation file path!")
}else if (!file.exists(opt$input)) {
  stop("SummarizedExperiment object after normalize and imputation does not exist!")
}else{
  input = opt$input
  se = readRDS(input)
}

if (is.null(opt$cluster)){
  stop("Please provide the cluster CSV file path!")
} else if (!file.exists(opt$cluster)) {
  stop("Cluster CSV file does not exist!")
}else{
  cluster_csv = opt$cluster
}

if (is.null(opt$seed)) {
  stop("Please provide the seed number!")
}else{
  seed = opt$seed
  set.seed(seed)
}

# Determine variables to plot
if (!is.null(opt$vars)) {
  var_cols <- trimws(strsplit(opt$vars, ",")[[1]])
} else if (!is.null(opt$var)) {
  var_cols <- opt$var
} else {
  stop("Please provide either --var or --vars with variable column name(s) to compare against k2!")
}

# Determine merge option
merge_plots <- if (!is.null(opt$merge)) {
  toupper(as.character(opt$merge)) %in% c("TRUE", "T", "YES", "Y", "1")
} else {
  FALSE
}

# If merge is TRUE, must have multiple variables
if (merge_plots && length(var_cols) == 1) {
  warning("Only one variable provided but merge=TRUE. Setting merge=FALSE.")
  merge_plots <- FALSE
}

# -------------------
# Main script
# -------------------
df <- make_df_for_plot(se, cluster_csv)

if (merge_plots && length(var_cols) > 1) {
  # Merge multiple plots into one page
  plot_list <- list()
  
  for (var_col in var_cols) {
    if (!var_col %in% names(df)) {
      warning(sprintf("Variable '%s' not found in data, skipping.", var_col))
      next
    }
    
    if (is.character(df[[var_col]]) || is.factor(df[[var_col]])) {
      result <- plot_cat_by_k2(df, k2_col = "k2",
                               var_col = var_col,
                               ha_colors = ha_colors, 
                               output_dir = output_dir,
                               return_plot_only = TRUE)
    } else {
      result <- plot_cont_by_k2(df,
                                k2_col = "k2",
                                var_col = var_col,
                                test_pref = "auto",
                                ha_colors = ha_colors,
                                output_dir = output_dir,
                                return_plot_only = TRUE)
    }
    plot_list[[var_col]] <- result$plot
  }
  
  if (length(plot_list) == 0) {
    stop("No valid variables found to plot!")
  }
  
  # Arrange plots in a grid
  n_plots <- length(plot_list)
  
  # Use user-specified rows/cols if provided, otherwise auto-calculate
  if (!is.null(opt$row) && !is.null(opt$col)) {
    nrow <- opt$row
    ncol <- opt$col
    if (nrow * ncol < n_plots) {
      warning(sprintf("Specified grid (%d x %d = %d) has fewer cells than plots (%d). Some plots may be omitted.", 
                     nrow, ncol, nrow * ncol, n_plots))
    }
  } else if (!is.null(opt$row)) {
    nrow <- opt$row
    ncol <- ceiling(n_plots / nrow)
  } else if (!is.null(opt$col)) {
    ncol <- opt$col
    nrow <- ceiling(n_plots / ncol)
  } else {
    # Auto-calculate grid layout
    if (n_plots <= 2) {
      ncol <- n_plots
      nrow <- 1
    } else if (n_plots <= 4) {
      ncol <- 2
      nrow <- 2
    } else if (n_plots <= 6) {
      ncol <- 3
      nrow <- 2
    } else {
      ncol <- 3
      nrow <- ceiling(n_plots / 3)
    }
  }
  
  # Create merged plot
  merged_plot <- do.call(gridExtra::arrangeGrob, 
                        c(plot_list, ncol = ncol, nrow = nrow))
  
  # Save merged plot
  file_name <- paste0("k2_vs_merged_", paste(var_cols, collapse = "_"), ".svg")
  save_path <- file.path(output_dir, file_name)
  ggplot2::ggsave(save_path, plot = merged_plot, 
                 width = 4.8 * ncol, 
                 height = 4.2 * nrow, 
                 dpi = 300)
  
  cat(sprintf("Merged plot saved to: %s\n", save_path))
  
} else {
  # Plot each variable separately
  for (var_col in var_cols) {
    if (!var_col %in% names(df)) {
      warning(sprintf("Variable '%s' not found in data, skipping.", var_col))
      next
    }
    
    if (is.character(df[[var_col]]) || is.factor(df[[var_col]])) {
      plot_cat_by_k2(df, k2_col = "k2",
                     var_col = var_col,
                     ha_colors = ha_colors, 
                     output_dir = output_dir)
    } else {
      plot_cont_by_k2(df,
                      k2_col = "k2",
                      var_col = var_col,
                      test_pref = "auto",
                      ha_colors = ha_colors,
                      output_dir = output_dir)
    }
  }
}



