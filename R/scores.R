# Patient-level score helpers shared by steps 07 and 07b.

# Mean of a vector; NA (not NaN, not an error) when there is nothing to average.
safe_mean <- function(x) if (length(x) == 0 || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# Score = mean over the available markers of each cell, then mean over the cells of each patient.
# `subset_df`: one row per cell with `patient_id` and one column per marker; `markers`: marker names
# (those absent from `subset_df` are ignored). Returns `patient_id`, `score`, `n_cells`.
score_by_patient <- function(subset_df, markers) {
  avail <- intersect(markers, colnames(subset_df))
  if (length(avail) == 0) return(tibble::tibble(patient_id = character(), score = numeric()))
  subset_df %>%
    dplyr::group_by(patient_id) %>%
    dplyr::summarise(score = safe_mean(rowMeans(dplyr::across(dplyr::all_of(avail)), na.rm = TRUE)),
                     n_cells = dplyr::n(), .groups = "drop")
}
