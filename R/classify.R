# Function to classify expression pattern with tunable parameters
classify_expression <- function(pseudotime, expression, threshold = 0.2, min_duration = 0.2, transient_threshold = 0.3, plateau_fraction = 0.2) {
  # Fit a smoothing spline to the data
  spline <- smooth.spline(pseudotime, expression, spar = 1)
  
  # Calculate the first derivative
  first_derivative <- predict(spline, pseudotime, deriv = 1)$y
  
  # Normalize pseudotime and set thresholds for sustained increase/decrease
  pseudotime_range <- max(pseudotime) - min(pseudotime)
  min_duration_points <- min_duration * pseudotime_range
  min_plateau_points <- plateau_fraction * length(pseudotime)
  
  # Check for multiple peaks (transient pattern)
  peaks <- which(diff(sign(diff(expression))) == -2) + 1
  if (length(peaks) > 0) {
    peak_idx <- which.max(expression)
    if (peak_idx > 1 && peak_idx < length(expression)) {
      # Identify baseline before peak and minimum after peak (if it exists)
      baseline_before_peak <- expression[1]
      if (peak_idx < length(expression)) {
        min_after_peak <- min(expression[(peak_idx + 1):length(expression)])
        
        # Check if the minimum after the peak is within the transient threshold of the baseline
        if (abs(min_after_peak - baseline_before_peak) <= transient_threshold * max(expression)) {
          return("Transient") 
        }
      } else {
        # If there is no clear minimum after the peak, consider it transient if it returns close to baseline
        if (abs(expression[length(expression)] - baseline_before_peak) <= transient_threshold * max(expression)) {
          return("Transient")
        }
      }
    }
  }
  
  # Check for plateauing: if the change after the local max is below the threshold and occupies a sufficient fraction of pseudotime
  local_max_idx <- which.max(expression)
  if (local_max_idx > 1 && local_max_idx < length(expression)) {
    after_max_change <- abs(expression[length(expression)] - expression[local_max_idx])
    if (after_max_change < threshold * max(expression) && (length(expression) - local_max_idx) >= min_plateau_points) {
      return("Plateau")
    }
  }
  
  # Determine the pattern based on the derivative and thresholds
  sustained_increase <- sum(first_derivative > threshold) >= min_duration_points
  sustained_decrease <- sum(first_derivative < -threshold) >= min_duration_points
  
  # Check if the gene is gradually increasing without local min/max
  if (length(peaks) == 0 && (expression[length(expression)] - expression[1]) > threshold * max(expression)) {
    return("Increasing")
  }
  
  # Check if the gene is gradually increasing
  if (sustained_increase) {
    return("Increasing")
  } else if (sustained_decrease) {
    return("Decreasing")
  }
  
  return("Unclassified")
}
