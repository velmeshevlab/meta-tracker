# Function to classify expression pattern with tunable parameters
classify_expression <- function(pseudotime, expression, threshold = 0.1, min_duration = 0.1, transient_threshold = 0.3) {
  # Fit a smoothing spline to the data
  spline <- smooth.spline(pseudotime, expression, spar = 1)
  
  # Calculate the first derivative
  first_derivative <- predict(spline, pseudotime, deriv = 1)$y
  
  # Normalize pseudotime and set thresholds for sustained increase/decrease
  pseudotime_range <- max(pseudotime) - min(pseudotime)
  min_duration_points <- min_duration * pseudotime_range
  
  # Determine the pattern based on the derivative and thresholds
  sustained_increase <- sum(first_derivative > threshold) >= min_duration_points
  sustained_decrease <- sum(first_derivative < -threshold) >= min_duration_points
  
  # Check if the gene is gradually increasing
  if (sustained_increase) {
    return("Increasing")
  } else if (sustained_decrease) {
    return("Decreasing")
  } else {
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
  }
  
  # Check for gradual increase without local min/max by comparing start and end
  start_expression <- expression[1]
  end_expression <- expression[length(expression)]
  if ((end_expression - start_expression) > threshold * max(expression)) {
    return("Increasing")
  }
  
  # Check for plateauing
  if (all(abs(first_derivative) < threshold) && !sustained_increase && !sustained_decrease) {
    return("Plateauing")
  }
  
  return("Unclassified")
}
