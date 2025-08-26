# Variable Standardization for Mixed Data Types

# Function to identify variable types 
identify_variable_types <- function(data, max_unique_for_discrete = 10) {
  var_types <- list(
    binary = character(),
    categorical = character(),
    continuous = character()
  )
  
  for (col in names(data)) {
    x <- data[[col]]
    unique_vals <- unique(x[!is.na(x)])
    n_unique <- length(unique_vals)
    
    # Check if character/factor
    if (is.character(x) || is.factor(x)) {
      if (n_unique == 2) {
        var_types$binary <- c(var_types$binary, col)
      } else {
        var_types$categorical <- c(var_types$categorical, col)
      }
    }
    # Check if numeric
    else if (is.numeric(x)) {
      # Binary numeric
      if (n_unique == 2) {
        var_types$binary <- c(var_types$binary, col)
      }
      # Likely categorical (few unique values)
      else if (n_unique <= max_unique_for_discrete && all(x %% 1 == 0, na.rm = TRUE)) {
        var_types$categorical <- c(var_types$categorical, col)
      }
      # Continuous
      else {
        var_types$continuous <- c(var_types$continuous, col)
      }
    }
  }
  
  return(var_types)
}

# Function to convert discrete variables to numeric format
convert_discrete_to_numeric <- function(data, binary_vars, categorical_vars) {
  data_converted <- data
  
  # Convert binary variables to 0/1
  for (var in binary_vars) {
    x <- data[[var]]
    
    if (is.character(x) || is.factor(x)) {
      # Convert to factor first if character
      if (is.character(x)) x <- as.factor(x)
      # Convert to 0/1
      data_converted[[var]] <- as.numeric(x) - 1
    } else if (is.numeric(x)) {
      # Already numeric, ensure it's 0/1
      unique_vals <- sort(unique(x[!is.na(x)]))
      if (!all(unique_vals %in% c(0, 1))) {
        # Map to 0/1
        data_converted[[var]] <- as.numeric(x == unique_vals[2])
      }
    }
  }
  
  # Convert categorical variables to numeric codes
  for (var in categorical_vars) {
    x <- data[[var]]
    
    if (is.character(x) || is.factor(x)) {
      # Convert to factor first if character
      if (is.character(x)) x <- as.factor(x)
      # Convert to numeric codes (1, 2, 3, ...)
      data_converted[[var]] <- as.numeric(x)
    }
    # If already numeric, keep as is
  }
  
  return(data_converted)
}

# Function to create a mapping dictionary for discrete variables
create_discrete_mappings <- function(data, binary_vars, categorical_vars) {
  mappings <- list()
  
  # Map binary variables
  for (var in binary_vars) {
    x <- data[[var]]
    unique_vals <- sort(unique(x[!is.na(x)]))
    
    if (is.character(x) || is.factor(x)) {
      mappings[[var]] <- data.frame(
        original = as.character(unique_vals),
        coded = c(0, 1),
        stringsAsFactors = FALSE
      )
    } else {
      mappings[[var]] <- data.frame(
        original = unique_vals,
        coded = c(0, 1),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Map categorical variables
  for (var in categorical_vars) {
    x <- data[[var]]
    
    if (is.character(x) || is.factor(x)) {
      if (is.character(x)) x <- as.factor(x)
      mappings[[var]] <- data.frame(
        original = levels(x),
        coded = seq_along(levels(x)),
        stringsAsFactors = FALSE
      )
    }
  }
  
  return(mappings)
}

# Main standardization function
standardize_mixed_data <- function(data, max_unique_for_discrete = 10, 
                                   center = TRUE, scale = TRUE) {
  
  # Identify variable types
  var_types <- identify_variable_types(data, max_unique_for_discrete)
  
  # Print variable type summary
  cat("Variable Type Summary:\n")
  cat("Binary variables (", length(var_types$binary), "):", 
      paste(var_types$binary, collapse = ", "), "\n")
  cat("Categorical variables (", length(var_types$categorical), "):", 
      paste(var_types$categorical, collapse = ", "), "\n")
  cat("Continuous variables (", length(var_types$continuous), "):", 
      paste(var_types$continuous, collapse = ", "), "\n\n")
  
  # Create mappings for discrete variables (before conversion)
  mappings <- create_discrete_mappings(data, var_types$binary, var_types$categorical)
  
  # Convert discrete variables to numeric
  data_processed <- convert_discrete_to_numeric(data, var_types$binary, var_types$categorical)
  
  # Standardize continuous variables
  if (length(var_types$continuous) > 0) {
    data_processed[var_types$continuous] <- scale(data_processed[var_types$continuous], 
                                                  center = center, 
                                                  scale = scale)
  }
  
  # Return results with metadata
  result <- list(
    data = data_processed,
    var_types = var_types,
    mappings = mappings,
    scaling_params = if(length(var_types$continuous) > 0) {
      list(
        means = attr(scale(data[var_types$continuous]), "scaled:center"),
        sds = attr(scale(data[var_types$continuous]), "scaled:scale")
      )
    } else NULL
  )
  
  class(result) <- "standardized_data"
  return(result)
}