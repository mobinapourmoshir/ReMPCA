################# Method for `hd` class #################
# Custom print method for the S3 object
print.hdClass <- function(x, ...) {
  cat("My Hybrid data Object:\n")
  print(unclass(x))  # Print the matrix without the class
  cat("\nCustom Attribute:\n")
  print(attr(x, "custom_attr"))
}
