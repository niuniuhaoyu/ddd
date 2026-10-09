options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("triplediff", quietly = TRUE)) {
  install.packages("triplediff", dependencies = TRUE)
}
cat("triplediff installed:", requireNamespace("triplediff", quietly = TRUE), "\n")
