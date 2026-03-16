# Helper: .is_trueish_env
# purpose: parse common truthy values in environment variables
# inputs:
#   x = character scalar
# outputs:
#   logical scalar
.is_trueish_env <- function(x) {
  tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y", "on")
}

# Helper: skip_if_not_integration
# purpose: skip tests unless integration mode is explicitly enabled
# inputs:
#   env_var = character scalar; environment variable name to inspect
# outputs:
#   NULL (invisibly) if integration tests should run; otherwise skips current test/file
skip_if_not_integration <- function(env_var = "RUN_INTEGRATION_TESTS") {
  if (!.is_trueish_env(Sys.getenv(env_var, unset = "false"))) {
    testthat::skip(sprintf(
      "Integration tests are disabled. Set %s=true to run.",
      env_var
    ))
  }
}
