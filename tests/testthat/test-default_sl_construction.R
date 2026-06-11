library(sl3)
# Test Script

# Create default Super Learner estimators used internally by EffectXshift
sls <- create_sls()

# create_sls() returns the two ensembles actually consumed by the pipeline:
#   - g_learner: a Super Learner for the (conditional density) exposure mechanism
#   - mu_learner: a stack of learners for the outcome regression
expect_equal(length(sls), 2)
expect_setequal(names(sls), c("g_learner", "mu_learner"))

# Check that the estimators are of the correct type
expect_true(inherits(sls$g_learner, "Lrnr_sl"))
expect_true(inherits(sls$mu_learner, "Stack"))
