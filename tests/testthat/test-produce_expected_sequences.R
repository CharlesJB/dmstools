library(dmstools)

reference <- "ACCTAG"

validate_test_results <- function(x) {
    expect_true(is(x, "data.frame"))
    expect_true(nrow(x) == (((nchar(reference)/3)*63)+1))
}

test_that("Valid params works", {
    validate_test_results(produce_expected_sequences(reference))
})
