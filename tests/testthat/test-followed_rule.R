context("Vectorised followed_rule()")

# The previous implementation, used as the reference
followed_rule_reference <- function(natural, shifted, A, mtp) {
  if (mtp) {
    return(rep(TRUE, nrow(natural)))
  }

  followed <- matrix(nrow = nrow(natural), ncol = length(A))

  for (i in seq_along(A)) {
    a <- A[i]
    followed[, i] <- mapply(function(x, y) isTRUE(all.equal(x, y)),
                            as.list(natural[, a]),
                            as.list(shifted[, a]))
  }

  apply(followed, 1, prod)
}

expect_same_as_reference <- function(natural, shifted, A, mtp = FALSE) {
  expect_identical(
    lmtp:::followed_rule(natural, shifted, A, mtp),
    followed_rule_reference(natural, shifted, A, mtp)
  )
}

test_that("binary and numeric treatments, including tolerance and missing values", {
  natural <- data.frame(
    A1 = c(0, 1, 1, 0, NA, NA, 1, 0),
    A2 = c(1.5, 2, 1e-10, 3, 4, NA, 1e5, 0),
    A3 = c(1L, 2L, 3L, NA, 5L, 6L, 7L, 8L)
  )
  shifted <- data.frame(
    A1 = c(0, 0, 1, 0, NA, 1, 1, 0),
    A2 = c(1.5 * (1 + 1e-9), 2.5, 2e-10, 3, NA, NA, 1e5 + 1, 1e-9),
    A3 = c(1, 2.1, 3, NA, 5, 6, 7, 8)
  )

  expect_same_as_reference(natural, shifted, "A1")
  expect_same_as_reference(natural, shifted, "A2")
  expect_same_as_reference(natural, shifted, "A3")
  expect_same_as_reference(natural, shifted, c("A1", "A2", "A3"))
})

test_that("infinite and NaN values", {
  natural <- data.frame(A = c(Inf, Inf, -Inf, NaN, NaN, 1, NA))
  shifted <- data.frame(A = c(Inf, 1, Inf, NaN, NA, NaN, NaN))
  expect_same_as_reference(natural, shifted, "A")
})

test_that("factor and character treatments", {
  natural <- data.frame(
    F1 = factor(c("a", "b", "c", NA), levels = c("a", "b", "c")),
    F2 = factor(c("a", "b", "c", "a"), levels = c("a", "b", "c")),
    C1 = c("x", "y", NA, "z"),
    stringsAsFactors = FALSE
  )
  shifted <- data.frame(
    F1 = factor(c("a", "c", "c", NA), levels = c("a", "b", "c")),
    F2 = factor(c("a", "b", "c", "a"), levels = c("c", "b", "a")),
    C1 = c("x", "x", NA, NA),
    stringsAsFactors = FALSE
  )

  expect_same_as_reference(natural, shifted, "F1")
  expect_same_as_reference(natural, shifted, "F2")
  expect_same_as_reference(natural, shifted, "C1")
  expect_same_as_reference(natural, shifted, c("F1", "C1"))
})

test_that("other types fall back to all.equal()", {
  natural <- data.frame(A = c(TRUE, FALSE, TRUE, NA))
  shifted <- data.frame(A = c(TRUE, TRUE, NA, NA))
  expect_same_as_reference(natural, shifted, "A")
})

test_that("modified treatment policies are always followed", {
  natural <- data.frame(A = c(0, 1))
  shifted <- data.frame(A = c(1, 1))
  expect_same_as_reference(natural, shifted, "A", mtp = TRUE)
})

test_that("random data", {
  set.seed(1)
  n <- 1000
  natural <- data.frame(A1 = rbinom(n, 1, 0.5), A2 = rnorm(n))
  natural$A1[sample(n, 50)] <- NA
  shifted <- natural
  change <- sample(n, 300)
  shifted$A1[change] <- 1 - shifted$A1[change]
  shifted$A2[change] <- shifted$A2[change] * (1 + sample(c(1e-10, 1e-6), 300, replace = TRUE))

  expect_same_as_reference(natural, shifted, c("A1", "A2"))
})
