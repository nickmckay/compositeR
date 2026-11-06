# Tests for input handling: interpretation directions, tibble input, and
# exactly dated archives.

test_that("interpDirectionSign reads current LiPD interpretation fields", {
  ts <- list(interpretation1_direction = "Positive", interpretation1_scope = "climate")
  expect_equal(interpDirectionSign(ts), 1)
  ts$interpretation1_direction <- "negative"
  expect_equal(interpDirectionSign(ts), -1)
  ts$interpretation1_direction <- "decrease"
  expect_equal(interpDirectionSign(ts), -1)
  ts$interpretation1_direction <- "increase"
  expect_equal(interpDirectionSign(ts), 1)
  ts$interpretation1_direction <- "decrease/increase"
  expect_true(is.na(interpDirectionSign(ts)))
  ts$interpretation1_direction <- -1
  expect_equal(interpDirectionSign(ts), -1)
  expect_true(is.na(interpDirectionSign(list(a = 1))))
})

test_that("interpDirectionSign picks the interpretation matching scope", {
  ts <- list(interpretation1_direction = "positive", interpretation1_scope = "isotope",
             interpretation2_direction = "negative", interpretation2_scope = "climate")
  expect_equal(interpDirectionSign(ts, scope = "climate"), -1)
  expect_equal(interpDirectionSign(ts, scope = "isotope"), 1)
  expect_equal(interpDirectionSign(ts, scope = NA), 1)
  # unscoped interpretations are used when none match
  expect_equal(interpDirectionSign(list(interpretation1_direction = "negative"),
                                   scope = "climate"), -1)
  # but not when another scope is explicitly set
  expect_true(is.na(interpDirectionSign(ts[1:2], scope = "climate")))
})

test_that("interpDirectionSign reads legacy interpDirection fields", {
  expect_equal(interpDirectionSign(list(climateInterpretation1_interpDirection = "negative")), -1)
  expect_true(is.na(interpDirectionSign(list(isotopeInterpretation1_interpDirection = "negative"),
                                        scope = "climate")))
})

test_that("simpleBinTs flips records with a negative direction", {
  ts <- list(age = 1:100, paleoData_values = as.numeric(1:100),
             interpretation1_direction = "negative", interpretation1_scope = "climate")
  bv <- seq(0, 100, by = 10)
  up <- simpleBinTs(ts, bv, alignInterpDirection = FALSE)
  down <- simpleBinTs(ts, bv, alignInterpDirection = TRUE)
  expect_equal(down, -up)
})

test_that("sampleEnsembleThenBinTs does not perturb ages of tree-ring records", {
  bv <- seq(1000, 2000, by = 50)
  tr <- list(year = 1001:2000, paleoData_values = sin(1:1000 / 20),
             archiveType = "Wood")
  set.seed(1)
  a <- sampleEnsembleThenBinTs(tr, bv, ageVar = "year", defaultUnc = 0)
  set.seed(2)
  b <- sampleEnsembleThenBinTs(tr, bv, ageVar = "year", defaultUnc = 0)
  expect_equal(a, b)
  expect_equal(a, simpleBinTs(tr, bv, ageVar = "year"))
  # with exactAgeArchives = character(0), ages are simulated with BAM
  set.seed(1)
  c1 <- sampleEnsembleThenBinTs(tr, bv, ageVar = "year", defaultUnc = 0,
                                exactAgeArchives = character(0))
  expect_false(isTRUE(all.equal(a, c1)))
})

test_that("compositeEnsembles accepts a tibble of timeseries", {
  usable <- which(sapply(Temp12k, function(x) {
    is.numeric(x$age) && is.null(dim(x$age)) &&
      length(x$paleoData_values) > 50 && !is.null(x$paleoData_TSid)
  }))[1:8]
  fTS <- Temp12k[usable]
  binvec <- seq(0, 10000, by = 500)
  tib <- tibble::tibble(
    paleoData_TSid = purrr::map_chr(fTS, "paleoData_TSid"),
    age = purrr::map(fTS, "age"),
    paleoData_values = purrr::map(fTS, "paleoData_values")
  )
  set.seed(5)
  fromList <- compositeEnsembles(fTS, binvec, nens = 3, binFun = simpleBinTs,
                                 stanFun = standardizeOverInterval,
                                 interval = c(2000, 6000),
                                 alignInterpDirection = FALSE, verbose = FALSE)
  set.seed(5)
  fromTib <- compositeEnsembles(tib, binvec, nens = 3, binFun = simpleBinTs,
                                stanFun = standardizeOverInterval,
                                interval = c(2000, 6000),
                                alignInterpDirection = FALSE, verbose = FALSE)
  expect_equal(fromTib$composite, fromList$composite)
})
