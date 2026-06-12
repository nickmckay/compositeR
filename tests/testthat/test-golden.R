# Golden-reference characterization tests. Reference outputs captured by
# tools/capture_reference.R from compositeR v0.2.0 (pre-stack-rewiring).
# These lock the numerical behavior of binning, spreading, standardization,
# and end-to-end compositing through the migration to the ens/lipdViz stack.

ref <- readRDS(test_path("reference_outputs.rds"))

# Temp12k ships with the package (LazyData)
usable <- which(sapply(Temp12k, function(x) {
  is.numeric(x$age) && is.null(dim(x$age)) &&
    length(x$paleoData_values) > 50 && !is.null(x$paleoData_TSid)
}))[1:8]
fTS <- Temp12k[usable]
binvec <- seq(0, 10000, by = 500)

makePdm <- function() {
  set.seed(9)
  ages <- seq(50, 9950, by = 100)
  pdm <- matrix(rnorm(length(ages) * 6,
                      mean = rep(c(10, 20, 5, 0, -5, 100), each = length(ages)),
                      sd = rep(c(1, 4, 2, 1, 0.5, 10), each = length(ages))),
                ncol = 6)
  pdm[3, 2] <- NA
  list(ages = ages, pdm = pdm)
}

test_that("simpleBinTs matches pre-migration reference", {
  expect_equal(simpleBinTs(fTS[[1]], binvec), ref$simpleBin)
})

test_that("sampleEnsembleThenBinTs matches pre-migration reference", {
  set.seed(7)
  expect_equal(sampleEnsembleThenBinTs(fTS[[1]], binvec), ref$sampleBin)
})

test_that("spreadPaleoData and removeConsecutiveDuplicates match reference", {
  sp <- spreadPaleoData(age = fTS[[1]]$age,
                        value = fTS[[1]]$paleoData_values,
                        spreadBy = 50, spreadMax = 500)
  expect_equal(sp, ref$spread)
  expect_equal(removeConsecutiveDuplicates(sp$spreadVal), ref$rcd)
})

test_that("standardization methods match pre-migration reference", {
  d <- makePdm()
  expect_equal(standardizeOverInterval(d$ages, d$pdm, interval = c(2000, 6000)),
               ref$stanInterval)
  set.seed(11)
  expect_equal(standardizeOverRandomInterval(d$ages, d$pdm, duration = 3000,
                                             searchRange = c(0, 9000)),
               ref$stanRandom)
  set.seed(12)
  expect_equal(standardizeMeanIteratively(d$ages, d$pdm, duration = 3000,
                                          searchRange = c(0, 9000)),
               ref$stanIterative)
  expect_equal(recordRMSE(td = ref$stanInterval, palMat = d$pdm), ref$rmse)
})

test_that("compositeEnsembles2 SCC-style matches pre-migration reference", {
  set.seed(20)
  out <- compositeEnsembles2(fTS, binvec, nens = 5,
                             stanFun = standardizeOverInterval,
                             interval = c(2000, 6000),
                             binFun = simpleBinTs,
                             verbose = FALSE)
  expect_equal(out, ref$compositeSCC)
  expect_s3_class(out, "paleoComposite")
})

test_that("compositeEnsembles2 DCC-style matches pre-migration reference", {
  set.seed(21)
  out <- compositeEnsembles2(fTS, binvec, nens = 5,
                             stanFun = standardizeOverRandomInterval,
                             duration = 3000,
                             searchRange = c(0, 9000),
                             binFun = sampleEnsembleThenBinTs,
                             verbose = FALSE)
  expect_equal(out, ref$compositeDCC)
  expect_s3_class(out, "paleoComposite")
})

test_that("compositeEnsembles2 exposes correctly-spelled parameters element", {
  set.seed(20)
  out <- compositeEnsembles2(fTS, binvec, nens = 5,
                             stanFun = standardizeOverInterval,
                             interval = c(2000, 6000),
                             binFun = simpleBinTs,
                             verbose = FALSE)
  expect_true("parameters" %in% names(out))
  # misspelled alias retained for one deprecation cycle, identical content
  expect_equal(out$paramaters, out$parameters)
})

test_that("summary.paleoComposite returns a per-bin quantile table", {
  set.seed(20)
  out <- compositeEnsembles2(fTS, binvec, nens = 5,
                             stanFun = standardizeOverInterval,
                             interval = c(2000, 6000),
                             binFun = simpleBinTs,
                             verbose = FALSE)
  s <- summary(out)
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), length(out$ages))
  expect_true(all(c("age", "q50", "nProxy") %in% names(s)))
})

test_that("compositeEnsembles warns once that it is superseded", {
  options(compositeR.compositeEnsembles.deprecated = NULL)
  expect_warning(
    compositeEnsembles(fTS, binvec, stanFun = standardizeOverInterval,
                       binFun = simpleBinTs, interval = c(2000, 6000)),
    "superseded by compositeEnsembles2"
  )
})
