# compositeEnsembles() puts records on one time scale using lipdR's time metadata

makeRecord <- function(id, ce = TRUE, ens = FALSE) {
  set.seed(1)
  yrs <- seq(1005, 1995, by = 10)
  vals <- sin(yrs / 50) + rnorm(length(yrs), sd = 0.1)
  r <- list(paleoData_TSid = id, dataSetName = id, archiveType = "Wood",
            paleoData_values = vals)
  if (ce) {
    r <- c(r, list(year = yrs, yearUnits = "CE", time = yrs, timeUnits = "yr AD",
                   timeDatum = 0, timeDirection = "prograde", timeExponent = 0))
  } else {
    r <- c(r, list(age = 1950 - yrs, ageUnits = "yr BP", time = 1950 - yrs, timeUnits = "yr BP",
                   timeDatum = 1950, timeDirection = "retrograde", timeExponent = 0))
  }
  if (ens) {
    r$ageEnsemble <- matrix(1950 - yrs, nrow = length(yrs), ncol = 3)
    r$ageEnsembleUnits <- "yr BP"
  }
  r
}

test_that("a CE record and the same record in yr BP give the same composite", {
  bv <- seq(0, 1000, by = 50)
  run <- function(fTS) {
    set.seed(3)
    compositeEnsembles(fTS, bv, nens = 2, binFun = simpleBinTs,
                       stanFun = standardizeOverInterval, interval = c(100, 900),
                       alignInterpDirection = FALSE, verbose = FALSE)
  }
  ce <- run(list(makeRecord("a", ce = TRUE), makeRecord("b", ce = TRUE)))
  bp <- run(list(makeRecord("a", ce = FALSE), makeRecord("b", ce = FALSE)))
  mixed <- run(list(makeRecord("a", ce = TRUE), makeRecord("b", ce = FALSE)))
  expect_equal(ce$composite, bp$composite)
  expect_equal(mixed$composite, bp$composite)
  expect_equal(ce$timeUnits, "yr BP")
  expect_true(any(is.finite(as.matrix(ce$composite))))
})

test_that("timeUnits sets the scale of binvec and the output", {
  bv <- seq(1000, 2000, by = 50)
  set.seed(3)
  out <- compositeEnsembles(list(makeRecord("a", ce = FALSE), makeRecord("b", ce = TRUE)),
                            bv, nens = 2, timeUnits = "CE", binFun = simpleBinTs,
                            stanFun = standardizeOverInterval, interval = c(1100, 1900),
                            alignInterpDirection = FALSE, verbose = FALSE)
  expect_equal(out$timeUnits, "CE")
  expect_equal(out$ages, seq(1025, 1975, by = 50))
  expect_true(any(is.finite(as.matrix(out$composite))))
  expect_output(printComposite(out), "CE")
  p <- plotComposite(out, combine = FALSE)
  expect_s3_class(p$compositeEns, "ggplot")
})

test_that("records without time metadata fall back to age, and missing time is named", {
  r <- makeRecord("a", ce = FALSE)
  r$time <- r$timeUnits <- r$timeDatum <- r$timeDirection <- r$timeExponent <- NULL
  expect_silent(fTS <- addCompositeTime(list(r), "yr BP"))
  expect_equal(fTS[[1]]$.compositeTime, r$age)
  r$age <- NULL
  expect_error(addCompositeTime(list(r), "yr BP"), "Can't put a on the yr BP scale")
})

test_that("an age ensemble is used, converted, unless binning with simpleBinTs", {
  r <- makeRecord("a", ce = TRUE, ens = TRUE)
  withEns <- addCompositeTime(list(r), "CE", ensemble = TRUE)[[1]]$.compositeTime
  expect_equal(dim(withEns), c(100, 3))
  expect_equal(withEns[, 1], r$year)
  noEns <- addCompositeTime(list(r), "CE", ensemble = FALSE)[[1]]$.compositeTime
  expect_equal(noEns, r$year)
})
