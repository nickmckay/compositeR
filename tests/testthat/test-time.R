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
    compositeEnsembles(fTS, bv, nens = 2, binFun = simpleBinTs, time.units = "yr BP",
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

test_that("time.units sets the scale of binvec and the output", {
  bv <- seq(1000, 2000, by = 50)
  set.seed(3)
  out <- compositeEnsembles(list(makeRecord("a", ce = FALSE), makeRecord("b", ce = TRUE)),
                            bv, nens = 2, time.units = "CE", binFun = simpleBinTs,
                            stanFun = standardizeOverInterval, interval = c(1100, 1900),
                            alignInterpDirection = FALSE, verbose = FALSE)
  expect_equal(out$timeUnits, "yr AD")
  expect_equal(out$ages, seq(1025, 1975, by = 50))
  expect_true(any(is.finite(as.matrix(out$composite))))
  expect_output(printComposite(out), "yr AD")
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


ceBp <- function() list(makeRecord("a", ce = TRUE), makeRecord("b", ce = FALSE))

test_that("compositing on CE and on ka gives the same composite", {
  bvCE <- seq(1000, 2000, by = 50)
  bvKa <- sort((1950 - bvCE) / 1000)
  run <- function(units, bv, interval, spread){
    set.seed(3)
    compositeEnsembles(ceBp(), bv, nens = 2, time.units = units, binFun = simpleBinTs, spread = spread,
                       stanFun = standardizeOverInterval, interval = interval,
                       alignInterpDirection = FALSE, verbose = FALSE)
  }
  flip <- function(m) as.matrix(m)[nrow(m):1, ]
  ce <- run("CE", bvCE, c(1100, 1900), FALSE)
  ka <- run("ka", bvKa, (1950 - c(1100, 1900)) / 1000, FALSE)
  expect_equal(ka$timeUnits, "ka")
  expect_equal(rev(ka$ages), (1950 - ce$ages) / 1000)
  expect_true(any(is.finite(as.matrix(ce$composite))))
  expect_equal(flip(ka$composite), as.matrix(ce$composite), ignore_attr = TRUE)

  # spreading is done on one orientation internally, so it gives the same
  # composite on every axis too
  ceS <- run("CE", bvCE, c(1100, 1900), TRUE)
  kaS <- run("ka", bvKa, (1950 - c(1100, 1900)) / 1000, TRUE)
  bpS <- run("yr BP", 1950 - rev(bvCE), 1950 - c(1900, 1100), TRUE)
  expect_true(any(is.finite(as.matrix(ceS$composite))))
  expect_equal(flip(kaS$composite), as.matrix(ceS$composite), ignore_attr = TRUE)
  expect_equal(flip(bpS$composite), as.matrix(ceS$composite), ignore_attr = TRUE)
  expect_equal(as.matrix(kaS$composite), as.matrix(bpS$composite), ignore_attr = TRUE)
})

test_that("spreadPaleoData gives the same result on yr BP, ka and CE axes", {
  set.seed(5)
  # uneven sampling, with grid points exactly halfway between samples (ties)
  bp <- c(12, 30, 48, 100, 140, 141, 300, 420, 610, 800)
  v <- rnorm(length(bp))
  sBP <- spreadPaleoData(bp, v, spreadBy = 1, spreadMax = 100)
  sKa <- spreadPaleoData(bp / 1000, v, spreadBy = 0.001, spreadMax = 0.1, time.units = "ka")
  sCE <- spreadPaleoData(1950 - bp, v, spreadBy = 1, spreadMax = 100, time.units = "CE")
  expect_equal(sKa$spreadAge * 1000, sBP$spreadAge)
  expect_equal(sKa$spreadVal, sBP$spreadVal)
  expect_equal(rev(1950 - sCE$spreadAge), sBP$spreadAge)
  expect_equal(rev(sCE$spreadVal), sBP$spreadVal)
  # output is in increasing order on its own axis
  expect_false(is.unsorted(stats::na.omit(sCE$spreadAge)))
  # the same with the default spreadBy and an explicit newAge
  expect_equal(rev(spreadPaleoData(1950 - bp, v, time.units = "CE")$spreadVal),
               spreadPaleoData(bp, v)$spreadVal)
  na <- seq(0, 800, by = 5)
  expect_equal(rev(spreadPaleoData(1950 - bp, v, newAge = 1950 - na, time.units = "CE")$spreadVal),
               spreadPaleoData(bp, v, newAge = na)$spreadVal)
})

test_that("spreadPaleoData honors maxAge and minAge on any axis", {
  bp <- seq(100, 1000, by = 100)
  v <- seq_along(bp)
  sBP <- spreadPaleoData(bp, v, spreadBy = 10, minAge = 200, maxAge = 500)
  expect_equal(range(sBP$spreadAge, na.rm = TRUE), c(200, 500))
  expect_equal(range(sBP$spreadVal, na.rm = TRUE), c(2, 5))
  sKa <- spreadPaleoData(bp / 1000, v, spreadBy = 0.01, minAge = 0.2, maxAge = 0.5, time.units = "ka")
  expect_equal(range(sKa$spreadAge, na.rm = TRUE), c(0.2, 0.5))
  # on a CE axis maxAge is the earliest year kept and minAge the latest
  sCE <- spreadPaleoData(1950 - bp, v, spreadBy = 10, minAge = 1750, maxAge = 1450, time.units = "CE")
  expect_equal(range(sCE$spreadAge, na.rm = TRUE), c(1450, 1750))
  expect_equal(rev(sCE$spreadVal), sBP$spreadVal)
  # a maxAge older than the record doesn't extend it
  expect_equal(spreadPaleoData(bp, v, spreadBy = 10, maxAge = 5000), spreadPaleoData(bp, v, spreadBy = 10))
})

test_that("time.variable picks a time column other than the primary one", {
  r <- makeRecord("a", ce = TRUE)
  r$calAge <- (1950 - r$year) / 1000
  r[c("calAgeUnits", "calAgeTimeUnits", "calAgeTimeDatum", "calAgeTimeDirection", "calAgeTimeExponent", "calAgeTimeSource")] <-
    list("ka", "ka", 1950, "retrograde", 3, "units")
  fTS <- addCompositeTime(list(r), "yr BP", time.variable = "calAge")
  expect_equal(fTS[[1]]$.compositeTime, 1950 - r$year)
})

test_that("ageVar and timeUnits still work, with a warning", {
  expect_warning(
    out <- compositeEnsembles(ceBp(), seq(1000, 2000, by = 50), nens = 1, timeUnits = "CE", ageVar = "time",
                              binFun = simpleBinTs, stanFun = standardizeOverInterval, interval = c(1100, 1900),
                              alignInterpDirection = FALSE, verbose = FALSE),
    "deprecated")
  expect_equal(out$timeUnits, "yr AD")
  expect_warning(simpleBinTs(makeRecord("a", ce = FALSE), seq(0, 1000, 50), ageVar = "age"), "deprecated")
})

test_that("without time.units, records sharing a scale stay on it, and mixed ones go to yr BP", {
  ce2 <- list(makeRecord("a", ce = TRUE), makeRecord("b", ce = TRUE))
  set.seed(3)
  expect_message(out <- compositeEnsembles(ce2, seq(1000, 2000, by = 50), nens = 1, binFun = simpleBinTs,
                                           stanFun = standardizeOverInterval, interval = c(1100, 1900),
                                           alignInterpDirection = FALSE, verbose = FALSE),
                 "Compositing on yr AD")
  expect_equal(out$timeUnits, "yr AD")
  expect_message(out2 <- compositeEnsembles(ceBp(), seq(0, 1000, by = 50), nens = 1, binFun = simpleBinTs,
                                            stanFun = standardizeOverInterval, interval = c(100, 900),
                                            alignInterpDirection = FALSE, verbose = FALSE),
                 "different time scales")
  expect_equal(out2$timeUnits, "yr BP")
})

test_that("a binvec on the wrong scale is caught", {
  expect_error(compositeEnsembles(ceBp(), seq(0, 1000, by = 50), nens = 1, time.units = "CE", binFun = simpleBinTs,
                                  stanFun = standardizeOverInterval, interval = c(100, 900),
                                  alignInterpDirection = FALSE, verbose = FALSE),
               "doesn't overlap")
})

test_that("assumed time scales give one warning", {
  r1 <- makeRecord("a", ce = FALSE); r2 <- makeRecord("b", ce = FALSE)
  for(nm in c("time", "timeUnits", "timeDatum", "timeDirection", "timeExponent", "ageUnits")){ r1[[nm]] <- NULL; r2[[nm]] <- NULL }
  w <- character(0)
  withCallingHandlers(addCompositeTime(list(r1, r2), "yr BP"),
                      warning = function(cnd){ w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_length(w, 1)
  expect_match(w, "2 record\\(s\\) was assumed")
})

test_that("spreading on a CE axis keeps early years, and stops at 2025 CE", {
  sp <- spreadPaleoData(age = c(-500, -400, -300), value = 1:3, spreadBy = 10, time.units = "CE")
  expect_true(any(is.finite(sp$spreadVal)))
  sp2 <- spreadPaleoData(age = c(2000, 2010, 2020, 2030), value = 1:4, spreadBy = 1, time.units = "CE")
  expect_true(max(sp2$spreadAge, na.rm = TRUE) <= 2025)
  sp3 <- spreadPaleoData(age = c(0.1, 0.2, 0.3), value = 1:3, spreadBy = 0.01, time.units = "ka")
  expect_gt(sum(is.finite(sp3$spreadVal)), 20)
})

test_that("plotComposite converts to another scale", {
  set.seed(3)
  out <- compositeEnsembles(ceBp(), seq(0, 1000, by = 50), nens = 2, time.units = "yr BP", binFun = simpleBinTs,
                            stanFun = standardizeOverInterval, interval = c(100, 900),
                            alignInterpDirection = FALSE, verbose = FALSE)
  p <- plotComposite(out, time.units = "CE", combine = FALSE)
  expect_equal(range(ggplot2::layer_data(p$dataUsed)$x), range(1950 - out$ages))
  expect_warning(plotComposite(out, ageUnits = "CE", combine = FALSE), "deprecated")
})

test_that("scaleComposite bins the scaling data the same way on CE and yr BP axes", {
  set.seed(8)
  yrs <- 1001:2000
  sd1 <- rnorm(length(yrs))
  bvCE <- seq(1000, 2000, by = 50)
  comp <- matrix(rnorm(20), ncol = 1)
  ce <- scaleComposite(comp, bvCE, yrs, sd1, time.units = "CE")
  bp <- scaleComposite(comp[20:1, ], rev(1950 - bvCE), 1950 - yrs, sd1, time.units = "yr BP")
  expect_equal(as.vector(ce), as.vector(bp)[20:1])
})
