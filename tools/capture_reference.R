# Capture golden reference outputs from the CURRENT implementation of
# compositeR before rewiring dependencies to the ens/lipdViz stack.
# Run from package root.
devtools::load_all(".", quiet = TRUE)

load("data/Temp12k.rda")

# small deterministic subset: first 8 records with plain numeric age vectors
usable <- which(sapply(Temp12k, function(x) {
  is.numeric(x$age) && is.null(dim(x$age)) &&
    length(x$paleoData_values) > 50 && !is.null(x$paleoData_TSid)
}))[1:8]
fTS <- Temp12k[usable]

binvec <- seq(0, 10000, by = 500)
ref <- list()

# 1. binning building blocks
ref$simpleBin <- simpleBinTs(fTS[[1]], binvec)
set.seed(7)
ref$sampleBin <- sampleEnsembleThenBinTs(fTS[[1]], binvec)

# 2. spreading and de-duplication
sp <- spreadPaleoData(age = fTS[[1]]$age,
                      value = fTS[[1]]$paleoData_values,
                      spreadBy = 50, spreadMax = 500)
ref$spread <- sp
ref$rcd <- removeConsecutiveDuplicates(sp$spreadVal)

# 3. standardization methods on a synthetic proxy matrix
set.seed(9)
ages <- seq(50, 9950, by = 100)
pdm <- matrix(rnorm(length(ages) * 6, mean = rep(c(10, 20, 5, 0, -5, 100), each = length(ages)),
                    sd = rep(c(1, 4, 2, 1, 0.5, 10), each = length(ages))),
              ncol = 6)
pdm[3, 2] <- NA
ref$stanInterval <- standardizeOverInterval(ages, pdm, interval = c(2000, 6000))
set.seed(11)
ref$stanRandom <- standardizeOverRandomInterval(ages, pdm, duration = 3000,
                                                searchRange = c(0, 9000))
set.seed(12)
ref$stanIterative <- standardizeMeanIteratively(ages, pdm, duration = 3000,
                                                searchRange = c(0, 9000))
ref$rmse <- recordRMSE(td = ref$stanInterval, palMat = pdm)

# 4. end-to-end composites (small nens for speed)
set.seed(20)
ref$compositeSCC <- compositeEnsembles2(fTS, binvec, nens = 5,
                                        stanFun = standardizeOverInterval,
                                        interval = c(2000, 6000),
                                        binFun = simpleBinTs,
                                        verbose = FALSE)
set.seed(21)
ref$compositeDCC <- compositeEnsembles2(fTS, binvec, nens = 5,
                                        stanFun = standardizeOverRandomInterval,
                                        duration = 3000,
                                        searchRange = c(0, 9000),
                                        binFun = sampleEnsembleThenBinTs,
                                        verbose = FALSE)

dir.create("tests/testthat", recursive = TRUE, showWarnings = FALSE)
saveRDS(ref, "tests/testthat/reference_outputs.rds")
cat("Saved", length(ref), "reference objects\n")
