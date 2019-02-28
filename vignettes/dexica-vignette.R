## ---- eval = FALSE-------------------------------------------------------
#  install.packages("devtools") # Skip this line if you aleady have devtools installed
#  library(devtools)
#  install_github("MPCary/DEXICA", build_vignettes = TRUE)
#  install_github("MPCary/DEXDATA.Celegans") # C. elegans data package
#  # Note: It is advisable to restart your R session after running install_github
#  # due to a known issue in which documentation files are incorrectly reported
#  # to be corrupt.

## ---- eval = TRUE--------------------------------------------------------
library(DEXICA)

# Load and preprocess microarray compendium
GPL200.50.pp = preprocessMatrix(GPL200.50.mini)

# Extract 10 modules from compendium
set.seed(1) # Set random seed for reproducibility
result = predictModules(GPL200.50.pp, n.comp = 10)
dim(result$S) # Dimensions of S matrix
result$S[1:3,1:3] # First few rows and columns of S matrix

## ---- fig.show='hold', eval = TRUE---------------------------------------
# Partition S matrix
S.fp = partition(result$S, method = "fixed", t = 3) # Fixed partitioning
S.vp = partition(result$S, method = "ann") # Variable partitioning

# Separate positive sets from negative sets
S.fp.split = splitMatrix(S.fp)
S.vp.split = splitMatrix(S.vp)

# Count the number of probesets in each set
fp.counts = apply(S.fp.split, 2, sum)
vp.counts = apply(S.vp.split, 2, sum)

# Plot counts
count.matrix = rbind(fp.counts, vp.counts)
count.matrix = count.matrix[, order(apply(count.matrix, 2, sum), decreasing = TRUE)]
barplot(count.matrix, beside = TRUE, xlab = "Hemi-modules", ylab = "Probesets",
        names.arg = rep("", ncol(count.matrix)), col = 1:2)
legend("topright", legend = c("Fixed", "ANN"), fill = 1:2, cex = 0.75)

## ---- fig.show='hold', eval = TRUE---------------------------------------
# Evaluate module predictions using GO annotations
fp.check.GO = checkModules(S.fp.split, GPL200.GO.mini)
vp.check.GO = checkModules(S.vp.split, GPL200.GO.mini)

data.frame(rbind(fp.check.GO, vp.check.GO))

## ---- eval = TRUE--------------------------------------------------------
# Create a simple DexBatch object
db = dexBatch(compen = list(compen1 = GPL200.50.mini),
              annmats = list(annmat1 = GPL200.GO.mini),
              params = parameterSet(n.comp = c(10,20)))
countJobs(db)

## ---- eval = FALSE-------------------------------------------------------
#  # Run DexBatch jobs
#  for(i in 1:countJobs(db)) {
#    runJob(db, i)
#  }

## ---- eval = TRUE--------------------------------------------------------
# Test a range of extracted components
n.comp.range = seq(from = 5, to = 100, by = 5)
p = parameterSet(n.comp = n.comp.range)
getParamSetCount(p)

## ---- eval = TRUE--------------------------------------------------------
# Test a range of extracted components and two preprocessing options
n.comp.range = seq(from = 5, to = 100, by = 5)
p = parameterSet(n.comp = n.comp.range, center.cols = c(T,F), scale.cols = c(T,F))
getParamSetCount(p)

## ---- eval = TRUE--------------------------------------------------------
# Test two preprocessing options 5 times each
my.seeds = c(1:5)
p = parameterSet(n.comp = 50, center.cols = c(T,F), scale.cols = c(T,F), w.init = my.seeds)
getParamSetCount(p)

## ---- eval = FALSE-------------------------------------------------------
#  # Create and save a DexBatch
#  my.seeds = c(1:3)
#  n.comp.range = seq(from = 5, to = 100, by = 5)
#  p = parameterSet(n.comp = n.comp.range, center.cols = c(T,F), scale.cols = c(T,F), w.init = my.seeds)
#  db = dexBatch(compen = list(compen1 = GPL200.50.mini),
#                annmats = list(annmat1 = GPL200.GO.mini),
#                params = p)
#  countJobs(db) # 240 jobs
#  
#  save(db, file = "db.Rdata")

## ---- eval = FALSE-------------------------------------------------------
#  # Call this script in parallel using cluster management software
#  load(file = "db.Rdata") # Load DexBatch object
#  j = getBatchJobID() # Get ID of current job
#  runJob(db, j)

## ---- eval = FALSE-------------------------------------------------------
#  maintainer("DEXICA")

## ---- eval = TRUE--------------------------------------------------------
sessionInfo()

