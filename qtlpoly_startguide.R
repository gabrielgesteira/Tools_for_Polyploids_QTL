## QTLpoly startguide
## More information: https://gabrielgesteira.github.io/QTLpoly/

## Installing package
install.package(qtlpoly)
## library(devtools)
## install_github("gabrielgesteira/QTLpoly@remim2")

## Loading packages
## library(doParallel)
library(qtlpoly) 

## Reading genotype probabilities
genoprob4x = readRDS("genoprob.err.rds")
(genoprob4x[[1]]$probs[1:10,1:2,1])

## Reading phenotypes
pheno4x = read.csv("pheno_RRD.csv", h=T, row.names = 1)
str(pheno4x)

## Putting geno + pheno together
data = read_data(ploidy = 4, geno.prob = genoprob4x, pheno = pheno4x, step = 1)

## Viewing dataset
print(data, detailed = TRUE) 

## Score-based resampling to assess genome-wide significance (don't run)
## data.sim = simulate_qtl(data = data, mu = 0, h2.qtl = NULL, var.error = 1, n.sim = 10, missing = TRUE, seed = 123)
## score.null = null_model(data = data.sim$results, n.clusters = 1)
## Reading score-statistics results
score.null = readRDS("score.null.rds")
min.pvl = unlist(lapply(score.null$results, function(x) return(x$pval[which.max(x$stat)])))
quantile(sort(min.pvl), c(0.2, 0.05))

## Adjusting null model 
null.mod = null_model(data = data, pheno.col = 1, n.clusters = 1)

## Viewing null model
print(null.mod)

## Performing first search round
search.mod = search_qtl(data = data, model = null.mod, w.size = 15, sig.fwd = 0.001, n.clusters = 1)

## Visualizing results
print(search.mod)

## Performing model optimization
optimize.mod = optimize_qtl(data = data, model = search.mod, sig.bwd = 0.0002, n.clusters = 1)

## Visualizing results
print(optimize.mod)

## Performing another search round
search.mod2 = search_qtl(data=data, model = optimize.mod, sig.fwd = 0.0002, n.clusters = 1)

## Visualizing results
print(search.mod2)

## Genomewide profiling with final set of QTL
profile.mod = profile_qtl(data = data, model = optimize.mod, d.sint = 1.5, polygenes = FALSE, n.clusters = 1)

## Visualizing results
print(profile.mod)
plot_profile(profile.mod)

## Checking lower SI
print(profile.mod, sint = "lower") 

## Checking upper SI
print(profile.mod, sint = "upper")

## Performing automatic search
remim.mod = remim(data = data, score.null = score.null, sig.fwd = 0.2, sig.bwd = 0.05, w.size = 15, d.sint = 1.5, n.clusters = 1)

## Visualizing results
print(remim.mod)

## Checking lower SI
print(remim.mod, sint = "lower")

## Checking upper SI
print(remim.mod, sint = "upper")

## Plotting profile
plot_profile(data = data, model = remim.mod, ylim = c(0,10))
plot_profile(data = data, model = profile.mod, ylim = c(0,10))

## Plotting profile
plot_profile(data = data, model = remim.mod, grid = TRUE) 
plot_profile(data = data, model = remim.mod, grid = FALSE) 

## Plotting QTL SI
plot_sint(data = data, model = remim.mod)

## Calculate metrics
fitted.mod = fit_model(data=data, model=remim.mod, probs="joint", polygenes="none")

## Viewing results
summary(fitted.mod)

## Plotting QTL info
plot_qtl(data = data, model = remim.mod, fitted = fitted.mod, drop.pheno = FALSE)

## Account for allele effects
est.effects = qtl_effects(ploidy = 4, fitted = fitted.mod)

## Plot allele effects
plot(est.effects) 

## Calculate breeding values
y.hat = breeding_values(data = data, fitted = fitted.mod)

## Plot BV
plot(y.hat) 

## FEIM
## Perform permutations
perm = permutations(data = data, pheno.col = 1, n.sim = 1000, n.clusters = 1)

## Check permutations
print(perm)

## Define significance threshold
(sig.lod = perm$sig.lod$`0.95`)

## Adjust FEIM model
feim.mod = feim(data = data, pheno.col = 1, w.size = 15, sig.lod = sig.lod)

## Viewing FEIM results
print(feim.mod)

## Plotting FEIM results
plot_profile(data = data, model = feim.mod, grid = TRUE)

## Exporting results to VIEWpoly
setwd("./viewpoly")
save(maps4x, file="mappoly.maps.RData")
save(data, file="qtlpoly.data.RData")
save(remim.mod, file="qtlpoly.remim.mod.RData")
save(fitted.mod, file="qtlpoly.fitted.mod.RData")
save(est.effects, file="qtlpoly.est.effects.RData")

## Loading VIEWpoly
install.packages("viewpoly")
library(viewpoly)
run_app()
