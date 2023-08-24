##########
##logRank.power: power analysis for log-rank test
##Compute sample size with specifying effect size (hazard ratio), significant level, and power 
##using the Freedman or Schoenfeld approaches. 

logRank.power <- function(S1, S2, hazard.ratio, censoring = TRUE, r = 1, alpha = 0.05, power = 0.8, alternative = c("one.sided", "two.sided"), method = c("Freedman", "Schoenfeld"), verbose = TRUE) 
{
alternative <- match.arg(alternative)
method <- match.arg(method)
stopifnot(alpha < 0.5, power >= 0.5)

if (censoring) {
	stopifnot(!missing(S1), (!missing(hazard.ratio))|(!missing(S2)))
	if (missing(S2)) {
		S2 <- S1^{hazard.ratio}
	} else  hazard.ratio <- log(S2)/log(S1)
} else {
	if (!missing(S1) & !missing(S2)) {
		hazard.ratio <- log(S2)/log(S1)
	} else stopifnot(!missing(hazard.ratio)) 
}

##quantiles of the standard normal
if (alternative == "two.sided"){
	za <- qnorm(1 - alpha/2)
} else za <- qnorm(1 - alpha)	
zb <- qnorm(power)

##total number of events 
if (method == "Freedman"){
	d <- (za + zb)^2*(1 + hazard.ratio*r)^2/(1 - hazard.ratio)^2/r
} else {
	d <- (za + zb)^2*((1 + r)/log(hazard.ratio))^2/r
}

##sample sizes
if (censoring) {
	prop.events <- (1 - S1 + (1 - S2)*r)/(1+r)
} else prop.events <- 1
n <- d/prop.events
n1 <- n/(1 + r) #number of samples in control group
n2 <- n*r/(1+r) #number of samples in treatment group

d <- round(d, 1)
n <- round(n, 1)
n1 <- round(n1, 1)
n2 <- round(n2, 1)

H1 <- ifelse(alternative == "one.sided", "HR < 1", "HR != 1")
if (verbose){
cat("Hazard ratio (HR):", hazard.ratio, "\n")
cat("Significant level: ", alpha, "\n")
cat("Power: ", power, "\n")
cat("Ratio of subjects in the two groups (N2/N1): ", r, "\n")
cat("Proportion of subjects in control group: ", 1/(1+r), "\n")
cat("Proportion of subjects in treatment group: ", r/(1+r), "\n\n")

cat("Estimated sample sizes for log-rank test\n")
cat("H0: HR = 1 versus H1: ", H1, "(", alternative, "alternative hypothesis)\n")
cat("method: ", method, "\nCensoring: ", censoring, "\n")
cat("Proportion of events: ", prop.events, "\n")
cat("Total number of events: ", d, "\n")
cat("Total number of samples: ", n, "\n")
cat("The number of samples in control group: ", n1, "\n")
cat("The number of samples in treatment group: ", n2, "\n\n")
}
 
c(N.events = d, N.samples = n, N.sample.control = n1, N.samples.treatment = n2)
}
