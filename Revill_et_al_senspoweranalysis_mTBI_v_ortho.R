library(pwr)

#
# Minimum detectable effect size sensitivity power analysis
#
# mTBI vs orthopaedic control comparison calculation
#

# Define parameters
N_mTBI <- 450
N_orthoinjured <- 1604

# Calculate total and proportion
total_n <- N_mTBI + N_orthoinjured
proportion_A <- N_mTBI / total_n  # Proportion of subjects in group A

# Set significance, level, required power and H0 / H1 probability
alpha <- 0.05  # Significance level
power <- 0.8   # Desired power
baseline_H0 <- 0.5  # Assume H0 equally as likely as H1

#
# Calculate for odds ratios
#

# Calculate the effect size (w) detectable with this sample size
# using Cohen's method for chi-square tests which can be converted to odds ratio
effect_size <- pwr.chisq.test(w = NULL, 
                              N = total_n,
                              df = 1, 
                              sig.level = alpha,
                              power = power)$w

# Function to calculate w from p1 and p2
w_from_p <- function(p1, p2) {
  return(abs(p2 - p1) / sqrt((p1 + p2) * (1 - p1) * (1 - p2) / 2))
}

# Function to find p2 given w and p1
find_p2 <- function(w, p1) {
  # Define function to optimize
  f <- function(p2) (w_from_p(p1, p2) - w)^2
  
  # Find p2 that gives desired w
  result <- optimize(f, c(0, 1))
  return(result$minimum)
}

# Calculate for fixed baseline prevalence of 0.5
p1 <- baseline_H0

# Find both possible p2 values (lower and higher than p1)
p2_lower <- find_p2(effect_size, p1)

# For the upper bound, we need to find a p2 > p1 that gives us the same effect size
f <- function(p2) (w_from_p(p1, p2) - effect_size)^2
p2_upper <- optimize(f, c(p1, min(1, p1 + 0.5)))$minimum

# Calculate odds ratios in both directions
or_lower <- (p2_lower / (1 - p2_lower)) / (p1 / (1 - p1))
or_upper <- (p2_upper / (1 - p2_upper)) / (p1 / (1 - p1))


#
# Calculate for Cohen's d
#

# Calculate minimum detectable effect size (Cohen's d) for two independent groups
d_result <- pwr.t2n.test(n1 = N_mTBI, 
                       n2 = N_orthoinjured, 
                       sig.level = alpha, 
                       power = power,
                       d = NULL)

cat("For mTBI (N = 450) vs orthopaedic injury control (N = 1604) comparison:\n")
cat("Minimum detectable odds ratio of above", round(or_upper, 3), "with 80% power at alpha=0.05\n")
cat("Minimum detectable Cohen's d of:", round(d_result$d, 3), "with 80% power at alpha=0.05\n")
