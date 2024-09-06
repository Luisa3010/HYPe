#### Helper script to find the likelihood of a binomial experiment

from scipy.stats import binom

# Given values
n = 998  # number of trials
p = 0.091  # probability of success in each trial

# Calculate the expected value (mean) and variance
mean = n * p
variance = n * p * (1 - p)

# Observed value (number of successes)
observed_value = 167

# Calculate the probability of observing exactly the observed_value successes
prob_exact = binom.pmf(observed_value, n, p)

# Calculate the cumulative probability of observing up to the observed_value successes
cumulative_prob = binom.cdf(observed_value, n, p)

# Probability of observing a value greater than the observed_value successes
prob_greater_than_observed = 1 - cumulative_prob

# Probability of observing a value less than or equal to the observed_value successes
prob_less_than_or_equal_observed = cumulative_prob

print(f"Probability of observing exactly {observed_value} errors: {prob_exact}")
print(f"Probability of observing <= {observed_value} errors: {prob_less_than_or_equal_observed + prob_exact}")
print(f"Probability of observing >= {observed_value} errors: {prob_greater_than_observed + prob_exact}")
