####### This script shows a regression plot of energies of blocks comparing
# the blocks present in the heads and in the tails


import math
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt

df = pd.read_csv('/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/block_counts_energy.csv', sep=',')
print(df.keys())
df['Change_Tails_to_heads'] = - df['Change_Tails_to_heads']
df['abs_Tails_to_Heads'] = df['Proportion_in_Heads'] - df['Proportion_in_Tails']


# Filter the DataFrame to only include blocks that occur at least once in 1000 seqs
filtered_df = df[(df['Proportion_in_Tails'] > 0.0001) | (df['Proportion_in_Tails'] > 0.0001)]


# Extract X and Y values
x =  filtered_df['Energy']
y = filtered_df['Change_Tails_to_heads']

# Add a constant to the independent variable (for the intercept)
# x = sm.add_constant(x)

# Compute weights as the sum of Proportion_in_Tails and Proportion_in_Heads
weights = filtered_df['Proportion_in_Tails'] + filtered_df['Proportion_in_Heads']
log_weights = [math.log(w) for w in weights]
alpha_by_log_weights = [(w - min(log_weights))/(max(log_weights)-min(log_weights)) for w in log_weights]
# Perform weighted linear regression
model = sm.WLS(y, x, weights=weights).fit()

# Print the regression results
# print(model.summary())

# Plot the data and the regression line
plt.figure(figsize=(10, 6))
plt.scatter(filtered_df['Energy'], y, label='Data',alpha = alpha_by_log_weights)
plt.plot(filtered_df['Energy'], model.predict(x), color='red', label='Weighted Fitted line')
plt.axhline(0, color='black', linewidth=1)
plt.xlabel('Energy')
plt.ylabel('Change Tails to Heads')
plt.title(f'Weighted Linear Regression of the Change in Number from Tails to Heads vs the Energy of the Block\n')
plt.legend()
plt.show()



# print(filtered_df[['Block','Change_Tails_to_heads']].sort_values(by = ['Change_Tails_to_heads']))
