import pandas as pd
import matplotlib.pyplot as plt

# Read the data from the text file
data = pd.read_csv('node_errors.csv', sep='\t')

# Display the data (optional)
print(data)

# Plot DoF vs. Error
plt.figure(figsize=(8, 6))
plt.plot(data['DoF'], data['Error'], marker='o', linestyle='-', color='b')

# Adding titles and labels
plt.title('Degrees of Freedom vs. Error')
plt.xlabel('Degrees of Freedom (DoF)')
plt.ylabel('Error')

# Optionally, add grid and logarithmic scale if appropriate
plt.grid(True)
# plt.xscale('log')
# plt.yscale('log')

# Save the plot as an image file
plt.savefig('DoF_vs_Error.png')

# Display the plot
plt.show()