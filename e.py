import pandas as pd
import matplotlib.pyplot as plt

# Read the data from the text file
try:
    # Replace '\t' with the appropriate delimiter if it's not tab-separated.
    data = pd.read_csv('node_errors_circle_HelmHoltz.txt', delimiter='\t')  # Assuming tab-separated values
except FileNotFoundError:
    print("Error: The file 'node_errors.txt' was not found.")
    exit()

# Ensure the required columns are present in the data
if 'NodeID' not in data.columns or 'Error' not in data.columns:
    print("Error: Required columns 'NodeID' and 'Error' are not present in the file.")
    exit()

# Display the first few rows of the data (optional)
print(data.head())

# Plot NodeID vs. Error (using NodeID as a proxy for DoF)
plt.figure(figsize=(8, 6))
plt.plot(data['NodeID'], data['Error'], marker='o', linestyle='-', color='b', label='Error')

# Adding titles and labels
plt.title('NodeID (DoF Proxy) vs. Error')
plt.xlabel('NodeID (Proxy for DoF)')
plt.ylabel('Error')

# Optionally, add grid and logarithmic scale if appropriate
plt.grid(True)
# Uncomment the following lines if logarithmic scale is needed
# plt.xscale('log')
# plt.yscale('log')

# Add legend
plt.legend()

# Save the plot as an image file
plt.savefig('NodeID_vs_Error.png')

# Display the plot
plt.show()