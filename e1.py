import pandas as pd
import matplotlib.pyplot as plt

# Define file names
file_names = ['node_errors_circle_poisson.txt',  # Replace with your actual file names
              'node_errors_circle_Laplace.txt',
              'node_errors_circle_HelmHoltz.txt']


# # Define file names
file_names = ['node_errors_square_poisson.txt',  # Replace with your actual file names
              'node_errors_square_Laplace.txt',
              'node_errors_square_HelmHoltz.txt']
errors = []
dofs = []
labels = []
for file_name in file_names:
    try:
        # Open and read the file
        with open(file_name, 'r') as file:
            lines = file.readlines()
            
            # Parse Error and Dof values
            for line in lines:
                if "Error :" in line:
                    error_value = float(line.split(":")[1].strip())
                    errors.append(error_value)
                elif "Dof :" in line:
                    dof_value = int(line.split(":")[1].strip())
                    dofs.append(dof_value)
            
            # Use file name as a label
            labels.append(file_name.replace('.txt', ''))

            print(f"Data parsed successfully from {file_name}!")

    except FileNotFoundError:
        print(f"Error: The file '{file_name}' was not found. Please check the file path.")
    except Exception as e:
        print(f"An error occurred while processing the file '{file_name}': {e}")

# Plot Error values for each file
try:
    plt.figure(figsize=(10, 7))
    plt.bar(labels, errors, color='skyblue', alpha=0.8)

    for i, error in enumerate(errors):
        plt.text(i, error, f"{error:.2f}", ha='center', va='bottom', fontsize=10)

    plt.title('Error Comparison Across Models', fontsize=16)
    plt.xlabel('Models', fontsize=14)
    plt.ylabel('Error', fontsize=14)
    plt.tight_layout()

    # Save the plot
    output_plot_file = 'error_comparison_square.png'
    plt.savefig(output_plot_file)
    print(f"Plot saved as '{output_plot_file}'.")

    # Show the plot
    plt.show()

except Exception as e:
    print(f"An error occurred while plotting: {e}")