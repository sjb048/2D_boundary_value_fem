import matplotlib.pyplot as plt

# Define file names
circle_files = ['node_errors_circle_poisson.txt', 
                'node_errors_circle_Laplace.txt',
                'node_errors_circle_HelmHoltz.txt']

square_files = ['node_errors_square_poisson.txt', 
                'node_errors_square_Laplace.txt',
                'node_errors_square_HelmHoltz.txt']

# Initialize lists for storing parsed data
circle_errors = []
circle_dofs = []
circle_labels = []

square_errors = []
square_dofs = []
square_labels = []

# Function to parse errors and DOFs
def parse_file_data(file_names, errors_list, dofs_list, labels_list):
    for file_name in file_names:
        try:
            with open(file_name, 'r') as file:
                lines = file.readlines()

                # Parse Error and Dof values
                for line in lines:
                    if "Error :" in line:
                        error_value = float(line.split(":")[1].strip())
                        errors_list.append(error_value)
                    elif "Dof :" in line:
                        dof_value = int(line.split(":")[1].strip())
                        dofs_list.append(dof_value)

                # Use file name as a label
                labels_list.append(file_name.replace('.txt', ''))

                print(f"Data parsed successfully from {file_name}!")

        except FileNotFoundError:
            print(f"Error: The file '{file_name}' was not found. Please check the file path.")
        except Exception as e:
            print(f"An error occurred while processing the file '{file_name}': {e}")

# Parse circle and square data
parse_file_data(circle_files, circle_errors, circle_dofs, circle_labels)
parse_file_data(square_files, square_errors, square_dofs, square_labels)

# Plot data
try:
    plt.figure(figsize=(12, 8))

    # Plot circle data
    plt.plot(circle_dofs, circle_errors, 'ro-', label='Circle Data')  # Red line, circle markers

    # Plot square data
    plt.plot(square_dofs, square_errors, 'gs-', label='Square Data')  # Green line, square markers

    # Annotate error values on the plot
    for dof, error in zip(circle_dofs, circle_errors):
        plt.text(dof, error, f"{error:.2f}", color='red', fontsize=10, ha='center', va='bottom')
    for dof, error in zip(square_dofs, square_errors):
        plt.text(dof, error, f"{error:.2f}", color='green', fontsize=10, ha='center', va='bottom')

    # Plot title and labels
    plt.title('Error Comparison: Circle vs Square Models', fontsize=18)
    plt.xlabel('Degrees of Freedom (DoF)', fontsize=14)
    plt.ylabel('Error', fontsize=14)
    plt.legend(fontsize=12)
    plt.grid(True)
    plt.tight_layout()

    # Save the plot
    output_plot_file = 'error_comparison_circle_square.png'
    plt.savefig(output_plot_file)
    print(f"Plot saved as '{output_plot_file}'.")

    # Show the plot
    plt.show()

except Exception as e:
    print(f"An error occurred while plotting: {e}")
