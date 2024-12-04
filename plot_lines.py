import matplotlib.pyplot as plt

# Data for circle and square problems
errors_circle = [0.00184879, 0.0158273, 0.59051]
errors_square = [0.264921, 0.375146, 0.00246654]

labels_circle = ["Poisson Circle", "Laplace Circle", "Helmholtz Circle"]
labels_square = ["Poisson Square", "Laplace Square", "Helmholtz Square"]

# Function to plot combined errors with connecting lines
def plot_combined_errors_with_lines(errors_circle, errors_square, labels_circle, labels_square, output_plot_file, title):
    try:
        plt.figure(figsize=(12, 8))

        # Plot Circle Errors
        plt.plot(labels_circle, errors_circle, marker='o', linestyle='-', color='blue', label='Circle Models')
        
        # Plot Square Errors
        plt.plot(labels_square, errors_square, marker='o', linestyle='-', color='green', label='Square Models')

        # Annotate each point with its error value
        for i, error in enumerate(errors_circle):
            plt.text(i, error, f"{error:.6f}", ha='center', va='bottom', fontsize=10)
        for i, error in enumerate(errors_square):
            plt.text(i + len(errors_circle), error, f"{error:.6f}", ha='center', va='bottom', fontsize=10)

        # Add legend
        plt.legend(loc="upper left", fontsize=12)

        plt.title(title, fontsize=16)
        plt.xlabel('Models', fontsize=14)
        plt.ylabel('Error', fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.tight_layout()

        # Save the plot
        plt.savefig(output_plot_file)
        print(f"Plot saved as '{output_plot_file}'.")

        # Show the plot
        plt.show()

    except Exception as e:
        print(f"An error occurred while plotting: {e}")

# Plot combined errors with connecting lines
plot_combined_errors_with_lines(
    errors_circle,
    errors_square,
    labels_circle,
    labels_square,
    'error_comparison_combined_line_with_lines.png',
    'Error Comparison Across Circle and Square Models'
)
