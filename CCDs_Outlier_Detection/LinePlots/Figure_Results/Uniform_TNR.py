import pandas as pd
import matplotlib.pyplot as plt

# Load the data from the Excel file
file_path = 'G:/code_working_folder/Simulation_Results/Figure_Results/Uniform_TPR_TNR.xlsx'
data = pd.read_excel(file_path)

# Extract and clean the data for TNR
plot_data = {
    'dimension': [],
    'algorithm': [],
    'data_size': [],
    'TNR': []
}

# Populating the dictionary correctly for TNR
for col in range(3, data.shape[1], 2):  # Start from the third column and step by 2 to skip TPR columns
    size = data.iloc[0, col - 1]  # Correctly retrieve the data size from the previous column
    for row in range(2, data.shape[0]):
        if pd.notna(data.iloc[row, 0]):
            current_dimension = data.iloc[row, 0].split('=')[1].strip()  # Update current dimension when a new one is found
        plot_data['dimension'].append(current_dimension)
        plot_data['algorithm'].append(data.iloc[row, 1].strip())
        plot_data['data_size'].append(size)
        plot_data['TNR'].append(float(data.iloc[row, col]))

# Convert to DataFrame
plot_df = pd.DataFrame(plot_data)

# Convert the 'dimension' column to categorical for plotting
dimension_order = sorted(plot_df['dimension'].unique(), key=lambda x: int(x))
plot_df['dimension'] = pd.Categorical(plot_df['dimension'], categories=dimension_order, ordered=True)

# Pivot the DataFrame for plotting again
pivot_df = plot_df.pivot_table(index=['dimension', 'algorithm'], columns='data_size', values='TNR', aggfunc='first').reset_index()

# Unique data sizes for plotting
unique_data_sizes = sorted(set(plot_df['data_size']))

# Ensure the order of algorithms is as specified: U-MCCDs, SU-MCCDs, UN-MCCDs, SUN-MCCDs
algorithm_order = ['U-MCCDs', 'SU-MCCDs', 'UN-MCCDs', 'SUN-MCCDs']

# Define line styles and colors
line_styles = ['-', '--', '-.', ':']
colors = ['black', 'dimgray', 'darkgray']

# Separate line plots for each algorithm, with each line representing different data sizes
fig, axes = plt.subplots(2, 2, figsize=(15, 7), sharex=True, sharey=True)

# Flatten the axes array for easy indexing
axes = axes.flatten()

# Plot each algorithm in the specified order
for ax, algorithm in zip(axes, algorithm_order):
    group = pivot_df[pivot_df['algorithm'] == algorithm]
    for i, size in enumerate(unique_data_sizes):  # Use unique data sizes correctly
        ax.plot(group['dimension'], group[size], linestyle=line_styles[i % len(line_styles)],
                color=colors[i % len(colors)], label=f"{size}")
    ax.set_title(algorithm, fontsize=14)
    ax.set_xlabel('Dimension $d$')
    ax.set_ylabel('TNR', fontsize=16)
    ax.legend(title='Data Size $n$')

plt.suptitle('TNRs with Uniform Clusters', fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 1])
plt.show()
