import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scripts.correlation.merged.AAL.AAL_Fuc import final_df

output_file = 'scripts/correlation/merged/AAL/Binding_vs_Flexibility_AAL.png'

# Set Seaborn style and consistent grid styling
sns.set(style="whitegrid")
plt.rcParams['grid.color'] = 'gray'
plt.rcParams['grid.linestyle'] = '--'

# Calculate correlation coefficients
corr_binding_flexibility = final_df['binding_score'].corr(final_df['weighted_mean_flexibility'])

# Print correlation coefficients
print(f"Correlation between Binding and Flexibility: {corr_binding_flexibility}")

# Create the scatter plot for Binding vs Flexibility
plt.figure(figsize=(8, 6))
sns.scatterplot(x='binding_score', y='weighted_mean_flexibility', data=final_df, alpha=0.7)
plt.title('Binding vs Flexibility (AAL)', fontsize=14)
plt.xlabel('Binding Score', fontsize=12)
plt.ylabel('Flexibility', fontsize=12)
plt.grid(True, color='gray', linestyle='--')

# Save the plot
plt.savefig(output_file, dpi=300)

# Show the plot
plt.show()
