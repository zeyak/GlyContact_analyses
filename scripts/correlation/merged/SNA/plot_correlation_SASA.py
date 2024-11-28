import matplotlib.pyplot as plt
import seaborn as sns
from scripts.correlation.merged.SNA.SNA_Neu5AcGc import final_df

output_file = 'scripts/correlation/merged/SNA/Binding_vs_SASA_SNA.png'

# Set Seaborn style
sns.set(style="whitegrid")
plt.rcParams['grid.color'] = 'gray'
plt.rcParams['grid.linestyle'] = '--'

# Calculate correlation coefficients
corr_binding_sasa_mean = final_df['binding_score'].corr(final_df['SASA_mean'])
corr_binding_sasa_median = final_df['binding_score'].corr(final_df['SASA_median'])
corr_binding_sasa_weighted = final_df['binding_score'].corr(final_df['SASA_weighted'])

# Print correlation coefficients
print(f"Correlation between Binding and SASA Mean: {corr_binding_sasa_mean}")
print(f"Correlation between Binding and SASA Median: {corr_binding_sasa_median}")
print(f"Correlation between Binding and SASA Weighted: {corr_binding_sasa_weighted}")

# Create subplots for each SASA metric
fig, axes = plt.subplots(3, 1, figsize=(6, 12), constrained_layout=True)

# Plot Binding vs SASA Mean
sns.scatterplot(ax=axes[0], x='binding_score', y='SASA_mean', data=final_df, alpha=0.7)
axes[0].set_title('Binding vs SASA Mean (SNA)', fontsize=12)
axes[0].set_xlabel('Binding Score', fontsize=10)
axes[0].set_ylabel('SASA Mean', fontsize=10)

# Plot Binding vs SASA Median
sns.scatterplot(ax=axes[1], x='binding_score', y='SASA_median', data=final_df, alpha=0.7, color='orange')
axes[1].set_title('Binding vs SASA Median (SNA)', fontsize=12)
axes[1].set_xlabel('Binding Score', fontsize=10)
axes[1].set_ylabel('SASA Median', fontsize=10)

# Plot Binding vs SASA Weighted
sns.scatterplot(ax=axes[2], x='binding_score', y='SASA_weighted', data=final_df, alpha=0.7, color='green')
axes[2].set_title('Binding vs SASA Weighted (SNA)', fontsize=12)
axes[2].set_xlabel('Binding Score', fontsize=10)
axes[2].set_ylabel('SASA Weighted', fontsize=10)

# Save the plot
plt.savefig(output_file, dpi=300)

# Show the plots
plt.show()
