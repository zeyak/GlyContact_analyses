import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scripts.correlation.merged.SNA.SNA_Neu5AcGc import final_df

output_file = 'scripts/correlation/merged/SNA/Binding_vs_Flexibility_SNA.png'

# Set Seaborn style and grid properties
sns.set(style="whitegrid")
plt.rcParams['grid.color'] = 'gray'
plt.rcParams['grid.linestyle'] = '--'

# Calculate correlation coefficients
corr_binding_flexibility = final_df['binding_score'].corr(final_df['weighted_mean_flexibility'])

# Print correlation coefficients
print(f"Correlation between Binding and Flexibility: {corr_binding_flexibility}")

# Create the scatter plot
plt.figure(figsize=(8, 6))
sns.scatterplot(x='binding_score', y='weighted_mean_flexibility', data=final_df, alpha=0.7)
plt.title('Binding vs Flexibility (SNA)')
plt.xlabel('Binding Score')
plt.ylabel('Flexibility')
plt.grid(True, color='gray', linestyle='--')
plt.tight_layout()
plt.savefig(output_file, dpi=300)
plt.show()

