import matplotlib.pyplot as plt
import seaborn as sns

out1 = 'scripts/correlation/merged/SNA/Binding_vs_Flexibility_SNA.png'
out2 = 'scripts/correlation/merged/SNA/Binding_vs_SASA_SNA.png'

# Set Seaborn style and grid properties
sns.set(style="whitegrid")
plt.rcParams['grid.color'] = 'gray'
plt.rcParams['grid.linestyle'] = '--'
plt.figure(figsize=(8, 6))

"""
Binding vs Flexibility
"""

def plot_corr_binding_flexibility(metric_df, lectin, binding_motif):

    # Calculate correlation coefficients
    corr_binding_flexibility = metric_df['binding_score'].corr(metric_df['weighted_mean_flexibility'])

    # Print correlation coefficients
    print(f"Correlation between Binding and Flexibility: {corr_binding_flexibility}")

    # Create the scatter plot
    sns.scatterplot(x='binding_score', y='weighted_mean_flexibility', data=metric_df, alpha=0.7)
    plt.set_title(f'Binding vs Flexibility Weighted {lectin} {binding_motif}', fontsize=12)
    plt.xlabel('Binding Score')
    plt.ylabel('Flexibility')
    plt.tight_layout()
    plt.savefig(out1, dpi=300)
    plt.show()

"""
Binding vs SASA
"""

def plot_corr_binding_SASA(metric_df, lectin, binding_motif):

    # Calculate correlation coefficients
    corr_binding_sasa_mean = metric_df['binding_score'].corr(metric_df['SASA_mean'])
    corr_binding_sasa_median = metric_df['binding_score'].corr(metric_df['SASA_median'])
    corr_binding_sasa_weighted = metric_df['binding_score'].corr(metric_df['SASA_weighted'])

    # Print correlation coefficients
    print(f"Correlation between Binding and SASA Mean: {corr_binding_sasa_mean}")
    print(f"Correlation between Binding and SASA Median: {corr_binding_sasa_median}")
    print(f"Correlation between Binding and SASA Weighted: {corr_binding_sasa_weighted}")

    # Create subplots for each SASA metric
    fig, axes = plt.subplots(3, 1, figsize=(6, 12), constrained_layout=True)

    # Plot Binding vs SASA Mean
    sns.scatterplot(ax=axes[0], x='binding_score', y='SASA_mean', data=metric_df, alpha=0.7)
    axes[0].set_title('Binding vs SASA Mean (SNA)', fontsize=12)
    axes[0].set_xlabel('Binding Score', fontsize=10)
    axes[0].set_ylabel('SASA Mean', fontsize=10)

    # Plot Binding vs SASA Median
    sns.scatterplot(ax=axes[1], x='binding_score', y='SASA_median', data=metric_df, alpha=0.7, color='orange')
    axes[1].set_title('Binding vs SASA Median (SNA)', fontsize=12)
    axes[1].set_xlabel('Binding Score', fontsize=10)
    axes[1].set_ylabel('SASA Median', fontsize=10)

    # Plot Binding vs SASA Weighted
    sns.scatterplot(ax=axes[2], x='binding_score', y='SASA_weighted', data=metric_df, alpha=0.7, color='green')
    axes[2].set_title(f'Binding vs SASA Weighted {lectin} {binding_motif}', fontsize=12)

    axes[2].set_xlabel('Binding Score', fontsize=10)
    axes[2].set_ylabel('SASA Weighted', fontsize=10)

    # Save the plot
    plt.savefig(out2, dpi=300)

    # Show the plots
    plt.show()
