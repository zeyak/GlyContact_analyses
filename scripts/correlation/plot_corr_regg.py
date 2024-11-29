import matplotlib.pyplot as plt
import seaborn as sns



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

    # Create the scatter plot with regression line
    sns.regplot(
        x='weighted_mean_flexibility',
        y='binding_score',
        data=metric_df,
        scatter_kws={'alpha': 0.7},
        line_kws={'color': 'red'}
    )

    # Set plot labels and title
    plt.title(f'Binding vs Flexibility Weighted {lectin} {binding_motif}', fontsize=12)
    plt.xlabel('Flexibility')
    plt.ylabel('Binding Score')
    plt.tight_layout()

    # Save and show the plot
    plt.savefig(f'scripts/correlation/plots/Binding_vs_Flexibility_{lectin}.png', dpi=300)
    plt.show()

"""
Binding vs SASA
"""


def plot_corr_binding_SASA(metric_df, lectin, binding_motif):
    # Calculate the correlation coefficient for SASA weighted score
    corr_binding_sasa_weighted = metric_df['binding_score'].corr(metric_df['SASA_weighted'])

    # Print the correlation coefficient
    print(f"Correlation between Binding and SASA Weighted: {corr_binding_sasa_weighted}")

    # Create a single plot for Binding vs SASA Weighted
    plt.figure(figsize=(8, 6))
    sns.regplot(
        x='SASA_weighted',
        y='binding_score',
        data=metric_df,
        scatter_kws={'alpha': 0.7},
        line_kws={'color': 'red'}
    )

    # Set plot title and axis labels
    plt.title(f'Binding vs SASA Weighted {lectin} {binding_motif}', fontsize=12)
    plt.xlabel('SASA Weighted', fontsize=10)
    plt.ylabel('Binding Score', fontsize=10)

    # Save the plot
    plt.savefig(f'scripts/correlation/plots/Binding_vs_Flexibility_{lectin}.png', dpi=300)

    # Show the plot
    plt.show()



def plot_combined(metric_df, lectin, binding_motif):
    # Calculate correlation coefficients
    corr_binding_flexibility = metric_df['binding_score'].corr(metric_df['weighted_mean_flexibility'])
    corr_binding_sasa_weighted = metric_df['binding_score'].corr(metric_df['SASA_weighted'])

    # Print correlation coefficients
    print(f"Correlation between Binding and Flexibility: {corr_binding_flexibility}")
    print(f"Correlation between Binding and SASA Weighted: {corr_binding_sasa_weighted}")

    # Create subplots for both plots side-by-side
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True, constrained_layout=True)

    # Plot Binding vs Flexibility
    sns.regplot(
        ax=axes[0],
        x='weighted_mean_flexibility',
        y='binding_score',
        data=metric_df,
        scatter_kws={'alpha': 0.7},
        line_kws={'color': 'red'}
    )
    axes[0].set_title(f'Binding vs Flexibility Weighted\n{lectin} {binding_motif}', fontsize=12)
    axes[0].set_xlabel('Flexibility')
    axes[0].set_ylabel('Binding Score')

    # Plot Binding vs SASA Weighted
    sns.regplot(
        ax=axes[1],
        x='SASA_weighted',
        y='binding_score',
        data=metric_df,
        scatter_kws={'alpha': 0.7},
        line_kws={'color': 'red'}
    )
    axes[1].set_title(f'Binding vs SASA Weighted\n{lectin} {binding_motif}', fontsize=12)
    axes[1].set_xlabel('SASA Weighted')

    # Save the combined plot
    plt.savefig(f'scripts/correlation/plots/Binding_vs_Flexibility_and_SASA_{lectin}.png', dpi=300)

    # Show the plots
    plt.show()


