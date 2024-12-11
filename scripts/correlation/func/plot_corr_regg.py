import matplotlib.pyplot as plt
import seaborn as sns


sns.set(style="whitegrid")
plt.rcParams['grid.color'] = 'gray'
plt.rcParams['grid.linestyle'] = '--'
plt.figure(figsize=(8, 6))

def plot_combined_colors(metric_df, lectin, binding_motif):
    """Plots Binding vs Flexibility and Binding vs SASA Weighted with colors for glycans."""

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))  # Create two subplots side by side

    # Plot Binding vs Flexibility without legend
    sns.scatterplot(
        ax=axes[0],
        x='weighted_mean_flexibility',
        y='binding_score',
        data=metric_df,
        hue="class",  # Color by glycan class
        palette="tab10",  # Use glycan column for coloring
        alpha=0.7
    )
    sns.regplot(
        ax=axes[0],
        x='weighted_mean_flexibility',
        y='binding_score',
        data=metric_df,
        scatter=False,  # Do not plot points again
        line_kws={'color': 'red'}
    )
    axes[0].set_title(f'Binding vs Flexibility Weighted\n{lectin} {binding_motif}', fontsize=12)
    axes[0].set_xlabel('Flexibility')
    axes[0].set_ylabel('Binding Score')
    axes[0].get_legend().remove()  # Remove legend from the first plot

    # Plot Binding vs SASA Weighted with the legend
    sns.scatterplot(
        ax=axes[1],
        x='SASA_weighted',
        y='binding_score',
        data=metric_df,
        hue='class',  # Use glycan column for coloring
        palette="tab10",
        alpha=0.7
    )
    sns.regplot(
        ax=axes[1],
        x='SASA_weighted',
        y='binding_score',
        data=metric_df,
        scatter=False,  # Do not plot points again
        line_kws={'color': 'red'}
    )
    axes[1].set_title(f'Binding vs SASA Weighted\n{lectin} {binding_motif}', fontsize=12)
    axes[1].set_xlabel('SASA Weighted')
    axes[1].set_ylabel('Binding Score')
    axes[1].legend(title="Glycan class", bbox_to_anchor=(1.05, 1), loc='upper left')  # Keep class legend

    # Adjust layout for better spacing
    plt.tight_layout()

    # Save the combined plot
    plt.savefig(f'scripts/correlation/plots/Binding_vs_Flexibility_and_SASA_{lectin}.png', dpi=300)

    # Show the plots
    plt.show()




