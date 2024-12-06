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




"""
Binding vs SASA
"""
def plot_corr_binding_SASA_(metric_df, lectin, binding_motif):
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

def plot_corr_binding_SASA(metric_df, lectin, binding_motif):
    """Plot correlation of binding scores with SASA-related metrics."""

    for metric_key in metric_df.columns:
        if metric_key.startswith("SASA"):
            # Calculate the correlation coefficient
            corr_value = metric_df['binding_score'].corr(metric_df[metric_key])

            # Print the correlation coefficient
            print(f"Correlation between Binding and {metric_key}: {corr_value}")

            # Create the plot
            plt.figure(figsize=(8, 6))
            sns.regplot(
                x=metric_key,
                y='binding_score',
                data=metric_df,
                scatter_kws={'alpha': 0.7},
                line_kws={'color': 'red'}
            )

            # Set plot title and axis labels
            plt.title(f'Binding vs {metric_key} ({lectin} - {binding_motif})', fontsize=12)
            plt.xlabel(metric_key, fontsize=10)
            plt.ylabel('Binding Score', fontsize=10)

            # Save the plot
            plt.savefig(f'scripts/correlation/plots/Binding_vs_{metric_key}_{lectin}.png', dpi=300)

            # Show the plot
            plt.show()

def plot_corr_binding_SASA_subplots(metric_df, lectin, binding_motif):
    """Plot correlation of binding scores with SASA-related metrics in subplots."""
    sasa_metrics = [metric for metric in metric_df.columns if metric.startswith("SASA_weighted")]

    # Create a figure with 3 rows and 3 columns of subplots
    fig, axes = plt.subplots(3, 3, figsize=(15, 15))
    axes = axes.flatten()  # Flatten the 2D array of axes for easier iteration

    # Iterate through SASA metrics and corresponding axes
    for i, metric_key in enumerate(sasa_metrics):
        if i < len(axes):  # Ensure we don't exceed the number of available subplots
            ax = axes[i]

            # Calculate correlation coefficient
            corr_value = metric_df['binding_score'].corr(metric_df[metric_key])
            print(f"Correlation between Binding and {metric_key}: {corr_value}")

            # Create the scatter plot with regression line
            sns.regplot(
                x=metric_key,
                y='binding_score',
                data=metric_df,
                ax=ax,
                scatter_kws={'alpha': 0.7},
                line_kws={'color': 'red'},
            )

            # Set subplot title and labels
            ax.set_title(f'{metric_key} (Corr: {corr_value:.2f})', fontsize=10)
            ax.set_xlabel(metric_key, fontsize=8)
            ax.set_ylabel('Binding Score', fontsize=8)
        else:
            break

    # Remove empty subplots if there are fewer metrics than subplots
    for j in range(len(sasa_metrics), len(axes)):
        fig.delaxes(axes[j])

    # Adjust layout for better spacing
    plt.tight_layout()

    # Save the plot
    plt.savefig(f'scripts/correlation/plots/set3/SASA/Binding_vs_SASA_Metrics_{lectin}.png', dpi=300)

    # Show the plot
    plt.show()

"""
Both plots combined
"""
def plot_combined(metric_df, lectin, binding_motif):
    # Calculate correlation coefficients
    corr_binding_flexibility = metric_df['binding_score'].corr(metric_df['weighted_mean_flexibility'])
    corr_binding_sasa_weighted = metric_df['binding_score'].corr(metric_df['SASA_weighted'])

    # Print correlation coefficients
    print("lectin: ", lectin)
    print("binding_motif: ", binding_motif)
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
        line_kws={'color': 'red'},
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
    plt.savefig(f'scripts/correlation/plots/set3/Binding_vs_Flexibility_and_SASA_{lectin}.png', dpi=300)

    # Show the plots
    plt.show()

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
    plt.savefig(f'scripts/correlation/plots/set3/Binding_vs_Flexibility_and_SASA_{lectin}.png', dpi=300)

    # Show the plots
    plt.show()




