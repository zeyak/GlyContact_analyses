import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import os

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
        hue="class",
        hue_order=['N', 'O', "free","lipid", 'NaN'],
        palette="tab10",
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
        hue='class',
        hue_order=['N', 'O', "free","lipid", 'NaN'],
        palette="tab10",
        # Use glycan column for coloring
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
    plt.savefig(f'scripts/correlation/plots/correlation/Binding_vs_Flexibility_and_SASA_{lectin}.png', dpi=300)

    # Show the plots
    plt.show()

def plot_separate_class(metric_df, lectin, binding_motif):
    """Plots Binding vs Flexibility and Binding vs SASA Weighted for N- and O-linked glycans in separate rows."""

    # Filter for O-glycans and N-glycans
    o_glycans = metric_df[metric_df['class'] == 'O']
    n_glycans = metric_df[metric_df['class'] == 'N']

    # Create a single figure with two rows and two columns
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))  # 2 rows for O and N, 2 columns for Flexibility and SASA

    # Define plotting function for individual glycan class
    def plot_single_metric(data, glycan_type, axis_flex, axis_sasa):
        # Binding vs Flexibility
        scatter_flex = sns.scatterplot(
            ax=axis_flex,
            x='weighted_mean_flexibility',
            y='binding_score',
            data=data,
            hue="class",
            hue_order=['N', 'O', "free","lipid", 'NaN'],  # Ensure consistent color scheme
            palette="tab10",
            alpha=0.7,
            legend=False  # Suppress legend for this plot
        )
        sns.regplot(
            ax=axis_flex,
            x='weighted_mean_flexibility',
            y='binding_score',
            data=data,
            scatter=False,
            line_kws={'color': 'red'}
        )
        axis_flex.set_title(f'{glycan_type} Glycans: Binding vs Flexibility Weighted\n{lectin} {binding_motif}', fontsize=12)
        axis_flex.set_xlabel('Flexibility')
        axis_flex.set_ylabel('Binding Score')

        # Binding vs SASA
        scatter_sasa = sns.scatterplot(
            ax=axis_sasa,
            x='SASA_weighted',
            y='binding_score',
            data=data,
            hue='class',
            hue_order=['N', 'O', "free","lipid", 'NaN'],  # Ensure consistent color scheme
            palette="tab10",
            alpha=0.7,
            legend=False  # Suppress legend for this plot
        )
        sns.regplot(
            ax=axis_sasa,
            x='SASA_weighted',
            y='binding_score',
            data=data,
            scatter=False,
            line_kws={'color': 'red'}
        )
        axis_sasa.set_title(f'{glycan_type} Glycans: Binding vs SASA Weighted\n{lectin} {binding_motif}', fontsize=12)
        axis_sasa.set_xlabel('SASA Weighted')
        axis_sasa.set_ylabel('Binding Score')

        # Return handles and labels for legend
        handles, labels = scatter_sasa.get_legend_handles_labels()
        return handles, labels

    # Row 1: O-linked glycans
    handles_o, labels_o = plot_single_metric(o_glycans, 'O-linked', axes[0, 0], axes[0, 1])

    # Row 2: N-linked glycans
    handles_n, labels_n = plot_single_metric(n_glycans, 'N-linked', axes[1, 0], axes[1, 1])

    # Create a single unified legend using one of the scatterplots
    fig.legend(
        handles=handles_o[:5],  # Use up to 5 classes for consistency
        labels=labels_o[:5],
        title="Glycan class",
        bbox_to_anchor=(1.02, 0.5),  # Position legend closer to the plot
        loc='center left',
        borderaxespad=0  # Reduce padding
    )

    # Adjust layout to make space for the legend
    plt.subplots_adjust(right=0.8)  # Adjust the right margin to make space for the legend

    # Adjust layout for better spacing
    plt.tight_layout()

    # Save the combined plot
    plt.savefig(f'scripts/correlation/plots/class/Binding_vs_Flexibility_and_SASA_{lectin}_O_N_separated.png', dpi=300)

    # Show the plots
    plt.show()

def visualize_mediation_results_with_class(metric_df,
                                           independent_var,
                                           class_var,
                                           dependent_var,
                                           effects,
                                           lectin,
                                           binding_motif):
    """
    Visualize mediation analysis results with a path diagram, scatterplots, and a bar chart.

    Parameters:
        metric_df (pd.DataFrame): The data containing variables for analysis.
        independent_var (str): Name of the independent variable.
        class_var (str): Name of the mediator variable (categorical, e.g., 'class').
        dependent_var (str): Name of the dependent variable.
        effects (dict): Dictionary containing total, direct, and indirect effects.
        lectin (str): Name of the lectin being analyzed.
        binding_motif (str): Name of the binding motif being analyzed.
        save_dir (str): Directory where plots will be saved. Defaults to the current directory.
    """
    save_dir = f"scripts/correlation/plots/mediation/{independent_var}/"
    os.makedirs(save_dir, exist_ok=True)

    # Unpack effects
    total_effect = effects['total_effect']
    direct_effect = effects['direct_effect']
    indirect_effect = effects['indirect_effect']

    # Step 1: Path Diagram
    def plot_mediation_path():
        G = nx.DiGraph()
        G.add_edge(independent_var, 'Class (Mediator)', weight=indirect_effect)
        G.add_edge('Class (Mediator)', dependent_var, weight=indirect_effect)
        G.add_edge(independent_var, dependent_var, weight=direct_effect)

        pos = {
            independent_var: (0, 1),
            'Class (Mediator)': (1, 0.5),
            dependent_var: (2, 1)
        }

        plt.figure(figsize=(8, 6))
        nx.draw_networkx_nodes(G, pos, node_size=2000, node_color='lightblue')
        edges = nx.get_edge_attributes(G, 'weight')
        nx.draw_networkx_edges(G, pos, arrowstyle='-|>', arrowsize=20)
        nx.draw_networkx_labels(G, pos, font_size=12, font_color='black')
        nx.draw_networkx_edge_labels(
            G, pos, edge_labels={k: f"{v:.3f}" for k, v in edges.items()}, font_color='red'
        )
        plt.title(f"Mediation Path Diagram\nLectin: {lectin}, Motif: {binding_motif}")
        plt.axis('off')
        plt.savefig(os.path.join(save_dir, f"{lectin}_{binding_motif}_mediation_path.png"))
        plt.show()

    # Step 2: Scatterplots
    def plot_scatterplots():
        plt.figure(figsize=(12, 6))

        # Define color palette for classes
        custom_palette = {"N": "blue", "O": "orange"}

        # Independent vs Class (Mediator)
        plt.subplot(1, 3, 1)
        sns.boxplot(x=class_var, y=independent_var, data=metric_df, palette=custom_palette)
        plt.title(f"{independent_var} vs. {class_var}\nLectin: {lectin}, Motif: {binding_motif}")
        plt.xlabel(class_var)
        plt.ylabel(independent_var)

        # Class (Mediator) vs Dependent
        plt.subplot(1, 3, 2)
        sns.boxplot(x=class_var, y=dependent_var, data=metric_df, palette=custom_palette)
        plt.title(f"{class_var} vs. {dependent_var}\nLectin: {lectin}, Motif: {binding_motif}")
        plt.xlabel(class_var)
        plt.ylabel(dependent_var)

        # Independent vs Dependent
        plt.subplot(1, 3, 3)
        sns.regplot(x=independent_var, y=dependent_var, data=metric_df, scatter_kws={'color': 'blue'})
        plt.title(f"{independent_var} vs. {dependent_var}\nLectin: {lectin}, Motif: {binding_motif}")
        plt.xlabel(independent_var)
        plt.ylabel(dependent_var)

        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, f"{lectin}_{binding_motif}_scatterplots.png"))
        plt.show()

    # Step 3: Bar Chart for Effects
    def plot_effects_bar_chart():
        effects_labels = ['Total Effect', 'Direct Effect', 'Indirect Effect']
        effects_values = [total_effect, direct_effect, indirect_effect]
        plt.figure(figsize=(8, 6))
        plt.bar(effects_labels, effects_values, color=['blue', 'green', 'orange'])
        plt.title(f"Mediation Effects\nLectin: {lectin}, Motif: {binding_motif}")
        plt.ylabel("Effect Value")
        plt.savefig(os.path.join(save_dir, f"{lectin}_{binding_motif}_effects_bar_chart.png"))
        plt.show()

    # Call visualizations
    print("\nStep 1: Mediation Path Diagram")
    plot_mediation_path()

    print("\nStep 2: Scatterplots")
    plot_scatterplots()

    print("\nStep 3: Mediation Effects Bar Chart")
    plot_effects_bar_chart()

