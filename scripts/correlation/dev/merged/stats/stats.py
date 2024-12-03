from scipy.stats import pearsonr, spearmanr, chi2_contingency
import numpy as np

#from scripts.correlation.func.metric_df import metric_df
#from scripts.correlation.func.plot_corr_regg import plot_combined

# Dictionary of lectins and their binding motifs
lectin_binding_motif = {
    "AAL": ["Fuc"],
    "SNA": ["Neu5Ac(a2-6)", "Neu5Gc(a2-6)"],
    "ConA": ["Man(a1-2)"],
    "MAL-II": ["Neu5Ac"],
    "PNA": ["Gal(b1-3)", "GalNAc"]
}


# Function to calculate correlation statistics
def calculate_correlation_stats(df, x_col, y_col):
    """
    Calculate Pearson and Spearman correlation and p-values.

    Args:
        df (pd.DataFrame): DataFrame containing the data.
        x_col (str): Column name for the independent variable (e.g., flexibility).
        y_col (str): Column name for the dependent variable (e.g., binding score).

    Returns:
        dict: Dictionary with correlation coefficients and p-values.
    """
    if df[x_col].isnull().any() or df[y_col].isnull().any():
        return {"pearson": (np.nan, np.nan), "spearman": (np.nan, np.nan)}

    pearson_corr, pearson_p = pearsonr(df[x_col], df[y_col])
    spearman_corr, spearman_p = spearmanr(df[x_col], df[y_col])

    return {
        "pearson": (pearson_corr, pearson_p),
        "spearman": (spearman_corr, spearman_p)
    }


# Function to calculate a contingency table and chi-square test
def calculate_chi_square(df, x_col, y_col, threshold=0.5):
    """
    Calculate chi-square test for contingency table.

    Args:
        df (pd.DataFrame): DataFrame containing the data.
        x_col (str): Column name for one categorical variable.
        y_col (str): Column name for the other categorical variable.
        threshold (float): Threshold to binarize numerical data into categories.

    Returns:
        float: p-value from the chi-square test.
    """
    df['x_bin'] = (df[x_col] > threshold).astype(int)
    df['y_bin'] = (df[y_col] > threshold).astype(int)

    contingency_table = pd.crosstab(df['x_bin'], df['y_bin'])
    _, p, _, _ = chi2_contingency(contingency_table)
    return p


# Main loop to process each lectin
for lectin, binding_motif in lectin_binding_motif.items():
    # Generate the metric DataFrame
    metric_df_instance = metric_df(lectin, binding_motif)

    # Perform statistical analysis for correlation
    x_col = "flexibility"  # Replace with the actual column name in your DataFrame
    y_col = "binding_score"  # Replace with the actual column name in your DataFrame
    correlation_stats = calculate_correlation_stats(metric_df_instance, x_col, y_col)

    # Perform chi-square test (optional)
    # Assuming you're comparing two binary conditions like "high flexibility" vs "high binding"
    chi_square_p = calculate_chi_square(metric_df_instance, x_col, y_col)

    # Print results for each lectin
    print(f"Lectin: {lectin}")
    print(f"Pearson Correlation: {correlation_stats['pearson'][0]:.2f} (p = {correlation_stats['pearson'][1]:.4f})")
    print(f"Spearman Correlation: {correlation_stats['spearman'][0]:.2f} (p = {correlation_stats['spearman'][1]:.4f})")
    print(f"Chi-Square Test p-value: {chi_square_p:.4f}")
    print("-" * 40)

    # Generate plots
    plot_combined(metric_df_instance, lectin, binding_motif)
