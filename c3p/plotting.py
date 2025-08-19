import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize


def plot_scatter(df, x, y, title="Scatterplot", results_dir=None, xlabel=None, ylabel=None):
    """
    Create a scatter plot with optional logarithmic x-axis.

    Args:
        df: DataFrame containing the data
        x: Name of x-axis column
        y: Name of y-axis column
        title: Plot title
        results_dir: Optional directory to save the plot
        xlabel: Optional x-axis label
        ylabel: Optional y-axis label
    """
    plt.clf()
    # check if x and y are in the dataframe
    if x not in df.columns or y not in df.columns:
        print(f"Columns {x} and {y} not found in DataFrame")
        return
    # drop na values
    df = df.dropna(subset=[x, y])
    if df.size == 0:
        return
    correlation = np.corrcoef(df[x], df[y])[0, 1]
    # print(f"Correlation between {x} and {y}: {correlation:.2f}")
    title = f"{title} (r = {correlation:.2f})"

    # Create scatter plot
    plt.scatter(df[x], df[y], s=5, alpha=0.5)

    # Draw horizontal line at mean y
    plt.axhline(df[y].mean(), color='red', linestyle='dashed', linewidth=1)

    # Fit trend line using transformed x values
    z = np.polyfit(df[x], df[y], 1)
    p = np.poly1d(z)

    # Generate x points for trend line
    plt.plot(df[x], p(df[x]), "r--", alpha=0.8)

    plt.xlabel(xlabel or x)
    plt.ylabel(ylabel or y)
    plt.title(title)
    plt.grid(True, alpha=0.3)

    if results_dir:
        plt.savefig(results_dir / f"scatterplot-{x}-{y}.png",
                    dpi=300,  # High resolution
                    bbox_inches='tight',  # Removes extra white spaces
                    pad_inches=0.1,  # Small padding around the figure
                    format='png',  # Format type
                    transparent=False,  # White background
                    facecolor='white',  # Figure face color
                    edgecolor='none',  # No edge color
                    )
    plt.close()
    return plt



from scipy import stats


def create_scatter_matrix(df: pd.DataFrame, compare='experiment_name', metric: str = 'f1', index='chemical_class') -> plt.Figure:
    """
    Create a scatter matrix to compare different experiments based on a specified metric.

    Args:
        df:
        compare:
        metric:
        index:

    Returns:

    """
    # Pivot the data to get experiments as columns and classes as rows
    df = df.dropna(subset=[metric])
    pivot_df = df.pivot(index=index, columns=compare, values=metric)

    # Calculate number of experiments
    n_experiments = len(pivot_df.columns)

    # Create figure and axis grid
    fig, axes = plt.subplots(n_experiments, n_experiments, figsize=(15, 15))

    # Add padding between subplots
    plt.subplots_adjust(hspace=0.6, wspace=0.3)

    # Iterate through each pair of experiments
    for i, exp1 in enumerate(pivot_df.columns):
        for j, exp2 in enumerate(pivot_df.columns):
            ax = axes[i, j]
            exp1_vals = pivot_df[exp1].fillna(0)
            exp2_vals = pivot_df[exp2].fillna(0)

            if i > j:  # Lower triangle: scatter plots
                # Create scatter plot
                ax.scatter(exp2_vals, exp1_vals, alpha=0.5)

                # Calculate correlation
                corr, _ = stats.pearsonr(exp2_vals, exp1_vals)

                # Add correlation coefficient text above the plot
                # ax.set_title(f'r = {corr:.3f}', pad=8, fontweight='bold', fontsize=10)

                ax.text(0.5, 1.15, f'r = {corr:.3f}',
                        ha='center', va='bottom',
                        transform=ax.transAxes,
                        fontweight='bold', fontsize=10)

                # Add correlation coefficient text
                #ax.text(0.05, 0.95, f'r = {corr:.3f}',
                #        transform=ax.transAxes,
                #        verticalalignment='top')

                # Set limits from 0 to 1 for F1 scores
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)

                # Add diagonal line
                ax.plot([0, 1], [0, 1], 'k--', alpha=0.3)

            elif i == j:  # Diagonal: experiment names
                ax.text(0.5, 0.5, exp1,
                        ha='center', va='center',
                        transform=ax.transAxes,
                        rotation=45)
                ax.set_xticks([])
                ax.set_yticks([])

            else:  # Upper triangle: keep empty
                ax.set_visible(False)

            # Only show labels on outer edges
            if i == n_experiments - 1:
                ax.set_xlabel(exp2)
            if j == 0:
                ax.set_ylabel(exp1)

    plt.suptitle('Model Comparison Matrix', size=16, y=1.02)
    return fig


def create_scatter_matrix_train_colored(df: pd.DataFrame, experiment_name: str, test_metric: str = 'f1', 
                                       train_metric: str = 'train_f1', index: str = 'chemical_class') -> plt.Figure:
    """
    Create a scatter matrix for a single experiment where points are colored by train F1 scores.
    
    Args:
        df: DataFrame containing the data
        experiment_name: Name of the experiment to filter for
        test_metric: Name of the test metric column (typically 'f1')
        train_metric: Name of the train metric column for coloring (typically 'train_f1')
        index: Column to use as index (typically 'chemical_class')
        
    Returns:
        matplotlib Figure object
    """
    # Filter for the specific experiment
    exp_df = df[df['experiment_name'] == experiment_name].copy()
    
    # Drop rows where either test or train metric is missing
    exp_df = exp_df.dropna(subset=[test_metric, train_metric])
    
    if len(exp_df) == 0:
        print(f"No data found for experiment '{experiment_name}' with both {test_metric} and {train_metric}")
        return None
    
    # Get the complexity metrics for scatter matrix
    complexity_metrics = ['complexity', 'max_indent', 'log_lines_of_code', 'methods_called_count', 
                         'returns_count', 'smarts_strings_count']
    
    # Filter to only include metrics that exist in the dataframe and have values
    available_metrics = [metric for metric in complexity_metrics if metric in exp_df.columns]
    available_metrics = [metric for metric in available_metrics if not exp_df[metric].isna().all()]
    
    # Add the test metric to the list
    plot_metrics = available_metrics + [test_metric]
    n_metrics = len(plot_metrics)
    
    if n_metrics < 2:
        print(f"Need at least 2 metrics for scatter matrix, found {n_metrics}")
        return None
    
    # Create figure and axis grid
    fig, axes = plt.subplots(n_metrics, n_metrics, figsize=(15, 15))
    
    # Handle case where we only have 2 metrics (axes won't be 2D)
    if n_metrics == 2:
        axes = axes.reshape(2, 2)
    
    # Add padding between subplots
    plt.subplots_adjust(hspace=0.6, wspace=0.3)
    
    # Prepare color mapping based on train metric
    train_values = exp_df[train_metric]
    norm = Normalize(vmin=train_values.min(), vmax=train_values.max())
    colormap = cm.viridis
    
    # Iterate through each pair of metrics
    for i, metric1 in enumerate(plot_metrics):
        for j, metric2 in enumerate(plot_metrics):
            ax = axes[i, j]
            
            if i > j:  # Lower triangle: scatter plots
                # Get values for both metrics, dropping NAs
                valid_data = exp_df[[metric1, metric2, train_metric]].dropna()
                
                if len(valid_data) == 0:
                    ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                    continue
                
                metric1_vals = valid_data[metric1]
                metric2_vals = valid_data[metric2]
                train_vals = valid_data[train_metric]
                
                # Create scatter plot with color based on train metric
                scatter = ax.scatter(metric2_vals, metric1_vals, c=train_vals, 
                                   cmap=colormap, norm=norm, alpha=0.7, s=30)
                
                # Calculate correlation
                if len(metric1_vals) > 1:
                    corr, _ = stats.pearsonr(metric2_vals, metric1_vals)
                    ax.text(0.5, 1.15, f'r = {corr:.3f}',
                            ha='center', va='bottom',
                            transform=ax.transAxes,
                            fontweight='bold', fontsize=10)
                
                # Set appropriate limits
                if metric1 == test_metric or metric2 == test_metric:
                    if metric1 == test_metric:
                        ax.set_ylim(0, 1)
                    if metric2 == test_metric:
                        ax.set_xlim(0, 1)
                
            elif i == j:  # Diagonal: metric names
                ax.text(0.5, 0.5, metric1,
                        ha='center', va='center',
                        transform=ax.transAxes,
                        rotation=45, fontweight='bold')
                ax.set_xticks([])
                ax.set_yticks([])
                
            else:  # Upper triangle: keep empty
                ax.set_visible(False)
            
            # Only show labels on outer edges
            if i == n_metrics - 1:
                ax.set_xlabel(metric2)
            if j == 0:
                ax.set_ylabel(metric1)
    
    # Add colorbar
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=colormap), ax=axes, shrink=0.8, aspect=30)
    cbar.set_label(f'{train_metric} (color scale)', rotation=270, labelpad=20)
    
    plt.suptitle(f'Scatter Matrix for {experiment_name}\n(colored by {train_metric})', size=16, y=1.02)
    return fig


def plot_method_comparison_scatter(df: pd.DataFrame, method1: str, method2: str, 
                                 test_metric: str = 'f1', train_metric: str = 'train_f1',
                                 index: str = 'chemical_class') -> plt.Figure:
    """
    Create a scatter plot comparing F1 scores of two methods, colored by train F1.
    
    Args:
        df: DataFrame containing the data
        method1: Name of first method (x-axis)
        method2: Name of second method (y-axis) 
        test_metric: Name of the test metric column (typically 'f1')
        train_metric: Name of the train metric column for coloring (typically 'train_f1')
        index: Column to use as index (typically 'chemical_class')
        
    Returns:
        matplotlib Figure object
    """
    # Pivot the data to get methods as columns and classes as rows
    pivot_df = df.pivot(index=index, columns='experiment_name', values=[test_metric, train_metric])
    
    # Extract F1 scores for both methods
    method1_f1 = pivot_df[test_metric][method1].dropna()
    method2_f1 = pivot_df[test_metric][method2].dropna()
    
    # Get train F1 for coloring - use whichever method has train data
    train_color_data = None
    color_method = None
    
    if train_metric in pivot_df.columns:
        if method1 in pivot_df[train_metric].columns and not pivot_df[train_metric][method1].isna().all():
            train_color_data = pivot_df[train_metric][method1]
            color_method = method1
        elif method2 in pivot_df[train_metric].columns and not pivot_df[train_metric][method2].isna().all():
            train_color_data = pivot_df[train_metric][method2] 
            color_method = method2
    
    # Find common classes (intersection of both methods)
    common_classes = method1_f1.index.intersection(method2_f1.index)
    
    if len(common_classes) == 0:
        print(f"No common classes found between {method1} and {method2}")
        return None
        
    # Filter to common classes
    method1_vals = method1_f1.loc[common_classes]
    method2_vals = method2_f1.loc[common_classes]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # If we have train data for coloring
    if train_color_data is not None:
        train_vals = train_color_data.loc[common_classes].dropna()
        # Only use classes that have train data for coloring
        color_classes = train_vals.index
        common_with_train = common_classes.intersection(color_classes)
        
        if len(common_with_train) > 0:
            # Plot points with color coding
            norm = Normalize(vmin=train_vals.min(), vmax=train_vals.max())
            colormap = cm.viridis
            
            scatter = ax.scatter(method1_vals.loc[common_with_train], 
                               method2_vals.loc[common_with_train],
                               c=train_vals.loc[common_with_train],
                               cmap=colormap, norm=norm, alpha=0.7, s=50)
            
            # Add colorbar
            cbar = fig.colorbar(scatter, ax=ax, shrink=0.8)
            cbar.set_label(f'{train_metric} ({color_method})', rotation=270, labelpad=20)
            
            # Plot remaining points (without train data) in gray
            remaining_classes = common_classes.difference(common_with_train)
            if len(remaining_classes) > 0:
                ax.scatter(method1_vals.loc[remaining_classes],
                          method2_vals.loc[remaining_classes], 
                          c='lightgray', alpha=0.5, s=50, label='No train data')
        else:
            # No train data available, plot all in single color
            ax.scatter(method1_vals, method2_vals, alpha=0.7, s=50)
    else:
        # No train data available, plot all in single color  
        ax.scatter(method1_vals, method2_vals, alpha=0.7, s=50)
    
    # Calculate correlation
    if len(method1_vals) > 1:
        corr, p_value = stats.pearsonr(method1_vals, method2_vals)
        ax.text(0.05, 0.95, f'r = {corr:.3f}, p = {p_value:.3f}',
                transform=ax.transAxes, fontsize=12,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Add diagonal line (perfect agreement)
    min_val = min(method1_vals.min(), method2_vals.min())
    max_val = max(method1_vals.max(), method2_vals.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, label='Perfect agreement')
    
    # Set labels and limits
    ax.set_xlabel(f'{method1} {test_metric}')
    ax.set_ylabel(f'{method2} {test_metric}')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Add legend if we have different point types
    if train_color_data is not None and len(common_classes.difference(train_vals.index)) > 0:
        ax.legend()
    
    # Title
    title = f'{method1} vs {method2} ({test_metric})'
    if train_color_data is not None:
        title += f'\nColored by {train_metric} ({color_method})'
    ax.set_title(title, fontsize=14)
    
    plt.tight_layout()
    return fig

