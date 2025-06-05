import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D  # Import Line2D for manual legend creation
import numpy as np  # For generating custom tick values

# Define the custom colors for each lineage (species class)
class_colors = {
    'Amphibians': "#228B22",  # Forest green
    'Birds': "#DAA520",  # Goldenrod
    'Fishes': "#4682B4",  # Steel blue
    'Invertebrates': "#DA70D6",  # Orchid
    'Mammals': "#B22222",  # Firebrick
    'Reptiles': "#556B2F"  # Dark olive green
}

def plotN50(path, title, output_trendlines, output_scatter, lookup_path, target):
    # Import data with your original logic
    vgp_assembly_dataset_nxContig = pd.read_csv(path, header=None, sep="\0")
    vgp_assembly_dataset_nxContig = vgp_assembly_dataset_nxContig[0].str.split('\t', expand=True, n=3)
    vgp_assembly_dataset_nxContig.columns = ['Accession', 'Tolid', 'Scientific_name', 'Data']

    # Split the 'Data' column into multiple columns
    values = vgp_assembly_dataset_nxContig['Data'].str.split('\t', expand=True)
    vgp_assembly_dataset_nxContig = pd.concat([vgp_assembly_dataset_nxContig, values], axis=1)
    vgp_assembly_dataset_nxContig.drop('Data', axis=1, inplace=True)

    # Melt the dataframe to long format
    vgp_assembly_dataset_nxContig = vgp_assembly_dataset_nxContig.melt(id_vars=['Accession', 'Tolid'],
                                                                       var_name="Data",
                                                                       value_vars=range(values.shape[1])).dropna()

    # Split the 'value' column into 'Size' and 'Percentage'
    vgp_assembly_dataset_nxContig[['Size', 'Percentage']] = vgp_assembly_dataset_nxContig['value'].str.split(",", expand=True).apply(pd.to_numeric)
    vgp_assembly_dataset_nxContig.drop('value', axis=1, inplace=True)

    # Convert 'Size' to Mbp
    vgp_assembly_dataset_nxContig['Size'] = vgp_assembly_dataset_nxContig['Size'] / 1000000

    # Read the lookup file for lineage info
    lookup = pd.read_csv(lookup_path, sep='\t')

    # Make sure to merge the lineage information with the assembly dataset
    lookup = lookup[['Accession # for main haplotype', 'Lineage']]
    vgp_assembly_dataset_nxContig = vgp_assembly_dataset_nxContig.merge(lookup, left_on='Accession',
                                                                        right_on='Accession # for main haplotype',
                                                                        how='left')

    # Create the first plot (trendlines)
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    ax1.set_yscale('log')
    plt.title(title + " - Trendlines")

    # List to store the handles and labels for the legend
    handles = []
    labels = []

    # Plot the trendlines (regression) without scatter points for each lineage in ax1
    sns.set(style="whitegrid")
    for lineage, grp in vgp_assembly_dataset_nxContig.groupby('Lineage'):
        # Get the color for this lineage from the class_colors dictionary
        color = class_colors.get(lineage, "#000000")  # Default to black if lineage not in class_colors

        # Plot the trendline for each lineage
        line = sns.regplot(data=grp, x='Percentage', y='Size', scatter=False, line_kws={'color': color}, ci=95, ax=ax1)

        # Manually add the legend entry using Line2D
        handles.append(Line2D([0], [0], color=color, lw=2))  # Use the color from the regression line
        labels.append(lineage)  # Use the lineage as the label

    ax1.set_xlabel('Nx (%)')
    ax1.set_ylabel('Mbp')

    # Set the y-axis tick format to show 1 decimal place
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.2f}'.format(x)))

    # Add target line
    ax1.axhline(y=target, color='black', linestyle='--', linewidth=2)

    # Custom x-axis ticks: every 0.1
    ax1.set_xticks(np.arange(min(vgp_assembly_dataset_nxContig['Percentage']),
                             max(vgp_assembly_dataset_nxContig['Percentage']),
                             0.1))

    # Set the x-tick format to show 1 decimal place
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x)))

    # Show legend for lineages in ax1, place it outside the plot
    ax1.legend(handles=handles, labels=labels, title='Lineage', loc='upper left', bbox_to_anchor=(1, 1), frameon=False)

    # Save the trendlines plot
    plt.tight_layout()
    plt.subplots_adjust(top=0.85, right=0.85)  # Adjust top to make room for the suptitle
    plt.savefig(output_trendlines, dpi=300)

    # Create the second plot (scatter plot)
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    ax2.set_yscale('log')
    plt.title(title + " - Scatter Plot")

    # Plot individual data points on ax2 (with smaller markers)
    for lineage, grp in vgp_assembly_dataset_nxContig.groupby('Lineage'):
        # Get the color for this lineage from the class_colors dictionary
        color = class_colors.get(lineage, "#000000")  # Default to black if lineage not in class_colors
        # Plot the scatter points with smaller size
        ax2.scatter(grp['Percentage'], grp['Size'], color=color, label=lineage, alpha=0.5, s=10)

    ax2.set_xlabel('Nx (%)')
    ax2.set_ylabel('Mbp')

    # Set the y-axis tick format to show 1 decimal place
    ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.2f}'.format(x)))

    # Add target line
    ax2.axhline(y=target, color='black', linestyle='--', linewidth=2)

    # Custom x-axis ticks: every 0.1
    ax2.set_xticks(np.arange(min(vgp_assembly_dataset_nxContig['Percentage']),
                             max(vgp_assembly_dataset_nxContig['Percentage']),
                             0.1))

    # Set the x-tick format to show 1 decimal place
    ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.1f}'.format(x)))

    # Show legend for lineages in ax2, place it outside the plot
    ax2.legend(title='Lineage', loc='upper left', bbox_to_anchor=(1, 1), frameon=False)

    # Save the scatter plot
    plt.tight_layout()
    plt.subplots_adjust(top=0.85, right=0.85)  # Adjust top to make room for the suptitle
    plt.savefig(output_scatter, dpi=300)


# Example usage
plotN50('gfastatsNxContig.tsv', 'vgp-assembly-manuscript Nx Contig', 'NxContig_Trendlines.png', 'NxContig_Scatter.png',
        'VGPPhase1-freeze-1.0.tsv', 1)
plotN50('gfastatsNxScaffold.tsv', 'vgp-assembly-manuscript Nx Scaffold', 'NxScaffold_Trendlines.png',
        'NxScaffold_Scatter.png', 'VGPPhase1-freeze-1.0.tsv', 10)
