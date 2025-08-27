# Visualization

Following three-dimensional disease network analysis, ***DiNetxify*** provides visualization tools to facilitate the interpretation of results. These tools enable the visualization of key findings from PheWAS, comorbidity networks, and disease trajectory analysis. The main visualization class, `DiNetxify.visualization.Plot`, accepts results in the form of a `pandas.DataFrame` and produces interactive plots.

## Initializing the Plot object

To create visualizations, you must initialize a `DiNetxify.visualization.Plot` object with your analysis results. You should pass in the PheWAS, comorbidity network, and disease trajectory results `pandas.DataFrame`, and optionally specify how the “exposure” (the primary risk factor and cohort definition) appear in the plots.

For example:

```python
from DiNetxify.visualization import Plot

# Suppose phewas_result, com_network_result, trajectory_result are obtained from three-dimensional disease network analysis
result_plot = Plot(
    phewas_result=phewas_result,           # Results of PheWAS analysis
    comorbidity_result=com_network_result, # Results of comorbidity network analysis
    trajectory_result=trajectory_result,   # Results of disease trajectory analysis
    exposure_name='Short LTL',             # Name of the exposure (for labeling the exposure node). Use None for exposed-only cohorts.
    exposure_size=15,                      # Relative size scaling for the exposure node (to make it prominent). None for exposed-only cohorts.
    exposure_location=(0, 0, 0)            # 3D coordinates for the exposure node. If None, it defaults to (0,0,0).
)

# If this were an exposed-only cohort (no explicit exposure variable), you would set:
result_plot = Plot(  
    phewas_result=phewas_result,  
    comorbidity_result=com_network_result,  
    trajectory_result=trajectory_result,  
    exposure_name=None,
    exposure_size=None,
    exposure_location=None                
)  
```

- **phewas_result** – Result (pandas.DataFrame) from the PheWAS analysis (must include columns for phecode identifier, disease system identifier, case counts identifier, significance identifier).
- **comorbidity_result** – Result (pandas.DataFrame) from the comorbidity network analysis. It must include columns for disease (D1, D2) identifiers, non-temporal pairs name identifier, the association metrics (beta or correlation) identifier, and significance identifier.
- **trajectory_result** – Result (pandas.DataFrame) from the disease trajectory analysis. It must have columns for the disease (D1, D2) identifiers, temporal pairs name identifier, the association metrics (beta or correlation), and significance identifier.
- **exposure_name** – Name for the exposure (the factor that defines exposed vs unexposed in the cohort). In our examples, the exposure was “Short LTL” (short leukocyte telomere length) in case study 1, hence the example. If you are analyzing something like “smoking” or “diabetes” as the exposure, you’d put that. For an exposed-only study, use None (because there isn’t a separate exposure node to highlight).
- **exposure_location** – The (x, y, z) coordinates where the exposure node should be placed in the plots. By default, if None, the exposure node will be placed at the origin (0,0,0). This is relevant only for 3D plotting; if exposure is None, this is ignored.
- **exposure_size** – A scaling factor for the exposure node’s size in the plots. Increase this to make the exposure node larger relative to disease nodes (to emphasize it). If None, in an exposed-only design, the exposure node is not present.

The `DiNetxify.visualization.Plot` class will internally verify that the required columns exist in the input `pandas.DataFrame` (for example, it expects certain default column names like `'phecode_d1'`, `'phecode_d2'` for pair identifiers, `'phecode'` for disease codes in results of PheWAS analysis, `'system'` for disease system category, `'name_disease_pair'` for a unique pair name, `'..._significance'` for significance identifier, etc. If you did not change column names, it uses these defaults). You can override these defaults by passing optional parameters if your DataFrame uses different column names (we’ll mention those as needed in each plot function below).

Now that `DiNetxify.visualization.Plot` (`result_plot`) is created, we can generate specific plots. Each plot either an interactive HTML file or a static image, as noted.

## PheWAS plot

The `DiNetxify.visualization.Plot.phewas_plot()` function generates a summary plot of the PheWAS results, essentially visualizing the diseases that occur significantly. In a standard or matched cohort study, this function will create a circular heatmap plot to show hazard ratios (HRs) for each significant disease outcome with the exposure. In an exposure-only group, this function will also generate a circular heatmap plot, but it is used to display the case count for each disease. Additionally, diseases are classified by category/system, and the shade of color in the plots represent the relative size of HRs or case count.

For example:

```python
# Generate a PheWAS plot  
result_plot.phewas_plot(  
    path="/your/project/path/phewas_plot.png",  # output file path (supports .png, .svg, .jpg)  
    is_exposure_only=False                      # False for cohort/matched designs; True if this is an exposed-only cohort  
)  
```

This function will save the plot to the specified file path. Supported formats include `.png` for static images (suitable for publications) and `.svg` for scalable vector graphics, among others

Parameters for `phewas_plot()` include:

- **path** – File path including filename and extension. Ensure the extension is one of the supported image types.
- **is_exposure_only** – Boolean flag; set to `True` if your study is an exposed-only cohort. For our example (standard/matched cohort), it’s `False`.

**Optional parameters:** (if your analysis results `pandas.DataFrame` columns differ from defaults or if you want to adjust aesthetics)

- **col_coef** – Name of the column in `phewas_result` that means effect size coefficients. *(Default: 'phewas_coef')*.
- **col_se** – Name of the column for standard error of the effect size. *(Default: 'phewas_se')*.
- **col_system** – Column name for the disease system/category. *(Default: 'system')*.
- **col_disease** – Column name for the disease name. *(Default: 'disease')*.
- **col_exposure** – Column name for number of cases in exposed group (used in an exposed-only cohort study). *(Default: 'N_cases_exposed')*.
- **disease_font_size** - Font size for disease labels in the plot. *(Default: 10)*.
- **system_font_size** - Font size for the disease system labels in the plot. *(Default: 17)*.
- **dpi** – Resolution of the output image (dots per inch). *(Default: 200)*.

## Comorbidity network plot

The `DiNetxify.visualization.Plot.comorbidity_network_plot()` function creates an interactive network visualization for the comorbidity network. It clusters diseases into communities based on their network connections (using the Louvain community detection algorithm), often in a circular layout where each community occupies a sector. The nodes (diseases) are colored by their disease system, and edges represent significant associations from the comorbidity network analysis. This plot is output as an HTML file, and you can hover to see details, zoom, etc.

For example:

```python
# Generate an interactive comorbidity network plot  
result_plot.comorbidity_network_plot(  
    path="/your/project/path/comorbidity_network.html"  # output file path (only supports .html)
)  
```

This function will generate an HTML file, which you can open with a web browser for viewing. In this plot, the nodes within the same community will cluster together, and you can view the relevant information of each node by hovering. This plot is helpful in identifying disease groups that occur frequently together beyond the expected frequency.

Parameters for `comorbidity_network_plot()` include:

- **path** – File path including filename and extension.

**Optional parameters:**

- **max_radius** – Maximum radial distance (in pixels) from the center that nodes can be placed. *(Default: 90.0)*.
- **min_radius** – Minimum radial distance (in pixels) from the center that nodes can be placed. *(Default: 35.0)*.
- **layer_distance** – Distance between concentric circles. *(Default: 40.0)*.
- **size_reduction** – A scaling factor (0 to 1) for node sizes to ensure they fit nicely (smaller values make nodes proportionally smaller). *(Default: 0.5)*.
- **line_width** – Width (pixels) of the lines (edges) connecting nodes. *(Default: 1.0)*.
- **line_color** – Color of the edges. Specify a named color (e.g., `'steelblue'`), hex code (`'#4682B4'`), or RGB tuple. *(Default: 'black')*.

These parameters allow fine-tuning the appearance if needed (for example, if node labels overlap, you might reduce node sizes or adjust radii). Usually, the defaults produce a clear separation of communities.

## Disease trajectory plot

The `DiNetxify.visualization.Plot.trajectory_plot()` function generates individual plots for each identified community, visualizing the disease trajectories contained within a single community. For each significant temporal disease pair, an arrow or directed edge is drawn from the antecedent disease to the consequent disease. The output is a set of static image files (one `.png` file per community) because it might be easier to print or inspect individually.

For example:

```python
# Generate disease trajectory plots (one per community)  
result_plot.trajectory_plot(
    path="/your/project/path/trajectory_plots/"  # output file path (just need a folder path and "/")
)
```

> **Important:** Here, `path` is a directory (you should include the trailing slash). The function will then save multiple files in this directory, named perhaps by community or numbered sequentially. Ensure the directory exists or the function might attempt to create it.

Parameters for `trajectory_plot()` include:

- **path** – File path (just need a folder path and "/").

**Optional parameters:**

- **source** – Column name in the DataFrames for the source disease (D1). *(Default: 'phecode_d1')*.
- **target** – Column name for the target disease (D2). *(Default: 'phecode_d2')*.
- **dpi** – Resolution of the output image (dots per inch). *(Default: 500 for these plots, to ensure high clarity since many arrows/labels might be present.)*
- **cluster_weight** – This parameter specifies which edge weight from the network to use when arranging the layout. *(Default: 'comorbidity_beta')*

The plots will illustrate the progression of diseases. Within each community, you might see something like a directed acyclic graph (hopefully acyclic if the data suggests a direction). The arrow from disease A to B indicates A tends to precede B. By generating one plot per community, it keeps the graphs manageable and interpretable. If a community has no significant trajectories, it is empty.

After running this, check the output directory for files. For instance, you might find files like `community_1.png`, `community_2.png`, etc., each showing the disease trajectory for that community of diseases.

## Three-dimensional plot

The `DiNetxify.visualization.Plot.three_dimension_plot()` function generates an integrated 3D interactive visualization combining both the comorbidity network and disease trajectory. In this representation, the comorbidity network is displayed on a horizontal plane, while disease trajectory are aligned along the vertical axis, providing a comprehensive three-dimensional perspective of both relationship types simultaneously. The exposure node (if present) is positioned at the center, with disease nodes arranged concentrically around it. The resulting plot is exported as an interactive HTML file, enabling user-driven rotation and zooming to facilitate detailed exploration of the 3D structure.

For example:

```python
# Generate a combined 3D network plot  
result_plot.three_dimension_plot(  
    path="/your/project/path/combined_network.html"  # output file path (only supports .html)
)  
```

This function will save an HTML file with the interactive 3D visualization. When you open it, you can drag to rotate the network in 3D space. One plane (e.g., viewed from top-down) will show the clusters and connections akin to the comorbidity network, and the other plane (viewed from the side) will show the directed trajectories. Where they intersect, you get a full picture of how diseases group and progress.

Parameters for `three_dimension_plot()` include:

- **path** – File path including filename and extension.

**Optional parameters:**

- **max_radius** – Maximum radius for node placement from the center (similar to the 2D network plot). *(Default: 180.0)*.
- **min_radius** – Minimum radius from center. *(Default: 35.0)*.
- **layer_distance** – Distance between layers in the radial direction (similar concept to above). *(Default: 40.0)*.
- **layout_width** – Width of the overall figure in pixels. *(Default: 900.0)*.
- **layout_height** – Height of the figure in pixels. *(Default: 900.0)*.
- **line_color** – Color for trajectory lines (since in 3D plot, maybe comorbidity edges are one style and trajectory edges another). *(Default: 'black')*.
- **line_width** – Width of trajectory lines. *(Default: 1.0)*.
- **size_reduction** – Node size scaling factor (0.1–1.0). *(Default: 0.5)*.
- **cluster_reduction_ratio** – A factor (0.1–1.0) to compress or spread out clusters in the 3D space. Lower means clusters are more tightly grouped. *(Default: 1.0)*.
- **cluster_weight** – Which weight to use for determining cluster layout (as in trajectory_plot, usually 'comorbidity_beta'). *(Default: 'comorbidity_beta')*.
- **font_style** – Font family for text elements (node labels, etc.). *(Default: 'Times New Roman')*.
- **font_size** – Base font size for text. *(Default: 15.0)*.

With all these visualization tools, even a beginner user can not only run the analysis with ***DiNetxify*** but also see the results in intuitive forms, which can greatly aid interpretation and presentation. Next, we provide a quick reference to the API of the main classes and functions for convenience.