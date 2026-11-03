# Tested R environment

The local execution used R 4.5.1 with Cairo PDF support on Linux. The following
versions were loaded or checked in that environment. `install_dependencies.R`
installs missing packages; it is a convenience installer, not a version lock.
Use the versions below when reproducing the checked environment.

| Package | Version |
|---|---|
| ape | 5.8.1 |
| circlize | 0.4.16 |
| ComplexHeatmap | 2.25.3 |
| cowplot | 1.2.0 |
| data.table | 1.17.8 |
| dplyr | 1.1.4 |
| ggnewscale | 0.5.2 |
| ggplot2 | 4.0.0 |
| ggrepel | 0.9.6 |
| ggsci | 5.2.0 |
| ggtree | 3.16.3 |
| jsonlite | 2.0.0 |
| patchwork | 1.3.2 |
| readr | 2.1.5 |
| scales | 1.4.0 |
| tibble | 3.3.0 |
| tidyr | 1.3.1 |
| vegan | 2.7.2 |

`grid`, `graphics`, `grDevices`, `stats`, `tools` and `utils` are supplied with R.
The default font is `sans`; no local font files or font-configuration paths are
required. No compiled analysis executables, Python environment or original
project scripts are required at run time. Input Newick trees, model results and
other upstream scientific products must already be available privately.

Some ggtree versions emit ggplot2 deprecation messages. In the tested environment
these messages did not prevent tree rendering. Changing package versions or fonts
can change layout; the scientific scopes and numerical checks remain explicit.
