# Figure-generating code

- `figures_12.py`: genetic evidence, cell atlases, complete association and correspondence panels.
- `figure_3.py`: fixed exclusions and normal-repair time course.
- `figures_45.py`: genetic candidates, TF context, embedding sensitivity and controlled-treatment results.
- `plot_supplementary_figures.py`: supporting panels from completed evidence tables.
- `.jsx` / `.jsxinc`: Adobe Illustrator text, arrows and panel-label composition. These are source code, not generated image files.

Only source code is supplied here. The plotting modules require a context providing `data`, `panel` and `save`, as documented in their docstrings; map input tables from the source index before use. Original illustration artwork is external and is not required to inspect scientific result tables. No generated images or manuscript files are uploaded.

## Portable numeric-panel rendering

After restoring the exact completed input tables listed in `config/figure_inputs.json`, run:

```sh
python figures/render_figures.py --registry config/figure_inputs.json --out outputs/figures --figures 1 2 3 4 5
```

This wrapper removes original host/path/resource orchestration while leaving the plotting modules and numeric inputs unchanged. It has passed syntax and CLI checks; it has not been certified by a new full render here. Figure 1A conceptual artwork and subsequent Adobe composition remain separate author-supplied display assets and are not required to inspect the numeric panels.
