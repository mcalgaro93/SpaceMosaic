# Methods figures

The figures in this directory are generated from the scripts in
`../scripts/`. Do not edit the exported image files manually.

Each script uses deterministic synthetic data and writes three formats:

- SVG for vector editing and web display;
- PDF for manuscript submission;
- high-resolution PNG for Markdown and Word previews.

Regenerate the multiscale embedding figure from the repository root with:

```r
Rscript manuscript/scripts/figure_multiscale_embedding.R
```

Regenerate the patch geometry figure with:

```r
Rscript manuscript/scripts/figure_patch_geometry.R
```

Regenerate the patch-iteration figure with:

```r
Rscript manuscript/scripts/figure_patch_iteration.R
```
