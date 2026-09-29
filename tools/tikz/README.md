# Vignette model diagrams

`restriction-models.tex` is the editable TikZ source. The vignette includes the
committed SVG for HTML and PDF for LaTeX output, so building the R package does
not require TikZ or a PDF-to-SVG converter.

From the repository root, rebuild with TeX Live (including standalone and TikZ)
and Poppler:

```sh
mkdir -p /tmp/semtests-vignette-tikz
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=/tmp/semtests-vignette-tikz tools/tikz/restriction-models.tex
cp /tmp/semtests-vignette-tikz/restriction-models.pdf vignettes/figures/restriction-models.pdf
pdftocairo -svg vignettes/figures/restriction-models.pdf vignettes/figures/restriction-models.svg
```

Both panels have identical scaling and geometry. Only the x2/x3 loadings are
highlighted. Intercepts and variances are omitted, as stated in the caption;
residual arrows are included to avoid suggesting error-free indicators.
