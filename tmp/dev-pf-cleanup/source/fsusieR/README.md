# fsusieR

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

The `fsusieR` package implements the Sum of Single Functions (fSuSiE) method for fine-mapping molecular QTLs from spatially structured traits like DNA methylation and histone modifications. fSuSiE extends the Sum of Single Effects (SuSiE) framework using wavelet-based functional regression to capture spatial correlations in molecular traits along the genome.

## System Requirements

### Software Dependencies

- R (≥ 3.5.0)
- Required R packages: wavethresh (≥ 4.6.8), ashr (≥ 2.2-47), mixsqp (≥ 0.3-43), matrixStats (≥ 0.62.0), smashr
- Optional: susieR and mvsusieR for comparation with fSuSiE

### Operating Systems

- macOS 12.0+
- Ubuntu 20.04+  
- Windows 10+

### Hardware

No special hardware required. Standard desktop computer with 8GB RAM sufficient for typical analyses.

## Installation

Install the latest version from GitHub (typical install time: < 1 minute):

```R
# install.packages("remotes")
remotes::install_github("stephenslab/smashr")
remotes::install_github("stephenslab/fsusieR")
```

## Demo

A brief demo with simulated methylation data can be found
[here](https://stephenslab.github.io/fsusieR/articles/methyl_demo.html). 

See the [pkgdown website](https://stephenslab.github.io/fsusieR) for
other examples.

## Citing this work

If you use `fsusieR` in your work, please cite:

> Denault, W.R.P., Sun, H., Carbonetto, P., Li, A.,  De Jager, L.P., Bennett, D, The Alzheimer’s Disease Functional Genomics Consortium, Wang, G. & Stephens, M. (2025). fSuSiE enables fine-mapping of QTLs from genome-scale molecular profiles. *bioRxiv* [DOI: 10.1101/2025.08.17.670732](https://doi.org/10.1101/2025.08.17.670732)

## License

This project is licensed under the BSD-3-Clause License - see the [LICENSE](LICENSE) file for details.

## Support

Please [post issues](https://github.com/stephenslab/fsusieR/issues) for questions or bug reports.
