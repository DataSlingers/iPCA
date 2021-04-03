# Integrated Principal Components Analysis (iPCA)
iPCA is a generalization of the classicial Principal Components Analysis to the data integration setting, where we observe multiple datasets. This repository provides the code associated with Tiffany M. Tang, Genevera I. Allen "Integrated Principal Components Analysis." arXiv preprint [arXiv:1810.00832](https://arxiv.org/abs/1810.00832) (2021+).

## Directory Structure

- **[data/](./data/)**: contains real genomics data for iPCA simulations
- **[functions/](./functions/)**: contains functions to perform iPCA and other existing data integration methods
- **[sims/](./sims/)**: contains files to run the simulations from the iPCA manuscript

## Usage

1. Download repository
2. Run example usage script: runiPCA.R

## Additional Notes

- As of 03/07/2020, the R package SpatioTemporal was removed from the CRAN repository. The R package QUIC has also since been removed as of 08/03/2020. To run functions that rely on these two packages (e.g., ``FFmleGlasso``, ``FFmleGlassoCorrelation``), we recommend install older versions of these packages from the archive:
```R
install.packages("https://cran.r-project.org/src/contrib/Archive/QUIC/QUIC_1.1.1.tar.gz", repos = NULL, type = "source”)
```


