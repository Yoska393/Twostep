# Twostep

<img src="GA.jpg" width="800">

The source codes for the Simulation and application in the "Integration of Proxy Intermediate Omics traits into a Nonlinear Two-Step model for accurate phenotypic prediction" are available here.



doi: [https://doi.org/10.1007/s00122-026-05171-3](https://doi.org/10.1007/s00122-026-05171-3)



## Repository structure

### Functions

| File | Model / Purpose |
|------|-----------------|
| `MyFunctions.R` | Custom helper functions |

### Simulation

| File | Model / Purpose |
|------|-----------------|
| `1_bb.Rmd` | BLUP–BLUP |
| `2_rr.Rmd` | RF–RF |
| `3_rb.Rmd` | RF–BLUP |
| `4_br.Rmd` | BLUP–RF |
| `5_compare.Rmd` | Model comparison |

### Application

| File | Description |
|------|-------------|
| `1_metpred.Rmd` | Metabolome prediction |
| `2_phenopred_GMet.Rmd` | Phenotype prediction using genomic + metabolomic data |
| `3_phenopred_GMicroMet.Rmd` | Phenotype prediction using genomic + micro + metabolomic data |


---

## Reference

Yoshioka, H., Mary-Huard, T., Aubert, J. et al. Integration of proxy intermediate omics traits into a nonlinear two-step model for accurate phenotypic prediction. Theor Appl Genet 139, 86 (2026). https://doi.org/10.1007/s00122-026-05171-3


