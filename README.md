# CoCoA: Conditional Correlation Models with Association Size

This R package implements the conditional correlation model in the paper, ["CoCoA: Conditional Correlation Models with Association Size"](https://doi.org/10.1093/biostatistics/kxac032) by Tu *et al*, now in print at *Biostatistics*. It can be installed via

```
devtools::install_github("danni-tu/rCoCoA")
```

The preprint for this paper is also available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.03.28.486098v1).

## Conditional Correlation Models

In the CoCoA model, we consider the coupling between two outcomes $X$ and $Y$ as a function of covariates $Z$ and $T$. In particular, we assume a bivariate distribution for $(X,Y)$ given $Z, T$ with the conditional correlation parameter

$$\mathrm{Corr}(X,Y|Z,T) = g^{-1}(\alpha + \beta Z + \gamma T),$$

where $\alpha, \beta$, and $\gamma$ are parameters to be estimated and $g^{-1}(\cdot)$ is a link function which ensures that correlations are between -1 and 1. The `rCoCoA` package implements three estimators of the the conditional correlation parameters:

1. Maximum likelihood (`get_params_mle`)
2. Restricted maximum likelihood (`get_params_reml`)
3. Second-order generalized estimating equations (`get_params_gee`).
