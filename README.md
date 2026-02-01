# ml-GAR

Multilevel Growth with Autoregressive Residuals — a Bayesian approach for modeling trends and intraindividual variability simultaneously in longitudinal panel data.

## What is this

Longitudinal data often contain two types of change: gradual developmental trends (e.g., growth in reading ability) and short-term state fluctuations (e.g., daily mood swings). Traditional two-stage approaches detrend first, then model the residuals — but this discards sampling uncertainty and performs poorly when time points are limited.

ml-GAR handles both in a single Bayesian framework. It uses a Gompertz curve for nonlinear growth trends and a multilevel AR(1) process for fluctuations around the trend. Simulations show that ml-GAR recovers parameters well even with as few as 5 time points per person, where two-stage methods break down.

For methodological details and simulation studies, see:

> Xiong, X., Li, Y., Hunter, M. D., & Chow, S.-M. (2025). Integrated trend and lagged modeling of multi-subject, multilevel, and short time series. *Multivariate Behavioral Research*. https://doi.org/10.1080/00273171.2025.2587286

## Dependencies

The model runs on JAGS:

- R (≥ 4.0)
- [JAGS](https://sourceforge.net/projects/mcmc-jags/) (≥ 4.3)
- R packages: `rjags`, `coda`, `ggplot2`, `MASS`

```r
install.packages(c("rjags", "coda", "ggplot2", "MASS"))
```

## Quick start

The tutorial `vignettes/mlGAR_Tutorial.Rmd` walks through data generation, model specification, MCMC sampling, and posterior diagnostics. Knit it to run the full workflow.

Basic steps:

1. Prepare data as a 3D array: `[time, variable, person]`.
2. Write the JAGS model. The core structure is a Level-1 observation model (Gompertz trend + AR state) and Level-2 random effect distributions.
3. Run MCMC. The tutorial uses 2 chains, 5000 adaptation, 10000 burn-in, 10000 sampling iterations with thinning = 2. Real data may need tuning.
4. Check convergence via RHAT (< 1.1) and ESS (higher is better).

## Model structure

Observation model (Equation 1 in the paper):

$$y_{t,i} = \tau_{t,i} + \eta_{t,i}$$

where $\tau_{t,i}$ is the Gompertz trend and $\eta_{t,i}$ is the AR(1) state fluctuation.

Trend component (Equation 2):

$$\tau_{t,i} = \theta_{1,i} \exp\left[-\theta_{2,i} \exp(-t \cdot \theta_{3,i})\right]$$

The three Gompertz parameters are asymptote ($\theta_1$), displacement ($\theta_2$), and growth rate ($\theta_3$).

State component:

$$\eta_{t,i} = \phi_i \cdot \eta_{t-1,i} + u_{t,i}, \quad u_{t,i} \sim N(0, \sigma^2_u)$$

All person-level parameters follow normal distributions at Level 2.

**Notation correspondence (code → paper):**

| Code | Paper |
|------|-------|
| `MU[id,t]` | $\tau_{t,i}$ |
| `MU_thetas[id,k]` | $\theta_{k,i}$ |
| `AR[id]` | $\phi_i$ |
| `IIV` | $\sigma^2_u$ |

## Priors

The tutorial uses weakly informative priors (see Table 3 in the paper). Tighten them if you have prior knowledge about parameter ranges. The AR parameter is truncated to (-1, 1) for stationarity.

## Output

`zcalc()` or `summarizePost()` returns posterior mean, median, mode, SD, 95% HDI, ESS, and RHAT for each parameter. Person-level parameters (`MU_thetas[id,k]`, `AR[id]`) are output by individual.

## Troubleshooting

If you run into convergence issues, try: widening priors, increasing adaptation iterations, changing initial values, or checking for outliers in your data.

## Citation

```bibtex
@article{xiong2025mlgar,
  title={Integrated trend and lagged modeling of multi-subject, multilevel, and short time series},
  author={Xiong, Xiaoyue and Li, Yanling and Hunter, Michael D. and Chow, Sy-Miin},
  journal={Multivariate Behavioral Research},
  year={2025},
  doi={10.1080/00273171.2025.2587286}
}
```

## License

MIT
