# SNAIL JTWPA parameter design

A design study for a DC-flux-biased, three-wave-mixing SNAIL travelling-wave parametric amplifier. The Python notebook is the entry point. The Julia nonlinear circuit simulation is still under development.

## Status

The Python workflow evaluates stable SNAIL minima, Taylor coefficients, junction parameters, RF line impedance, exact discrete Bloch dispersion, cold phase mismatch and loop-inductance sensitivity. Its validation record is in [VALIDATION.md](VALIDATION.md).

The supplied Julia run stopped with `ArgumentError: No matching component found in allowedcomponents.` The notebook retains that limitation in a prominent status note, even though old outputs have been cleared. Pumped gain, phase matching, noise performance and layout readiness have **not** been validated. Earlier trial gain curves and failed solver experiments are not included as performance evidence.

## Files

- [Python parameter notebook](notebooks/SNAIL_JTWPA_parameter_design.ipynb)
- [Provisional JSON handoff](notebooks/snail_jtwpa_design_parameters.json)
- [Julia nonlinear workflow, incomplete](notebooks/SNAIL_JTWPA_JosephsonCircuits_simulation.ipynb)
- `figures/`: plots from the executed Python workflow

## Baseline candidate

| Parameter | Model value |
| --- | --- |
| Large junctions per SNAIL | 3 |
| Small/large critical-current ratio | 0.29 |
| Flux bias | approximately 0.4089 Φ₀ |
| Cell inductance | 145 pH |
| RF shunt capacitance | 58 fF |
| Nominal line impedance | 50 Ω |
| Pump frequency | 12 GHz |
| Signal band under study | 4–8 GHz |

These are starting values. Junction capacitance, loop inductance, bias coupling and process assumptions must be checked against actual fabrication and EM extraction. A near-zero bare quartic coefficient does not establish a Kerr-free pumped device. Cold phase mismatch is not the final nonlinear phase-matching criterion.

The revised baseline does not include a resonant phase-matching network by default. It differs from earlier periodically loaded trials; do not combine their results.

## Run the Python workflow

```shell
python -m pip install -r requirements.txt
jupyter lab
```

Open the Python notebook and run all cells. Its JSON handoff is written beside the notebook. See `VALIDATION.md` for the tested versions.

For the unfinished Julia workflow, activate this directory with `julia --project=.` and run `using Pkg; Pkg.instantiate()` before using an IJulia notebook environment. `Project.toml` declares direct dependencies; no `Manifest.toml` or verified solver-version lock is supplied. This does not resolve the component error.

## Validation still needed

1. Confirm a single SNAIL’s DC operating point and small-signal inductance.
2. Match unpumped 5–20-cell dispersion to the Python model.
3. Establish converged signal–idler conversion at `fi = fp - fs`.
4. Check nonlinear phase mismatch, pump depletion and phase excursions.
5. Scale cell count only after those checks pass.

[Project explanation](https://mukul-quantum-researcher.github.io/projects/snail-jtwpa/) · [Portfolio](https://mukul-quantum-researcher.github.io)
