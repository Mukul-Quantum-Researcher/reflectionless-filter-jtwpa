# JTWPA gain with and without reflectionless filtering

I simulated how external termination impedance changes the gain of a Josephson travelling-wave parametric amplifier. The three cases use the same amplifier: ideal 50-ohm loads, nonideal circulators, and a reflectionless-filter/circulator cascade.

The JTWPA circuit and parameters follow the official JosephsonCircuits example based on Macklin et al. The comparison is motivated by [Kedziora et al., Reflection-less filter for superconducting quantum circuits](https://arxiv.org/abs/2506.20645). The circulator and filter S-parameters are explicitly defined surrogate models; they are not measured Touchstone data from that paper.

## Implementation

- A 2,048-cell JTWPA with a 7.12 GHz, 1.85 microamp pump.
- Tuple-based circuit construction and harmonic-balance/small-signal analysis in JosephsonCircuits.jl.
- Frequency-dependent complex termination impedances derived from network reflection coefficients.
- A comparison of gain and the reflection presented to the amplifier.

## Files

- `notebooks/kedziora_jtwpa_fig1_josephsoncircuits.ipynb`: reviewed copy of the supplied notebook, preserving recorded plots and solver warnings. Large printed solver objects were removed for readability.
- `src/model.jl`: the same circuit and surrogate termination definitions, extracted into an executable model.
- `scripts/verify_terminations.py`: sampled checks of the surrogate-network impedances and filter scattering matrix.
- `figures/`: recorded output from the supplied notebook.
- `results/`: numerical exports and validation records.

## Reading the result

The recorded comparison illustrates the effect being studied: the filtered curve is closer to the ideal-load baseline than the circulator-only curve. However, the supplied circulator-only run reports non-convergence, so that original gain trace is provisional. No fresh Julia harmonic-balance run is claimed. The separate checks cover the surrogate networks only.

These are model comparisons, not measured amplifier performance or an exact numerical reproduction of the paper. Filter noise/quantum efficiency and gain measured after external insertion loss are not established here. Frequency-dependent port impedances also require care when comparing the normalisation of S-parameters across cases.

## Run

Use Julia 1.12 with the declared project dependencies and an IJulia environment to run the notebook. The supplied sweep has 1,301 signal frequencies (0.01 GHz spacing). Check convergence before using any quantitative gain result.

For the lightweight external-model check, run `python scripts/verify_terminations.py` with NumPy installed. This does not validate harmonic-balance gain or physical causality of the surrogate networks.

[Project page](https://mukul-quantum-researcher.github.io/projects/reflectionless-filter/) · [JosephsonCircuits documentation](https://josephsoncircuits.org/stable/)
