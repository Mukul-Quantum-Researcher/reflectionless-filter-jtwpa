using JosephsonCircuits

# JTWPA parameters
@variables Rleft Rright Cg Lj Cj Cc Cr Lr
@variables w

Nj = 2048
pmrpitch = 4

R0 = 50.0

Cg_value = 45.0e-15
Lj_value = IctoLj(3.4e-6)
Cj_value = 55.0e-15
Cc_value = 30.0e-15
Cr_value = 2.8153e-12
Lr_value = 1.70e-10

println("Lj = ", Lj_value*1e12, " pH")
println("Ic = 3.4 uA")


function build_jtwpa()
    circuit = Tuple{String,String,String,Num}[]

    # Input port and its reference impedance
    push!(circuit, ("P1_0", "1", "0", 1))
    push!(circuit, ("R1_0", "1", "0", Rleft))

    # First half shunt capacitance
    push!(circuit, ("C1_0", "1", "0", Cg/2))

    # First junction cell
    push!(circuit, ("Lj1_2", "1", "2", Lj))
    push!(circuit, ("C1_2", "1", "2", Cj))

    j = 2

    for i = 2:Nj-1
        if mod(i, pmrpitch) == pmrpitch ÷ 2
            # Modified junction cell at the phase-matching resonator position
            push!(circuit, ("C$(j)_0", "$(j)", "0", Cg - Cc))
            push!(circuit, ("Lj$(j)_$(j+2)", "$(j)", "$(j+2)", Lj))
            push!(circuit, ("C$(j)_$(j+2)", "$(j)", "$(j+2)", Cj))

            # Phase-matching resonator
            push!(circuit, ("C$(j)_$(j+1)", "$(j)", "$(j+1)", Cc))
            push!(circuit, ("C$(j+1)_0", "$(j+1)", "0", Cr))
            push!(circuit, ("L$(j+1)_0", "$(j+1)", "0", Lr))

            j += 1
        else
            push!(circuit, ("C$(j)_0", "$(j)", "0", Cg))
            push!(circuit, ("Lj$(j)_$(j+1)", "$(j)", "$(j+1)", Lj))
            push!(circuit, ("C$(j)_$(j+1)", "$(j)", "$(j+1)", Cj))
        end

        j += 1
    end

    # Last half shunt capacitance, output termination and port
    push!(circuit, ("C$(j)_0", "$(j)", "0", Cg/2))
    push!(circuit, ("R$(j)_0", "$(j)", "0", Rright))
    push!(circuit, ("P$(j)_0", "$(j)", "0", 2))

    return circuit
end

circuit_template = build_jtwpa()
println("Number of netlist elements: ", length(circuit_template))


# Helper: convert a one-port reflection coefficient into input impedance.
function gamma_to_z(gamma, z0)
    return z0 * (1 + gamma) / (1 - gamma)
end

# --- Surrogate circulator model ---
#
# This is NOT the measured circulator from the paper.
# It is deliberately simple: good match in-band and increasingly reflective
# outside the nominal 4-10 GHz operating band.

function gamma_circulator_positive(wrad)
    f = wrad / (2*pi)

    f_lo = 4.0e9
    f_hi = 10.0e9
    fc = sqrt(f_lo*f_hi)
    bw = f_hi - f_lo

    x = (f - fc)/(bw/2)

    # In-band return loss ≈ 20 dB; out-of-band reflection rises smoothly.
    window = abs(x) <= 1 ? 1.0 : exp(-((abs(x)-1)*2.2)^2)
    RL_dB = 20.0*window + 1.0*(1.0-window)

    # A modest electrical delay gives the reflection a frequency-dependent phase.
    phase_delay = 0.35e-9

    return 10.0^(-RL_dB/20) * exp(1im*2*pi*f*phase_delay)
end

function Z_circulator(wrad)
    if wrad < 0
        return conj(Z_circulator(-wrad))
    end
    gamma = gamma_circulator_positive(wrad)
    return gamma_to_z(gamma, R0)
end

@register_symbolic Z_circulator(w)


# --- Surrogate reflectionless-filter model ---
#
# The actual paper uses a measured/nonideal reflectionless filter.
# Its exact measured Touchstone data are not contained in the manuscript.
#
# We model the key behavior needed for Fig. 1:
#   * low reflection across the useful band
#   * low insertion loss in the passband
#   * strong attenuation in the out-of-band region
#
# The filter is represented here by its 2-port S parameters, then reduced
# to the input reflection coefficient of "filter + circulator".

function reflectionless_filter_S_positive(wrad)
    f = wrad / (2*pi)

    f_low = 3.2e9
    f_high = 11.8e9

    # Smooth passband envelope.
    n = 8.0
    low_edge = 1.0 / (1.0 + (f_low/max(f, 1.0))^n)
    high_edge = 1.0 / (1.0 + (f/max(f_high, 1.0))^n)
    transmission_envelope = sqrt(low_edge * high_edge)

    insertion_loss_dB = 0.75
    delay = 80e-12

    S21 = 10.0^(-insertion_loss_dB/20) *
          transmission_envelope *
          exp(-1im*2*pi*f*delay)

    # Reflectionless behavior: keep S11/S22 small over the full band.
    reflection_dB = -28.0
    S11 = 10.0^(reflection_dB/20) *
          exp(1im*(2*pi*f*35e-12 + 0.3))

    S22 = S11 * exp(1im*0.2)

    return ComplexF64[
        S11  S21;
        S21  S22
    ]
end

function reflectionless_filter_S(wrad)
    if wrad < 0
        return conj.(reflectionless_filter_S(-wrad))
    end
    return reflectionless_filter_S_positive(wrad)
end

# Input reflection of a 2-port network A terminated by a one-port load ΓL.
function input_gamma_from_S(SA, gammaL)
    return SA[1,1] +
           SA[1,2]*SA[2,1]*gammaL /
           (1 - SA[2,2]*gammaL)
end

# Configuration (c): reflectionless filter followed by the circulator.
function gamma_filter_plus_circulator_positive(wrad)
    Sf = reflectionless_filter_S_positive(wrad)
    gammaL = gamma_circulator_positive(wrad)
    return input_gamma_from_S(Sf, gammaL)
end

function Z_filter(wrad)
    if wrad < 0
        return conj(Z_filter(-wrad))
    end
    gamma = reflectionless_filter_S_positive(wrad)[1,1]
    return gamma_to_z(gamma, R0)
end

function Z_filter_plus_circulator(wrad)
    if wrad < 0
        return conj(Z_filter_plus_circulator(-wrad))
    end
    gamma = gamma_filter_plus_circulator_positive(wrad)
    return gamma_to_z(gamma, R0)
end

@register_symbolic Z_filter(w)
@register_symbolic Z_filter_plus_circulator(w)


circuit_a = build_jtwpa()
circuit_b = build_jtwpa()
circuit_c = build_jtwpa()

# Numeric component values for the common JTWPA.
common_defs = Dict(
    Lj => Lj_value,
    Cg => Cg_value,
    Cc => Cc_value,
    Cr => Cr_value,
    Lr => Lr_value,
    Cj => Cj_value,
)

defs_a = merge(common_defs, Dict(
    Rleft => R0,
    Rright => R0,
))

defs_b = merge(common_defs, Dict(
    Rleft => Z_circulator(w),
    Rright => Z_circulator(w),
))

defs_c = merge(common_defs, Dict(
    Rleft => Z_filter(w),
    Rright => Z_filter_plus_circulator(w),
))

println("Configuration (a) elements: ", length(circuit_a))
println("Configuration (b) elements: ", length(circuit_b))
println("Configuration (c) elements: ", length(circuit_c))

