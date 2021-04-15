"""
Functions to calculate different line broadenings.
"""

# For now
const mass_H = 1.008 * m_u
const mass_He = 4.003 * m_u
const abund_He = 10^10.99 / 10^12  # From RH


#=----------------------------------------------------------------------------
                    Linear Stark broadening: C_2 / r^2
----------------------------------------------------------------------------=#
"""
    γ_linear_stark(
        electron_density::NumberDensity{T},
        n_upper::Integer,
        n_lower::Integer
    ) where T <: AbstractFloat

Calculate linear Stark broadening according to the approximation of
[Sutton (1978)](https://ui.adsabs.harvard.edu/abs/1978JQSRT..20..333S/), eq (24),
so that it can be added into a Voigt profile and avoid the computation of a Holtsmark
profile. Valid up to electron densities of 1e19 m^-3 in the chromosphere.

# Arguments
- `electron_density`: electron density per volume
- `n_upper`: principal quantum number of upper level
- `n_lower`: principal quantum number of lower level
"""
function γ_linear_stark(
    electron_density::NumberDensity{T},
    n_upper::Integer,
    n_lower::Integer
) where T <: AbstractFloat
    @assert n_upper > n_lower
    @assert n_upper > 1
    if n_upper - n_lower == 1
        a1 = convert(T, 0.642)
    else
        a1 = convert(T, 1)
    end
    γβ = convert(T, 0.425)
    power = convert(T, 2/3)
    zcoeff = convert(T, 6e-5)u"m^2/s"
    # Factor of 2 to convert from half half-width to half-width
    return 2 * γβ * a1 * zcoeff * (n_upper^2 - n_lower^2) * electron_density^power
end


#=----------------------------------------------------------------------------
                    Quadratic Stark broadening: C_4 / r^4
----------------------------------------------------------------------------=#
"""
    c4_traving(line::AtomicLine)

Calculate the \$C_4\$ interaction constant (quadratic Stark effect) using the
recipe of Traving (1960), "Uber die Theorie der Druckverbreiterung von Spektrallinien",
p 93.
"""
function c4_traving(line::AtomicLine)
    n_eff_u = n_eff(line.χ∞, line.χj, line.Z)
    n_eff_l = n_eff(line.χ∞, line.χi, line.Z)
    C4 = (e^2 * inv_4πε0 * a_0^3 * 2 * π / (h * 18 * line.Z^4) *
        ((n_eff_u * (5 * n_eff_u^2 + 1))^2 - (n_eff_l * (5 * n_eff_l^2 + 1))^2))
    return C4 |> u"m^4 / s"
end


"""
    const_quadratic_stark(line::AtomicLine;
                          mean_atomic_weight::Unitful.Mass=28 * m_u,
                          scaling::Real=1)

Calculate height-independent constant to use in `γ_quadratic_stark`, using the recipe
from RH, which is based on the following estimate:

\$\$
\\gamma = 11.37 \\cdot vrel^{1/3} * C_4^{2/3} * (n_e + n_{ion}),
\$\$

Using the estimate for \$C_4\$ from Traving (1960), "Uber die Theorie der
Druckverbreiterung von Spektrallinien", p 93., and \$n_{ion}\\approx n_e\$
(following Gray).
"""
function const_quadratic_stark(line::AtomicLine;
                               mean_atomic_weight::Unitful.Mass=28 * m_u,
                               scaling::Real=1)
    C = ustrip(8 * k_B / (π * line.atom_weight) |> u"J/(K * kg)")
    Cm = ((1 + line.atom_weight / m_e)^(1/6) +
          (1 + line.atom_weight / mean_atomic_weight)^(1/6))
    C4 = ustrip(c4_traving(line) |> u"m^4 / s")
    cStark23 = 11.37u"m^3 / s" * (scaling * C4)^(2/3)
    return C^(1/6) * cStark23 * Cm
end


"""
    γ_quadratic_stark(
        electron_density::NumberDensity,
        temperature::Unitful.Temperature;
        stark_constant::Unitful.VolumeFlow=1.0u"m^3 / s",
    )

Compute quadratic Stark broadening for a given `electron_density`.  If `temperature` is
nonzero, then it will apply the standard recipe of RH (using \$C_4\$ from Traving 1960).
The `stark_constant` can be obtained either from atomic data sources, or, if using the RH
recipe, using the function `const_quadratic_stark`.
"""
function γ_quadratic_stark(
    electron_density::NumberDensity,
    temperature::Unitful.Temperature;
    stark_constant::Unitful.VolumeFlow=1.0u"m^3 / s",
)
    if temperature > 0u"K"
        t_factor = ustrip(temperature |> u"K")^(1/6)
    else
        t_factor = 1.0
    end
    return stark_constant * t_factor * electron_density
end


"""
    γ_quadratic_stark_gray(
        electron_density::NumberDensity,
        temperature::Unitful.Temperature,
        c4::Quantity{<:AbstractFloat, Unitful.𝐋^4 / Unitful.𝐓},
    )

Compute quadratic Stark broadening using the recipe of Gray (2005), page 244, eq 11.27.
The interaction constant `c4` should be provided, either from atomic data or from the
estimate of Traving (1960) using `c4_traving`.
"""
function γ_quadratic_stark_gray(
    electron_density::NumberDensity,
    temperature::Unitful.Temperature,
    c4::Quantity{<:AbstractFloat, Unitful.𝐋^4 / Unitful.𝐓},
)
    # Formula assumes CGS
    ne_term = ustrip(electron_density * k_B |> u"erg / (K * cm^3)")
    c4_term = ustrip(c4 |> u"cm^4 / s")
    t_term = ustrip(temperature |> u"K")
    log10γ = 19 + (2/3) * log10(c4_term) + log10(ne_term) + (1/6) * log10(t_term)
    return (10^log10γ) * u"s^-1"
end


#=----------------------------------------------------------------------------
                 van der Waals broadening, C_6 / r^6
----------------------------------------------------------------------------=#
"""
    function const_unsold(line::AtomicLine; H_scaling::Real=1, He_scaling::Real=1)

Compute height-independent constant for γ_unsold, to be used in function `γ_unsold`.
Based on expressions from RH broad.c, which uses formula in Mihalas (1978),
pp 282, 286-287, eq. (9-50) for v_rel, table 9-1 and eq. (9-76) for the interaction
coefficient C6.  Input is an `AtomicLine`, which contains the necessary energies and
element weight. The van der Waals broadening can be scaled for both H and He perturbers
using `H_scaling` and `He_scaling`.
"""
function const_unsold(line::AtomicLine; H_scaling::Real=1, He_scaling::Real=1)
    Δr = (Ry^2 * (1 / (line.χ∞ - line.χj)^2 - 1 / (line.χ∞ - line.χi)^2)) |> u"J/J"
    C6 = ustrip((2.5 * e^2 * αp * inv_4πε0^2 * 2 * π *
                 (line.Z * a_0)^2 / h * Δr) |> u"C^2 * m^6 / (F * J * s)")
    v_rel_const = ustrip(8 * k_B / (π * line.atom_weight) |> u"J/(K * kg)")
    v_rel_H = v_rel_const * (1 + line.atom_weight / mass_H)
    v_rel_He = v_rel_const * (1 + line.atom_weight / mass_He)
    return 8.08 * (H_scaling * v_rel_H^0.3 + He_scaling * abund_He * v_rel_He^0.3) * C6^0.4
end


"""
    function γ_unsold(unsold_const::AbstractFloat, temperature::Unitful.Temperature,
                      h_ground_pop::NumberDensity)

Compute van der Waals broadening in Lindholm theory using Unsöld's approximation
for the interaction coefficient \$C_6\$. Based on Mihalas (1978), pp 282, 286-287.
Takes the atmosphere-indepenent `unsold_const` from `γ_unsold_const`, temperature,
and ground-level populations of neutral hydrogen, and returns broadening in units of s^-1.
"""
function γ_unsold(
    unsold_const::AbstractFloat,
    temperature::Unitful.Temperature,
    h_ground_pop::NumberDensity
)
    return (unsold_const * ustrip(temperature |> u"K")^0.3 *
            ustrip(h_ground_pop |> u"m^-3") * u"s^-1")
end


#=----------------------------------------------------------------------------
                    Radiation utilities
----------------------------------------------------------------------------=#
"""
    function calc_Aji(λ0::Unitful.Length, g_ratio::Real, f_value::AbstractFloat)

Compute the spontaneous deexcitation rate \$A_{ul}\$ (natural broadening)
for a bound-bound transition, using the SI expression:

\$\$
A_{ul} = \\frac{2\\pi e^2}{\\varepsilon_0 m_e c} \\frac{g_l}{g_u} \\frac{f_{lu}}{\\lambda^2}
\$\$

for a given rest wavelength `λ0`, ration between statistical weights of lower and
upper levels (`g_ratio` = gl / gu), and `f_value` .
"""
function calc_Aji(λ0::Unitful.Length, g_ratio::Real, f_value::AbstractFloat)
    (2π * e^2 / (ε_0 * m_e * c_0) * g_ratio * f_value / λ0^2) |> u"s^-1"
end


"""
    function calc_Bji(λ0::Unitful.Length, Aji::Unitful.Frequency)

Compute the induced deexcitation rate \$B_{ul}\$ for a bound-bound transition,
using the SI expression *per wavelength*:

\$\$
B_{ul} = \\frac{\\lambda^5}{2 h c} A_{ul}
\$\$

for a given rest wavelength `λ0`, and spontaneous deexcitation rate `Aji.`
"""
calc_Bji(λ0::Unitful.Length, Aji::Unitful.Frequency) = (λ0^5 * Aji / (2h * c_0^2)) |> u"m^3 / J"


"""
Compute damping parameter.
"""
function damping(γ::Unitful.Frequency, λ::Unitful.Length, ΔλD::Unitful.Length)
    return (γ * λ^2 / (4 * π * c_0 * ΔλD)) |> u"m/m"
end
