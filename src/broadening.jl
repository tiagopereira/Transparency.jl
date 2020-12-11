"""
Functions to calculate different line broadenings.
"""

const E_∞ = R_∞ * c_0 * h
const αp = 4.5 * 4*π * ε_0 * a_0^3   # Polarisability of hydrogen [F m^2]
const inv_4πε0 = 1. / (4 * π * ε_0)

# For now
const mass_H = 1.008 * m_u
const mass_He = 4.003 * m_u
const abund_He = 10^10.99 / 10^12  # From RH


"""
    function γ_unsold_const(line::AtomicLine; H_scaling::Real=1, He_scaling::Real=1)

Compute height-independent constant of γ_unsold, to be used in function `γ_unsold`.
Based on expressions from RH broad.c, which uses formula in Mihalas (1978), 
pp 282, 286-287, eq. (9-50) for v_rel, table 9-1 and eq. (9-76) for the interaction
coefficient C6.  Input is an `AtomicLine`, which contains the necessary energies and 
element weight. The van der Waals broadening can be scaled for both H and He perturbers 
using `H_scaling` and `He_scaling`.
"""
function γ_unsold_const(line::AtomicLine; H_scaling::Real=1, He_scaling::Real=1)
    Δr = (E_∞^2 * (1 / (line.χ∞ - line.χj)^2 - 1 / (line.χ∞ - line.χi)^2)) |> u"J/J"
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
function γ_unsold(unsold_const::AbstractFloat, temperature::Unitful.Temperature, 
                  h_ground_pop::NumberDensity)
    unsold_const * ustrip(temperature |> u"K")^0.3 * ustrip(h_ground_pop |> u"m^-3") * u"s^-1"
end


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
    (γ * λ^2 / (4 * π * c_0 * ΔλD)) |> u"m/m"
end