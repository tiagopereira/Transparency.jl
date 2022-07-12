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
    zcoeff = convert(T, 6e-5)u"s^-1"  # should be m^2/s, dropping m to allow ne=ustrip(...)
    ne = ustrip(electron_density |> u"m^-3")
    # Factor of 2 to convert from half half-width to half-width
    return 2 * γβ * a1 * zcoeff * (n_upper^2 - n_lower^2) * ne^power
end


#=----------------------------------------------------------------------------
                    Quadratic Stark broadening: C_4 / r^4
----------------------------------------------------------------------------=#
"""
    c4_traving(χup, χlo, χ∞, Z)

Calculate the \$C_4\$ interaction constant (quadratic Stark effect) using the
recipe of Traving (1960), "Uber die Theorie der Druckverbreiterung von Spektrallinien",
p 93.

# Arguments
- χup, χlo, χ∞: energies of upper, lower, and ionisation
- Z: effective nuclear charge of the ionised level
"""
function c4_traving(χup, χlo, χ∞, Z)
    n_eff_u = n_eff(χ∞, χup, Z)
    n_eff_l = n_eff(χ∞, χlo, Z)
    C4 = (e^2 * inv_4πε0 * a_0^3 * 2 * π / (h * 18 * Z^4) *
        ((n_eff_u * (5 * n_eff_u^2 + 1))^2 - (n_eff_l * (5 * n_eff_l^2 + 1))^2))
    return C4 |> u"m^4 / s"
end


"""
    const_quadratic_stark(atomic_mass::Unitful.Mass, χup::Unitful.Energy,
                          χlo::Unitful.Energy, χ∞::Unitful.Energy, Z::Real;
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
function const_quadratic_stark(atomic_mass::Unitful.Mass, χup::Unitful.Energy,
                               χlo::Unitful.Energy, χ∞::Unitful.Energy, Z::Real;
                               mean_atomic_weight::Unitful.Mass=28 * m_u,
                               scaling::Real=1)
    C = ustrip(8 * k_B / (π * atomic_mass) |> u"J/(K * kg)")
    Cm = ((1 + atomic_mass / m_e)^(1/6) +
          (1 + atomic_mass / mean_atomic_weight)^(1/6))
    C4 = ustrip(c4_traving(χup, χlo, χ∞, Z) |> u"m^4 / s")
    cStark23 = 11.37u"m^3 / s" * (scaling * C4)^(2/3)
    return C^(1/6) * cStark23 * Cm
end

# Deprecated
function const_quadratic_stark(line::AtomicLine;
                               mean_atomic_weight::Unitful.Mass=28 * m_u,
                               scaling::Real=1)
    @warn "Calling const_quadratic_stark with AtomicLine is deprecated"
    const_quadratic_stark(line.atom_weight, line.χj, line.χi, line.χ∞, line.Z;
                          mean_atomic_weight, scaling)
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
    function const_unsold(atom_mass::Unitful.Mass, χup::Unitful.Energy, χlo::Unitful.Energy,
                          χ∞::Unitful.Energy, Z::Real; H_scaling=1, He_scaling=1)

Compute atmosphere-independent constant for γ_unsold, to be used in function `γ_unsold`.
Based on expressions from RH broad.c, which uses formula in Mihalas (1978),
pp 282, 286-287, eq. (9-50) for v_rel, table 9-1 and eq. (9-76) for the interaction
coefficient C6. Arguments are line parameters, where Z is the nuclear charge of the
upper level plus one (e.g. 1 for neutral, 2 for singly ionised).

The van der Waals broadening can be scaled for both H and He perturbers
using `H_scaling` and `He_scaling`.
"""
function const_unsold(atomic_mass::Unitful.Mass, χup::Unitful.Energy, χlo::Unitful.Energy,
                      χ∞::Unitful.Energy, Z::Real; H_scaling=1, He_scaling=1)
    Δr = (Ry^2 * (1 / (χ∞ - χup)^2 - 1 / (χ∞ - χlo)^2)) |> u"J/J"
    C6 = ustrip((2.5 * e^2 * αp * inv_4πε0^2 * 2 * π *
                 (Z * a_0)^2 / h * Δr) |> u"C^2 * m^6 / (F * J * s)")
    v_rel_const = ustrip(8 * k_B / (π * atomic_mass) |> u"J/(K * kg)")
    v_rel_H = v_rel_const * (1 + atomic_mass / mass_H)
    v_rel_He = v_rel_const * (1 + atomic_mass / mass_He)
    return 8.08 * (H_scaling * v_rel_H^0.3 + He_scaling * abund_He * v_rel_He^0.3) * C6^0.4
end

# For compatibility, now deprecated.
function const_unsold(line::AtomicLine; H_scaling::Real=1, He_scaling::Real=1)
    @warn "Calling const_unsold with AtomicLine is deprecated and will be removed soon"
    const_unsold(line.atom_weight, line.χj, line.χi, line.χ∞, line.Z;
                 H_scaling= H_scaling, He_scaling=He_scaling)
end


"""
    function γ_unsold(unsold_const::AbstractFloat, temperature::Unitful.Temperature,
                      h_neutral_density::NumberDensity)

Compute van der Waals broadening in Lindholm theory using Unsöld's approximation
for the interaction coefficient \$C_6\$. Based on Mihalas (1978), pp 282, 286-287.
Takes the atmosphere-indepenent `unsold_const` from `γ_unsold_const`, temperature,
and populations of neutral hydrogen, and returns broadening in units of s^-1.
"""
function γ_unsold(
    unsold_const::AbstractFloat,
    temperature::Unitful.Temperature,
    h_neutral_density::NumberDensity
)
    return (unsold_const * ustrip(temperature |> u"K")^0.3 *
            ustrip(h_neutral_density |> u"m^-3") * u"s^-1")
end


"""
    function const_barklem(atomic_weight::Unitful.Mass, α::Real, σ::Real)

Compute the atmosphere-independent constant used to calculate broadening from collisions
with neutral hydrogen using the recipes of Barklem/O'Mara/Anstee, in the function
`γ_barklem`. The calculation performed here follows eq (3) of
[Anstee & O'Mara (1995)](https://ui.adsabs.harvard.edu/abs/1995MNRAS.276..859A).

# Arguments
- `atomic_weight::Unitful.Mass`: atomic weight of element
- `α::Real`: velocity exponent from Barklem/O'Mara/Anstee tables
- `σ::Real`: line broadening cross section in atomic units (a_0^2) for a collision
   velocity of 10 km/s, from Barklem/O'Mara/Anstee tables.

# Returns
- `Unitful.VolumeFlow`: line broadening width per neutral hydrogen atom. Needs
   to be multiplied by temperature ^ ((1 - α)/2) to give proper temperature dependence.
"""
function const_barklem(atomic_mass::Unitful.Mass, α::Real, σ::Real)
    α < 0 && error("α must be non-negative")
    σ < 0 && error("σ must be non-negative")
    μ = m_u / (1 / Ar_H + 1 / (atomic_mass / m_u))
    # Using 1 K to keep units right for later multiplication by correct temperature
    v_bar = sqrt(8 * k_B * u"K"/ (π * μ)) |> u"m/s"
    v_ratio = (1e4u"m/s" / v_bar) |> u"m/m"
    # Squared Bohr radius is to convert from atomic units to m^2, factor of 2 from HW to FW
    return (a_0^2 * 2 * (4 / π)^(α / 2) * gamma((4 - α) / 2) * v_bar * σ *
            v_ratio^α) |> u"m^3 / s"
end


"""
    function γ_barklem(
        α::AbstractFloat,
        barklem_const::Unitful.VolumeFlow,
        temperature::Unitful.Temperature,
        h_neutral_density::NumberDensity,
    )

Compute van der Waals broadening from collisions with neutral hydrogen atoms following the
theory from Barklem/O'Mara/Anstee.

# Arguments
- `α::AbstractFloat`: velocity exponent from Barklem/O'Mara/Anstee tables
- `barklem_const::Unitful.VolumeFlow`: atmosphere-independent constant computed from
   `const_barklem()`.
- `temperature::Unitful.Temperature`
- `h_neutral_density::NumberDensity`: number density of neutral hydrogen atoms

# Returns
- `γ::Unitful.Frequency`: broadening in units of s^-1.

# Notes
This computes broadening only from hydrogen atoms. For collisions with helium atoms,
it is recommended to add van der Waals broadening using Unsöld's approximation (see example).

# Examples
Example for Ca II 854.2 nm line:
```
julia> Ca8542 = AtomicLine(25414.400u"cm^-1", 13710.880u"cm^-1", 95785.470u"cm^-1",
4, 6, 7.242e-02, 40.08 * m_u, 20);

julia> temp = 6000u"K";

julia> h_density = 1e23u"m^-3";

julia> bconst = const_barklem(Ca8542.atom_weight, 0.275, 291)
7.495208174533257e-16 m³ s⁻¹

julia> γ = γ_barklem(0.275, bconst, temp, h_density)
1.3596876505340942e11 s⁻¹
```

Now adding van der Waals broadening for helium as well:
```
julia> uconst = const_unsold(Ca8542; H_scaling=0, He_scaling=1)
2.482484115415461e-16

julia> γ = γ_barklem(0.275, bconst, temp, h_density) + γ_unsold(uconst, temp, h_density)
1.3630631018876224e11 s⁻¹
```

# References
- [Anstee & O'Mara (1995)](https://ui.adsabs.harvard.edu/abs/1995MNRAS.276..859A)
- [Barklem & O'Omara (1997)](https://ui.adsabs.harvard.edu/abs/1997MNRAS.290..102B)
- [Barklem, O'Mara & Ross (1998)](https://ui.adsabs.harvard.edu/abs/1998MNRAS.296.1057B)
- [Barklem & O'Mara (1998)](https://ui.adsabs.harvard.edu/abs/1998MNRAS.300..863B)
"""
function γ_barklem(
    α::AbstractFloat,
    barklem_const::Unitful.VolumeFlow,
    temperature::Unitful.Temperature,
    h_neutral_density::NumberDensity,
)
    return barklem_const * ustrip(temperature |> u"K")^((1 - α)/2) * h_neutral_density
end


"""
    function const_deridder_rensbergen(
        atomic_mass_pert::Unitful.Mass,
        atomic_mass_rad::Unitful.Mass
        α::Real,
        β::Real
    )

Compute the atmosphere-independent constant used to calculate broadening using
the recipes of [Deridder & Rensbergen (1976)](https://ui.adsabs.harvard.edu/abs/1976A%26AS...23..147D),
in the function `γ_deridder_rensbergen`.

# Arguments
- `atomic_weight_pert::Unitful.Mass`: atomic mass of perturbing element (hydrogen or helium)
- `atomic_weight_rad::Unitful.Mass`: atomic mass of line-producing species
- `α::Real`: α parameter as taken from the tables of Deridder & Rensbergen (1976),
   in units of 10^-8 cm^3/s (ignoring the dimensions of the temperature exponent)
- `β::Real`: β parameter as taken from the tables of Deridder & Rensbergen (1976),
   dimensionless.

# Returns
- `Unitful.VolumeFlow`: line broadening width per perturber atom. Needs
   to be multiplied by temperature ^ β to give proper temperature dependence.
"""
function const_deridder_rensbergen(
    atomic_mass_pert::Unitful.Mass,
    atomic_mass_rad::Unitful.Mass,
    α::Real,
    β::Real,
)
    α < 0 && error("α must be non-negative")
    α = (α * 1e-8u"cm^3/s") |> u"m^3/s"  # Convert from paper's 1e-8 units to SI
    mass_corr = (1 + atomic_mass_pert / atomic_mass_rad) ^ β
    return α * mass_corr
end


"""
    function γ_deridder_rensbergen(
        β::Real,
        deridder_const::Unitful.VolumeFlow,
        temperature::Unitful.Temperature,
        h_neutral_density::NumberDensity,
    )

Compute van der Waals broadening from collisions with neutral hydrogen or helium
atoms following [Deridder & Rensbergen (1976)](https://ui.adsabs.harvard.edu/abs/1976A%26AS...23..147D).

# Arguments
- `β::Real`: β parameter as taken from the tables of Deridder & Rensbergen (1976),
   dimensionless.
- `deridder_const::Unitful.VolumeFlow`: atmosphere-independent constant computed from
   `const_deridder_rensbergen()`.
- `temperature::Unitful.Temperature`
- `perturber_density::NumberDensity`: number density of perturber atoms, either
   neutral hydrogen or helium.

# Returns
- `γ::Unitful.Frequency`: broadening in units of s^-1.
"""
function γ_deridder_rensbergen(
    β::Real,
    deridder_const::Unitful.VolumeFlow,
    temperature::Unitful.Temperature,
    perturber_density::NumberDensity,
)
    return deridder_const * ustrip(temperature |> u"K")^β * perturber_density
end

#=----------------------------------------------------------------------------
                    Radiation utilities
----------------------------------------------------------------------------=#
"""
    function calc_Aji(λ0::Unitful.Length, g_ratio::Real, f_value::AbstractFloat)

Compute the spontaneous deexcitation rate \$A_{ul}\$ (natural broadening)
for a bound-bound transition, using the SI expression *per wavelength*:

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
