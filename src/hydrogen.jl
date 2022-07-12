"""
Computes the extinction for hydrogen-related bound-free and free-free transitions.
Functions are organised in two types:

1. `σ_*` functions compute the cross section (in m^2) or "cross section
   coefficient" (in m^5)
2. `α_*` functions compute the linear extinction coefficient (in m^-1),
   by multiplying

extinction coefficient \$\\alpha_\\nu\$ for hydrogen-related
bound-free and free-free transitions.

Includes:
* Neutral Hydrogen bound-free and free-free.
* H\$^-\$ bound-free and free-free.
* H\$_2^-\$ free-free.
* H\$_2^+\$ free-free.
* Rayleigh scattering by molecular H\$_2\$
"""

#=----------------------------------------------------------------------------
                            Catch-all functions
----------------------------------------------------------------------------=#
"""
    σ_hminus_ff(
        λ::Unitful.Length,
        temperature::Unitful.Temperature,
        recipe::String="stilley"
    )

Compute free-free extinction from H minus ion. Recipe can be one of:

- `stilley` (default): Interpolates table from
  [Stilley & Callaway (1970)](https://ui.adsabs.harvard.edu/abs/1970ApJ...160..245S/abstract),
  which is valid for λ up to 9113 nm.

- `john`: Follows
  [John (1988)](https://ui.adsabs.harvard.edu/abs/1988A%26A...193..189J/abstract),
  which is valid beyond 9113 nm but may not be good below 364.5 nm.
"""
function σ_hminus_ff(λ::Unitful.Length, temperature::Unitful.Temperature;
                     recipe::String="stilley"
)
    if recipe == "stilley"
        σ = σ_hminus_ff_stilley(λ, temperature)
    elseif recipe == "john"
        σ = σ_hminus_ff_john(λ, temperature)
    else
        throw("NotImplemented recipe $recipe")
    end
    return σ
end

"""
    α_hminus_ff(
        λ::Unitful.Length,
        temperature::Unitful.Temperature,
        h_neutral_density::NumberDensity,
        electron_density::NumberDensity;
        recipe::String="stilley"
    )

Compute free-free extinction from H minus ion. Recipe can be one of:

- `stilley` (default): Interpolates table from
  [Stilley & Callaway (1970)](https://ui.adsabs.harvard.edu/abs/1970ApJ...160..245S/abstract),
  which is valid for λ up to 9113 nm.

- `john`: Follows
  [John (1988)](https://ui.adsabs.harvard.edu/abs/1988A%26A...193..189J/abstract),
  which is valid beyond 9113 nm but may not be good below 364.5 nm.
"""
function α_hminus_ff(
    λ::Unitful.Length,
    temperature::Unitful.Temperature,
    h_neutral_density::NumberDensity,
    electron_density::NumberDensity;
    recipe::String="stilley"
)
    σ = σ_hminus_ff(λ, temperature; recipe=recipe)
    return σ * h_neutral_density * electron_density
end

"""
    σ_hminus_bf(
        λ::Unitful.Length,
        temperature::Unitful.Temperature,
        h_neutral_density::NumberDensity,
        electron_density::NumberDensity;
        recipe::String="wbr"
    )

Compute bound-free extinction from H minus ion. Recipe can be one of:

- `geltman` (default): Uses recipe from
  [Geltman (1962)](https://ui.adsabs.harvard.edu/abs/1962ApJ...136..935G/abstract)

- `john`: Follows
  [John (1988)](https://ui.adsabs.harvard.edu/abs/1988A%26A...193..189J/abstract),
  which is valid beyond 9113 nm but may not be good below 364.5 nm.

- `wbr`: Follows
  [Wishart (1979)](https://ui.adsabs.harvard.edu/abs/1979MNRAS.187P..59W) for λ > 175 nm,
  and [Broad and Reinhardt (1976)](https://ui.adsabs.harvard.edu/abs/1976PhRvA..14.2159B)
  for λ <= 164 nm.
"""
function σ_hminus_bf(λ::Unitful.Length, temperature::Unitful.Temperature;
                     recipe::String="wbr"
)
    if recipe == "geltman"
        σ = σ_hminus_bf_geltman(λ, temperature)
    elseif recipe == "john"
        σ = σ_hminus_bf_john(λ, temperature)
    elseif recipe == "wbr"
        σ = σ_hminus_bf_wbr(λ, temperature)
    else
        throw("NotImplemented recipe $recipe")
    end
    return σ
end

"""
    α_hminus_bf(
        λ::Unitful.Length,
        temperature::Unitful.Temperature,
        h_neutral_density::NumberDensity,
        electron_density::NumberDensity;
        recipe::String="wbr"
    )

Compute bound-free extinction from H minus ion. Recipe can be one of:

- `geltman`: Uses recipe from
  [Geltman (1962)](https://ui.adsabs.harvard.edu/abs/1962ApJ...136..935G/abstract)

- `john`: Follows
  [John (1988)](https://ui.adsabs.harvard.edu/abs/1988A%26A...193..189J/abstract),
  which is valid beyond 9113 nm but may not be good below 364.5 nm.

- `wbr` (default): Follows
  [Wishart (1979)](https://ui.adsabs.harvard.edu/abs/1979MNRAS.187P..59W) for λ > 175 nm,
  and [Broad and Reinhardt (1976)](https://ui.adsabs.harvard.edu/abs/1976PhRvA..14.2159B)
  for λ <= 164 nm.
"""
function α_hminus_bf(
    λ::Unitful.Length,
    temperature::Unitful.Temperature,
    h_neutral_density::NumberDensity,
    electron_density::NumberDensity;
    recipe::String="wbr"
)
    σ = σ_hminus_bf(λ, temperature; recipe=recipe)
    return σ * h_neutral_density * electron_density
end


#=----------------------------------------------------------------------------
                   Recipes from Karzas and Latter / Kurucz
----------------------------------------------------------------------------=#
#= Tabulated Gaunt factors from Kurucz (1970, SAO Special Report no. 309), page 77,
 which is a fit to  figures 3-5 of Karzas and Latter (1961, ApJ Suppl 6, 167)

- Note from Mats Carlsson in RH:
    There is extrapolation outside the range of lg(gamma2) (T outside [1570,1.57e8K])
    or outside range of lg(U) (for lambda=3.2mm for T>49 kK). This should be OK for
    reasonable extrapolation distances (and better than setting to constant end value
    or zero). Interpolation tested against table in Gustafsson (1973) with results
    within 1%.
=#
const kurucz_ff_table = [5.53  5.49  5.46  5.43  5.40  5.25  5.00  4.69  4.48  4.16  3.85
                         4.91  4.87  4.84  4.80  4.77  4.63  4.40  4.13  3.87  3.52  3.27
                         4.29  4.25  4.22  4.18  4.15  4.02  3.80  3.57  3.27  2.98  2.70
                         3.64  3.61  3.59  3.56  3.54  3.41  3.22  2.97  2.70  2.45  2.20
                         3.00  2.98  2.97  2.95  2.94  2.81  2.65  2.44  2.21  2.01  1.81
                         2.41  2.41  2.41  2.41  2.41  2.32  2.19  2.02  1.84  1.67  1.50
                         1.87  1.89  1.91  1.93  1.95  1.90  1.80  1.68  1.52  1.41  1.30
                         1.33  1.39  1.44  1.49  1.55  1.56  1.51  1.42  1.33  1.25  1.17
                         0.90  0.95  1.00  1.08  1.17  1.30  1.32  1.30  1.20  1.15  1.11
                         0.55  0.58  0.62  0.70  0.85  1.01  1.15  1.18  1.15  1.11  1.08
                         0.33  0.36  0.39  0.46  0.59  0.76  0.97  1.09  1.13  1.10  1.08
                         0.19  0.21  0.24  0.28  0.38  0.53  0.76  0.96  1.08  1.09  1.09]
# table_y is log10(hν / kT)
const kurucz_ff_table_y = -4:0.5:1.5
# table_x is log10(3.28805e15 * (Z^2 h) / (kT))
const kurucz_ff_table_x = -3:0.5:2.0
const h_k = h / k_B
const hc_k = h * c_0 / k_B
const hminusχ = 0.754u"eV"
const saha_const = h^2 / (2 * π * m_e * k_B)  # Constant for Saha equation
const αff_const = (4 / (3 * h * c_0) * (e^2 / (4 * π * ε_0))^3 *
                   sqrt(2 * π / (3 * m_e^3 * k_B))) |> u"K^(1/2) * m^5 / s^3"
const αbf_const = (4 * e^2 / (3 * π * sqrt(3) * ε_0 * m_e * c_0^2 * R_∞)) |> u"m^2"
const kurucz_ff_interp = LinearInterpolation((kurucz_ff_table_y, kurucz_ff_table_x),
                                             kurucz_ff_table, extrapolation_bc=Line())

"""
    gaunt_ff(ν::Unitful.Frequency, temperature::Unitful.Temperature, charge::Int)
    gaunt_ff(λ::Unitful.Length, temperature::Unitful.Temperature, charge::Int)

Compute Gaunt factor for free-free based on [Karzas and Latter (1961, ApJ Suppl 6, 167)]
(https://ui.adsabs.harvard.edu/abs/1961ApJS....6..167K/abstract)
fit in
[Kurucz (1970, SAO Special Report no. 309), page 77](https://ui.adsabs.harvard.edu/abs/1970SAOSR.309.....K/abstract)
"""
function gaunt_ff(ν::Unitful.Frequency, temperature::Unitful.Temperature, charge::Int)
    lookup_y = log10(h_k * ν / temperature)
    lookup_x = log10(3.28805e15u"Hz" * charge^2 * h_k / temperature)
    return kurucz_ff_interp(lookup_y, lookup_x)::Float64
end

function gaunt_ff(λ::Unitful.Length, temperature::Unitful.Temperature, charge::Int)
    return gaunt_ff(c_0 / λ, temperature, charge)
end


#=----------------------------------------------------------------------------
                            Recipes from Seaton
----------------------------------------------------------------------------=#
"""
    gaunt_bf(charge::Int, n_eff::Number, λ::Unitful.Length)::Float64

Compute bound-free Gaunt factor for a given nuclear charge Z, effective principal
quantum number and wavelength λ. Taken from RH. Formula from
[Seaton (1960), Rep. Prog. Phys. 23, 313](https://ui.adsabs.harvard.edu/abs/1960RPPh...23..313S/abstract),
page 316.
"""
function gaunt_bf(λ::Unitful.Length, Z::Real, n_eff::Real)::Float64
    x = ustrip(1 / (λ * R_∞ * Z^2) |> u"m/m")
    x3 = x^(1/3)
    nsqx = 1 / (n_eff^2 * x)
    g_bf = 1 + 0.1728 * x3 * (1 - 2 * nsqx) - 0.0496 * x3^2 * (1 - (1 - nsqx) * 0.66666667 * nsqx)
    @assert g_bf >= 0 "gaunt_bf negative, calculation will not be reliable"
    return g_bf
end

function gaunt_bf(ν::Unitful.Frequency, Z::Real, n_eff::Real)::Float64
    return gaunt_bf(c_0 / ν, Z, n_eff)
end


"""
    n_eff(energy_upper::Unitful.Energy, energy_lower::Unitful.Energy, Z::Integer)

Compute the effective principal quantum number for a given energy difference
and nuclear charge Z.
"""
function n_eff(energy_upper::Unitful.Energy, energy_lower::Unitful.Energy, Z::Real)
    return Z * sqrt(Ryh / (energy_upper - energy_lower))
end

#=----------------------------------------------------------------------------
                            Recipes from Stilley
----------------------------------------------------------------------------=#
const stilley_ff_λ = [0.0, 303.8, 455.6, 506.3, 569.5, 650.9, 759.4, 911.3,
                      1013.0, 1139.0, 1302.0, 1519.0, 1823.0, 2278.0, 3038.0,
                      4556.0, 9113.0]  # in nm
const stilley_ff_t = 5040.0 ./ [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
                                1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]  # in K
const stilley_ff_table =
   [[0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00  #=    0.0 nm
=#   0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00];
    [3.44e-02 4.18e-02 4.91e-02 5.65e-02 6.39e-02 7.13e-02 7.87e-02 8.62e-02  #=  303.8 nm
=#   9.36e-02 1.01e-01 1.08e-01 1.16e-01 1.23e-01 1.30e-01 1.38e-01 1.45e-01];
    [7.80e-02 9.41e-02 1.10e-01 1.25e-01 1.40e-01 1.56e-01 1.71e-01 1.86e-01  #=  455.6 nm
=#   2.01e-01 2.16e-01 2.31e-01 2.45e-01 2.60e-01 2.75e-01 2.89e-01 3.03e-01];
    [9.59e-02 1.16e-01 1.35e-01 1.53e-01 1.72e-01 1.90e-01 2.08e-01 2.25e-01  #=  506.3 nm
=#   2.43e-01 2.61e-01 2.78e-01 2.96e-01 3.13e-01 3.30e-01 3.47e-01 3.64e-01];
    [1.21e-01 1.45e-01 1.69e-01 1.92e-01 2.14e-01 2.36e-01 2.58e-01 2.80e-01  #=  569.5 nm
=#   3.01e-01 3.22e-01 3.43e-01 3.64e-01 3.85e-01 4.06e-01 4.26e-01 4.46e-01];
    [1.56e-01 1.88e-01 2.18e-01 2.47e-01 2.76e-01 3.03e-01 3.31e-01 3.57e-01  #=  650.9 nm
=#   3.84e-01 4.10e-01 4.36e-01 4.62e-01 4.87e-01 5.12e-01 5.37e-01 5.62e-01];
    [2.10e-01 2.53e-01 2.93e-01 3.32e-01 3.69e-01 4.06e-01 4.41e-01 4.75e-01  #=  759.4 nm
=#   5.09e-01 5.43e-01 5.76e-01 6.08e-01 6.40e-01 6.72e-01 7.03e-01 7.34e-01];
    [2.98e-01 3.59e-01 4.16e-01 4.70e-01 5.22e-01 5.73e-01 6.21e-01 6.68e-01  #=  911.3 nm
=#   7.15e-01 7.60e-01 8.04e-01 8.47e-01 8.90e-01 9.32e-01 9.73e-01 1.01e+00];
    [3.65e-01 4.39e-01 5.09e-01 5.75e-01 6.39e-01 7.00e-01 7.58e-01 8.15e-01  #= 1013.0 nm
=#   8.71e-01 9.25e-01 9.77e-01 1.03e+00 1.08e+00 1.13e+00 1.18e+00 1.23e+00];
    [4.58e-01 5.50e-01 6.37e-01 7.21e-01 8.00e-01 8.76e-01 9.49e-01 1.02e+00  #= 1139.0 nm
=#   1.09e+00 1.15e+00 1.22e+00 1.28e+00 1.34e+00 1.40e+00 1.46e+00 1.52e+00];
    [5.92e-01 7.11e-01 8.24e-01 9.31e-01 1.03e+00 1.13e+00 1.23e+00 1.32e+00  #= 1302.0 nm
=#   1.40e+00 1.49e+00 1.57e+00 1.65e+00 1.73e+00 1.80e+00 1.88e+00 1.95e+00];
    [7.98e-01 9.58e-01 1.11e+00 1.25e+00 1.39e+00 1.52e+00 1.65e+00 1.77e+00  #= 1519.0 nm
=#   1.89e+00 2.00e+00 2.11e+00 2.21e+00 2.32e+00 2.42e+00 2.51e+00 2.61e+00];
    [1.14e+00 1.36e+00 1.58e+00 1.78e+00 1.98e+00 2.17e+00 2.34e+00 2.52e+00  #= 1823.0 nm
=#   2.68e+00 2.84e+00 3.00e+00 3.15e+00 3.29e+00 3.43e+00 3.57e+00 3.70e+00];
    [1.77e+00 2.11e+00 2.44e+00 2.75e+00 3.05e+00 3.34e+00 3.62e+00 3.89e+00  #= 2278.0 nm
=#   4.14e+00 4.39e+00 4.63e+00 4.86e+00 5.08e+00 5.30e+00 5.51e+00 5.71e+00];
    [3.10e+00 3.71e+00 4.29e+00 4.84e+00 5.37e+00 5.87e+00 6.36e+00 6.83e+00  #= 3038.0 nm
=#   7.28e+00 7.72e+00 8.14e+00 8.55e+00 8.95e+00 9.33e+00 9.71e+00 1.01e+01];
    [6.92e+00 8.27e+00 9.56e+00 1.08e+01 1.19e+01 1.31e+01 1.42e+01 1.52e+01  #= 4556.0 nm
=#   1.62e+01 1.72e+01 1.82e+01 1.91e+01 2.00e+01 2.09e+01 2.17e+01 2.25e+01];
    [2.75e+01 3.29e+01 3.80e+01 4.28e+01 4.75e+01 5.19e+01 5.62e+01 6.04e+01  #= 9133.0 nm
=#   6.45e+01 6.84e+01 7.23e+01 7.60e+01 7.97e+01 8.32e+01 8.67e+01 9.01e+01]]
const stilley_ff_interp = LinearInterpolation((stilley_ff_λ, stilley_ff_t[end:-1:1]),
                             stilley_ff_table[:, end:-1:1], extrapolation_bc=Line())

"""
    σ_hminus_ff_stilley(λ::Unitful.Length, temperature::Unitful.Temperature)

Compute free-free cross section coefficient from H minus ion, for a given wavelength
and temperature. Units are m^5, needs to be multiplied by electron density and
density of neutral hydrogen atoms to obtain linear extinction. Interpolates table from
[Stilley & Callaway (1970)](https://ui.adsabs.harvard.edu/abs/1970ApJ...160..245S/abstract),
page 255, which is valid for λ up to 9113 nm.
"""
function σ_hminus_ff_stilley(λ::Unitful.Length, temperature::Unitful.Temperature)
    λi = ustrip(λ |> u"nm")   # convert to units of table
    temp = ustrip(temperature |> u"K")
    kappa = stilley_ff_interp(λi, temp)::Float64 * 1e-29u"m^4/N"
    return  k_B * temperature * kappa |> u"m^5"
end


"""
    α_hminus_ff_stilley(
        λ::Unitful.Length,
        temperature::Unitful.Temperature,
        h_neutral_density::NumberDensity,
        electron_density::NumberDensity
    )

Compute free-free extinction from H minus ion, for a given wavelength, temperature,
density of neutral hydrogen atoms `h_neutral_density`, and electron density.
Based on `σ_hminus_ff_stilley`.
"""
function α_hminus_ff_stilley(
    λ::Unitful.Length,
    temperature::Unitful.Temperature,
    h_neutral_density::NumberDensity,
    electron_density::NumberDensity
)
    σ = σ_hminus_ff_stilley(λ, temperature)
    return σ * h_neutral_density * electron_density |> u"m^-1"
end

#=----------------------------------------------------------------------------
                            Recipes from Geltman
----------------------------------------------------------------------------=#
const geltman_bf_λ = [   0.0,   50.0,  100.0,  150.0,  200.0,  250.0,  300.0,
                       350.0,  400.0,  450.0,  500.0,  550.0,  600.0,  650.0,
                       700.0,  750.0,  800.0,  850.0,  900.0,  950.0, 1000.0,
                      1050.0, 1100.0, 1150.0, 1200.0, 1250.0, 1300.0, 1350.0,
                      1400.0, 1450.0, 1500.0, 1550.0, 1600.0, 1641.9]  # in nm
const geltman_bf_σ = [0.00, 0.15, 0.33, 0.57, 0.85, 1.17, 1.52, 1.89, 2.23,
                      2.55, 2.84, 3.11, 3.35, 3.56, 3.71, 3.83, 3.92, 3.95,
                      3.93, 3.85, 3.73, 3.58, 3.38, 3.14, 2.85, 2.54, 2.20,
                      1.83, 1.46, 1.06, 0.71, 0.40, 0.17, 0.0]  # in 1e-21 m^2

const geltman_bf_interp = LinearInterpolation(geltman_bf_λ, geltman_bf_σ, extrapolation_bc=0)


"""
    σ_hminus_bf_geltman(λ::Unitful.Length, temperature::Unitful.Temperature)

Compute bound-free cross section from H minus ion. Uses recipe from
[Geltman (1962)](https://ui.adsabs.harvard.edu/abs/1962ApJ...136..935G/abstract)
Units are m^5, needs to be multiplied by density of neutral H atoms and electron
density to obtain linear extinction.
"""
function σ_hminus_bf_geltman(λ::Unitful.Length, temperature::Unitful.Temperature)
    λi = ustrip(λ |> u"nm")   # convert to units of table
    stimulated_emission = exp(-hc_k / (λ * temperature))
    # Get H- fraction to convert from σ per H- atom to σ per H atom per electron
    hminus_frac = calc_hminus_density(1.0u"m^-3", temperature, 1.0u"m^-3") * u"m^6"
    σ = geltman_bf_interp(λi)::Float64 * 1e-21u"m^2" *
        (1 - stimulated_emission) * hminus_frac
    return σ
end


"""
    α_hminus_bf_geltman(
        λ::Unitful.Length,
        temperature::Unitful.Temperature,
        h_minus_density::NumberDensity
    )
    α_hminus_bf_geltman(
        λ::Unitful.Length,
        temperature::Unitful.Temperature,
        h_neutral_density::NumberDensity,
        electron_density::NumberDensity
    )

Compute extinction from H minus ion, from input H minus populations.
Uses recipe from
[Geltman (1962)](https://ui.adsabs.harvard.edu/abs/1962ApJ...136..935G/abstract)
"""

function α_hminus_bf_geltman(
    λ::Unitful.Length,
    temperature::Unitful.Temperature,
    h_minus_density::NumberDensity
)
    σ = σ_hminus_bf_geltman(λ, temperature)
    λi = ustrip(λ |> u"nm")   # convert to units of table
    stimulated_emission = exp(-hc_k / (λ * temperature))
    σ = geltman_bf_interp(λi)::Float64 * 1e-21u"m^2" * (1 - stimulated_emission)
    return σ * h_minus_density
end


function α_hminus_bf_geltman(
    λ::Unitful.Length,
    temperature::Unitful.Temperature,
    h_neutral_density::NumberDensity,
    electron_density::NumberDensity
)
    h_minus_density = calc_hminus_density(h_neutral_density, temperature, electron_density)
    return α_hminus_bf_geltman(λ, temperature, h_minus_density)
end


"""
    calc_hminus_density(
        h_neutral_density::NumberDensity,
        temperature::Unitful.Temperature,
        electron_density::NumberDensity
    )

Compute H minus populations based on electron density, temperature, and
density of neutral hydrogen atoms.
"""
function calc_hminus_density(
    h_neutral_density::NumberDensity,
    temperature::Unitful.Temperature,
    electron_density::NumberDensity
)
    tmp = (ustrip(saha_const) / ustrip(temperature |> u"K"))^(3/2) * u"m^3"
    ϕ = tmp * exp(hminusχ / (k_B * temperature)) / 4
    return h_neutral_density * electron_density * ϕ
end


#=----------------------------------------------------------------------------
                            Recipes from John
----------------------------------------------------------------------------=#
"""
    σ_hminus_ff_john(λ::Unitful.Length, temperature::Unitful.Temperature)

Compute free-free cross section coefficient from H minus ion, for a given wavelength
and temperature. Units are m^5, needs to be multiplied by electron density and
density of neutral hydrogen atoms to obtain linear extinction. Uses recipe from
[John (1988)](https://ui.adsabs.harvard.edu/abs/1988A%26A...193..189J/abstract),
which is valid beyond 9113 nm but may not be good below 364.5 nm.

Includes stimulated emission.
"""
function σ_hminus_ff_john(λ::Unitful.Length, temperature::Unitful.Temperature)
    λμ = ustrip(λ |> u"μm")
    if λμ > 0.3645
        table =
          SA[    0.0000     0.0000      0.0000      0.0000     0.0000    0.0000
              2483.3460   285.8270  -2054.2910   2827.7760 -1341.5370  208.9520
             -3449.8890 -1158.3820   8746.5230 -11485.6320  5303.6090 -812.9390
              2200.0400  2427.7190 -13651.1050  16755.5240 -7510.4940 1132.7380
              -696.2710 -1841.4000   8624.9700 -10051.5300  4400.0670 -655.0200
                88.2830   444.5170  -1863.8640   2095.2880  -901.7880  132.9850]
    else
        table =
          SA[  518.1021  -734.8666   1021.1775   -479.0721    93.1373   -6.4285
               473.2636  1443.4137  -1977.3395    922.3575  -178.9275   12.3600
              -482.2089  -737.1616   1096.8827   -521.1341   101.7963   -7.0571
               115.5291   169.6374   -245.6490    114.2430   -21.9972    1.5097
                 0.0000     0.0000      0.0000      0.0000     0.0000    0.0000
                 0.0000     0.0000      0.0000      0.0000     0.0000    0.0000]
    end
    sqrtθ = sqrt((5040.0u"K" / temperature) |> u"K/K")
    λinv = 1.0 / λμ
    κ = 0.0
    for i in 1:6
        κ += sqrtθ^(1 + i) * (λμ^2 * table[i, 1] + table[i, 2] +
                              λinv * (table[i, 3] + λinv * (table[i, 4] +
                              λinv * (table[i, 5] + λinv * table[i, 6]))))
    end
    return κ * 1e-32u"m^4/N" * k_B * temperature |> u"m^5"
end

"""
    α_hminus_ff_john(
        λ::Unitful.Length,
        temperature::Unitful.Temperature,
        h_neutral_density::NumberDensity,
        electron_density::NumberDensity
    )

Compute free-free extinction from H minus ion. Based on `σ_hminus_ff_john`.
"""
function α_hminus_ff_john(
    λ::Unitful.Length,
    temperature::Unitful.Temperature,
    h_neutral_density::NumberDensity,
    electron_density::NumberDensity
)
    σ = σ_hminus_ff_john(λ, temperature)
    # NOTE: in RH temperature from electron pressure is set to 5040 K!
    # Also, RH uses sometimes nHminus = nH * pe, other times atmos.nHmin...
    return σ * h_neutral_density * electron_density
end


"""
    σ_hminus_bf_john(λ::Unitful.Length, temperature::Unitful.Temperature)

Compute free-free cross section coefficient from H minus ion, for a given wavelength
and temperature. Units are m^5, needs to be multiplied by electron density and
density of neutral hydrogen atoms to obtain linear extinction. Uses recipe from
[John (1988)](https://ui.adsabs.harvard.edu/abs/1988A%26A...193..189J/abstract),
which is valid for  0.125 <= λ (μm) <= 1.6419. Seems to already include stimulated
emission.
"""
function σ_hminus_bf_john(λ::Unitful.Length, temperature::Unitful.Temperature)
    table = SA[152.519, 49.534, -118.858, 92.536, -34.194, 4.982]
    λμ = ustrip(λ |> u"μm")
    λ0 = ustrip(1.6419u"μm")
    λ1 = ustrip(0.125u"μm")  # edge wavelength when approximation for fλ no longer valid
    temp = ustrip(temperature |> u"K")
    λidiff = max(0.0, 1.0 / λμ - 1.0 / λ0)  # cases beyond λ0 set to zero
    σλ = 1e-18 * λμ^3 * λidiff^1.5
    fλ = 0.0
    if λμ < λ1
        λidiff = 1.0 / λ1 - 1.0 / λ0
    end
    for n in 1:6
        fλ += table[n] * λidiff^((n-1)/2)
    end
    σλ *= fλ
    α = h_k * c_0
    κ = 0.750 * sqrt(temp)^-5 * exp(α / (λ0*u"μm" * temperature)) *
       (1 - exp(-α / (λ * temperature))) * σλ * 1e-3u"m^4/N"
    return κ * k_B * temperature |> u"m^5"
end


"""
    α_hminus_bf_john(
        λ::Unitful.Length,
        temperature::Unitful.Temperature,
        h_neutral_density::NumberDensity,
        electron_density::NumberDensity
    )

Compute extinction from H minus ion. Based on `σ_hminus_bf_john`.
"""
function α_hminus_bf_john(
    λ::Unitful.Length,
    temperature::Unitful.Temperature,
    h_neutral_density::NumberDensity,
    electron_density::NumberDensity
)
    σ = σ_hminus_bf_john(λ, temperature)
    return σ * h_neutral_density * electron_density
end

#=----------------------------------------------------------------------------
        Recipes from Wishart (1979) and Broad and Reinhardt (1976)
----------------------------------------------------------------------------=#
const wbr_λ = [     18, 19.6, 21.4, 23.6, 26.4, 29.8, 34.3, 40.4, 49.1, 62.6,  121, 139,
                   164,  175,  200,  225,  250,  275,  300,  325,  350,  375,  400, 425,
                   450,  475,  500,  525,  550,  575,  600,  625,  650,  675,  700, 725,
                   750,  775,  800,  825,  850,  875,  900,  925,  950,  975, 1000, 1025,
                  1050, 1075, 1100, 1125, 1150, 1175, 1200, 1225, 1250, 1275, 1300, 1325,
                  1350, 1375, 1400, 1425, 1450, 1475, 1500, 1525, 1550, 1575, 1600, 1610,
                  1620, 1630]   # in nm
const wbr_σ = [0.067, 0.088, 0.117, 0.155, 0.206, 0.283, 0.414, 0.703,  1.24,  2.33,
                5.43,  5.91,  7.29, 7.918, 9.453, 11.08, 12.75, 14.46, 16.19, 17.92,
               19.65, 21.35, 23.02, 24.65, 26.24, 27.77, 29.23, 30.62, 31.94, 33.17,
               34.32, 35.37, 36.32, 37.17, 37.91, 38.54, 39.07, 39.48, 39.77, 39.95,
               40.01, 39.95, 39.77, 39.48, 39.06, 38.53, 37.89, 37.13, 36.25, 35.28,
               34.19, 33.01, 31.72, 30.34, 28.87, 27.33, 25.71, 24.02, 22.26, 20.46,
               18.62, 16.74, 14.85, 12.95, 11.07, 9.211, 7.407, 5.677, 4.052, 2.575,
               1.302, 0.8697, 0.4974, 0.1989]  # in 1e-22 m^2
const wbr_bf_interp = LinearInterpolation(wbr_λ, wbr_σ, extrapolation_bc=0)


"""
    σ_hminus_bf_wbr(λ::Unitful.Length, temperature::Unitful.Temperature)

Compute cross section coefficient for bound-free from H minus ion (units m^5).
Needs to be multiplied by density of neutral hydrogen atoms and electron density
to obtain linear extinction. Uses recipe from
[Wishart (1979)](https://ui.adsabs.harvard.edu/abs/1979MNRAS.187P..59W) for λ down to 175 nm,
and recipe from [Broad and Reinhardt (1976)](https://ui.adsabs.harvard.edu/abs/1976PhRvA..14.2159B)
for λ=164 nm and below, following the recommendation from Mathisen (1984, MSC thesis).
"""
function σ_hminus_bf_wbr(λ::Unitful.Length, temperature::Unitful.Temperature)
    λi = ustrip(λ |> u"nm")   # convert to units of table
    κ = wbr_bf_interp(λi)::Float64 * 1e-22u"m^2"
    stimulated_emission = exp(-hc_k / (λ * temperature))
    # Get H- fraction to convert from σ per H- atom to σ per H atom per electron
    hminus_frac = calc_hminus_density(1.0u"m^-3", temperature, 1.0u"m^-3") * u"m^6"
    return κ * (1 - stimulated_emission) * hminus_frac
end

"""
    α_hminus_bf_wbr(
        λ::Unitful.Length,
        temperature::Unitful.Temperature,
        h_minus_density::NumberDensity
    )
    α_hminus_bf_wbr(
        λ::Unitful.Length,
        temperature::Unitful.Temperature,
        h_neutral_density::NumberDensity,
        electron_density::NumberDensity
    )

Compute extinction from H minus ion, from input H minus populations. Uses recipe from
[Wishart (1979)](https://ui.adsabs.harvard.edu/abs/1979MNRAS.187P..59W) for λ down to 175 nm,
and recipe from [Broad and Reinhardt (1976)](https://ui.adsabs.harvard.edu/abs/1976PhRvA..14.2159B)
for λ=164 nm and below, following the recommendation from Mathisen (1984, MSC thesis).
"""
function α_hminus_bf_wbr(
    λ::Unitful.Length,
    temperature::Unitful.Temperature,
    h_minus_density::NumberDensity
)
    λi = ustrip(λ |> u"nm")   # convert to units of table
    κ = wbr_bf_interp(λi)::Float64 * 1e-22u"m^2"
    stimulated_emission = exp(-hc_k / (λ * temperature))
    return κ * h_minus_density * (1 - stimulated_emission)
end

function α_hminus_bf_wbr(
    λ::Unitful.Length,
    temperature::Unitful.Temperature,
    h_neutral_density::NumberDensity,
    electron_density::NumberDensity
)
    h_minus_density = calc_hminus_density(h_neutral_density, temperature, electron_density)
    return α_hminus_bf_wbr(λ, temperature, h_minus_density)
end

#=----------------------------------------------------------------------------
                            Recipes from Mihalas
----------------------------------------------------------------------------=#
"""
    σ_hydrogenic_ff_σ(
        ν::Unitful.Frequency,
        charge::Real
    )

Compute free-free cross section (units m^5) for a hydrogen-like species.
Following Mihalas (1978) p. 101 and
[Rutten's IART](https://www.uio.no/studier/emner/matnat/astro/AST4310/h20/pensumliste/iart.pdf)
p 69. For linear extinction, needs to be multiplied by electron density
and species density. Includes stimulated emission.

Will not work in single-precision.
"""
function σ_hydrogenic_ff(
    ν::Unitful.Frequency,
    temperature::Unitful.Temperature,
    charge::Real
)
    ν = ν |> u"s^-1"
    stimulated_emission = exp(-h_k * ν / temperature)
    return (αff_const * charge^2 / sqrt(temperature) * ν^-3 *
            gaunt_ff(ν, temperature, charge)) * (1 - stimulated_emission)
end

"""
    α_hydrogenic_ff(
        ν::Unitful.Frequency,
        temperature::Unitful.Temperature,
        electron_density::NumberDensity,
        species_density::NumberDensity,
        charge::Int
    )

Compute free-free linear extinction using `hydrogenic_ff_σ`.
For the hydrogen case, `ion_density` is the proton density (H II).
"""
function α_hydrogenic_ff(
    ν::Unitful.Frequency,
    temperature::Unitful.Temperature,
    electron_density::NumberDensity,
    species_density::NumberDensity,
    charge::Real
)
    α = σ_hydrogenic_ff(ν, temperature, charge) * electron_density * species_density
    return α
end

"""
    σ_hydrogenic_bf(
        ν::Unitful.Frequency,
        ν_edge::Unitful.Frequency,
        temperature::Unitful.Temperature,
        species_density::NumberDensity,
        charge::Real,
        n_eff::AbstractFloat
    )

Compute bound-free cross section for a hydrogen-like species in m^2. Multiply
by number density of species to obtain linear extinction. Following Mihalas (1978) p. 99 and
[Rutten's IART](https://www.uio.no/studier/emner/matnat/astro/AST4310/h20/pensumliste/iart.pdf)
p. 70. Using simplified constant and expression for threshold cross section.
"""
function σ_hydrogenic_bf(
    ν::Unitful.Frequency,
    ν_edge::Unitful.Frequency,
    temperature::Unitful.Temperature,
    charge::Real,
    n_eff::AbstractFloat
)
    if ν < ν_edge
        return 0 * αbf_const
    else
        ν3_ratio = (ν_edge / ν)^3
        stimulated_emission = exp(-h_k * ν / temperature)
        return (αbf_const * charge^4 * ν3_ratio * n_eff *
                (1 - stimulated_emission) * gaunt_bf(ν, charge, n_eff))
    end
end


"""
    α_hydrogenic_bf(
        ν::Unitful.Frequency,
        ν_edge::Unitful.Frequency,
        temperature::Unitful.Temperature,
        species_density::NumberDensity,
        charge::Real,
        n_eff::AbstractFloat
    )

Compute bound-free extinction for a hydrogen-like species.
Based on `σ_hydrogenic_bf`.
"""
function α_hydrogenic_bf(
    ν::Unitful.Frequency,
    ν_edge::Unitful.Frequency,
    temperature::Unitful.Temperature,
    species_density::NumberDensity,
    charge::Real,
    n_eff::AbstractFloat
)
    σ = σ_hydrogenic_bf(ν, ν_edge, temperature, charge, n_eff)
    return σ * species_density
end


"""
    σ_hydrogenic_bf_scaled(
        σ0::Unitful.Area,
        ν::Unitful.Frequency,
        ν_edge::Unitful.Frequency,
        charge::Real,
        n_eff::AbstractFloat
    )

    σ_hydrogenic_bf_scaled(
        σ0::Unitful.Area,
        λ::Unitful.Length,
        λ_edge::Unitful.Length,
        charge::Real,
        n_eff::AbstractFloat
)

Compute bound-free cross section for a hydrogen-like species by scaling
a peak cross section σ0 with frequency and the appropriate Gaunt factor.
No stimulated emission is added.
"""
function σ_hydrogenic_bf_scaled(
    σ0::Unitful.Area,
    ν::Unitful.Frequency,
    ν_edge::Unitful.Frequency,
    charge::Real,
    n_eff::AbstractFloat
)
    if ν < ν_edge
        σ = 0 * σ0
    else
        σ = σ0 * (ν_edge / ν)^3 * (
            gaunt_bf(ν, charge, n_eff) / gaunt_bf(ν_edge, charge, n_eff))
    end
    return σ
end

function σ_hydrogenic_bf_scaled(
    σ0::Unitful.Area,
    λ::Unitful.Length,
    λ_edge::Unitful.Length,
    charge::Real,
    n_eff::AbstractFloat
)
    return σ_hydrogenic_bf_scaled(σ0, c_0 / λ, c_0 / λ_edge, charge, n_eff)
end


#=----------------------------------------------------------------------------
                            Recipes from Bell (1980)
----------------------------------------------------------------------------=#
const bell_ff_λ = [   0.0,  350.5,  414.2,   506.3,   569.6,  650.9,  759.4,
                    911.3, 1139.1, 1518.8,  1822.6,  2278.3, 3037.7, 3645.2,
                   4556.5, 6075.3, 9113.0, 11391.3, 15188.3]   # in nm
const bell_ff_t = 5040.0 ./ [0.5, 0.8, 1.0, 1.2, 1.6, 2.0, 2.8, 3.6]  # in K
const bell_ff_κ =
    [0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00
     4.17e-02  6.10e-02  7.34e-02  8.59e-02  1.11e-01  1.37e-01  1.87e-01  2.40e-01
     5.84e-02  8.43e-02  1.01e-01  1.17e-01  1.49e-01  1.82e-01  2.49e-01  3.16e-01
     8.70e-02  1.24e-01  1.46e-01  1.67e-01  2.10e-01  2.53e-01  3.39e-01  4.27e-01
     1.10e-01  1.54e-01  1.80e-01  2.06e-01  2.55e-01  3.05e-01  4.06e-01  5.07e-01
     1.43e-01  1.98e-01  2.30e-01  2.59e-01  3.17e-01  3.75e-01  4.92e-01  6.09e-01
     1.92e-01  2.64e-01  3.03e-01  3.39e-01  4.08e-01  4.76e-01  6.13e-01  7.51e-01
     2.73e-01  3.71e-01  4.22e-01  4.67e-01  5.52e-01  6.33e-01  7.97e-01  9.63e-01
     4.20e-01  5.64e-01  6.35e-01  6.97e-01  8.06e-01  9.09e-01  1.11e+00  1.32e+00
     7.36e-01  9.75e-01  1.09e+00  1.18e+00  1.34e+00  1.48e+00  1.74e+00  2.01e+00
     1.05e+00  1.39e+00  1.54e+00  1.66e+00  1.87e+00  2.04e+00  2.36e+00  2.68e+00
     1.63e+00  2.14e+00  2.36e+00  2.55e+00  2.84e+00  3.07e+00  3.49e+00  3.90e+00
     2.89e+00  3.76e+00  4.14e+00  4.44e+00  4.91e+00  5.28e+00  5.90e+00  6.44e+00
     4.15e+00  5.38e+00  5.92e+00  6.35e+00  6.99e+00  7.50e+00  8.32e+00  9.02e+00
     6.47e+00  8.37e+00  9.20e+00  9.84e+00  1.08e+01  1.16e+01  1.28e+01  1.38e+01
     1.15e+01  1.48e+01  1.63e+01  1.74e+01  1.91e+01  2.04e+01  2.24e+01  2.40e+01
     2.58e+01  3.33e+01  3.65e+01  3.90e+01  4.27e+01  4.54e+01  4.98e+01  5.33e+01
     4.03e+01  5.20e+01  5.70e+01  6.08e+01  6.65e+01  7.08e+01  7.76e+01  8.30e+01
     7.16e+01  9.23e+01  1.01e+02  1.08e+02  1.18e+02  1.26e+02  1.38e+02  1.47e+02]
# Linear extrapolation in λ, flat extrapolation in θ
const bell_ff_interp = LinearInterpolation((bell_ff_λ, bell_ff_t[end:-1:1]),
                   bell_ff_κ[:, end:-1:1], extrapolation_bc=(Line(), Flat()))

"""
    σ_h2minus_ff(λ::Unitful.Length, temperature::Unitful.Temperature)

Compute free-free cross section coefficient from H2^- molecule.  Units are m^5,
needs to be multiplied by electron density anddensity of neutral hydrogen atoms
to obtain linear extinction. Follows recipe from
[Bell (1980)](https://ui.adsabs.harvard.edu/abs/1980JPhB...13.1859B/abstract),
page 1863. Stimulated emission is included.
"""
function σ_h2minus_ff(λ::Unitful.Length, temperature::Unitful.Temperature)
    λi = ustrip(λ |> u"nm")   # convert to units of table
    temp = ustrip(temperature |> u"K")
    κ = bell_ff_interp(λi, temp)::Float64 * 1e-29u"m^4/N"
    return κ * k_B * temperature |> u"m^5"
end


"""
    α_h2minus_ff(
        λ::Unitful.Length,
        temperature::Unitful.Temperature,
        h2_density::NumberDensity,
        electron_density::NumberDensity
    )

Compute extinction from H2^- molecule. Based on `σ_h2minus_ff`.
"""
function α_h2minus_ff(
    λ::Unitful.Length,
    temperature::Unitful.Temperature,
    h2_density::NumberDensity,
    electron_density::NumberDensity
)
    σ = σ_h2minus_ff(λ, temperature)
    return σ * h2_density * electron_density |> u"m^-1"
end


#=----------------------------------------------------------------------------
                            Recipes from Bates (1952)
----------------------------------------------------------------------------=#
const bates_λ = 1f7 ./ [   500,   1000,   1500,   2000,   2500,   3000,
                          3500,   4000,   5000,   6000,   7000,   8000,
                          9000, 10_000, 12_000, 14_000, 16_000, 18_000,
                        20_000, 22_000, 24_000, 26_000,  Inf32]  # in nm
const bates_t = [2.5f+03, 3.0f+03, 3.5f+03, 4.0f+03, 5.0f+03,
                    6.0f+03, 7.0f+03, 8.0f+03, 1.0f+04, 1.2f+04]  # in K
# H2plus bf + ff extinction coefficient in 1e-49 m^5
const bates_κ = convert(Array{Float32, 2},
                [  1.14  0.94  0.80  0.69  0.55  0.46  0.39  0.34  0.27  0.226
                   1.61  1.32  1.12  0.97  0.77  0.63  0.54  0.47  0.38  0.31
                   1.97  1.60  1.35  1.17  0.92  0.76  0.64  0.56  0.44  0.37
                   2.28  1.83  1.53  1.32  1.03  0.85  0.72  0.63  0.50  0.41
                   2.56  2.04  1.69  1.45  1.13  0.92  0.78  0.68  0.54  0.44
                   2.84  2.23  1.84  1.56  1.21  0.99  0.83  0.72  0.57  0.47
                   3.1   2.42  1.97  1.67  1.28  1.04  0.88  0.76  0.60  0.49
                   3.4   2.60  2.10  1.77  1.35  1.09  0.92  0.79  0.62  0.51
                   4.0   2.98  2.36  1.96  1.47  1.17  0.98  0.84  0.66  0.54
                   4.8   3.4   2.63  2.15  1.57  1.25  1.04  0.89  0.69  0.57
                   5.6   3.9   2.91  2.33  1.67  1.31  1.08  0.92  0.71  0.58
                   6.7   4.4   3.2   2.53  1.77  1.37  1.12  0.95  0.73  0.59
                   7.9   5.0   3.5   2.74  1.87  1.43  1.16  0.97  0.74  0.60
                   9.3   5.6   3.9   2.95  1.97  1.48  1.19  0.99  0.75  0.61
                  13.0   7.2   4.7   3.4   2.18  1.58  1.25  1.03  0.77  0.62
                  18.1   9.3   5.7   4.0   2.40  1.69  1.30  1.06  0.78  0.62
                  25.1  11.9   7.0   4.7   2.64  1.80  1.36  1.09  0.78  0.62
                  35.   15.2   8.4   5.4   2.91  1.91  1.41  1.11  0.79  0.61
                  47.   19.3  10.2   6.3   3.2   2.03  1.46  1.14  0.79  0.61
                  64.   24.3  12.2   7.3   3.5   2.16  1.52  1.16  0.79  0.60
                  86.   31.   14.6   8.4   3.8   2.29  1.57  1.18  0.79  0.59
                 114.   38.   17.3   9.6   4.2   2.42  1.63  1.21  0.79  0.58
                   0.    0.    0.    0.    0.    0.    0.    0.    0.    0.  ])
const bates_bf_fraction = convert(Array{Float32, 2},
        [0.059  0.046  0.037  0.031  0.022  0.017  0.014  0.011  0.008  0.006
         0.135  0.107  0.087  0.072  0.053  0.041  0.033  0.027  0.020  0.015
         0.214  0.171  0.141  0.118  0.088  0.069  0.056  0.046  0.034  0.026
         0.291  0.236  0.196  0.166  0.125  0.098  0.080  0.067  0.049  0.038
         0.363  0.298  0.250  0.214  0.162  0.129  0.105  0.088  0.065  0.050
         0.430  0.357  0.303  0.260  0.200  0.159  0.131  0.110  0.082  0.064
         0.490  0.413  0.353  0.305  0.237  0.190  0.157  0.132  0.099  0.077
         0.546  0.464  0.400  0.349  0.273  0.221  0.183  0.155  0.116  0.091
         0.640  0.556  0.486  0.429  0.342  0.280  0.234  0.200  0.151  0.120
         0.715  0.632  0.561  0.501  0.406  0.336  0.284  0.243  0.186  0.148
         0.775  0.696  0.625  0.564  0.464  0.388  0.331  0.285  0.220  0.176
         0.822  0.748  0.680  0.619  0.517  0.437  0.375  0.326  0.254  0.204
         0.859  0.792  0.727  0.667  0.564  0.482  0.417  0.364  0.286  0.232
         0.888  0.827  0.767  0.709  0.607  0.524  0.456  0.400  0.317  0.258
         0.929  0.881  0.829  0.777  0.680  0.597  0.526  0.467  0.375  0.309
         0.954  0.917  0.874  0.829  0.739  0.658  0.586  0.525  0.428  0.356
         0.970  0.942  0.906  0.867  0.786  0.708  0.638  0.577  0.476  0.400
         0.980  0.959  0.930  0.897  0.824  0.751  0.683  0.622  0.519  0.440
         0.987  0.970  0.947  0.919  0.854  0.786  0.721  0.661  0.558  0.476
         0.991  0.978  0.960  0.936  0.879  0.816  0.754  0.695  0.593  0.510
         0.994  0.984  0.969  0.949  0.898  0.841  0.782  0.726  0.625  0.541
         0.996  0.988  0.976  0.959  0.914  0.862  0.806  0.752  0.653  0.569
         0.     0.     0.     0.     0.     0.     0.     0.     0.     0.   ])
const bates_ff_fraction = 1 .- bates_bf_fraction
const bates_ff_κ = bates_κ .* bates_ff_fraction
const bates_bf_κ = bates_κ .* bates_bf_fraction
const bates_ff_interp = LinearInterpolation((bates_λ[end:-1:1], bates_t),
                          bates_ff_κ[end:-1:1, :], extrapolation_bc=(Line(), Flat()))
const bates_bf_interp = LinearInterpolation((bates_λ[end:-1:1], bates_t),
                          bates_bf_κ[end:-1:1, :], extrapolation_bc=(Line(), Flat()))


"""
    σ_h2plus_ff(λ::Unitful.Length, temperature::Unitful.Temperature)

Compute free-free cross section coefficient from H2plus molecule, according to recipe
from [Bates (1952)](https://ui.adsabs.harvard.edu/abs/1952MNRAS.112...40B/abstract),
page 43.
"""
function σ_h2plus_ff(λ::Unitful.Length, temperature::Unitful.Temperature)
    λi = convert(Float32, ustrip(λ |> u"nm"))   # convert to units of table
    temp = convert(Float32, ustrip(temperature |> u"K"))
    σ = bates_ff_interp(λi, temp)::Float32 * u"m^5" * 1e-49
    return σ
end

"""
    α_h2plus_ff(λ::Unitful.Length, temperature::Unitful.Temperature,
              h_neutral_density::NumberDensity, proton_density::NumberDensity)

Compute free-free extinction from H2plus molecule. Based on `σ_h2plus_ff`.
"""
function α_h2plus_ff(
    λ::Unitful.Length,
    temperature::Unitful.Temperature,
    h_neutral_density::NumberDensity,
    proton_density::NumberDensity
)
    σ = σ_h2plus_ff(λ, temperature)
    return σ * h_neutral_density * proton_density
end


"""
    σ_h2plus_bf(λ::Unitful.Length, temperature::Unitful.Temperature)

Compute bound-free extinction from H2plus molecule, according to recipe from
[Bates (1952)](https://ui.adsabs.harvard.edu/abs/1952MNRAS.112...40B/abstract),
page 43.
"""
function σ_h2plus_bf(λ::Unitful.Length, temperature::Unitful.Temperature)
    λi = convert(Float32, ustrip(λ |> u"nm"))  # convert to units of table
    temp = convert(Float32, ustrip(temperature |> u"K"))
    # Table in 1e-49 m^5
    return bates_bf_interp(λi, temp)::Float32 * u"m^5" * 1e-49
end

"""
    h2plus_bf(λ::Unitful.Length, temperature::Unitful.Temperature,
              h_neutral_density::NumberDensity, proton_density::NumberDensity)

Compute bound-free extinction from H2plus molecule, based on `σ_h2plus_bf`.

"""
function α_h2plus_bf(
    λ::Unitful.Length,
    temperature::Unitful.Temperature,
    h_neutral_density::NumberDensity,
    proton_density::NumberDensity
)
    σ = σ_h2plus_bf(λ, temperature)
    return σ * h_neutral_density * proton_density
end

#=----------------------------------------------------------------------------
                  Recipes from Victor & Dalgarno (1969)
----------------------------------------------------------------------------=#
const victor_h2_λ = SA[121.57, 130.00, 140.00, 150.00, 160.00, 170.00, 185.46,
                       186.27, 193.58, 199.05, 230.29, 237.91, 253.56, 275.36,
                       296.81, 334.24, 404.77, 407.90, 435.96, 546.23, 632.80]  # in nm
const victor_h2_σ = SA[2.35E-28, 1.22E-28, 6.80E-29, 4.24E-29, 2.84E-29, 2.00E-29,
                       1.25E-29, 1.22E-29, 1.00E-29, 8.70E-30, 4.29E-30, 3.68E-30,
                       2.75E-30, 1.89E-30, 1.36E-30, 8.11E-31, 3.60E-31, 3.48E-31,
                       2.64E-31, 1.04E-31, 5.69E-32]  # in m2
const victor_h2_interp = LinearInterpolation(victor_h2_λ, victor_h2_σ, extrapolation_bc=Line())

"""
    function σ_rayleigh_h2(λ::Unitful.Length)

Compute cross section from Rayleigh scattering from H2 molecules. Uses recipe from
[Victor and Dalgarno (1969)](https://aip.scitation.org/doi/pdf/10.1063/1.1671412),
J. Chem. Phys. 50, 2535, page 2538 for λ <= 632.80, and the recipe from
[Tarafdar & Vardya (1973)](https://ui.adsabs.harvard.edu/abs/1973MNRAS.163..261T/abstract)
page 272, for λ > 632.80.
"""
function σ_rayleigh_h2(λ::Unitful.Length)
    λi = ustrip(λ |> u"nm")
    if λi >= victor_h2_λ[1]
        if λi <= victor_h2_λ[end]
            σ_h2 = victor_h2_interp(λi)::Float64 * u"m^2"
        else
            λ2 =1 / λi^2
            # Tarafdar coeffs converted from λ^4, λ^6, λ^8 in Å to nm and cm^2 to m^2:
            σ_h2 = (8.779e-21 + (1.323e-16 + 2.245e-12 * λ2) * λ2) * λ2^2 * u"m^2"
        end
    else
        σ_h2 = 0.0u"m^2"
    end
    return σ_h2
end

"""
    function α_rayleigh_h2(λ::Unitful.Length, h2_density::NumberDensity)

Compute extinction from Rayleigh scattering from H2 molecules. Based on `σ_rayleigh_h2`.
"""
function α_rayleigh_h2(λ::Unitful.Length, h2_density::NumberDensity)
    σ_h2 = σ_rayleigh_h2(λ)
    return σ_h2 * h2_density
end


#=----------------------------------------------------------------------------
                  Recipes from Dalgarno (1962)
----------------------------------------------------------------------------=#
"""
    σ_rayleigh_h(λ::Unitful.Length, h_ground_density::NumberDensity)

Compute cross section from Rayleigh scattering from H atoms. To obtain extinction,
multiply by the number density of neutral hydrogen atoms in the ground state.
Uses recipe from Dalgarno (1962), Geophysics Corp. of America, Technical Report No. 62-28-A
(unavailable), which is accurate to 1% for λ > 125.0 nm.
"""
function σ_rayleigh_h(λ::Unitful.Length)
    λi = ustrip(λ |> u"Å")
    if λi >= 1215.7
        λ2 = 1 / λi^2
        # First coefficient has conversion from Mbarn to m^2. From RH:
        σ_h = (5.81e-17 * λ2^2 * (1 + 2.452e6 * λ2 +  4.801e12 * λ2^2)) * u"m^2"
    else
        σ_h = 0.0u"m^2"
    end
    return σ_h
end


"""
    α_rayleigh_h(λ::Unitful.Length, h_ground_density::NumberDensity)

Compute extinction from Rayleigh scattering from H atoms. `h_ground_density`
is the number density of neutral hydrogen atoms in the ground state.
Based on `σ_rayleigh_h`.
"""
function α_rayleigh_h(λ::Unitful.Length, h_ground_density::NumberDensity)
    σ_h = σ_rayleigh_h(λ)
    return σ_h * h_ground_density
end
