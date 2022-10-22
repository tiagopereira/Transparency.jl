"""
Functions to calculate atomic collisional rates.
"""

# Dimension for CE and CI collisional tables, in SI units m^3 K^-1/2 s-1
@derived_dimension CEI_dimension Unitful.𝐋^3 / (Unitful.𝚯^(1/2)  * Unitful.𝐓)

const Ω_c0 = Ry / sqrt(m_e) * π * a_0^2 * sqrt(8 / (π * k_B))

"""
    coll_CE(
        rate_interpolant::Interpolations.AbstractInterpolation{<:CEI_dimension, 1},
        g_ratio::Real,
        electron_density::NumberDensity,
        temperature::Unitful.Temperature
    )

Calculate collisional de-excitation by electrons of a bound-bound transition, using
the CE expression from RH / MULTI.

# Arguments
- `rate_interpolant`: a generic interpolant that takes as argument a value of
temperature, and returns in units similar to m^3 K^-1/2 s^-1
- `g_ratio`: ratio g_l / g_u, statistical weights of lower and upper level
- `electron_density`: electron density
- `temperature`: gas temperature

# Returns
- `coll_deexc`: collisional de-excitations per second (from upper level to lower level)
"""
function coll_CE(
    rate_interpolant::Interpolations.AbstractInterpolation{<:CEI_dimension, 1},
    g_ratio::Real,
    electron_density::NumberDensity,
    temperature::Unitful.Temperature
)
    return (rate_interpolant(temperature) * g_ratio *
            electron_density * sqrt(temperature)) |> u"s^-1"
end


"""
    coll_CI(
        rate_interpolant::Interpolations.AbstractInterpolation{<:CEI_dimension, 1},
        dE::Unitful.Energy,
        electron_density::NumberDensity,
        temperature::Unitful.Temperature
    )

Calculate collisional ionisation by electrons of a given level, using
the CI expression from RH / MULTI.

# Arguments
- `rate_interpolant`: a generic interpolant that takes as argument a value of
temperature, and returns in units similar to m^3 K^-1/2 s^-1
- `dE`: energy difference between continuum and level, dE = E_cont - E_level.
- `electron_density`: electron density
- `temperature`: gas temperature

# Returns
- `coll_ion`: collisional ionisations per second (from level to continuum)
"""
function coll_CI(
    rate_interpolant::Interpolations.AbstractInterpolation{<:CEI_dimension, 1},
    dE::Unitful.Energy,
    electron_density::NumberDensity,
    temperature::Unitful.Temperature
)
    @assert dE > 0u"J" "dE should be positive"
    return (rate_interpolant(temperature) * exp(-dE / (k_B * temperature)) *
            electron_density * sqrt(temperature)) |> u"s^-1"
end


"""
    coll_Ω(
        rate_interpolant::Interpolations.AbstractInterpolation{<:CEI_dimension, 1},
        g_u::Integer,
        electron_density::NumberDensity,
        temperature::Unitful.Temperature
    )

Calculate collisional de-excitation by electrons of a bound-bound transition, using
the OMEGA expression from RH (OHMEGA from MULTI), from a tabulated dimensionless Ω
(collision strength).

# Arguments
- `rate_interpolant`: a generic interpolant that takes as argument a value of
temperature, and returns a dimensionless Ω
- `g_u`: statistical weight of upper level
- `electron_density`: electron density
- `temperature`: gas temperature

# Returns
- `coll_deexc`: collisional de-excitations per second (from upper level to lower level)
"""
function coll_Ω(
    rate_interpolant::Interpolations.AbstractInterpolation{<:Real, 1},
    g_u::Integer,
    electron_density::NumberDensity,
    temperature::Unitful.Temperature
)
    return (rate_interpolant(temperature) * Ω_c0 * electron_density /
            (g_u * sqrt(temperature))) |> u"s^-1"
end


#=----------------------------------------------------------------------------
Collisional rates for hydrogen from Johnson (1972), based on make_h.c from RH
----------------------------------------------------------------------------=#
const johnson_c0 = sqrt(8 * k_B / (π * m_e)) * 2 * π * a_0^2
const johnson_c1 = 32 / (3 * sqrt(3) * π)

"""
    coll_exc_hydrogen_johnson(
        n_l::Integer,
        n_u::Integer,
        electron_density::NumberDensity,
        temperature::Unitful.Temperature,
    )

Calculates the rates for collisional excitation for a hydrogen line, resulting
from collisions with electrons, using the recipe of
[Johnson (1972)](https://ui.adsabs.harvard.edu/abs/1972ApJ...174..227J/abstract)

# Arguments
- `n_l`: lower level (principal quantum number) of transition
- `n_u`: upper level (principal quantum number) of transition
- `electron_density`: electron density
- `temperature`: gas temperature

# Returns
- `coll_exc`: collisional excitations per second from n_l -> n_u
"""
function coll_exc_hydrogen_johnson(
    n_l::Integer,
    n_u::Integer,
    electron_density::NumberDensity,
    temperature::Unitful.Temperature,
)
    @assert n_l > 0 "n_l must be > 0"
    @assert n_u > n_l "n_u must be > n_l"
    rn = _rn(n_l)
    bn = _bn(n_l)
    x = 1 - (n_l / n_u)^2
    rnn = rn * x
    fnn = johnson_c1 * n_l / (n_u * x)^3 * (g0(n_l) + (g1(n_l) + g2(n_l) / x) / x)
    Ann = 2 * n_l^2 / x * fnn
    Bnn = 4 * n_l * (n_l / n_u)^3 * (1 + 4 / (3 * x) + bn / x^2) / x^2
    y = x * Ryh  / (n_l^2 * k_B * temperature)
    z = rnn + y
    coll_exc = (johnson_c0 * (n_l * y)^2 / x * sqrt(temperature) * electron_density *
                (Ann * ((1 / y + 0.5) * expint(1, y) - (1 / z + 0.5) * expint(1, z)) +
                 (Bnn - Ann * log(2 * n_l^2 / x)) * (expint(2, y) / y - expint(2, z) / z)))
    return coll_exc |> u"s^-1"
end


"""
    coll_ion_hydrogen_johnson(
        n::Integer,
        electron_density::NumberDensity,
        temperature::Unitful.Temperature,
    )

Calculates the rates for collisional ionisation for a hydrogen level, resulting
from collisions with electrons, using the recipe of
[Johnson (1972)](https://ui.adsabs.harvard.edu/abs/1972ApJ...174..227J/abstract).

# Arguments
- `n`: principal quantum number of hydrogen level
- `electron_density`: electron density
- `temperature`: gas temperature

# Returns
- `coll_ion`: collisional ionisations per second from level n
"""
function coll_ion_hydrogen_johnson(
    n::Integer,
    electron_density::NumberDensity,
    temperature::Unitful.Temperature,
)
    @assert n > 0 "n must be > 0"
    An = johnson_c1 * n * (g0(n) / 3 + g1(n) / 4 + g2(n) / 5)
    Bn = 2 * n^2 / 3 * (5 + _bn(n))
    yn = Ryh / (n^2 * k_B * temperature)
    zn = _rn(n) + yn
    coll_ion = (johnson_c0 * (n * yn)^2 * sqrt(temperature) * electron_density *
                (An * (expint(1, yn) / yn - expint(1, zn) / zn) +
                 (Bn - An * log(2 * n^2)) * (ξ(yn) - ξ(zn))))
    return coll_ion |> u"s^-1"
end


"""
    CI_RH_hydrogen(n::Integer, temperature::Unitful.Temperature)

Calculate CI coefficients for collisional ionisation for a hydrogen level,
in the SI units of RH, using the recipe of
[Johnson (1972)](https://ui.adsabs.harvard.edu/abs/1972ApJ...174..227J/abstract).

# Arguments
- `n`: principal quantum number of hydrogen level
- `temperature`: gas temperature
"""
function CI_RH_hydrogen(n::Integer, temperature::Unitful.Temperature)
    tmp = coll_ion_hydrogen_johnson(n, 1.0u"m^-3", temperature)
    yn = Ryh / (n^2 * k_B * temperature)
    return tmp * exp(yn) / sqrt(temperature) / 1.0u"m^-3"
end


"""
    CE_RH_hydrogen(n_l::Integer, n_u::Integer, temperature::Unitful.Temperature)

Calculate CE coefficients for collisional deexcitation for a hydrogen transition,
in the SI units of RH, using the recipe of
[Johnson (1972)](https://ui.adsabs.harvard.edu/abs/1972ApJ...174..227J/abstract).

# Arguments
- `n_l`: lower level (principal quantum number) of transition
- `n_u`: upper level (principal quantum number) of transition
- `temperature`: gas temperature
"""
function CE_RH_hydrogen(n_l::Integer, n_u::Integer, temperature::Unitful.Temperature)
    tmp = coll_exc_hydrogen_johnson(n_l, n_u, 1.0u"m^-3", temperature)
    x = 1 - (n_l / n_u)^2
    y = x * Ryh  / (n_l^2 * k_B * temperature)
    return tmp * exp(y) / sqrt(temperature) / 1.0u"m^-3"
end


#=----------------------------------------------------------------------------
                        Utility functions from Johnson (1972)
----------------------------------------------------------------------------=#
# bn and rn expressions, from eqs (24) and (32) of Johnson (1972)
_bn(n::Integer) = (n == 1) ? -0.603 : (4.0+ (-18.63 + (36.24 - 28.09 / n) / n) / n) / n
_rn(n::Integer) = (n == 1) ? 0.45 : 1.94 * n^-1.57

# ξ(t), eq (42)
ξ(t::AbstractFloat) = expint(0, t) - 2 * expint(1, t) + expint(2, t)


# Gaunt factor coefficients, Table 1 of Johnson (1972)
function g0(n::Integer)
    if n == 1
        return 1.1330f0
    elseif n == 2
        return 1.0785f0
    else
        return 0.9935f0 + (0.2328f0 - 0.1296f0 / n) / n
    end
end

function g1(n::Integer)
    if n == 1
        return -0.4059f0
    elseif n == 2
        return -0.2319f0
    else
        return -(0.6282 - (0.5598 - 0.5299 / n) / n) / n
    end
end

function g2(n::Integer)
    if n == 1
        return 0.07014f0
    elseif n == 2
        return 0.02947f0
    else
        return (0.3887f0 - (1.181f0 - 1.4700f0 / n) / n) / n^2
    end
end


#=----------------------------------------------------------------------------
        Collisional rates for hydrogen from Przybilla & Butler (2004)
----------------------------------------------------------------------------=#

const PB04_temp = [0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3,
                    4, 5, 6, 8, 10, 15, 20, 25] * 10000u"K"
const PB04_Ω =
    [[6.40f-1 6.98f-1 7.57f-1 8.09f-1 8.97f-1 9.78f-1 1.06f+0 1.15f+0   #= 1->2
=#    1.32f+0 1.51f+0 1.68f+0 2.02f+0 2.33f+0 2.97f+0 3.50f+0 3.95f+0];
     [2.20f-1 2.40f-1 2.50f-1 2.61f-1 2.88f-1 3.22f-1 3.59f-1 3.96f-1   #= 1->3
=#    4.64f-1 5.26f-1 5.79f-1 6.70f-1 7.43f-1 8.80f-1 9.79f-1 1.06f+0];
     [9.93f-2 1.02f-1 1.10f-1 1.22f-1 1.51f-1 1.80f-1 2.06f-1 2.28f-1   #= 1->4
=#    2.66f-1 2.95f-1 3.18f-1 3.55f-1 3.83f-1 4.30f-1 4.63f-1 4.88f-1];
     [4.92f-2 5.84f-2 7.17f-2 8.58f-2 1.12f-1 1.33f-1 1.50f-1 1.64f-1   #= 1->5
=#    1.85f-1 2.01f-1 2.12f-1 2.29f-1 2.39f-1 2.59f-1 2.71f-1 2.81f-1];
     [2.97f-2 4.66f-2 6.28f-2 7.68f-2 9.82f-2 1.14f-1 1.25f-1 1.33f-1   #= 1->6
=#    1.45f-1 1.53f-1 1.58f-1 1.65f-1 1.70f-1 1.77f-1 1.82f-1 1.85f-1];
     [5.03f-2 6.72f-2 7.86f-2 8.74f-2 1.00f-1 1.10f-1 1.16f-1 1.21f-1   #= 1->7
=#    1.27f-1 1.31f-1 1.34f-1 1.36f-1 1.37f-1 1.39f-1 1.39f-1 1.40f-1];
     [2.35f+1 2.78f+1 3.09f+1 3.38f+1 4.01f+1 4.71f+1 5.45f+1 6.20f+1   #= 2->3
=#    7.71f+1 9.14f+1 1.05f+2 1.29f+2 1.51f+2 1.93f+2 2.26f+2 2.52f+2];
     [1.07f+1 1.15f+1 1.23f+1 1.34f+1 1.62f+1 1.90f+1 2.18f+1 2.44f+1   #= 2->4
=#    2.89f+1 3.27f+1 3.60f+1 4.14f+1 4.56f+1 5.31f+1 5.83f+1 6.23f+1];
     [5.22f+0 5.90f+0 6.96f+0 8.15f+0 1.04f+1 1.23f+1 1.39f+1 1.52f+1   #= 2->5
=#    1.74f+1 1.90f+1 2.03f+1 2.23f+1 2.37f+1 2.61f+1 2.78f+1 2.89f+1];
     [2.91f+0 4.53f+0 6.06f+0 7.32f+0 9.17f+0 1.05f+1 1.14f+1 1.21f+1   #= 2->6
=#    1.31f+1 1.38f+1 1.44f+1 1.51f+1 1.56f+1 1.63f+1 1.68f+1 1.71f+1];
     [5.25f+0 7.26f+0 8.47f+0 9.27f+0 1.03f+1 1.08f+1 1.12f+1 1.14f+1   #= 2->7
=#    1.17f+1 1.18f+1 1.19f+1 1.19f+1 1.20f+1 1.19f+1 1.19f+1 1.19f+1];
     [1.50f+2 1.90f+2 2.28f+2 2.70f+2 3.64f+2 4.66f+2 5.70f+2 6.72f+2   #= 3->4
=#    8.66f+2 1.04f+3 1.19f+3 1.46f+3 1.67f+3 2.08f+3 2.39f+3 2.62f+3];
     [7.89f+1 9.01f+1 1.07f+2 1.26f+2 1.66f+2 2.03f+2 2.37f+2 2.68f+2   #= 3->5
=#    3.19f+2 3.62f+2 3.98f+2 4.53f+2 4.95f+2 5.68f+2 6.16f+2 6.51f+2];
     [4.13f+1 6.11f+1 8.21f+1 1.01f+2 1.31f+2 1.54f+2 1.72f+2 1.86f+2   #= 3->6
=#    2.08f+2 2.24f+2 2.36f+2 2.53f+2 2.65f+2 2.83f+2 2.94f+2 3.02f+2];
     [7.60f+1 1.07f+2 1.25f+2 1.37f+2 1.52f+2 1.61f+2 1.68f+2 1.72f+2   #= 3->7
=#    1.78f+2 1.81f+2 1.83f+2 1.85f+2 1.86f+2 1.87f+2 1.86f+2 1.87f+2];
     [5.90f+2 8.17f+2 1.07f+3 1.35f+3 1.93f+3 2.47f+3 2.96f+3 3.40f+3   #= 4->5
=#    4.14f+3 4.75f+3 5.25f+3 6.08f+3 6.76f+3 8.08f+3 9.13f+3 1.00f+4];
     [2.94f+2 4.21f+2 5.78f+2 7.36f+2 1.02f+3 1.26f+3 1.46f+3 1.64f+3   #= 4->6
=#    1.92f+3 2.15f+3 2.33f+3 2.61f+3 2.81f+3 3.15f+3 3.36f+3 3.51f+3];
     [4.79f+2 7.06f+2 8.56f+2 9.66f+2 1.11f+3 1.21f+3 1.29f+3 1.34f+3   #= 4->7
=#    1.41f+3 1.46f+3 1.50f+3 1.55f+3 1.57f+3 1.61f+3 1.62f+3 1.63f+3];
     [1.93f+3 2.91f+3 4.00f+3 5.04f+3 6.81f+3 8.20f+3 9.29f+3 1.02f+4   #= 5->6
=#    1.15f+4 1.26f+4 1.34f+4 1.49f+4 1.63f+4 1.97f+4 2.27f+4 2.54f+4];
     [1.95f+3 3.24f+3 4.20f+3 4.95f+3 6.02f+3 6.76f+3 7.29f+3 7.70f+3   #= 5->7
=#    8.26f+3 8.63f+3 8.88f+3 9.21f+3 9.43f+3 9.78f+3 1.00f+4 1.02f+4];
     [6.81f+1 1.17f+4 1.50f+4 1.73f+4 2.03f+4 2.21f+4 2.33f+4 2.41f+4   #= 6->7
=#    2.52f+4 2.60f+4 2.69f+4 2.90f+4 3.17f+4 3.94f+4 4.73f+4 5.50f+4]]
# Index of transition data in matrix form: PB04_index[n_l, n_u] = i,
# where i is row index in PB04_Ω
const PB04_index = [0  1  2  3  4  5  6
                    0  0  7  8  9 10 11
                    0  0  0 12 13 14 15
                    0  0  0  0 16 17 18
                    0  0  0  0  0 19 20
                    0  0  0  0  0  0 21]
const PB04_interp = [linear_interpolation(PB04_temp, PB04_Ω[i, :],
                                          extrapolation_bc=Line()) for i in 1:16]

"""
    coll_deexc_hydrogen_PB04(
        n_l::Integer,
        n_u::Integer,
        g_u::Integer,
        electron_density::NumberDensity,
        temperature::Unitful.Temperature,
    )

Calculate collisional de-excitation by electrons of a bound-bound hydrogen transition,
using the Ω data from Table 3 of
[Przybilla & Butler (2004)](https://ui.adsabs.harvard.edu/abs/2004ApJ...610L..61P/abstract).
Data available only up to n=7.

# Arguments
- `n_l`: lower level (principal quantum number) of transition
- `n_u`: upper level (principal quantum number) of transition
- `g_u`: statistical weight of upper level
- `electron_density`: electron density
- `temperature`: gas temperature

# Returns
- `coll_deexc`: collisional de-excitations per second (from upper level to lower level)
"""
function coll_deexc_hydrogen_PB04(
    n_l::Integer,
    n_u::Integer,
    g_u::Integer,
    electron_density::NumberDensity,
    temperature::Unitful.Temperature,
)
    @assert 7 > n_l > 0 "Must have 7 > n_l > 0"
    @assert 8 > n_u > n_l "Must have 8 > n_u > n_l"
    return coll_Ω(PB04_interp[PB04_index[n_l, n_u]], g_u, electron_density, temperature)
end
