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