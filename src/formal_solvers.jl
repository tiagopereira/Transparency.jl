"""
Tools for computing the formal solution of the radiative transfer equation.
"""


"""
Compute intensity by "brute force" trapezoidal integration. Doesn't work well in many
cases, its use is discouraged.
"""
function calc_intensity_brute_force(
    distance::Array{<:Unitful.Length, 1},
    extinction::Array{<:PerLength, 1},
    source_function::Array{<:UnitsIntensity_λ, 1}
)
    @assert distance[2] > distance[1] "Distance must be monotonically increasing"
    # Since integration functions don't work with units,
    # need to ensure quantities are in compatible units
    dist = ustrip.(distance .|> u"m")
    ext = ustrip.(extinction .|> u"m^-1")
    source = ustrip.(source_function .|> u"kW / (m^2 * nm)")
    τ = cumul_integrate(dist, ext, TrapezoidalFast())
    return integrate(τ, source .* exp.(-τ))u"kW / (m^2 * nm)"
end


"""
    piecewise_1D_nn(
        z::Array{<:Unitful.Length, 1},
        α::Array{<:PerLength, 1},
        source_function::Array{<:Unitful.Quantity, 1};
        to_end::Bool=false,
        initial_condition=:source,
    )

Compute piecewise integration of the radiative transfer equation,
assuming nearest-neighbour integration of the source function, for a given
height `z`, extinction `α` and `source_function`. The optional
keyword argument `to_end` defines the direction of the integration:
if `false` (default) will start to integrate intensity from the last
element to the first, and if `true` will integrate from the first
element to the last. `initial_condition` can take two values: `:zero` for
no radiation, or `:source` (default) to take the source function at the
starting point.
"""
function piecewise_1D_nn(
    z::Array{<:Unitful.Length, 1},
    α::Array{<:PerLength, 1},
    source_function::Array{<:Unitful.Quantity, 1};
    to_end::Bool=false,
    initial_condition=:source
)
    ndep = length(z)
    if to_end
        start = 1
        incr = 1
        depth_range = 2:ndep
    else
        start = ndep
        incr = -1
        depth_range = ndep-1:-1:1
    end
    intensity = similar(source_function)
    if initial_condition == :source
        intensity[start] = source_function[start]
    elseif initial_condition == :zero
        intensity[start] *= 0
    else
        error("NotImplemented initial condition $initial_condition")
    end
    for i in depth_range
        Δτ = abs(z[i] - z[i-incr]) * (α[i] + α[i-incr]) / 2
        w, _ = _w2(Δτ)
        intensity[i] = ((1 - w)*intensity[i-incr] +
                         w * (source_function[i] + source_function[i-incr]) / 2)
    end
    return intensity
end


"""
    piecewise_1D_linear(
        z::Array{<:Unitful.Length, 1},
        α::Array{<:PerLength, 1},
        source_function::Array{<:Unitful.Quantity, 1};
        to_end::Bool=false,
        initial_condition=:source,
    )

Compute piecewise integration of the radiative transfer equation,
assuming linear integration of the source function, for a given
height `z`, extinction `α` and `source_function`. The optional
keyword argument `to_end` defines the direction of the integration:
if `false` (default) will start to integrate intensity from the last
element to the first, and if `true` will integrate from the first
element to the last. `initial_condition` can take two values: `:zero` for
no radiation, or `:source` (default) to take the source function at the
starting point.
"""
function piecewise_1D_linear(
    z::Array{<:Unitful.Length, 1},
    α::Array{<:PerLength, 1},
    source_function::Array{<:Unitful.Quantity, 1};
    to_end::Bool=false,
    initial_condition=:source,
)
    ndep = length(z)
    if to_end
        start = 1
        incr = 1
        depth_range = 2:ndep
    else
        start = ndep
        incr = -1
        depth_range = ndep-1:-1:1
    end
    intensity = similar(source_function)
    if initial_condition == :source
        intensity[start] = source_function[start]
    elseif initial_condition == :zero
        intensity[start] *= 0
    else
        error("NotImplemented initial condition $initial_condition")
    end
    for i in depth_range
        Δτ = abs(z[i] - z[i-incr]) * (α[i] + α[i-incr]) / 2
        ΔS = (source_function[i-incr] - source_function[i]) / Δτ
        w1, w2 = _w2(Δτ)
        intensity[i] = (1 - w1)*intensity[i-incr] + w1*source_function[i] + w2*ΔS
    end
    return intensity
end


"""
Computes weights for linear integration of source function,
approximating `exp(-Δτ)` for very small and very large values of `Δτ`.
"""
function _w2(Δτ::T) where T <: AbstractFloat
    if Δτ < 5e-4
        w1 = Δτ * (1 - Δτ / 2)
        w2 = Δτ^2 * (0.5f0 - Δτ / 3)
    elseif Δτ > 50
        w1 = w2 = one(T)
    else
        expΔτ = exp(-Δτ)
        w1 = 1 - expΔτ
        w2 = w1 - Δτ * expΔτ
    end
    return w1, w2
end


"""
    function feautrier(
        z::Array{<:Unitful.Length{T}, 1},
        α::Array{<:PerLength{T}, 1},
        source_function::Array{<:Unitful.Quantity, 1}
    ) where T <: AbstractFloat

Calculate solution to radiative transfer equation using the Feautrier method.
Uses algorithm from [Rybicki & Hummer, 1991, A&A 245](https://ui.adsabs.harvard.edu/abs/1991A&A...245..171R).
Returns the height-dependent Feautrier variable `P`:

```math
P \\equiv 1/2 (I^+ + I^-)
```

Currently operates under the following assumptions:
* The first index of the variables is the top of the atmosphere
* The boundary conditions are zero radiation at the top and source function at the bottom

Therefore, the emergent intensity is `2 * P[1]`, since ``I^-[1]=0``.

Not properly tested, use with care!
"""
function feautrier(
    z::Array{<:Unitful.Length{T}, 1},
    α::Array{<:PerLength{T}, 1},
    source_function::Array{<:Unitful.Quantity, 1}
) where T <: AbstractFloat
    ndep = length(z)
    F = Array{T}(undef, ndep)
    Z = similar(source_function)
    P = similar(source_function)
    Δτ = Array{T}(undef, ndep)
    Δτ[end] = zero(T)
    for i in 1:ndep-1
        Δτ[i] = (z[i] - z[i+1]) * (α[i] + α[i+1]) / 2
    end
    H1 = one(T) + 2 / Δτ[1]
    C1 = 2 / Δτ[1]^2
    Hn = one(T) + 2 / Δτ[end - 1]
    An = 2 / Δτ[end - 1]^2
    # Start elimination
    F[1] = H1 / C1
    Z[1] = source_function[1] / (H1 + C1)
    for i in 2:ndep-1
        Δτ_mid = (Δτ[i] + Δτ[i - 1]) / 2
        A = 1 / (Δτ_mid * Δτ[i - 1])
        C = 1 / (Δτ_mid * Δτ[i])
        F[i] = (one(T) + A * F[i - 1] / (one(T) + F[i - 1])) / C
        Z[i] = (source_function[i] + A * Z[i - 1]) / (C * (one(T) + F[i]))
    end
    # Now backsubstitution
    P[end] = (source_function[end] + An * Z[end - 1]) /
        (Hn + An * (F[end - 1] / (one(T) + F[end - 1])))
    for i in ndep-1:-1:1
        P[i] = P[i + 1] / (one(T) + F[i]) + Z[i]
    end
    return P
end
