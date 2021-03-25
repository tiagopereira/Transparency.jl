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
        to_end=false
    )

Compute piecewise integration of the radiative transfer equation,
assuming nearest-neighbour integration of the source function, for a given
height `z`, extinction `α` and `source_function`. The optional
keyword argument `to_end` defines the direction of the integration:
if `false` (default) will start to integrate intensity from the last
element to the first, and if `true` will integrate from the first
element to the last.

Currently, the initial condition for intensity at the start point
is set to the source function.
"""
function piecewise_1D_nn(
    z::Array{<:Unitful.Length, 1},
    α::Array{<:PerLength, 1},
    source_function::Array{<:Unitful.Quantity, 1};
    to_end=false
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
    intensity[start] = source_function[start]
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
        to_end=false
    )

Compute piecewise integration of the radiative transfer equation,
assuming linear integration of the source function, for a given
height `z`, extinction `α` and `source_function`. The optional
keyword argument `to_end` defines the direction of the integration:
if `false` (default) will start to integrate intensity from the last
element to the first, and if `true` will integrate from the first
element to the last.

Currently, the initial condition for intensity at the start point
is set to the source function.
"""
function piecewise_1D_linear(
    z::Array{<:Unitful.Length, 1},
    α::Array{<:PerLength, 1},
    source_function::Array{<:Unitful.Quantity, 1};
    to_end=false
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
    intensity[start] = source_function[start]
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
