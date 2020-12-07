"""
Set of recipes to compute opacities.
"""
module Transparency
    export hminus_ff, hminus_bf, hydrogenic_ff, hydrogenic_bf, h2minus_ff, 
           h2plus_ff, h2plus_bf, rayleigh_h2, rayleigh_h, thomson, humlicek,
           voigt_profile, dispersion_profile, doppler_width
    using Unitful
    using Interpolations
    import PhysicalConstants.CODATA2018: h, k_B, R_∞, c_0, m_e, e, ε_0
    @derived_dimension NumberDensity Unitful.𝐋^-3
    @derived_dimension PerLength Unitful.𝐋^-1

    include("hydrogen.jl")
    include("thomson.jl")
    include("voigt.jl")
end
