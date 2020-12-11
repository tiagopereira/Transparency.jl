"""
Set of recipes to compute opacities.
"""
module Transparency
    export hminus_ff, hminus_bf, hydrogenic_ff, hydrogenic_bf, h2minus_ff, 
           h2plus_ff, h2plus_bf, rayleigh_h2, rayleigh_h, thomson, humlicek,
           voigt_profile, dispersion_profile, doppler_width, 
           γ_unsold_const, γ_unsold, calc_Aji, calc_Bji, damping,
           AtomicLine, αline_λ, jline_λ, blackbody_λ, blackbody_ν, calc_intensity
    using Unitful
    using Interpolations
    import NumericalIntegration: integrate, cumul_integrate, TrapezoidalFast
    import PhysicalConstants.CODATA2018: h, k_B, R_∞, c_0, m_e, m_u, e, ε_0, a_0
    @derived_dimension NumberDensity Unitful.𝐋^-3
    @derived_dimension PerLength Unitful.𝐋^-1
    @derived_dimension UnitsIntensity_λ Unitful.𝐋^-1 * Unitful.𝐌 * Unitful.𝐓^-3

    include("hydrogen.jl")
    include("thomson.jl")
    include("voigt.jl")
    include("line.jl")
    include("broadening.jl")
end
