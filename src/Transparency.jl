"""
Set of recipes to compute opacities.
"""
module Transparency
    export σ_hminus_ff, σ_hminus_bf, α_hminus_ff, α_hminus_bf
    export σ_hydrogenic_ff, σ_hydrogenic_bf, α_hydrogenic_ff, α_hydrogenic_bf
    export σ_hydrogenic_bf_scaled
    export σ_h2minus_ff, σ_h2plus_ff, σ_h2plus_bf, α_h2minus_ff, α_h2plus_ff, α_h2plus_bf
    export σ_rayleigh_h2, σ_rayleigh_h, σ_thomson, α_rayleigh_h2, α_rayleigh_h, α_thomson
    export humlicek, voigt_profile, dispersion_profile
    export calc_Aul, calc_Bul, damping, doppler_width
    export const_unsold, γ_unsold, γ_stark_linear
    export const_barklem, γ_barklem
    export const_deridder_rensbergen, γ_deridder_rensbergen
    export const_stark_quadratic, γ_stark_quadratic, γ_stark_quadratic_gray
    export blackbody_λ, blackbody_ν
    export piecewise_1D_nn, piecewise_1D_linear, calc_intensity_brute_force, feautrier
    export coll_CE, coll_CI, coll_Ω
    export coll_deexc_hydrogen_PB04, coll_exc_hydrogen_johnson, coll_ion_hydrogen_johnson
    export CE_RH_hydrogen, CI_RH_hydrogen

    using Interpolations
    using StaticArrays
    using Unitful
    import NumericalIntegration: integrate, cumul_integrate, TrapezoidalFast
    import PhysicalConstants.CODATA2018: h, k_B, R_∞, c_0, m_e, m_u, e, ε_0, a_0
    import SpecialFunctions: expint, gamma

    @derived_dimension NumberDensity Unitful.𝐋^-3
    @derived_dimension PerLength Unitful.𝐋^-1
    @derived_dimension UnitsIntensity_λ Unitful.𝐋^-1 * Unitful.𝐌 * Unitful.𝐓^-3

    const Ar_H = 1.007975  # Atomic weight of hydrogen
    const Ry = R_∞ * c_0 * h  # Rydberg energy
    const Ryh = Ry / (1 + m_e / (Ar_H * m_u))  # Hydrogen ionisation energy
    const αp = 4.5 * 4*π * ε_0 * a_0^3   # Polarisability of hydrogen [F m^2]
    const inv_4πε0 = 1. / (4 * π * ε_0)

    include("line.jl")
    include("broadening.jl")
    include("collisions.jl")
    include("formal_solvers.jl")
    include("hydrogen.jl")
    include("thomson.jl")
    include("voigt.jl")
end
