"""
Set of recipes to compute opacities.
"""
module Transparency
    export hminus_ff, hminus_bf, hydrogenic_ff, hydrogenic_bf
    export h2minus_ff, h2plus_ff, h2plus_bf
    export rayleigh_h2, rayleigh_h, thomson
    export humlicek, voigt_profile, dispersion_profile
    export calc_Aji, calc_Bji, damping, doppler_width, γ_unsold_const, γ_unsold
    export AtomicLine, αline_λ, jline_λ, blackbody_λ, blackbody_ν, calc_intensity
    export coll_CE, coll_CI, coll_Ω
    export coll_exc_hydrogen_johnson, coll_ion_hydrogen_johnson
    export CE_RH_hydrogen, CI_RH_hydrogen

    using Unitful
    using Interpolations
    import SpecialFunctions: expint
    import NumericalIntegration: integrate, cumul_integrate, TrapezoidalFast
    import PhysicalConstants.CODATA2018: h, k_B, R_∞, c_0, m_e, m_u, e, ε_0, a_0

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
    include("hydrogen.jl")
    include("thomson.jl")
    include("voigt.jl")
end
