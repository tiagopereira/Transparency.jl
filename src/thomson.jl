const σ_thomson = (e^4 / (6π * ε_0^2 * m_e^2 * c_0^4)) |> u"m^2"

"""
    α_thomson(electron_density::NumberDensity)

Compute the Thomson extinction as a function of electron density.
"""
α_thomson(electron_density::NumberDensity) = σ_thomson * electron_density
