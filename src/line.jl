"""
Computes line extinction and associated quantities.
"""

"""
    function AtomicLine(χu::Unitful.Energy{T}, χl::Unitful.Energy{T},
                        χ∞::Unitful.Energy{T}, gu::Int, gl::Int, f_value::T,
                        atom_weight::Unitful.Mass{T}, Z::Int)  where T <: AbstractFloat

Structure for atomic line.
"""
struct AtomicLine{T <: AbstractFloat}
    Aji::Unitful.Frequency{T}
    # Units of Bij/Bji defined for J_lambda
    Bji::Unitful.Quantity{T, Unitful.𝐋 * Unitful.𝐓^2 / Unitful.𝐌}
    Bij::Unitful.Quantity{T, Unitful.𝐋 * Unitful.𝐓^2 / Unitful.𝐌}
    λ0::Unitful.Length{T}
    χi::Unitful.Energy{T}
    χj::Unitful.Energy{T}
    # Properties of atom, not line, but keeping here for now
    χ∞::Unitful.Energy{T}
    atom_weight::Unitful.Mass{T}
    # Z was used incorrectly before. It should be the effective nuclear charge of
    # the upper level plus one. E.g. 1 for neutral, 2 for singly ionised.
    Z::Int
    function AtomicLine(χu::Quantity{T}, χl::Quantity{T}, χ∞::Quantity{T},
                        gu::Int, gl::Int, f_value::T, atom_weight::Unitful.Mass{T},
                        Z::Int)  where T <: AbstractFloat
        χu = wavenumber_to_energy(χu)
        χl = wavenumber_to_energy(χl)
        χ∞ = wavenumber_to_energy(χ∞)
        # Add conversion from cm^-1 to aJ, if type of χu is L^-1
        @assert χ∞ > χu
        @assert χu > χl
        @assert gu > 0
        @assert gl > 0
        @assert f_value > 0
        @assert atom_weight > 0u"kg"
        @assert Z >= 1
        λ0 = convert(Quantity{T, Unitful.𝐋}, ((h * c_0) / (χu - χl)) |> u"nm")
        Aul = convert(Quantity{T, Unitful.𝐓^-1}, calc_Aji(λ0, gl / gu, f_value))
        Bul = calc_Bji(λ0, Aul)
        Blu = gu / gl * Bul
        new{T}(Aul, Bul, Blu, λ0, χl, χu, χ∞, atom_weight, Z)
    end
end


"""
If input is in wavenumber, convert to energy. Otherwise keep as energy.
"""
function wavenumber_to_energy(a::Quantity{T}) where T <: AbstractFloat
    if typeof(a) <: PerLength
        a = convert(Unitful.Quantity{T, Unitful.𝐋^2 * Unitful.𝐓^-2 * Unitful.𝐌},
                    (h * c_0 * a) |> u"aJ")
    end
    @assert typeof(a) <: Unitful.Energy{T} "Input units must either be wavenumber or energy"
    return a
end


"""
Compute line extinction given an `AtomicLine` struct, `profile` defined per wavelength,
and upper and lower population densities `n_u` and `n_l`.
"""
function αline_λ(line::AtomicLine, profile::PerLength, n_u::NumberDensity, n_l::NumberDensity)
    (h * c_0 / (4 * π * line.λ0) * profile * (n_l * line.Bij - n_u * line.Bji)) |> u"m^-1"
end


"""
Compute line emissivity given an `AtomicLine` struct, `profile` defined per wavelength,
and upper population density `n_u`.
"""
function jline_λ(line::AtomicLine, profile::PerLength, n_u::NumberDensity)
    (h * c_0 / (4 * π * line.λ0) * n_u * line.Aji * profile) |> u"W / (m^3 * nm)"
end


"""
Calculates the Blackbody (Planck) function per wavelength, for given
wavelength and temperature.
"""
function blackbody_λ(λ::Unitful.Length, temperature::Unitful.Temperature)
    (2h * c_0^2) / (λ^2 * λ^3 * (exp((h * c_0 / k_B) / (λ * temperature)) - 1))
end


"""
Calculates the Blackbody (Planck) function per frequency, for a given
frequency and temperature.
"""
function blackbody_ν(ν::Unitful.Frequency, temperature::Unitful.Temperature)
    (2h / c_0^2) * (ν^3 / (exp((h / k_B) * (ν / temperature)) - 1))
end
