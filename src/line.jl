"""
Computes line extinction and associated quantities.
"""

"""
    function calc_Aul(λ0::Unitful.Length, g_ratio::Real, f_value::AbstractFloat)

Compute the spontaneous deexcitation rate \$A_{ul}\$ (natural broadening)
for a bound-bound transition, using the SI expression *per wavelength*:

\$\$
A_{ul} = \\frac{2\\pi e^2}{\\varepsilon_0 m_e c} \\frac{g_l}{g_u} \\frac{f_{lu}}{\\lambda^2}
\$\$

for a given rest wavelength `λ0`, ration between statistical weights of lower and
upper levels (`g_ratio` = gl / gu), and `f_value` .
"""
function calc_Aul(λ0::Unitful.Length, g_ratio::Real, f_value::AbstractFloat)
    (2π * e^2 / (ε_0 * m_e * c_0) * g_ratio * f_value / λ0^2) |> u"s^-1"
end


"""
    function calc_Bul(λ0::Unitful.Length, Aul::Unitful.Frequency)

Compute the induced deexcitation rate \$B_{ul}\$ for a bound-bound transition,
using the SI expression *per wavelength*:

\$\$
B_{ul} = \\frac{\\lambda^5}{2 h c} A_{ul}
\$\$

for a given rest wavelength `λ0`, and spontaneous deexcitation rate `Aul.`
"""
calc_Bul(λ0::Unitful.Length, Aul::Unitful.Frequency) = (λ0^5 * Aul / (2h * c_0^2)) |> u"m^3 / J"


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
