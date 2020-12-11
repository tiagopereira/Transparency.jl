"""
Computes line profiles and associated quantities.
"""

const invSqrtPi = 1. / sqrt(π)

""" 
Compute scaled complex complementary error function using 
[Humlicek (1982, JQSRT 27, 437)](https://doi.org/10.1016/0022-4073(82)90078-4) 
W4 rational approximations. 
Here, z is defined as z = v + i * a, and returns w(z) = H(a,v) + i * L(a, v).
"""
function humlicek(z::Complex)
    s = abs(real(z)) + imag(z)
    if s > 15.0
        # region I
        w = im * invSqrtPi * z / (z * z - 0.5)
    elseif s > 5.5
        # region II
        zz = z * z
        w = im * (z * (zz * invSqrtPi - 1.4104739589)) / (0.75 + zz * (zz - 3.0))
    else
        x, y = real(z), imag(z)
        t = y - im * x
        if y >= 0.195 * abs(x) - 0.176
            # region III
            w = ((16.4955 + t * (20.20933 + t * (11.96482 + t * (3.778987 + 0.5642236 * t))))
               / (16.4955 + t * (38.82363 + t * (39.27121 + t * (21.69274 + t * (6.699398 + t))))))
        else
            # region IV
            u = t * t
            nom = t * (36183.31 - u * (3321.99 - u * (1540.787 -  u * 
                   (219.031 - u * (35.7668 - u * (1.320522 - u * .56419))))))
            den = 32066.6 - u * (24322.8 - u * (9022.23 - u * (2186.18 - 
                    u * (364.219 - u * (61.5704 - u * (1.84144 - u))))))
            w = exp(u) - nom / den
        end
    end
    return w
end


"""
Compute the normalised Voigt profile for wavelength, given a damping constant a, 
dimensionless wavelength v = (λ - λ0) / ΔλD, and Doppler width ΔλD. 
Uses Humlicek's W4 approximation.
"""
function voigt_profile(a::AbstractFloat, v::AbstractFloat, ΔλD::Unitful.Length)
    z = v + a * im
    profile = real(humlicek(z))
    return profile * invSqrtPi / ΔλD
end


"""
Compute the normalised dispersion (or Faraday) profile for wavelength, given a 
damping constant a, dimensionless wavelength v = (λ - λ0) / ΔλD, 
and Doppler width ΔλD.
Uses Humlicek's W4 approximation.
"""
function dispersion_profile(a::AbstractFloat, v::AbstractFloat, ΔλD::Unitful.Length)
    z = v + a * im
    profile = imag(humlicek(z))
    return profile * invSqrtPi / ΔλD
end


"""
Compute Doppler width for wavelength ΔλD, given a wavelength λ0, 
mass of atom/ion, and temperature.
"""
function doppler_width(λ0::Unitful.Length, mass::Unitful.Mass, temperature::Unitful.Temperature)
    (λ0 / c_0 * sqrt(2 * k_B * temperature / mass)) |> u"nm"
end