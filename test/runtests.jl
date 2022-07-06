using Interpolations
using SpecialFunctions: erfcx
using Test
using Transparency
using Unitful
import PhysicalConstants.CODATA2018: h, k_B, R_∞, c_0, m_e, e, ε_0, m_u

@testset "Thomson" begin
    @test thomson(1.e29u"m^-3") ≈ 6.652458732173518u"m^-1"
    @test thomson(0u"m^-3") == 0.0u"m^-1"
    @test_throws MethodError thomson(1)
    @test_throws MethodError thomson(1u"m")
end

@testset "Hydrogen" begin
    @testset "hminus_ff" begin
        θ = [0.5, 1, 2]
        temp = 5040u"K" ./ θ
        λ = [91130, 11390, 3038]u"Å"
        # Testing against values in the original table
        @test all(Transparency.hminus_ff_stilley.(
            λ, temp, 1u"m^-3", 1u"Pa" ./ (k_B .* temp)) ≈ [2.75e-28, 8.76e-30, 1.45e-30]u"m^-1")
        # Only testing against implementation (no table, only coefficients)
        @test all(Transparency.hminus_ff_john.(
            λ, temp, 1u"m^-3", 1u"Pa" ./ (k_B .* temp)) ≈ [2.7053123415611454e-28,
                8.7520478601176e-30, 1.6977098265784962e-30]u"m^-1")
        @test hminus_ff.(λ, temp, 1u"m^-3", 1u"m^-3", recipe="stilley") ≈
            Transparency.hminus_ff_stilley.(λ, temp, 1u"m^-3", 1u"m^-3")
        @test hminus_ff.(λ, temp, 1u"m^-3", 1u"m^-3", recipe="john") ≈
            Transparency.hminus_ff_john.(λ, temp, 1u"m^-3", 1u"m^-3")
        @test_throws "ErrorException" hminus_ff.(λ, temp, 1u"m^-3", 1u"m^-3", recipe="aaa")
    end
    @testset "hminus_bf" begin
        # Values from the table, ensuring no stimulated emission
        λ = [18, 121, 164, 500, 850, 1600]u"nm"  # three values from each table
        @test all(Transparency.hminus_bf_wbr.(λ, 0u"K", 1u"m^-3") ≈
            [6.7e-24, 5.43e-22, 7.29e-22, 2.923e-21, 4.001e-21, 1.302e-22]u"m^-1")
        @test Transparency.hminus_bf_wbr(0u"nm", 0u"K", 1u"m^-3") == 0u"m^-1"
        λ = [500, 8000, 16000]u"Å"
        @test all(Transparency.hminus_bf_geltman.(λ, 0u"K", 1u"m^-3") ≈ [1.5e-22,
            3.92e-21, 1.7e-22]u"m^-1")
        @test Transparency.hminus_bf_geltman(0u"nm", 0u"K", 1u"m^-3") == 0u"m^-1"
        # Only testing against implementation (no table, only coefficients)
        @test all(Transparency.hminus_bf_geltman.(λ, 5000u"K", 1e24u"m^-3", 1e25u"m^-3") ≈ [
            2.5276368595110785, 64.24516492134173, 2.3904062765360763]u"m^-1")
        @test all(Transparency.hminus_bf_john.(λ, 5000u"K", 1e24u"m^-3", 1e25u"m^-3") ≈ [
            2.5257481813577, 65.22804371400161, 1.8270928643449478]u"m^-1")
        temp = [2000, 5000, 10000]u"K"
        @test Transparency.calc_hminus_density.(1e13u"m^-3", temp, 1e14u"m^-3") ≈
            [91.94566524525405, 1.6850912396740527, 0.24835857088939245]u"m^-3"
        @test hminus_bf.(λ, temp, 1u"m^-3", 1u"m^-3", recipe="geltman") ≈
            Transparency.hminus_bf_geltman.(λ, temp, 1u"m^-3", 1u"m^-3")
        @test hminus_bf.(λ, temp, 1u"m^-3", 1u"m^-3", recipe="john") ≈
            Transparency.hminus_bf_john.(λ, temp, 1u"m^-3", 1u"m^-3")
        @test_throws "ErrorException" hminus_bf.(λ, temp, 1u"m^-3", 1u"m^-3", recipe="aaa")
    end
    @testset "hydrogenic" begin
        t1 = (3.28805e15u"Hz" * h / k_B) |> u"K"
        temp = [t1 * 1e3, t1, t1 * 1e-2]
        ν = k_B / h .* [temp[1] * 1e-4, temp[2] * 1e-1, temp[3] * 1e1]
        # Test a few points from Kurucz's table
        @test all(Transparency.gaunt_ff.(ν, temp, 1) .≈ [5.53, 1.8, 1.08])
        @test Transparency.gaunt_ff(500u"nm", 5000u"K", 1) == Transparency.gaunt_ff(
            c_0 / 500u"nm", 5000u"K", 1)
        @test Transparency.gaunt_bf(500u"nm", 2, 5) == Transparency.gaunt_bf(c_0 / 500u"nm", 2, 5)
        # Note: gaunt_bf only tested against RH version, which is < 0 in some cases!
        @test Transparency.gaunt_bf(500u"nm", 1, 1.1) ≈ 0.023132239149173173
        @test Transparency.gaunt_bf(100u"nm", 2, 5) ≈ 1.0517607988011628
        @test_throws AssertionError Transparency.gaunt_bf(500u"nm", 1, 1)
        ν = [0.15, 0.6, 300]u"PHz"
        # Only testing against implementation (no table, only expression)
        @test all(hydrogenic_ff.(ν, 5000u"K", 1e24u"m^-3", 1e24u"m^-3", 1) ≈ [
            134.2760926064222, 2.662541955122332, 2.0317078183858428e-8]u"m^-1")
        @test all(hydrogenic_ff_σ.(ν, 5000u"K", 1) ≈ [
            1.759801422878951e-46, 2.6709661157204164e-48, 2.0317078183858426e-56]u"m^5")
        @test all(hydrogenic_bf.(ν, ν / 1.1, 5000u"K", 1e22u"m^-3", 1., 2.) ≈ [
            2.4903105889694794, 9.569685175346825, 9.154320813938323]u"m^-1")
        @test hydrogenic_bf(ν[1], ν[1] * 1.1, 1u"K", 1u"m^-3", 1, 2.) == 0u"m^-1"
        @test hydrogenic_bf(ν[1], ν[1] * 1.1, 1u"K", 1u"m^-3", 1, 2.) == 0u"m^-1"
        @test hydrogenic_bf_σ_scaled(1u"m^2", ν[3], ν[3] * 1.1, 1., 1.) == 0u"m^2"
        @test hydrogenic_bf_σ_scaled(42u"m^2", ν[3], ν[3], 1., 1.) == 42u"m^2"
        @test hydrogenic_bf_σ_scaled(1u"m^2", c_0 / ν[3], c_0 /ν[3], 2., 5.) ≈
                 hydrogenic_bf_σ_scaled(1u"m^2", ν[3], ν[3], 2., 5.)
    end
    @testset "h2minus" begin
        # Test a few points from the table
        λ = [350.5, 1139.1, 15_188.3]u"nm"
        temp = 5040u"K" ./ [0.5, 1.6, 3.6]
        @test all(h2minus_ff.(λ, temp, 1e31u"m^-3", 1u"Pa" ./ (k_B .* temp)) ≈ [
            4.17, 80.60, 14700]u"m^-1")
    end
    @testset "h2plus" begin
        λ = 1f7u"nm" ./ [500, 4000, 26_000]
        temp = [12_000, 5000, 2500]u"K"
        # Test opacity table directly, sum of bf and ff
        κ_total = h2plus_ff.(λ, temp, 1f25u"m^-3", 1f24u"m^-3") .+
                  h2plus_bf.(λ, temp, 1f25u"m^-3", 1f24u"m^-3")
        @test all(κ_total ≈ [0.226, 1.35, 114.]u"m^-1")
        # Test bf fraction table
        @test all(h2plus_bf.(λ, temp, 1f25u"m^-3", 1f24u"m^-3") ./
                        κ_total ≈ [0.006, 0.273, 0.996])
    end
    @testset "rayleigh" begin
        # Testing against implementation
        @test rayleigh_h2(100u"nm", 1u"m^-3") == 0.0u"m^-1"
        λ = [200, 600, 2000]u"nm"
        @test all(rayleigh_h2.(λ, 1e30u"m^-3") ≈ [
            8.565893085787453, 0.0747454429941088, 0.0005507634570312499]u"m^-1")
        @test rayleigh_h(100u"nm", 1u"m^-3") == 0.0u"m^-1"
        @test all(rayleigh_h.(λ, 1e30u"m^-3") ≈ [
            6.946808203124999, 0.04804975738502133, 0.00036536185226953123]u"m^-1")
    end
end

@testset "Voigt" begin
    a = im .* 10 .^LinRange(-4, -0.3, 20)'
    v = 10 .^LinRange(-1, 2, 20)'
    z = a' .+ hcat(-v, v)
    @testset "humlicek" begin
        # Test against more precise erfcx function
        @test isapprox(humlicek.(z), erfcx.(-im .* z), rtol=1e-4)
        @test humlicek(im * 0) == 1.
        @test_throws MethodError humlicek(0.)
    end
    @testset "profiles" begin
        wave = 1u"nm"
        @test ustrip(voigt_profile(0.1, 0.5, wave)) ≈ real(humlicek(0.5 + 0.1 *im)) / sqrt(π)
        @test ustrip(dispersion_profile(0.1, 0.5, wave)) ≈ imag(humlicek(0.5 + 0.1 *im)) / sqrt(π)
        # Symmetry / anti-symmetry
        @test voigt_profile.(0.1, v, wave) == voigt_profile.(0.1, -v, wave)
        @test dispersion_profile.(0.1, v, wave) == -dispersion_profile.(0.1, -v, wave)
    end
    @test doppler_width(ustrip(c_0) * u"m", 1u"kg",  ustrip(0.5 / k_B) * u"K") ≈ 1u"m"
    @test doppler_width(ustrip(c_0) / u"s", 1u"kg",  ustrip(0.5 / k_B) * u"K") ≈ 1u"Hz"
end

@testset "Broadening" begin
    line = AtomicLine(1.1u"aJ", 1.0u"aJ", 1.5u"aJ", 18, 8, 1.0, 1.0u"kg", 1)
    @testset "van der Waals" begin
        # Testing against implementation
        @test const_unsold(line) ≈ 1.1131993895644783e-15 rtol=1e-10
        @test γ_unsold.(1.0, 1u"K", [1, 2]u"m^-3") ≈ [1, 2]u"s^-1"
        @test γ_unsold(1.0, 1000u"K", 1u"m^-3") ≈ (1000^0.3)u"s^-1"
        @test const_barklem(m_u * 1, 0.3, 300) ≈ 1.0853795252714703e-15u"m^3 / s"
        @test const_barklem(m_u * 1, 1, 2) == 2 * const_barklem(m_u * 1, 1, 1)
        @test_throws ErrorException const_barklem(m_u * 1, 1, -1)
        @test_throws ErrorException const_barklem(m_u * 1, -1, 1)
        @test γ_barklem(0.3, 123u"m^3/s", 1u"K", 1u"m^-3") == 123.0u"s^-1"
        @test γ_barklem(0.3, 1e-16u"m^3/s", 6000u"K", 1e23u"m^-3") ≈ 2.100646320154e8u"s^-1"
        @test const_deridder_rensbergen(m_u * 1, m_u * 1, 3, 1) ≈ 3*2e-14u"m^3/s"
        @test_throws ErrorException const_deridder_rensbergen(m_u * 1, m_u * 1, -1, 1)
        @test γ_deridder_rensbergen(0.5, 2e-14u"m^3/s", 6000u"K", 1e22u"m^-3") ≈
            1.5491933384829668e10u"s^-1"
        @test γ_deridder_rensbergen(0.5, 1u"m^3/s", 1e3u"K", 1u"m^-3") ≈ sqrt(1e3)u"s^-1"
    end
    @testset "Linear Stark" begin
        # Test against Sutton formula
        @test γ_linear_stark.([0., 1.]u"m^-3", 3, 1) ≈ [0, 0.000408]u"s^-1"
        tmp = 0.00016371u"s^-1"
        @test γ_linear_stark.([1, 1e20]u"m^-3", 3, 2) ≈ [tmp, tmp * (1e20)^(2/3)]
        @test_throws AssertionError γ_linear_stark(1.0u"m^-3", 1, 1)
        @test_throws AssertionError γ_linear_stark(1.0u"m^-3", 1, 0)
    end
    @testset "Quadratic Stark" begin
        # Testing against implementation
        @test Transparency.c4_traving(line) ≈ 3.744741607310466e-23u"m^4 / s"
        @test const_quadratic_stark(line) ≈ 2.7236711602117037e-13u"m^3 / s"
        @test const_quadratic_stark(line; scaling=2) ≈ const_quadratic_stark(line) * 2^(2/3)
        @test γ_quadratic_stark(1.2345u"m^-3", 0u"K") ≈ 1.2345u"s^-1"
        @test γ_quadratic_stark(1e10u"m^-3", 10000u"K") ≈ 1e10u"s^-1" * 10000^(1/6)
        # Testing against implementation
        temp = [5000, 10000]u"K"
        @test (γ_quadratic_stark_gray.(1e22u"m^-3", temp, 1e-20u"m^4/s") ≈
                [5.709239783376956e11, 6.40840498153864e11]u"s^-1")
    end
    @testset "Radiation quantities" begin
        @test calc_Aji(1u"m", ustrip(ε_0 * m_e * c_0), 1 / ustrip(2π * e^2)) ≈ 1u"s^-1"
        @test calc_Bji(1u"m", 0u"Hz") ≈ 0u"m^3 / J"
        @test calc_Bji(1000u"nm", 1e9u"Hz") ≈ 8.396002689872053e-6u"m^3 / J"
        @test damping(1u"Hz", 1u"m", 1u"m") ≈ ustrip(1 / (4π * c_0))
    end
end

@testset "Line" begin
    line = AtomicLine(1.39728917u"aJ", 1.0u"aJ", 1.5u"aJ", 1, 1, 1.0, 1.0u"kg", 1)
    @testset "AtomicLine" begin
        @test line.λ0 ≈ 500u"nm"
        @test_throws MethodError AtomicLine(1.1u"aJ", 1u"aJ", 1.5u"aJ", 1, 1, 1.0, 1.0u"kg", 1)
        @test_throws AssertionError AtomicLine(.0u"aJ", .0u"aJ", .5u"aJ", 1, 1, 1.0, .1u"kg", 1)
    end
    @test αline_λ(line, 1u"m^-1", 1u"m^-3", 1u"m^-3") == 0u"m^-1"
    # Testing against implementation
    @test αline_λ(line, 1u"nm^-1", 1e17u"m^-3", 1e18u"m^-3") ≈ 1.9918846554254643u"m^-1"
    @test jline_λ(line, 1u"m^-1", 1e21u"m^-3") ≈ 8.435271743930054u"W / (m^3 * nm)"
    @testset "Blackbody" begin
        λ = 500.0u"nm"
        @test blackbody_λ(λ, 5000u"K") ≈ 12.107190590398108u"kW / (m^2 * nm)"
        @test blackbody_ν(c_0 / λ, 5000u"K") ≈ 10.096310186694323u"nW / (m^2 * Hz)"
        @test blackbody_λ(λ, 5000u"K") ≈ blackbody_ν(c_0 / λ, 5000u"K") * c_0 / λ^2
    end
    @test Transparency.wavenumber_to_energy(ustrip(1 / (h * c_0)) * u"m^-1") ≈ 1u"J"
end

@testset "Collisions" begin
    @testset "Johnson" begin
        ne = 1e20u"m^-3"
        @test_throws AssertionError coll_exc_hydrogen_johnson(2, 1, ne, 1u"K")
        @test_throws AssertionError coll_exc_hydrogen_johnson(-1, 1, ne, 1u"K")
        @test_throws AssertionError coll_ion_hydrogen_johnson(-1, ne, 1u"K")
        @test coll_exc_hydrogen_johnson(1, 2, ne, 1u"K") ≈ 0.0u"s^-1"
        @test coll_ion_hydrogen_johnson(1, ne, 1u"K") ≈ 0.0u"s^-1"
        @test all(Transparency._bn.([1, 2, 3]) ≈ [-0.603, 0.116875, 0.25876543209876557])
        @test all(Transparency._rn.([1, 2, 3]) ≈ [0.45, 1.94 * 2^-1.57, 1.94 * 3^-1.57])
        @test Transparency.ξ(1e20) ≈ 0.0
        @test Transparency.ξ(0.1) ≈ 6.125071285714834
        # Testing against implementation
        temp = [3000, 5000, 10000]u"K"
        @test all(coll_exc_hydrogen_johnson.(1, 2, ne, temp) ≈
            [3.930707156378747e-11, 0.0002263838629287018, 24.38373538697928]u"s^-1")
        @test all(coll_exc_hydrogen_johnson.(1, 4, ne, temp) ≈
            [1.3966268525286798e-16, 4.2003964667891764e-8, 0.08902629232613525]u"s^-1")
        @test all(coll_exc_hydrogen_johnson.(2, 3, ne, temp) ≈
            [41449.22586174524, 712944.3194152612, 6.358527475844925e6]u"s^-1")
        @test all(coll_ion_hydrogen_johnson.(1, ne, temp) ≈
            [2.0663760818067978e-18, 3.980869050632006e-9, 0.04717372440180093]u"s^-1")
        @test all(coll_ion_hydrogen_johnson.(2, ne, temp) ≈
            [5.6897483527787776, 1746.3754923299948, 166085.00320954874]u"s^-1")
        @test all(coll_ion_hydrogen_johnson.(3, ne, temp) ≈
            [35134.766878609924, 592137.0862432675, 6.093655265672546e6]u"s^-1")
        @test all(CE_RH_hydrogen.(1, [2, 3, 4], 5000u"K") ≈ [6.098416316523097e-16,
            1.151527001714162e-16, 4.20364505409004e-17]u"m^3 / (K^(1/2) * s)")
        @test all(CI_RH_hydrogen.([1, 2, 3], 5000u"K") ≈ [2.86396932776755e-17,
            6.595860989488697e-16, 2.791765778377678e-15]u"m^3 / (K^(1/2) * s)")
    end
    @testset "RH rates" begin
        temp = [3000, 5000, 10000]u"K"
        data1 = [1, 1, 1]u"m^3 / (K^(1/2) * s)"
        data2 = [1, 1, 1]
        interp1 = LinearInterpolation(temp, data1)
        interp2 = LinearInterpolation(temp, data2)
        # Using wrong units of data:
        @test_throws MethodError coll_CE(interp2, 1, 1u"m^-3", 5000u"K")
        @test_throws MethodError coll_CI(interp2, 1u"J", 1u"m^-3", 5000u"K")
        @test_throws MethodError coll_Ω(interp1, 1, 1u"m^-3", 5000u"K")
        @test_throws AssertionError coll_CI(interp1, -1u"J", 1u"m^-3", 5000u"K")
        @test coll_CE(interp1, 1, 1u"m^-3", 10000u"K") ≈ 100.0u"s^-1"
        @test coll_CI(interp1, 1e-50u"J", 1u"m^-3", 10000u"K") ≈ 100.0u"s^-1"
        @test coll_Ω(interp2, 1, 1u"m^-3"/8.629132180819955e-14, 10000u"K") ≈ 1.0u"s^-1"
    end
    @testset "Przybilla & Butter" begin
        temp = [5000, 10000, 30000, 250000]u"K"
        ne = 1u"m^-3"
        c0 = sqrt.(temp) ./ Transparency.Ω_c0
        @test_throws AssertionError coll_deexc_hydrogen_PB04(1, 8, 1, ne, temp[1])
        @test_throws AssertionError coll_deexc_hydrogen_PB04(3, 1, 1, ne, temp[1])
        # Compare with a few random values from the original table
        @test isapprox(ustrip.(coll_deexc_hydrogen_PB04.(1, 2, 1, ne, temp) .* c0),
                       [0.698, 0.809, 1.15, 3.95], atol=1e-3)
        @test isapprox(ustrip.(coll_deexc_hydrogen_PB04.(1, 4, 1, ne, temp) .* c0),
                       [0.102, 0.122, 0.228, 0.488], atol=1e-3)
        @test isapprox(ustrip.(coll_deexc_hydrogen_PB04.(2, 3, 1, ne, temp) .* c0),
                       [27.8, 33.8, 62.0, 252], atol=1e-3)
        @test isapprox(ustrip.(coll_deexc_hydrogen_PB04.(2, 7, 1, ne, temp) .* c0),
                       [7.26, 9.27, 11.4, 11.9], atol=1e-3)
        @test isapprox(ustrip.(coll_deexc_hydrogen_PB04.(4, 5, 1, ne, temp) .* c0),
                       [817, 1350, 3400, 10000], atol=1e-3)
    end
end

@testset "Formal solvers" begin
    @testset "Weights" begin
        @test all(Transparency._w2(60.) .== (1, 1))
        @test all(Transparency._w2(1.) .≈ (1-exp(-1), 1 - 2*exp(-1)))
        @test all(Transparency._w2(1f-6) .≈ (9.999995f-7, 4.9999967f-13))
    end
    @testset "Piecewise" begin
        # Constant source function
        z = collect(LinRange(1, 1e6, 20))u"m"
        alpha = ones(20) * 1e-20u"m^-1"
        S = ones(20)*u"kW / (nm * m^2)"
        @test piecewise_1D_linear(z, alpha, S) ≈ S
        @test piecewise_1D_nn(z, alpha, S) ≈ S
        alpha = ones(20) * 1e-1u"m^-1"
        @test (calc_intensity_brute_force(z, alpha, S) * 2 ≈
               calc_intensity_brute_force(z, alpha, S*2))
        @test piecewise_1D_linear(z, alpha, S) ≈ S
        @test piecewise_1D_linear(z, alpha, S;
                                  initial_condition=:zero)[[1, end]] ≈ [S[1], S[1]*0]
        @test piecewise_1D_nn(z, alpha, S) ≈ S
        @test piecewise_1D_nn(z, alpha, S;
                              initial_condition=:zero)[[1, end]] ≈ [S[1], S[1]*0]
        # Linear extinction and source function, test reversibility
        alpha = collect(LinRange(1e-3, 1e-5, 20)u"m^-1")
        S = collect(LinRange(1, 100, 20)u"kW / (nm * m^2)")
        @test (piecewise_1D_linear(z, reverse(alpha), reverse(S); to_end=true) ≈
               reverse(piecewise_1D_linear(z, alpha, S)))
        @test (piecewise_1D_nn(z, reverse(alpha), reverse(S); to_end=true) ≈
               reverse(piecewise_1D_nn(z, alpha, S)))
        # Exceptions
        @test_throws ErrorException piecewise_1D_linear(z, alpha, S; initial_condition=:aaa)
        @test_throws ErrorException piecewise_1D_nn(z, alpha, S; initial_condition=:aaa)
    end
    @testset "Feautrier" begin
        z = collect(LinRange(2e6, -1e5, 20))u"m"
        alpha = 1e-5 * ones(20)u"m^-1"
        S = zeros(20)u"kW / (nm * m^2)"
        # Simple tests
        @test feautrier(z, alpha, S) ≈ S
        S = 100 * ones(20)u"kW / (nm * m^2)"
        # Against implementation
        @test feautrier(z, alpha, S)[1] * 2 ≈ 106.65292755967045u"kW / (nm * m^2)"
        # Test reversibility
        @test feautrier(z, alpha, S)[1] ≈ feautrier(z, reverse(alpha), reverse(S))[end]
    end
end
