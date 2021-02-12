using Interpolations
using SpecialFunctions: erfcx
using Test
using Transparency
using Unitful
import PhysicalConstants.CODATA2018: h, k_B, R_∞, c_0, m_e, e, ε_0

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
    end
    @testset "hminus_bf" begin
        # Values from the table, ensuring no stimulated emission
        λ = [500, 8000, 16000]u"Å"
        @test all(Transparency.hminus_bf_geltman.(λ, 0u"K", 1u"m^-3") ≈ [1.5e-22,
            3.92e-21, 1.7e-22]u"m^-1")
        @test Transparency.hminus_bf_geltman(0u"nm", 0u"K", 1u"m^-3") == 0u"m^-1"
        # Only testing against implementation (no table, only coefficients)
        @test all(Transparency.hminus_bf_john.(λ, 5000u"K", 1e24u"m^-3", 1e25u"m^-3") ≈ [
            2.5257481813577, 65.22804371400161, 1.8270928643449478]u"m^-1")
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
        @test all(hydrogenic_bf.(ν, ν / 1.1, 5000u"K", 1e22u"m^-3", 1., 2.) ≈ [
            2.4903105889694794, 9.569685175346825, 9.154320813938323]u"m^-1")
        @test hydrogenic_bf(ν[1], ν[1] * 1.1, 1u"K", 1u"m^-3", 1, 2.) == 0u"m^-1"
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
end

@testset "Broadening" begin
    line = AtomicLine(1.1u"aJ", 1.0u"aJ", 1.5u"aJ", 18, 8, 1.0, 1.0u"kg", 1)
    # Testing against implementation
    @test isapprox(γ_unsold_const(line), 1.1131993895644783e-15, rtol=1e-15)
    @test all(γ_unsold.(1.0, 1u"K", [1, 2]u"m^-3") ≈ [1, 2]u"s^-1")
    @test γ_unsold(1.0, 1000u"K", 1u"m^-3") ≈ (1000^0.3)u"s^-1"
    @test calc_Aji(1u"m", ustrip(ε_0 * m_e * c_0), 1 / ustrip(2π * e^2)) ≈ 1u"s^-1"
    @test calc_Bji(1u"m", 0u"Hz") ≈ 0u"m^3 / J"
    @test calc_Bji(1000u"nm", 1e9u"Hz") ≈ 8.396002689872053e-6u"m^3 / J"
    @test damping(1u"Hz", 1u"m", 1u"m") ≈ ustrip(1 / (4π * c_0))
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
    @testset "Intensity" begin
        dist = collect(LinRange(0, 100, 1001))u"m"
        ext = ones(1001)u"m^-1"
        S = ones(1001)u"kW / (m^2 * nm)"
        @test isapprox(calc_intensity(dist, ext, S), 0.9055913235757069u"kW / (m^2 * nm)",
                       atol=1e-3u"kW / (m^2 * nm)")
        @test_throws AssertionError calc_intensity(-dist, ext, S)
    end
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
end