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
        λ = 1e7u"nm" ./ [500, 4000, 26_000]
        temp = [12_000, 5000, 2500]u"K"
        # Test opacity table directly, sum of bf and ff
        κ_total = h2plus_ff.(λ, temp, 1e49u"m^-3", 1u"m^-3") .+
                  h2plus_bf.(λ, temp, 1e49u"m^-3", 1u"m^-3")
        @test all(κ_total ≈ [0.226, 1.35, 114.]u"m^-1")
        # Test bf fraction table
        @test all(h2plus_bf.(λ, temp, 1e49u"m^-3", 1u"m^-3") ./ 
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