# Koen van der Heijden, November - December 2023

module BornModel

using Parameters
using Statistics
using Random

function heaviside(x, k = 1e5) # was 1e5)
    return (1/2) * (1 + tanh(k*x))
end 

#= 
--------------------------------------------------------------------------
Dimensional Model
--------------------------------------------------------------------------
=# 

function density(T, S)
    # following the definition from Born (2014)
    α = 0.11  # kg/(m3 K)
    β = 0.77  # kg/(m3 psu)
    return -α*T + β*S
end

function calculate_Ud(T3, S3, p)
    @unpack α, β, g, w, d, h, f, ρ0, c1, c2, τ = p.cons
    @unpack F, Ubtp, T2, S2, T4, S4, T0 = p

    return Ubtp .- (g*d)/(2*f*ρ0*w) .* (-α.*(T4 .- T3) + β.*(S4 .- S3))
end

function calculate_Us(Ud, T1, S1, p)
    @unpack α, β, g, w, d, h, f, ρ0, c1, c2, τ = p.cons
    @unpack F, Ubtp, T2, S2, T4, S4, T0 = p

    return Ud  .- (g*h)/(2*f*ρ0*w) .* (-α.*(T2 .- T1) + β.*(S2 .- S1))
end

# function calculate_M(Ud, Us, p)
#     @unpack α, β, g, w, d, h, f, ρ0, c1, c2, τ = p.cons
#     @unpack F, Ubtp, T2, S2, T4, S4, T0 = p

#     return Us*w*h + Ud*w*d
# end

function born!(dx, x, p, t = nothing)
    @unpack α, β, g, w, d, h, f, S0, ρ0, c1, c2, τ = p.cons
    @unpack F, Ubtp, T2, S2, T4, S4, T0 = p
    
    T1, S1, T3, S3 = x

    σ1 = density(T1, S1)
    σ2 = density(T2, S2)
    σ3 = density(T3, S3)
    σ4 = density(T4, S4)

    Ud = Ubtp - (g*d)/(2*f*ρ0*w) * (σ4 - σ3)
    Us = Ud   - (g*h)/(2*f*ρ0*w) * (σ2 - σ1)
 
    dx[1] = c1*Us*(T2 - T1) + (d/(h+d))*c2*heaviside(σ1 - σ3)*(T3 - T1) + (1/τ)*(T0 - T1)
    dx[2] = c1*Us*(S2 - S1) + (d/(h+d))*c2*heaviside(σ1 - σ3)*(S3 - S1) - F*(S0/h)
    dx[3] = c1*Ud*(T4 - T3) - (h/(h+d))*c2*heaviside(σ1 - σ3)*(T3 - T1)
    dx[4] = c1*Ud*(S4 - S3) - (h/(h+d))*c2*heaviside(σ1 - σ3)*(S3 - S1)

    dx
end

function seasonal_born!(dx, x, p_seas, t = nothing)
    @unpack α, β, g, w, d, h, f, S0, ρ0, c1, c2, τ, T = p_seas.cons
    @unpack F, Ubtp, T2, S2, T4, S4, T0, Tamp = p_seas
    
    T1, S1, T3, S3 = x

    σ1 = density(T1, S1)
    σ2 = density(T2, S2)
    σ3 = density(T3, S3)
    σ4 = density(T4, S4)

    Ud = Ubtp - (g*d)/(2*f*ρ0*w) * (σ4 - σ3)
    Us = Ud   - (g*h)/(2*f*ρ0*w) * (σ2 - σ1)

    dx[1] = c1*Us*(T2 - T1) + (d/(h+d))*c2*heaviside(σ1 - σ3)*(T3 - T1) + (1/τ)*(T0 - Tamp*cos((2*π*t)/T) - T1)
    dx[2] = c1*Us*(S2 - S1) + (d/(h+d))*c2*heaviside(σ1 - σ3)*(S3 - S1) - F*(S0/h)
    dx[3] = c1*Ud*(T4 - T3) - (h/(h+d))*c2*heaviside(σ1 - σ3)*(T3 - T1)
    dx[4] = c1*Ud*(S4 - S3) - (h/(h+d))*c2*heaviside(σ1 - σ3)*(S3 - S1) 

    dx
end

#= 
--------------------------------------------------------------------------
Non-Dimensionalized Model
--------------------------------------------------------------------------
=# 

function nondimensional_born!(dx, x, p, t = nothing)
    @unpack T2, S2, T4, S4, T0, Tamp, η, μ1, μ2, μ3, μ4, r = p

    T1, S1, T3, S3 = x

    Ud = 1 - η*((S4 - 1) - (S3 - T3))
    Us = 1 - η*((S4 - 1) - (S3 - T3)) - η*r*((S2 - T2) - (S1 - T1))

    dx[1] = μ1*Us*(T2 - T1) + μ2*(T3 - T1)*heaviside((S1 - T1) - (S3 - T3)) + μ3*(T0 - T1)
    dx[2] = μ1*Us*(S2 - S1) + μ2*(S3 - S1)*heaviside((S1 - T1) - (S3 - T3)) - μ4
    dx[3] = μ1*Ud*(T4 - T3) - μ2*r*(T3 - T1)*heaviside((S1 - T1) - (S3 - T3))
    dx[4] = μ1*Ud*(S4 - S3) - μ2*r*(S3 - S1)*heaviside((S1 - T1) - (S3 - T3))

    dx    
end 

function nondimensional_seasonal_born!(dx, x, p, t = nothing) 
    @unpack T2, S2, T4, S4, T0, Tamp, η, μ1, μ2, μ3, μ4, r  = p

    T1, S1, T3, S3 = x

    Ud = 1 - η*((S4 - 1) - (S3 - T3))
    Us = 1 - η*((S4 - 1) - (S3 - T3)) - η*r*((S2 - T2) - (S1 - T1))

    dx[1] = μ1*Us*(T2 - T1) + μ2*(T3 - T1)*heaviside((S1 - T1) - (S3 - T3)) + μ3*(T0 - (Tamp * cos(2 * π * t)) - T1)
    dx[2] = μ1*Us*(S2 - S1) + μ2*(S3 - S1)*heaviside((S1 - T1) - (S3 - T3)) - μ4
    dx[3] = μ1*Ud*(T4 - T3) - μ2*r*(T3 - T1)*heaviside((S1 - T1) - (S3 - T3))
    dx[4] = μ1*Ud*(S4 - S3) - μ2*r*(S3 - S1)*heaviside((S1 - T1) - (S3 - T3))

    dx
end 

### conversions to nondimensional form

function nondimensional_T(T)
    """
    Calculates the nondimensional temperature from a given temperature (in degrees Celsius)
    """
    return (T .+ 273.15)/(4 .+ 273.15)
end

function nondimensional_S(S)
    """
    Calculates the nondimensional salinity from a given salinity (in unit 1 / psu)
    """
    return (0.77 .* S)/(0.11 .* 277.15)
end

function nondimensional_M(M)
    """
    Calculates the nondimensional transport from a given transport (in Sv)
    """
    return M*(1e6)/(0.133*1400*100e3)

end

function nondimensional_U(U)
    """
    Calculates the nondimensional velocity from a given velocity (in m/s)
    """
    return U / 0.133
end

function nondimensional_F(F)
    """
    Calculates μ4 from F
    """
    return (F*35*0.77) / (100*0.11*277.15)
end

function nondimensional_σ(σ)
    """
    Calculates nondimensional σ from dimensional σ (in kg/m3)
    """
    return σ/(0.11*277.15)
end

### conversions to dimensional form 

function dimensional_T(T)
    """
    Calculates the salinity (in unit 1 / psu) from a given nondimensional salinity
    """
    return T.*(4 .+ 273.15) .- 273.15
end

function dimensional_S(S)
    """
    Calculates the temperature (in degrees Celsius) from a given nondimensional temperature
    """
    return (0.11 .* 277.15 .* S)/0.77
end

function dimensional_M(M)
    """
    Calculates the transport (in Sv) from a given nondimensional transport M
    """
    return M*(0.133*1400*100e3)/(1e6)

end

function dimensional_U(U)
    """
    Calculates the velocity (in m/s) from a given nondimensional velocity U
    """
    return U * 0.133
end

function dimensional_F(μ4)
    """
    Calculates F from μ4
    """
    return (μ4*100*0.11*277.15) / (35*0.77)
end

function dimensional_σ(σ)
    """
    Calculates dimensional σ (in kg/m3) from nondimensional σ
    """
    return σ*((0.11*277.15))
end

### calculate M from integrated solution

function calculate_M(odesol, pars, timesteps_per_year = nothing)
    """
    takes odesol and pars, returns M at the final time step
    """
    @unpack T2, S2, T4, S4, T0, Tamp, η, μ1, μ2, μ3, μ4, r  = pars

    if timesteps_per_year == nothing
        T1, S1, T3, S3 = odesol[1, :], odesol[2, :], odesol[3, :], odesol[4, :]

        Ud = 1 .- η.*((S4 .- 1) .- (S3 .- T3)) 
        Us = 1 .- η.*((S4 .- 1) .- (S3 .- T3)) .- η.*r.*((S2 .- T2) .- (S1 .- T1))
        Mreturn  = r.*Us .+ Ud
    else 
        T1, S1, T3, S3 = odesol[1, end - timesteps_per_year:end], mean(odesol[2, end - timesteps_per_year:end]), mean(odesol[3, end - timesteps_per_year:end]), mean(odesol[4, end - timesteps_per_year:end])

        Ud = 1 .- η.*((S4 .- 1) .- (S3 .- T3)) 
        Us = 1 .- η.*((S4 .- 1) .- (S3 .- T3)) .- η.*r.*((S2 .- T2) .- (S1 .- T1))
        M  = r.*Us .+ Ud
        Mreturn = mean(M) 
    end

    return Mreturn
end

#= 
--------------------------------------------------------------------------
Noise 
--------------------------------------------------------------------------
=# 

# function stochastic_nondimensional_born!(dx, x, p, t = nothing)
#     @unpack T2, S2, T4, S4, T0, Tamp, η, μ1, μ2, μ3, μ4, μ5, r = p

#     T1, S1, T3, S3 = x

#     Ud = 1 - η*((S4 - 1) - (S3 - T3))
#     Us = 1 - η*((S4 - 1) - (S3 - T3)) - η*r*((S2 - T2) - (S1 - T1))

#     dx[1] = μ1*Us*(T2 - T1) + μ2*(T3 - T1)*heaviside((S1 - T1) - (S3 - T3)) + μ3*(T0 - (Tamp * cos(2 * π * t)) - T1)
#     dx[2] = μ1*Us*(S2 - S1) + μ2*(S3 - S1)*heaviside((S1 - T1) - (S3 - T3)) - μ4 + μ5 * randn()
#     dx[3] = μ1*Ud*(T4 - T3) - μ2*r*(T3 - T1)*heaviside((S1 - T1) - (S3 - T3))
#     dx[4] = μ1*Ud*(S4 - S3) - μ2*r*(S3 - S1)*heaviside((S1 - T1) - (S3 - T3))

#     dx    
# end 

# function stochastic_S_nondimensional_born!(dx, x, p, t = nothing)
#     @unpack T2, S2, T4, S4, T0, Tamp, η, μ1, μ2, μ3, μ4, μ5_S, r = p

#     T1, S1, T3, S3 = x

#     Ud = 1 - η*((S4 - 1) - (S3 - T3))
#     Us = 1 - η*((S4 - 1) - (S3 - T3)) - η*r*((S2 - T2) - (S1 - T1))

#     dx[1] = μ1*Us*(T2 - T1) + μ2*(T3 - T1)*heaviside((S1 - T1) - (S3 - T3)) + μ3*(T0 - (Tamp * cos(2 * π * t)) - T1)
#     dx[2] = μ1*Us*(S2 + μ5_S*randn() - S1) + μ2*(S3 - S1)*heaviside((S1 - T1) - (S3 - T3)) - μ4 
#     dx[3] = μ1*Ud*(T4 - T3) - μ2*r*(T3 - T1)*heaviside((S1 - T1) - (S3 - T3))
#     dx[4] = μ1*Ud*(S4 - S3) - μ2*r*(S3 - S1)*heaviside((S1 - T1) - (S3 - T3))

#     dx    
# end 
 
function stochastic_nondimensional_born!(dx, x, p, t = nothing)
    @unpack T2, S2, T4, S4, T0, Tamp, η, μ1, μ2, μ3, μ4, r, τT0, τS2, τF = p

    T1, S1, T3, S3 = x

    Ud = 1 - η*((S4 - 1) - (S3 - T3))
    Us = 1 - η*((S4 - 1) - (S3 - T3)) - η*r*((S2 + x[6] - T2) - (S1 - T1))

    # fully deterministic part
    dx[1] = μ1*Us*(T2 - T1) + μ2*(T3 - T1)*heaviside((S1 - T1) - (S3 - T3)) + μ3*(T0 + x[5] - (Tamp * cos(2 * π * t)) - T1)
    dx[2] = μ1*Us*(S2 + x[6] - S1) + μ2*(S3 - S1)*heaviside((S1 - T1) - (S3 - T3)) - μ4 + x[7] 
    dx[3] = μ1*Ud*(T4 - T3) - μ2*r*(T3 - T1)*heaviside((S1 - T1) - (S3 - T3))
    dx[4] = μ1*Ud*(S4 - S3) - μ2*r*(S3 - S1)*heaviside((S1 - T1) - (S3 - T3))

    # deterministic part of Ornstein-Uhlenbeck process; if τXX = 0 this is set to zero
    if τT0 == 0
        dx[5] = 0
    else
        dx[5] = -x[5]/τT0
    end

    if τS2 == 0
        dx[6] = 0
    else
        dx[6] = -x[6]/τS2
    end

    if τF == 0
        dx[7] = 0
    else
        dx[7] = -x[7]/τF
    end

    dx
end 

function noisefunction(noisevec)
    function σ(dx, x, p, t = nothing)  
        dx[1] = 0
        dx[2] = 0
        dx[3] = 0
        dx[4] = 0
        dx[5] = noisevec[1]
        dx[6] = noisevec[2]
        dx[7] = noisevec[3]
    end

    return σ
end

end


 