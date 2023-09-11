using Plots, Statistics, Printf, NetCDF, GibbsSeaWater
using FileIO, JLD2
using Dates: AbstractTime

# Mukund's plotting example

################################################################################
# User choices
run_name = "wind5_size256_conc58_rand_thick"

################################################################################
# Functions

function FloeShape(h,k,r)
    θ = LinRange(0,2*π,200)
    h .+ r*sin.(θ), k .+ r*cos.(θ)
end

function FloeOrient(x,y,r,θ)
    l = LinRange(0,r,100)
    l*cos(θ) .+ x, l*sin(θ) .+ y
end

function get_curl(fldx,fldy,dx,dy)
    # fldx must be on the u-grid point and fldy on the v-grid
    # returns field on the omega grid
    nx,ny = size(fldx)
    dvdx = zeros((ny,nx))
    dudy = zeros((ny,nx))
    dvdx[1:end-1,:] = diff(fldy,dims=1)/dx
    dudy[:,1:end-1] = diff(fldx,dims=2)/dy
    return (dvdx - dudy)
end

function get_shear_strain(fldx,fldy,dx,dy)
    # fldx must be on the u-grid point and fldy on the v-grid
    # returns field on the omega grid
    nx,ny = size(fldx)
    dvdx = zeros((ny,nx))
    dudy = zeros((ny,nx))
    dvdx[1:end-1,:] = diff(fldy,dims=1)/dx
    dudy[:,1:end-1] = diff(fldx,dims=2)/dy
    return (dvdx + dudy)
end

function get_normal_strain(fldx,fldy,dx,dy)
    # fldx must be on the u-grid point and fldy on the v-grid
    # returns field on the omega grid
    nx,ny = size(fldx)
    dudx = zeros((ny,nx))
    dvdy = zeros((ny,nx))
    dudx[1:end-1,:] = diff(fldx,dims=1)/dx
    dvdy[:,1:end-1] = diff(fldy,dims=2)/dy
    return (dudx - dvdy)
end

function get_div(fldx,fldy,dx,dy)
    # fldx must be on the u-grid point and fldy on the v-grid
    # returns field on the omega grid
    nx,ny = size(fldx)
    dudx = zeros((ny,nx))
    dvdy = zeros((ny,nx))
    dudx[1:end-1,:] = diff(fldx,dims=1)/dx
    dvdy[:,1:end-1] = diff(fldy,dims=2)/dy
    return (dudx + dvdy)
end

maybe_int(t) = isinteger(t) ? Int(t) : t
minute = 60
hour = 3600
day = 24*3600
year = 360*day

function prettytime(t)
    # Modified from: https://github.com/JuliaCI/BenchmarkTools.jl/blob/master/src/trials.jl
    iszero(t) && return "0 seconds"

    t = maybe_int(t)

    if t < 1e-9
        # No point going to picoseconds, just return something readable.
        return @sprintf("%.3e seconds", t)
    elseif t < 1e-6
        value, units = t * 1e9, "ns"
    elseif t < 1e-3
        value, units = t * 1e6, "μs"
    elseif t < 1
        value, units = t * 1e3, "ms"
    elseif t < minute
        value = t
        units = value == 1 ? "second" : "seconds"
    elseif t < hour
        value = maybe_int(t / minute)
        units = value == 1 ? "minute" : "minutes"
    elseif t < day
        value = maybe_int(t / hour)
        units = value == 1 ? "hour" : "hours"
    elseif t < year
        value = maybe_int(t / day)
        units = value == 1 ? "day" : "days"
    else
        value = maybe_int(t / year)
        units = value == 1 ? "year" : "years"
    end

    if isinteger(value)
        return @sprintf("%d %s", value, units)
    else
        return @sprintf("%.3f %s", value, units)
    end
end

function plot_floes(xfloes,yfloes,rfloes,hfloes,theta_floes,zeta_floes,Lx,Ly,i)
    nfloes = length(rfloes)
    for ifloe in 1:nfloes
        if hfloes[i,ifloe] != 0.0
            fillalpha = 0.05
            #fillalpha = 0.8*hfloes[i,ifloe]/max_floe_h
            plot!(FloeShape(xfloes[i,ifloe]/1000,yfloes[i,ifloe]/1000,rfloes[ifloe]/1000),seriestype = [:shape,], lw = 0.5, c=:blue,linecolor=:black,fillalpha=fillalpha,xlims=(0,Lx/1000),ylims=(0,Ly/1000),aspect_ratio=1,xlabel="x [km]",ylabel="y [km]")
            clr = :red
            if zeta_floes[i,ifloe]<0
                clr=:green
            end
            plot!(FloeOrient(xfloes[i,ifloe]/1000,yfloes[i,ifloe]/1000,rfloes[ifloe]/1000,theta_floes[i,ifloe]),seriestype=[:shape,],lw = 6,linecolor=clr,fillalpha=fillalpha,xlims=(0,Lx/1000),ylims=(0,Ly/1000),aspect_ratio=1,xlabel="x [km]",ylabel="y [km]")
        end
    end
end

################################################################################

# Constants
Q = 300 # Solar constant for Arctic summer (W/m2)
A = 90 # Upward flux constant A (W/m2)
B = 15 # Upward flux constant B (W/m2/K)
alpha_I = 0.7 # Ice albedo
Cd_IO = 0.0055
Cd_IA = 0.00125
Cd_AO = 0.00125
rho_O = 1027.1719
rho_A = 1.25
omega = 2*π/(3600*24)
f = 2*omega*sin(70*π/180)
alpha = 4.9467e-05
beta = 7.8137e-04
Tref = 0
Sref = 34
Tf = -1.8
rho_I = 1000 # Ice floe density [kg/m3]
Lf = 334*1000  # Latent heat of freezing (J/kg)

# Loading surface data
fname = "../data/" * run_name * "/surface.nc"
xc = ncread(fname, "xC")
yc = ncread(fname, "yC")
Nx = length(xc)
Ny = length(yc)
usurf = ncread(fname, "u")[1:Nx,1:Ny,:,:]
vsurf = ncread(fname, "v")[1:Nx,1:Ny,:,:]
wsurf = ncread(fname, "w")
Tsurf = ncread(fname, "T")
Ssurf = ncread(fname, "S")
t = ncread(fname, "time")
sigma = gsw_sigma0.(Ssurf, Tsurf)
nit = length(t)
dx = xc[2] - xc[1]
dy = yc[2] - yc[1]
Lx = xc[end] + dx/2
Ly = yc[end] + dy/2

# Loading floe data
fname = "../data/" * run_name * "/floes.nc"
xfloes = ncread(fname, "x")[1:nit,:]
yfloes = ncread(fname, "y")[1:nit,:]
hfloes = ncread(fname, "h")[1:nit,:]
rfloes = ncread(fname, "r")
ufloes = ncread(fname, "u")[1:nit,:]
vfloes = ncread(fname, "v")[1:nit,:]
theta_floes = ncread(fname, "theta")[1:nit,:]
zeta_floes = ncread(fname, "zeta")[1:nit,:]
nfloes = length(rfloes)
max_floe_h = maximum(hfloes)

# # # Loading BCs data
# fname = "../data/" * run_name * "/surface_bcs.nc"
# hflx = ncread(fname, "hflx")
# taux = ncread(fname, "taux")
# taux = ncread(fname, "taux")
# # tauy = ncread(fname, "tauy")
# # # SIfract = ncread(fname, "SIfract")

i = nit
anim = @animate for i=2:nit
    vort = get_curl(usurf[:,:,1,i],vsurf[:,:,1,i],dx,dy) # curl
    # strain = get_normal_strain(usurf[:,:,1,i],vsurf[:,:,1,i],dx,dy) + get_shear_strain(usurf[:,:,1,i],vsurf[:,:,1,i],dx,dy) # Total strain rate
    # div = get_div(usurf[:,:,1,i],vsurf[:,:,1,i],dx,dy) # Divergence
    #OW = (strain.^2 - vort.^2)/f^2 # Okubo-Weiss parameter
    Ro = vort/f

    #wek = get_curl(taux[:,:,i],tauy[:,:,i],dx,dy)/f/rho_O # Ekman pumping
    conv = -get_div(usurf[:,:,1,i],vsurf[:,:,1,i],dx,dy)
    #
    title_plt = plot(title = prettytime(t[i]), grid = false, showaxis = false, bottom_margin = -30Plots.px,ticks=nothing)
    #
    Ro_plt = contour(xc/1000, yc/1000, Ro', xlabel="x [km]",ylabel="y [km]",title="(c) Ro",linewidth = 0, fill=true,color=:balance,clims=(-0.6, 0.6),legend=false,colorbar=true)
    plot_floes(xfloes,yfloes,rfloes,hfloes,theta_floes,zeta_floes,Lx,Ly,i)
    #
    T_plt = contour(xc/1000, yc/1000, Tsurf[:,:,1,i]', xlabel="x [km]", ylabel="y [km]",title="(b) T surf [°C]",linewidth = 0, fill=true,color=:thermal,legend=false,colorbar=true)#,clims=(Tf, 1))
    plot_floes(xfloes,yfloes,rfloes,hfloes,theta_floes,zeta_floes,Lx,Ly,i)
    #
    S_plt = contour(xc/1000, yc/1000, Ssurf[:,:,1,i]', xlabel="x [km]", ylabel="y [km]",title="(a) S surf [psu]",linewidth = 0, fill=true,color=:haline,legend=false,colorbar=true,clims=(27.5, 29.8))
    plot_floes(xfloes,yfloes,rfloes,hfloes,theta_floes,zeta_floes,Lx,Ly,i)
    #
    wsurf_plt = contour(xc/1000, yc/1000, 3600*24*wsurf[:,:,1,i]', xlabel="x [km]",ylabel="y [km]",title="(d) w surf [m/day]",linewidth = 0, fill=true,color=:balance,clims=(-5, 5),legend=false,colorbar=true)
    plot_floes(xfloes,yfloes,rfloes,hfloes,theta_floes,zeta_floes,Lx,Ly,i)
    #
    plot(title_plt,S_plt,T_plt,Ro_plt,wsurf_plt,size=(1500, 1200),layout = @layout([A{0.01h};[B C]; [D E] ]))
end
mp4(anim, "../videos/" * run_name * ".mp4", fps = 5) # hide
