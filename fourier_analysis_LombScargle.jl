using NetCDF, PyPlot, Statistics, FFTW, StatsBase, LombScargle, ProgressMeter

vx = ncread("SVD/2013-01-01_2023-02-17_vx_vy_5dmean.nc", "vx")
vy = ncread("SVD/2013-01-01_2023-02-17_vx_vy_5dmean.nc", "vy")
t  = ncread("SVD/2013-01-01_2023-02-17_vx_vy_5dmean.nc", "time")
v  = sqrt.(vx.^2 .+ vy.^2)
nx, ny, nt = size(v)

# calculate the Fourier transform
v_flat = reshape(v, nx*ny, nt)
Vhats = []
Ids = []
m1s = []; m2s = []; m3s = []
@showprogress for i in axes(v_flat, 1)
    if !all(isnan.(v_flat[i,:]) .|| v_flat[i,:] .== 0)
        Ind = findall(.!isnan.(v_flat[i,:]))
        Vhat = lombscargle(t[Ind], v_flat[i,Ind])
        t1 = findall(300 .< period(Vhat) .< 450)
        t2 = findall(160 .< period(Vhat) .< 200)
        t3 = findall(450 .< period(Vhat) .< 600)
        m1 = maximum(power(Vhat)[t1])
        m2 = maximum(power(Vhat)[t2])
        m3 = maximum(power(Vhat)[t3])
        push!(Vhats,Vhat)
        push!(m1s,m1)
        push!(m2s,m2)
        push!(m3s,m3)
        push!(Ids,i)
    end
end

# plot peak amplitude of certain period
p1 = zeros(nx*ny); p1[Ids] = m1s
p2 = zeros(nx*ny); p2[Ids] = m2s
p3 = zeros(nx*ny); p3[Ids] = m3s

p1_ = reshape(p1,nx,ny)
p2_ = reshape(p2,nx,ny)
p3_ = reshape(p3,nx,ny)

# figure(); pcolormesh((p2_./p1_)',clim=(0,0.5)); colorbar()

# plot freq and time domain side by side
a = zeros(Int,size(v_flat,1))
a[Ids] = range(1,length(Ids))
b = reshape(a,nx,ny)
is = [(70, 242), # channel-like feature
      (93, 226), # upstream, single peak
      (30, 316), # at the snout, very bright yellow
    #   (58, 334), # snout more to the east
    #   (58, 258),
      (152,330)] # Valery

for (i,j) in is
    figure()
    subplot(1,2,1)
    vv = Vhats[b[i,j]]
    plot(period(vv),power(vv)); title("$i,$j"); xlim(0,900)
    subplot(1,2,2)
    plot(t,v[i,j,:])
end
