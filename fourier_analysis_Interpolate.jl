using NetCDF, PyPlot, Statistics, FFTW, StatsBase, Interpolations, ProgressMeter

vx = ncread("SVD/2013-01-01_2023-02-17_vx_vy_5dmean.nc", "vx")
vy = ncread("SVD/2013-01-01_2023-02-17_vx_vy_5dmean.nc", "vy")
t  = Float64.(ncread("SVD/2013-01-01_2023-02-17_vx_vy_5dmean.nc", "time"))
v  = sqrt.(vx.^2 .+ vy.^2)
nx, ny, nt = size(v)

dt = 5
t_new = t[1]:dt:t[end]
n = length(t_new)
freqs = 1 ./(dt*n)*(0:n)
perds = 1 ./ freqs         # period in days

pows    = zeros(nx,ny,n)
pows_rel= zeros(nx,ny,n)
# vhat_clean = zeros(ComplexF64,nx,ny,n)
# v_clean = zeros(ComplexF64,nx,ny,n)
@showprogress for i in axes(v,1)
    for j in axes(v,2)
        if !all(isnan.(v[i,j,:]) .|| v[i,j,:] .== 0)
            v_interp          = LinearInterpolation(t,v[i,j,:]).(t_new)
            vhat              = fft(v_interp)
            pow               = real.(vhat .* conj.(vhat))
            pows[i,j,:]       = pow
            pows_rel[i,j,:]   = pow ./ maximum(pow[2:end-1])  # normalize   (ends of the spectrum go very high)
            # ind_clean         = findall(pow[2:end-1] .> 0.08) .+ 1
            # vhat_clean[i,j,ind_clean] = vhat[ind_clean]
            # v_clean[i,j,:]    = ifft(vhat_clean[i,j,:])  # not good...
        end
    end
end

figure(figsize=(20,20))
pcolormesh(pows[1:250,150:450,19]',clim=(0,5e9)); colorbar(); axis("equal"); title("absolute power of 2-peak frequency")
savefig("double-peak-power.jpg")
figure(figsize=(20,20)) # this one does not seem to show the intersting stuff
pcolormesh(pows_rel[:,:,19]',clim=(0,0.3)); colorbar(); axis("equal"); title("power of 2-peak freq. relative to 1-peak freq.")
savefig("double-peak-rel-power.jpg")

figure(figsize=(12,10))
plot(perds[4:50], pows[30, 316,4:50],label="double peak, terminus")
plot(perds[4:50], pows[54, 294,4:50],label="'hole' at the terminus")
plot(perds[4:50], pows[93, 226,4:50],label="single peak")
plot(perds[4:50], pows[111,287,4:50],label="valery")
xlabel("Period (d)")
ylabel("power")
legend()
savefig("power-spectrum.jpg")

figure(figsize=(20,20))
plot(v[30, 316,:],label="double peak, terminus")
plot(v[54, 294,:],label="'hole' at the terminus")
plot(v[93, 226,:],label="single peak")
plot(v[111,287,:],label="valery")
legend()
savefig("vel-time-series.jpg")
