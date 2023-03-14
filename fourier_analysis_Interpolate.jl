using NetCDF, PyPlot, Statistics, FFTW, StatsBase, Interpolations, ProgressMeter

vx = ncread("SVD/2013-01-01_2023-02-17_vx_vy_5dmean.nc", "vx")
vy = ncread("SVD/2013-01-01_2023-02-17_vx_vy_5dmean.nc", "vy")
t  = Float64.(ncread("SVD/2013-01-01_2023-02-17_vx_vy_5dmean.nc", "time"))
tvec = Dates.DateTime(2013,6,30)+Dates.Day.(ncread(nf,"time"))
v  = sqrt.(vx.^2 .+ vy.^2)
nx, ny, nt = size(v)

dt = 1
t_new = t[1]:dt:t[end]
n = length(t_new)
freqs = 1 ./(dt*n)*(0:n)
perds = 1 ./ freqs         # period in days

pows    = zeros(nx,ny,n)
phase   = zeros(nx,ny,n)
# vhat_clean = zeros(ComplexF64,nx,ny,n)
# v_clean = zeros(ComplexF64,nx,ny,n)
@showprogress for i in axes(v,1)
    for j in axes(v,2)
        if !all(isnan.(v[i,j,:]) .|| v[i,j,:] .== 0)
            v_interp          = LinearInterpolation(t,v[i,j,:]).(t_new)
            vhat              = fft(v_interp)
            pow               = real.(vhat .* conj.(vhat))
            pows[i,j,:]       = pow
            phase[i,j,:]      = atan.(real.(vhat), imag.(vhat))
            # ind_clean         = findall(pow[2:end-1] .> 0.08) .+ 1
            # vhat_clean[i,j,ind_clean] = vhat[ind_clean]
            # v_clean[i,j,:]    = ifft(vhat_clean[i,j,:])  # not good...
        end
    end
end

# double peak power
figure(figsize=(20,20))
pows_2 = sum(pows[:,:,12:19],dims=3)[:,:,1]
pcolormesh(pows_2[1:250,150:450]',clim=(0, 3e11)); colorbar(); axis("equal"); title("summed power of higher frequencies")
savefig("double-peak-power.jpg")
# single peak power
figure(figsize=(20,20))
pcolormesh(pows[1:250,150:450,10]',clim=(0, 6e11)); colorbar(); axis("equal"); title("absolute power of 1/yr frequency")
savefig("single-peak-power.jpg")
# relative power
figure(figsize=(20,20)) # this one does not seem to show the intersting stuff
pcolormesh(pows_2' ./ pows[:,:,10]',clim=(0,0.65)); colorbar(); axis("equal"); title("power of higher freqs. relative to 1/yr freq.")
savefig("double-peak-rel-power.jpg")
# double peak power
figure(figsize=(20,20))
pcolormesh(phase[1:250,150:450,10]' ./ pi,clim=(0, 0.5)); colorbar(); axis("equal"); title("Phase of main velocity peak (*pi)")
savefig("phase.jpg")

figure(figsize=(12,10))
plot(perds[4:50], pows[30, 316,4:50],label="double peak, terminus")
plot(perds[4:50], pows[54, 294,4:50],label="'hole' at the terminus")
plot(perds[4:50], pows[93, 226,4:50],label="single peak")
plot(perds[4:50], pows[111,287,4:50],label="valerie")
xlabel("Period (d)")
ylabel("power")
legend()
savefig("power-spectrum.jpg")

figure(figsize=(20,20))
plot(v[30, 316,:],label="double peak, terminus")
plot(v[54, 294,:],label="'hole' at the terminus")
plot(v[93, 226,:],label="single peak")
plot(v[150,326,:],label="valerie")
legend()
savefig("vel-time-series.jpg")

figure(figsize=(20,20))
# plot(v[26, 329,:],label="0")
plot(v[30, 316,:],label="1")
plot(v[40, 307,:],label="2")
plot(v[60, 260,:],label="3")
plot(v[49, 331,:],label="side")
legend()
savefig("vel-time-series_phase_shift.jpg")


# compare to terminus positions
vtime   = DateTime(2013,6,30) + Day.(t)

vals, head  = readdlm("terminus_data/terminus_positions.csv", ',',header=true)
tertime = DateTime.(vals[:,1])
ter_pos = Float64.(vals[:,2])
tertime_sorted = tertime[sortperm(tertime)]
ter_pos_sorted = ter_pos[sortperm(tertime)]

# plot together
fig, ax1 = subplots()
color = "tab:red"
ax1.set_xlabel("Year")
ax1.set_ylabel("Relative terminus position (m)"; color)
# ax1.plot(tertime_sorted, ter_pos_sorted; color)
ax1.plot(tertime_sorted[2:end], diff(ter_pos_sorted); color)
ax1.tick_params(axis="y", labelcolor=color)

ax2 = ax1.twinx()
color = "tab:blue"
ax2.set_ylabel("velocity"; color)
ax2.plot(tvec, v[26,327,:]; color)
ax2.tick_params(axis="y", labelcolor=color)
# time_idx = Int64.(indexin(tertime_sorted, vtime))

fig.tight_layout() 
