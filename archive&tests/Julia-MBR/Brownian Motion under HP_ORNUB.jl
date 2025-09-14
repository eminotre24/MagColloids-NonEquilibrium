using Random, Distributions, Statistics, Plots

#Constantes
k = 0.5 
g = 1.0        
x0 = 0.0
v0 = 0.0
γ = 0.75      
T = 10.0       
n = 1000       
dt = T / n
sims = 5000  

t = range(0, stop=T, length=n)

function approx(x,val,ϵ)
    if x > val - ϵ && x < val + ϵ
        return true
    else
        return false
    end
end

function ornub(x0, n, k, dt, g)
    x = zeros(n) 
    x[1] = x0    
    for i in 2:n
        dW = sqrt(dt) * randn()        
        dx = -k * x[i-1] * dt + g * dW 
        x[i] = x[i-1] + dx             
    end
    return x
end

function harm_brown(ic, n::Int, k::Float64, γ::Float64, g::Float64, dt::Float64)
    x = zeros(n)
    v = zeros(n)
    x[1] = ic[1]
    v[1] = ic[2]    
    for i in 2:n
        dW = sqrt(dt) * randn()
        dv = -k * x[i-1] * dt - γ * v[i-1] * dt
        v[i] = v[i-1] + dv 
        dx = v[i] * dt + g * dW
        x[i] = x[i-1] + dx         
    end
    return x
end

#Guardar las runs
est = zeros(Float64, sims, n)
plot()

for i ∈ 1:sims
    est[i, :] = harm_brown([x0 v0], n, k, γ, g, dt)
    plot!(t, est[i,:], xlabel="Tiempo", ylabel="x(t)", title="Brownian Motion under Harmonic Oscilator", label = false)
end
savefig("Prueba_Brownian_M_HO.png")

#Filtracion de Datos
ϵ = 0.25
p1 = 1
p2 = 0
t1 = 1
t2 = 2

#Conversion de tiempo a frame
pos1 = Int(t1/dt)
pos2 = Int(t2/dt)
accepted_trays = Matrix{Float64}(undef, 0, n)

plot()
for i ∈ 1:sims
    trj_i = est[i, :]
    if approx(trj_i[pos1], p1, ϵ) && approx(trj_i[pos2], p2, ϵ)
        global accepted_trays = vcat(accepted_trays, trj_i')
        plot!(t, trj_i, color = "red", alpha=0.1, label="")
    end
end
mean_traj = vec(mean(accepted_trays, dims = 1))
plot!(t, mean_traj, title="Prueba MBR", label = "Promedio", color="black")

#scatter!([pos1 pos2], [p1 p2], yerr = ϵ, label = "Cond. Points", legend=:topright)

savefig("Prueba_MBR_Brownian_M_HO.png")