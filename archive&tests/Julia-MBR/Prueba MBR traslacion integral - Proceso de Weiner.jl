using Random, Distributions, Statistics, Plots

#Constantes
k = 1.0 
g = 1.0        
x0 = 0.0
γ = 1.5     
T = 20.0       
n = 1000       
dt = T / n
sims = 6000  

t = range(0, stop=T, length=n)

function approx(x,val,ϵ)
    if x > val - ϵ && x < val + ϵ
        return true
    else
        return false
    end
end

function weiner(x0, n, dt)
    x = zeros(n)
    x[1] = x0
    for i in 2:n
        dW = sqrt(dt) * randn()
        x[i] = x[i-1] + dW  
    end
    return x
end

function ornub(x0, n, k, γ, dt, g)
    x = zeros(n) 
    x[1] = x0    
    for i in 2:n
        dW = sqrt(dt) * randn()        
        dx = (-k / γ) * x[i-1] * dt + g * dW 
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

#Filtracion de Datos
t1 = 1
t2 = 2

#Conversion de tiempo a frame
pos1 = Int(t1/dt)
pos2 = Int(t2/dt)
accepted_trays = Matrix{Float64}(undef, 0, n)

for i in 1:sims
    trj_i = weiner(x0, n, dt)
    x_0 = trj_i[pos2]
    d = trj_i[pos2] - trj_i[pos1] #trajectoria x_τ = x_0 - d -> d = x_0 - x_τ
    if d > 0.25
        #Reajuste de la trayectoria
        mbrtrj = -(trj_i .- x_0) ./ d
        global accepted_trays = vcat(accepted_trays, mbrtrj')
        plot!(t, mbrtrj, color = "red", alpha=0.1, label="")
    end
end

f(x) = 0.5

mean_traj = vec(mean(accepted_trays, dims = 1))
hline!(f.(t), color="blue", label = "0.5")
plot!(t, mean_traj, title="Prueba MBR", label = "Promedio", color="black", ylims=(-4, 4))
#plot!(t, ones(n).*0.5, label="y = 0.5", linestyle=:dash, color=:red)

#scatter!([pos1 pos2], [p1 p2], yerr = ϵ, label = "Cond. Points", legend=:topright)

savefig("Prueba_MBR_Weiner.png")

