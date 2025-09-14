using Random
using Distributions
using Plots
using Statistics

k = 0.5 
g = 1.0        
x0 = 0.0       
T = 10.0       
n = 1000       
dt = T / n
sims = 300  

t = range(0, stop=T, length=n)

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
# Simular el proceso de Ornstein-Uhlenbeck
est = zeros(Float64, sims, n)
plot()

for i ∈ 1:sims
    est[i, :] = ornub(x0, n, k, dt, g)
    plot!(t, est[i,:], xlabel="Tiempo", ylabel="x(t)", title="Proceso de Ornstein-Uhlenbeck", label = false)
end
savefig("Procesos de Ornstein Uhlenbeck.png")

means = zeros(Float64, n)
vars = zeros(Float64, n)
for i ∈ 1:n
    means[i] = mean(est[:,i])
    vars[i] = varm(est[:,i], means[i])
end

plot(t, means, xlabel="Tiempo", ylabel="x(t)", title="Promedio de Proceso de Ornstein-Uhlenbeck", label = "Promedio")
savefig("Promedios Procesos de Ornstein Uhlenbeck.png")

plot(t, vars, xlabel="Tiempo", ylabel="var(x)", title="Varianza de Proceso de Ornstein-Uhlenbeck", label = "Varianza")
savefig("Varianza Procesos de Ornstein Uhlenbeck.png")

println(g^2/(2k))
