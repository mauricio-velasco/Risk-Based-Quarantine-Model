#Mortality (and the split into high-risk and low risk taken from ICL group)
using DifferentialEquations
using Plots

# 1. I_1
# 2. I_2
# 3. S_1
# 4. S_2
# 5. D_1
# 6. D_2
# 7. R
# 8. D

#Each scenario has a name and mortalities mu_jov and mu_old
name = "Imperial_"

function simulacion( name; num_dias_simulacion = 350, mu_jov = 1.2921/1000, mu_old= 3.24/100.0, NumC = 2.28, Dur=5.8 )
    #Common parameters
    c = NumC/Dur
    #Death and recovery rates:
    mu_1=mu_jov/Dur
    rec_1=(1-mu_jov)/Dur
    mu_2=mu_old/Dur
    rec_2=(1-mu_old)/Dur

    function SIRD(du,u,p,t)
        al_1 = alpha_1(t)#These functions are defined outside simulacion
        al_2 = alpha_2(t)#These functions are defined outside simulacion
        du[1] =  c*al_1*(1/N)*(al_1*u[1] + al_2*u[2])*u[3]-(rec_1+mu_1)*u[1]
        du[2] =  c*al_2*(1/N)*(al_1*u[1] + al_2*u[2])*u[4]-(rec_2+mu_2)*u[2]
        du[3] = (-1)*c*al_1*(1/N)*(al_1*u[1] + al_2*u[2])*u[3]
        du[4] = (-1)*c*al_2*(1/N)*(al_1*u[1] + al_2*u[2])*u[4]
        du[5] = mu_1*u[1]
        du[6] = mu_2*u[2]
        du[7] = rec_1*u[1] + rec_2*u[2]
        du[8] = mu_1*u[1]+mu_2*u[2]
    end
    #Initial conditions (Bogota):
    N = 8281.0
    u0 = [0.1; 0.01; 7183.035-0.1; 1097.995-0.01; 0.0 ; 0.0 ; 0.0; 0.0]
    #

    tspan = (0.0,num_dias_simulacion)
    prob = ODEProblem(SIRD,u0,tspan)
    sol = solve(prob)
    return(sol)
end

function print_res(sol,nombre,tofile=false)
    numDeads_1 = string(round(last(sol)[5],digits=2))
    numDeads_2 = string(round(last(sol)[6],digits=2))
    numDeads_Tot = string(round(last(sol)[8],digits=2))

    if (tofile)
        write(file, nombre)
        write(file,"\n")
        write(file,string([numDeads_1,numDeads_2, numDeads_Tot]))
        write(file,"\n")
    else
        print(nombre)
        print("\n")
        print([numDeads_1,numDeads_2, numDeads_Tot])
        print("\n")
    end
end

#For each simulation we define the quarantine functions and then run simulate
#Cumulative death counts:
font = Plots.font("Helvetica", 7)
deathsPlot1 = plot(title= "Cumulative number of deaths per group:", legend=:topleft, legendfont= font, size =(800,600), lw=15)

font = Plots.font("Helvetica", 10)
deathsPlot2 = plot(title= "Cumulative number of deaths:", legend=:topleft, legendfont= font, size =(800,600), lw=15)

font = Plots.font("Helvetica", 10)
InfectedPlot = plot(title= "Number of Infected individuals each group", legend=:topright, legendfont= font, size =(800,600), lw=15)

#Quarentena maxima
max_quar = 100
# Escenario 1: NQ Nada de cuarentena
function alpha_1(t)
    return 1.0
end
function alpha_2(t)
    return 1.0
end
nombre = "NQ"
sol = simulacion(nombre)
plot!(deathsPlot1,sol, linecolor= :blue, vars = [(0,5)], linestyle = :dot, label="NQ:D_1",yformatter = :plain)
plot!(deathsPlot1,sol, linecolor= :red, vars = [(0,6)],linestyle = :dot, label="NQ:D_2",yformatter = :plain)

plot!(InfectedPlot,sol, linecolor= :blue, vars = [(0,1)], linestyle = :dot, label="NQ:I_1",yformatter = :plain)
plot!(InfectedPlot,sol, linecolor= :red, vars = [(0,2)],linestyle = :dot, label="NQ:I_2",yformatter = :plain)


plot!(deathsPlot2,sol, vars = [(0,8)], linestyle = :dot, label="NQ",yformatter = :plain)
print_res(sol, nombre)


#LINESTYLES :auto, :dash, :dashdot, :dot, :solid


# Escenario 2a: RBQ short serp.
function alpha_1(t)
    t<90 ? 0.4 : 1.0
end
function alpha_2(t)
    t<max_quar ? 0.4 : 1.0
end
nombre = "RBQ_short"
sol = simulacion(nombre)

plot!(deathsPlot1, sol, linecolor=:blue, vars = (0,5), linestyle = :dash, label="RBQ_S:D_1")
plot!(deathsPlot1, sol, linecolor=:red , vars = (0,6), linestyle = :dash, label="RBQ_S:D_2")

plot!(InfectedPlot, sol, linecolor=:blue, vars = (0,1), linestyle = :dash, label="RBQ_S:I_1")
plot!(InfectedPlot, sol, linecolor=:red , vars = (0,2), linestyle = :dash, label="RBQ_S:I_2")


plot!(deathsPlot2, sol, vars = (0,8), linestyle = :dash, label="RBQ_S",yformatter = :plain)
print_res(sol, nombre)

# Escenario 2b: RBQ medium sep.
function alpha_1(t)
    t<40 ? 0.4 : 1.0
end
function alpha_2(t)
    t<max_quar ? 0.4 : 1.0
end
nombre = "RBQ_medium"
sol = simulacion(nombre)
plot!(deathsPlot1, sol, linecolor=:blue, vars = (0,5), linestyle = :dashdot, label="RBQ_M:D_1")
plot!(deathsPlot1, sol, linecolor=:red , vars = (0,6), linestyle = :dashdot, label="RBQ_M:D_2")

plot!(InfectedPlot, sol, linecolor=:blue, vars = (0,1), linestyle = :dashdot, label="RBQ_M:I_1")
plot!(InfectedPlot, sol, linecolor=:red , vars = (0,2), linestyle = :dashdot, label="RBQ_M:I_2")

plot!(deathsPlot2, sol, vars = (0,8), linestyle = :dashdot, label="RBQ_M")
print_res(sol, nombre)

# Escenario 2c: RBQ long sep
function alpha_1(t)
    t<10 ? 0.4 : 1.0
end
function alpha_2(t)
    t<max_quar ? 0.4 : 1.0
end
nombre = "RBQ_long"
sol = simulacion(nombre)
plot!(deathsPlot1, sol, linecolor=:blue, vars = (0,5), linestyle = :solid, label="RBQ_L:D_1")
plot!(deathsPlot1, sol, linecolor=:red , vars = (0,6), linestyle = :solid, label="RBQ_L:D_2")

plot!(InfectedPlot, sol, linecolor=:blue, vars = (0,1), linestyle = :solid, label="RBQ_L:I_1")
plot!(InfectedPlot, sol, linecolor=:red , vars = (0,2), linestyle = :solid, label="RBQ_L:I_2")


plot!(deathsPlot2, sol, vars = (0,8), linestyle = :solid, label="RBQ_L",yformatter = :plain)
print_res(sol, nombre)


# Escenario 3: EQ quarantine
function alpha_1(t)
    t<max_quar ? 0.4 : 1.0
end
function alpha_2(t)
    t<max_quar ? 0.4 : 1.0
end

nombre = "EQ"
sol = simulacion(nombre)

plot!(deathsPlot1, sol, linecolor=:blue, vars = (0,5), linestyle = :dashdotdot, label="EQ:D_1")
plot!(deathsPlot1, sol, linecolor=:red , vars = (0,6), linestyle = :dashdotdot, label="EQ:D_2")

plot!(InfectedPlot, sol, linecolor=:blue, vars = (0,1), linestyle = :dashdotdot, label="EQ:I_1")
plot!(InfectedPlot, sol, linecolor=:blue, vars = (0,2), linestyle = :dashdotdot, label="EQ:I_2")


plot!(deathsPlot2, sol, vars = (0,8), linestyle = :dashdotdot, label="EQ",yformatter = :plain)
print_res(sol, nombre)

savefig(deathsPlot1, name*"_deathsPlot_v1.png")
savefig(deathsPlot2, name*"_deathsPlot_v2.png")
savefig(InfectedPlot,name*"_InfectedPlot.png")
