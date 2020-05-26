#Mortality (and the split into high-risk and low risk taken from Spain Data)
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
name = "Spain_"
#Valores originales NumC = 2.28, Dur=5.8
function simulacion( name; num_dias_simulacion = 350, mu_jov = 3.32664/1000, mu_old= 9.35358/100.0, NumC = 2.28, Dur=5.8, fun_1 = alpha_1, fun_2=alpha_2)
    #Common parameters
    c = NumC/Dur
    #Death and recovery rates:
    mu_1=mu_jov/Dur
    rec_1=(1-mu_jov)/Dur
    mu_2=mu_old/Dur
    rec_2=(1-mu_old)/Dur

    function SIRD(du,u,p,t; fun_1=fun_1, fun_2=fun_2)
        al_1 = fun_1(t)#These functions are defined outside simulacion
        al_2 = fun_2(t)#These functions are defined outside simulacion
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

function print_res(sol,nombre, tofile=false)
    numDeads_1 = string(round(last(sol)[5],digits=2))
    numDeads_2 = string(round(last(sol)[6],digits=2))
    numDeads_Tot = string(round(last(sol)[8],digits=2))

    if (tofile)
        write(file, nombre)
        write(file,"\n")
        write(file,string([numDeads_1,numDeads_2, numDeads_Tot]))
        write(file,"\n")
    else
        print("\n")
        print(name*nombre)
        print("\n")
        print([numDeads_1,numDeads_2, numDeads_Tot])
        print("\n")
    end
end



for NumC=[2.2,2.5,2.8]
    for dur=[4,6,8]
        for Avail=[0.4^(1/2)]
            #For each simulation we define the quarantine functions and then run simulate
            #Cumulative death counts:
            font = Plots.font("Helvetica", 9)
            deathsPlot1 = plot(title= "Cumulative number of deaths per group:",xaxis=("t (in days)"), yaxis = ("number of deaths (in thousands)"), legend=:topleft, legendfont= font, size =(800,600), lw=17)

            font = Plots.font("Helvetica", 11)
            deathsPlot2 = plot(title= "Cumulative number of deaths:",xaxis=("t (in days)"), yaxis = ("number of deaths (in thousands)"), legend=:topleft, legendfont= font, size =(800,600), lw=17)

            font = Plots.font("Helvetica", 11)
            InfectedPlot = plot(title= "Number of Infected individuals in each group", xaxis=("t (in days)"), yaxis = ("number of deaths (in thousands)"), legend=:topright, legendfont= font, size =(800,600), lw=17)

            avail = Avail
            #Quarentena maxima
            max_quar = 100
            #NumC = 2.28
            #avail=0.4

            # Escenario 1: NQ Nada de cuarentena
            function alpha_1(t)
                return 1.0
            end
            function alpha_2(t)
                return 1.0
            end
            paramsName = "_R0_"*string(NumC)*"_tau_"*string(dur)*"_avail_"*string(avail)*":"


            nombre = "NQ"
            sol = simulacion(nombre, NumC=NumC, Dur=dur, fun_1 = alpha_1, fun_2=alpha_2)
            plot!(deathsPlot1,sol, linecolor= :blue, vars = [(0,5)], linestyle = :dot, label="NQ:D_l",yformatter = :plain)
            plot!(deathsPlot1,sol, linecolor= :red, vars = [(0,6)],linestyle = :dot, label="NQ:D_h",yformatter = :plain)

            plot!(InfectedPlot,sol, linecolor= :blue, vars = [(0,1)], linestyle = :dot, label="NQ:I_l",yformatter = :plain)
            plot!(InfectedPlot,sol, linecolor= :red, vars = [(0,2)],linestyle = :dot, label="NQ:I_h",yformatter = :plain)


            plot!(deathsPlot2,sol, vars = [(0,8)], linestyle = :dot, label="NQ",yformatter = :plain)
            print_res(sol, nombre*paramsName)


            #LINESTYLES :auto, :dash, :dashdot, :dot, :solid
            # Escenario 2a: RBQ short sep.
            function beta_1(t)
                t<90 ? avail : 1.0
            end
            function beta_2(t)
                t<max_quar ? avail : 1.0
            end
            nombre = "RBQ_short"
            sol = simulacion(nombre,  NumC=NumC, Dur=dur, fun_1 = beta_1, fun_2=beta_2)

            plot!(deathsPlot1, sol, linecolor=:blue, vars = (0,5), linestyle = :dash, label="RBQ_S:D_l")
            plot!(deathsPlot1, sol, linecolor=:red , vars = (0,6), linestyle = :dash, label="RBQ_S:D_h")

            plot!(InfectedPlot, sol, linecolor=:blue, vars = (0,1), linestyle = :dash, label="RBQ_S:I_l")
            plot!(InfectedPlot, sol, linecolor=:red , vars = (0,2), linestyle = :dash, label="RBQ_S:I_h")


            plot!(deathsPlot2, sol, vars = (0,8), linestyle = :dash, label="RBQ_S",yformatter = :plain)
            print_res(sol, nombre*paramsName)

            # Escenario 2b: RBQ medium sep.
            function delta_1(t)
                t<40 ? avail : 1.0
            end
            function delta_2(t)
                t<max_quar ? avail : 1.0
            end
            nombre = "RBQ_medium"

            sol = simulacion(nombre,NumC=NumC, Dur=dur, fun_1 = delta_1, fun_2=delta_2)
            plot!(deathsPlot1, sol, linecolor=:blue, vars = (0,5), linestyle = :dashdot, label="RBQ_M:D_l")
            plot!(deathsPlot1, sol, linecolor=:red , vars = (0,6), linestyle = :dashdot, label="RBQ_M:D_h")

            plot!(InfectedPlot, sol, linecolor=:blue, vars = (0,1), linestyle = :dashdot, label="RBQ_M:I_l")
            plot!(InfectedPlot, sol, linecolor=:red , vars = (0,2), linestyle = :dashdot, label="RBQ_M:I_h")


            plot!(deathsPlot2, sol, vars = (0,8), linestyle = :dash, label="RBQ_M",yformatter = :plain)
            print_res(sol, nombre*paramsName)



            # Escenario 2c: RBQ long sep
            function gamma_1(t)
                t<10.0 ? avail : 1.0
            end
            function gamma_2(t)
                t<max_quar ? avail : 1.0
            end
            nombre = "RBQ_long"
            sol = simulacion(nombre,  NumC=NumC, Dur=dur, fun_1 = gamma_1, fun_2=gamma_2)
            plot!(deathsPlot1, sol, linecolor=:blue, vars = (0,5), linestyle = :solid, label="RBQ_L:D_l")
            plot!(deathsPlot1, sol, linecolor=:red , vars = (0,6), linestyle = :solid, label="RBQ_L:D_h")

            plot!(InfectedPlot, sol, linecolor=:blue, vars = (0,1), linestyle = :solid, label="RBQ_L:I_l")
            plot!(InfectedPlot, sol, linecolor=:red , vars = (0,2), linestyle = :solid, label="RBQ_L:I_h")


            plot!(deathsPlot2, sol, vars = (0,8), linestyle = :solid, label="RBQ_L",yformatter = :plain)
            print_res(sol, nombre*paramsName)


            # Escenario 3: EQ quarantine
            function eta_1(t)
                t<max_quar ? avail : 1.0
            end
            function eta_2(t)
                t<max_quar ? avail : 1.0
            end

            nombre = "EQ"
            sol = simulacion(nombre, NumC=NumC, Dur=dur, fun_1 = eta_1, fun_2=eta_2)

            plot!(deathsPlot1, sol, linecolor=:blue, vars = (0,5), linestyle = :dashdotdot, label="EQ:D_l")
            plot!(deathsPlot1, sol, linecolor=:red , vars = (0,6), linestyle = :dashdotdot, label="EQ:D_h")

            plot!(InfectedPlot, sol, linecolor=:blue, vars = (0,1), linestyle = :dashdotdot, label="EQ:I_l")
            plot!(InfectedPlot, sol, linecolor=:blue, vars = (0,2), linestyle = :dashdotdot, label="EQ:I_h")


            plot!(deathsPlot2, sol, vars = (0,8), linestyle = :dashdotdot, label="EQ",yformatter = :plain)
            print_res(sol, nombre*paramsName)

            savefig(deathsPlot1, name*paramsName*"_deathsPlot_v1.png")
            savefig(deathsPlot2, name*paramsName*"_deathsPlot_v2.png")
            savefig(InfectedPlot,name*paramsName*"_InfectedPlot.png")
        end
    end
end
