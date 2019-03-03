#####################################
########### P4 Data Clean ###########
#####################################
####### ADMIN AND PACKAGES ######

# Set path to solvers
# Replace with appropriate path where I saved solvers in the beginning of the course.
PATH_TO_SOLVERS = ENV["ERE291_SOLVERS"]

# Boot up packages
using JuMP
using AmplNLWriter
using Pkg
using Plots
using DataFrames


#Now use knitro
ampl_solver="knitro"


using CSV

data = CSV.read("P4_test_data_wind_price.csv")
#test test test
#Let's try this

display(data)
#Wind Data
WINDENERGY=data[:,4]
WINDPRICE=data[:,3]


#Parameters
NUMBERPASSIVE=4
HEIGHTDROP=200 #[m]
STORAGEACTIVEAREA=pi*40^2 #[m^2]
HEIGHTPASSIVEMAX=40 #[m]
HEIGHTPASSIVEMIN=0.5 #[m]
STORAGEPASSIVEMAX=6057 #[m^3]
STORAGEACTIVEMAX=153943 #[m^3]
HEIGHTACTIVEMAX=30.626 #[m]
HEIGHTACTIVEMIN=0.5 #[m]
OVERALLTURBINEEFFICIENCY=0.8 #This includes the Penstock, Maniforld, Nozzle, Runner (a bit low), drive, and generator (also lowball)
OVERALLPUMPEFFICIENCY=0.91 #This reflects waterways losses, pump losses, motor losses, transformer losses.
STORAGEPASSIVEAREA=151.425 #[m^2]
HEADLOSSFRACTION=0.95 #Effective height pressure lost from friction in penstock.
#LAGTIME=30 #[s]
DISCHARGERATEMAX=2*3600 #[m^3/hr]
CHARGERATEMAX=2*3600
DENSITYWATER=1000 #[kg/m^3]
DISCHARGECOEFFICIENT=0.6
GRAVITY=9.81 #m/s^2

#Solved from maximum Discharge Rate and maximum passive height.
TRANSFERRADIUS=0.0973 #m^3
<<<<<<< HEAD
#T=length(WINDENERGY)
T = 30
=======
T=30 #length(WINDENERGY)
>>>>>>> 892f6338bbec0f56dc7aa5bd10b8c1351ff8b335
TIMESTEP=1 #[hr]

m=Model(solver=AmplNLSolver(joinpath(PATH_TO_SOLVERS,ampl_solver),["outlev=2"]))

#Note that changing expressions use LowerCase while unchanging parameters are all UpperCase

<<<<<<< HEAD
@variable(m,HeightPassive[t=1:T]>=0,start=HEIGHTPASSIVEMIN)
@variable(m,HeightActive[t=1:T]>=0,start=HEIGHTACTIVEMAX)
#@variable(m,GridEnergy[t=1:T]>=0,start=1)
#@variable(m, SellWindEnergy[t=1:T] >= 0,start=1)
#@variable(m, SaveWindEnergy[t=1:T] >= 0,start=1)
=======
@variable(m,HeightPassive[t=1:T]>=0,start=HEIGHTPASSIVEMAX)
@variable(m,HeightActive[t=1:T]>=0,start=HEIGHTACTIVEMAX)
@variable(m,GridEnergy[t=1:T]>=0,start=20)
#@variable(m, SellWindEnergy[t=1:T] >= 0,start=1)
@variable(m, SaveWindEnergy[t=1:T] >= 0,start=0)
>>>>>>> 892f6338bbec0f56dc7aa5bd10b8c1351ff8b335
#@variable(m, StoreWindEnergyOut[t=1:T] >= 0,start=1)
#@variable(m,InStorage[t=1:T+1] >=0,start=1)
@variable(m,PumpEnergy[t=1:T]>=0,start=1)

#Radius of Penstock Pipe (transfer radius)


#Discharge Rate Equation
<<<<<<< HEAD

@NLexpression(m,DischargeRate[t=1:T],NUMBERPASSIVE*pi*(TRANSFERRADIUS^2)*DISCHARGECOEFFICIENT*sqrt(2*GRAVITY*HeightPassive[t])*3600) #[m^3/hr]
=======
@NLexpression(m,DischargeRate[t=1:T],3600*NUMBERPASSIVE*pi*(TRANSFERRADIUS^2)*DISCHARGECOEFFICIENT*sqrt(2*GRAVITY*HeightPassive[t])) #[m^3/hr]
>>>>>>> 892f6338bbec0f56dc7aa5bd10b8c1351ff8b335

#Charge Rate
#Need an expression for PumpEnergy
#@NLexpression(m,PumpEnergy[t=1:T],GridEnergy[t]+SaveWindEnergy[t])

<<<<<<< HEAD
@NLexpression(m,ChargeRate[t=1:T],3600E6*PumpEnergy[t]/(DENSITYWATER*OVERALLPUMPEFFICIENCY*GRAVITY*HEADLOSSFRACTION*HEIGHTDROP*TIMESTEP)) #[m^3/hr]
=======
@NLexpression(m,ChargeRate[t=1:T],3600*(10^6)*PumpEnergy[t]/(DENSITYWATER*OVERALLPUMPEFFICIENCY*GRAVITY*HEADLOSSFRACTION*HEIGHTDROP)) #[m^3/hr]
>>>>>>> 892f6338bbec0f56dc7aa5bd10b8c1351ff8b335

#Electricity Generation
@NLexpression(m,HydroElectricityGen[t=1:T],OVERALLTURBINEEFFICIENCY*DischargeRate[t]*DENSITYWATER*GRAVITY*HEIGHTDROP*TIMESTEP/(3.6*10^9))

#This fulfills degrees of freedom and the wind energy balance.
#@NLexpression(m,SellWindEnergy[t=1:T],WINDENERGY[t]-SaveWindEnergy[t])

#energy balance constraint
#@NLconstraint(m, [t = 1:T], (SellWindEnergy[t] + SaveWindEnergy[t] + HydroElectricityGen[t]) >= (GridEnergy[t] + WINDENERGY[t]) )

#Initial Water Depths
@NLconstraint(m,HeightPassive[1]==HEIGHTPASSIVEMAX)
@NLconstraint(m,HeightActive[1]==HEIGHTACTIVEMIN)

#Water Depth States
@NLconstraint(m,[t=1:T-1],HeightPassive[t+1]==HeightPassive[t]+DischargeRate[t]*TIMESTEP/(NUMBERPASSIVE*STORAGEPASSIVEAREA)-ChargeRate[t]*TIMESTEP/(NUMBERPASSIVE*STORAGEPASSIVEAREA))
@NLconstraint(m,[t=1:T-1],HeightActive[t+1]==HeightActive[t]-DischargeRate[t]*TIMESTEP/(STORAGEACTIVEAREA)-ChargeRate[t]*TIMESTEP/(STORAGEACTIVEAREA))

#Water Depth Limits
@NLconstraint(m,[t=1:T],HeightPassive[t]<=HEIGHTPASSIVEMIN)
@NLconstraint(m,[t=1:T],HeightPassive[t]>=HEIGHTPASSIVEMAX)

@NLconstraint(m,[t=1:T],HeightActive[t]<=HEIGHTACTIVEMAX)
@NLconstraint(m,[t=1:T],HeightActive[t]>=HEIGHTACTIVEMIN)

# Charge and Discharge Rate Limits
@NLconstraint(m,[t=1:T],DischargeRate[t]<=DISCHARGERATEMAX)
@NLconstraint(m,[t=1:T],ChargeRate[t]<=CHARGERATEMAX)
#Energy Balance Constraints
#@NLconstraint(m,[t=1:T],WINDENERGY[t]==SellWindEnergy[t]+SaveWindEnergy[t])

@NLobjective(m,Max,sum(((WINDENERGY[t]+HydroElectricityGen[t]-PumpEnergy[t])*WINDPRICE[t]) for t=1:T))
solve(m)

###################################################################
println("Maximized Revenue ",getobjectivevalue(m))
#println("Wind Energy We Sell ",getvalue(SellWindEnergy[1:5]))
println("Hydro Electricity Sold ",getvalue(HydroElectricityGen[1:5]))
#println("Grid Energy Purchased ",getvalue(GridEnergy[1:5]))
println("Passive Storage Height ",getvalue(HeightPassive[1:5]))
println("Active Storage Height ",getvalue(HeightActive[1:5]))

####################################################################

<<<<<<< HEAD
GridEnergyPlot = plot(1:T,
    getvalue(HydroElectricityGen[1:T]),
  ylabel = "Hydro electricity generation (MWh)",
 xlabel = "Time")
display(GridEnergyPlot)

WindSoldPlot = plot(1:T,
    getvalue(PumpEnergy[1:T]),
   seriestype= :bar,
  ylabel = "Pump energy used (MWh)",
 xlabel = "Time")
display(WindSoldPlot)

HeightPassivePlot = plot(1:T,
        getvalue(HeightPassive[1:T]),
        seriestype= :bar,
        ylabel = "Height Passive (m)",
        xlabel = "Time")

display(HeightPassivePlot)

HeightActivePlot = plot(1:T,
        getvalue(HeightActive[1:T]),
        seriestype= :bar,
        ylabel = "Height Active (m)",
        xlabel = "Time")

display(HeightActivePlot)
=======

PricePlot = plot(1:T,
data.USD_per_MWh[1:T],
ylabel = "Price (USD/MWh)",
 xlabel = "Time",
 title="Locational Marginal Price",
 lw = 2)

 HeightPlot = plot(1:T,
         [getvalue(HeightActive[1:T]), getvalue(HeightPassive[1:T])],
         ylabel = "Height (m)",
         xlabel = "Time",
         label=["Active Reservoir Height" "Passive Reservoir Height"],
         lw = 2)

TimeEnergyPlot = plot(1:T,
    [ getvalue(GridEnergy[1:T]), getvalue(SaveWindEnergy[1:T]), getvalue(HydroElectricityGen[1:T])],
  ylabel = "Energy (MWh)",
 xlabel = "Time",
 title="Energy Sale and Purchase",
 label=["Energy Purchased from Grid" "Energy Stored" "Hydro Energy Sold"],
 lw = 2)

#getvalue(SellWindEnergy[1:T]),
#"Wind Energy Sold"

display(plot(TimeEnergyPlot, HeightPlot))





#display(HeightPassivePlot)

#HeightActivePlot = plot(1:T,
#        getvalue(HeightActive[1:T]),
#        seriestype= :bar,
#        ylabel = "Height Active (m)",
#        xlabel = "Time")


#WindSoldPlot = plot(1:T,
#    getvalue(SellWindEnergy[1:T]),
#   seriestype= :bar,
#  ylabel = "Wind Energy Sold (MWh)",
# xlabel = "Time")
#display(WindSoldPlot)
>>>>>>> 892f6338bbec0f56dc7aa5bd10b8c1351ff8b335
