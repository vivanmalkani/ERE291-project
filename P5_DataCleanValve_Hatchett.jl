#####################################
########### P4 Data Clean ###########
#####################################
####### ADMIN AND PACKAGES ######

# Set path to solvers
# Replace with appropriate path where I saved solvers in the beginning of the course.
PATH_TO_SOLVERS = ENV["ERE291_SOLVERS"]

# Boot up packages
using JuMP
using Statistics
using AmplNLWriter
using Pkg
using Plots
using DataFrames


#Now use knitro
ampl_solver="knitro"


using CSV

data = CSV.read("P5_Hatchett.csv")

#display(data)
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
LAGTIME=30 #[s]
DISCHARGERATEMAX=2*3600 #[m^3/hr]
CHARGERATEMAX=2*3600
DENSITYWATER=1000 #[kg/m^3]
DISCHARGECOEFFICIENT=0.6
GRAVITY=9.81 #m/s^2

#Solved from maximum Discharge Rate and maximum passive height.
TRANSFERRADIUS=0.0973 #m^3
T=length(WINDENERGY)
TIMESTEP=1 #[hr]

m=Model(solver=AmplNLSolver(joinpath(PATH_TO_SOLVERS,ampl_solver),["outlev=2"]))

#Note that changing expressions use LowerCase while unchanging parameters are all UpperCase

@variable(m,HeightPassive[t=1:T]>=0,start=HEIGHTPASSIVEMAX)
#@variable(m,HeightActive[t=1:T]>=0,start=HEIGHTACTIVEMAX)
@variable(m,GridEnergy[t=1:T]>=0,start=20)
#@variable(m, SellWindEnergy[t=1:T] >= 0,start=1)
@variable(m, SaveWindEnergy[t=1:T] >= 0,start=15)
#@variable(m, StoreWindEnergyOut[t=1:T] >= 0,start=1)
#@variable(m,InStorage[t=1:T+1] >=0,start=1)
@variable(m,DischargeFlow[t=1:T]>=0, start=0)
#Radius of Penstock Pipe (transfer radius)

#Active Reservoir Height
@NLexpression(m,HeightActive[t=1:T],HEIGHTACTIVEMAX-(NUMBERPASSIVE*HeightPassive[t]*STORAGEPASSIVEAREA)/STORAGEACTIVEAREA)

#Discharge Rate Equation
@NLexpression(m,DischargeRate[t=1:T],3600*NUMBERPASSIVE*pi*(TRANSFERRADIUS^2)*DISCHARGECOEFFICIENT*sqrt(2*GRAVITY*HeightPassive[t])) #[m^3/hr]

#Charge Rate
#Need an expression for PumpEnergy
@NLexpression(m,PumpEnergy[t=1:T],GridEnergy[t]+SaveWindEnergy[t])

@NLexpression(m,ChargeRate[t=1:T],3600*10^6*OVERALLPUMPEFFICIENCY*HEADLOSSFRACTION*PumpEnergy[t]/(DENSITYWATER*GRAVITY*HEIGHTDROP)) #[m^3/hr]

#Electricity Generation
@NLexpression(m,HydroElectricityGen[t=1:T],OVERALLTURBINEEFFICIENCY*DischargeFlow[t]*DENSITYWATER*GRAVITY*HEIGHTDROP/(3.6*10^9))

#This fulfills degrees of freedom and the wind energy balance.
@NLexpression(m,SellWindEnergy[t=1:T],WINDENERGY[t]-SaveWindEnergy[t])

#energy balance constraint
#@NLconstraint(m, [t = 1:T], (SellWindEnergy[t] + SaveWindEnergy[t] + HydroElectricityGen[t]) >= (GridEnergy[t] + WINDENERGY[t]) )

#Initial Water Depths
@NLconstraint(m,HeightPassive[1]==HEIGHTPASSIVEMAX)
#@NLconstraint(m,HeightActive[1]==HEIGHTACTIVEMIN)

#Water Depth States
@NLconstraint(m,[t=1:T-1],HeightPassive[t+1]==HeightPassive[t]+DischargeFlow[t]/(NUMBERPASSIVE*STORAGEPASSIVEAREA)-ChargeRate[t]*TIMESTEP/(NUMBERPASSIVE*STORAGEPASSIVEAREA))
#@NLconstraint(m,[t=1:T-1],HeightActive[t+1]==HeightActive[t]-DischargeRate[t]*TIMESTEP/(STORAGEACTIVEAREA)-ChargeRate[t]*TIMESTEP/(STORAGEACTIVEAREA))

#Water Depth Limits
@NLconstraint(m,[t=1:T],HeightPassive[t]<=HEIGHTPASSIVEMAX)
@NLconstraint(m,[t=1:T],HeightPassive[t]>=HEIGHTPASSIVEMIN)

@NLconstraint(m,[t=1:T],HeightActive[t]<=HEIGHTACTIVEMAX)
@NLconstraint(m,[t=1:T],HeightActive[t]>=HEIGHTACTIVEMIN)

# Charge and Discharge Rate Limits
#@NLconstraint(m,[t=1:T],DischargeRate[t]<=DISCHARGERATEMAX)
@NLconstraint(m,[t=1:T],DischargeFlow[t]<=DischargeRate[t]*TIMESTEP)
@NLconstraint(m,[t=1:T],ChargeRate[t]<=CHARGERATEMAX)
#Energy Balance Constraints
#@NLconstraint(m,[t=1:T],WINDENERGY[t]==SellWindEnergy[t]+SaveWindEnergy[t])

@NLobjective(m,Max,sum(((SellWindEnergy[t]+HydroElectricityGen[t]-GridEnergy[t])*WINDPRICE[t]) for t=1:T))
solve(m)
### Number of results to print
j=10
###################################################################
println("Maximized Revenue ",getobjectivevalue(m))
println("Wind Energy We Sell ",getvalue(SellWindEnergy[1:j]))
println("Hydro Electricity Sold ",getvalue(HydroElectricityGen[1:j]))
println("Grid Energy Purchased ",getvalue(GridEnergy[1:j]))
println("Passive Storage Height ",getvalue(HeightPassive[1:j]))
println("Active Storage Height ",getvalue(HeightActive[1:j]))
println("Revenue without storage ", sum(WINDPRICE[t]*WINDENERGY[t] for t=1:T))
####################################################################

daily_data = hcat(getvalue(GridEnergy), getvalue(HydroElectricityGen), getvalue(SellWindEnergy), getvalue(HeightActive), getvalue(HeightPassive))
daily_data = convert(DataFrame, daily_data)
daily_data = hcat(daily_data, data)
rename!(daily_data, :x1 => :Grid_Energy)
rename!(daily_data, :x2 => :Hydro_Elec)
rename!(daily_data, :x3 => :Wind_Sold)
rename!(daily_data, :x4 => :Height_Active)
rename!(daily_data, :x5 => :Height_Passive)

daily_data = by(daily_data, :OPR_DT,
Grid_Energy = :Grid_Energy => sum,
Hydro_Elec = :Hydro_Elec => sum,
Wind_Sold = :Wind_Sold => sum,
 Height_Active = :Height_Active => mean,
  Height_Passive = :Height_Passive => mean)

PricePlot = plot(1:T,
data.USD_per_MWh[1:T],
ylabel = "Price (USD/MWh)",
 xlabel = "Time (Hours)",
 label = ["Price varying by hour"],
 title="Locational Marginal Price",
 lw = 2)

 HeightPlot = plot(1:length(daily_data.OPR_DT),
         [daily_data.Height_Active, daily_data.Height_Passive],
         ylabel = "Height (m)",
         xlabel = "Day of Month (January)",
         title="Daily Average Height Trend",
         label = ["Height (Active)" "Height (Passive)"],
         lw = 2)

EnergyPlot = plot(1:length(daily_data.OPR_DT),
                 [daily_data.Grid_Energy, daily_data.Hydro_Elec, daily_data.Wind_Sold],
                 ylabel = "Energy (MWh)",
                 xlabel = "Day of Month (January)",
                 title="Daily Energy Transactions",
                 label = ["Grid Energy" "Hydro Electricity" "Wind Energy Sold"],
                 lw = 2)

png(EnergyPlot, "EnergyPlot")
png(HeightPlot, "HeightPlot")
png(PricePlot, "PricePlot")





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
