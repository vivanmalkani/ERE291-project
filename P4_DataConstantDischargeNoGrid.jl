#####################################
########### P4 Data Clean ###########
#####################################
####### ADMIN AND PACKAGES ######

# Set path to solvers
# Replace with appropriate path where I saved solvers in the beginning of the course.
PATH_TO_SOLVERS = ENV["ERE291_SOLVERS"]

# Boot up packages
using JuMP
using Cbc
using Pkg
using DataFrames


#Now use knitro
ampl_solver="snopt"


using CSV

data = CSV.read("P4_test_data_wind_price.csv")
#test test test
#Let's try this

#display(data)
#Wind Data
WINDENERGY=data[:,4] #[MWh]
WINDPRICE=data[:,3] #[$/MWh]


#Parameters
NUMBERPASSIVE=4
HEIGHTDROP=200 #[m]
STORAGEACTIVEAREA=pi*1600 #[m^2]
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
#T=100
TIMESTEP=1 #[hr]

m = Model(solver=CbcSolver())

#Note that changing expressions use LowerCase while unchanging parameters are all UpperCase

@variable(m,HeightPassive[t=1:T]>=0)
#@variable(m,HeightActive[t=1:T]>=0,start=HEIGHTACTIVEMAX)

#@variable(m,DischargeRate[t=1:T]>=0,start=DISCHARGERATEMAX)
#@variable(m,GridEnergy[t=1:T]>=0,start=0.1) #[MWH]
#@variable(m, SellWindEnergy[t=1:T] >= 0,start=1)
@variable(m, SaveWindEnergy[t=1:T] >= 0)
@variable(m,ChooseDischarge[t=1:T],Bin)
#@variable(m, StoreWindEnergyOut[t=1:T] >= 0,start=1)
#@variable(m,InStorage[t=1:T+1] >=0,start=1)

#Radius of Penstock Pipe (transfer radius)


#Discharge Rate Equation

#@NLexpression(m,HeightPassive[t=1:T],((DischargeRate[t]/(3600*NUMBERPASSIVE*pi*(TRANSFERRADIUS^2)*DISCHARGECOEFFICIENT))^2)/(2*GRAVITY))
#@NLexpression(m,DischargeRate[t=1:T],3600*NUMBERPASSIVE*pi*(TRANSFERRADIUS^2)*DISCHARGECOEFFICIENT*sqrt(2*GRAVITY*HeightPassive[t])) #[m^3/hr]
@expression(m,HeightActive[t=1:T],HEIGHTACTIVEMAX-(NUMBERPASSIVE*HeightPassive[t]*STORAGEPASSIVEAREA)/STORAGEACTIVEAREA)
#Charge Rate
#Need an expression for PumpEnergy
#@NLexpression(m,PumpEnergy[t=1:T],GridEnergy[t]+SaveWindEnergy[t]) #[MWh]

@expression(m,ChargeRate[t=1:T],3600*10^6*SaveWindEnergy[t]*OVERALLPUMPEFFICIENCY*HEADLOSSFRACTION/(DENSITYWATER*GRAVITY*HEIGHTDROP)) #[m^3/hr]

#Electricity Generation
@expression(m,HydroElectricityGen[t=1:T],ChooseDischarge[t]*OVERALLTURBINEEFFICIENCY*DISCHARGERATEMAX*DENSITYWATER*GRAVITY*HEIGHTDROP/(3.6*10^9))

#This fulfills degrees of freedom and the wind energy balance.
@expression(m,SellWindEnergy[t=1:T],WINDENERGY[t]-SaveWindEnergy[t])

#energy balance constraint
#@NLconstraint(m, [t = 1:T], (SellWindEnergy[t] + SaveWindEnergy[t] + HydroElectricityGen[t]) >= (GridEnergy[t] + WINDENERGY[t]) )

#Initial Water Depths
@constraint(m,HeightPassive[1]==HEIGHTPASSIVEMIN)

#Should be unnecessary from NLexpression.
#@NLconstraint(m,HeightActive[1]==HEIGHTACTIVEMAX-(NUMBERPASSIVE*HEIGHTPASSIVEMIN*STORAGEPASSIVEAREA)/STORAGEACTIVEAREA)

#Water Depth States
@constraint(m,[t=1:T-1],HeightPassive[t+1]==HeightPassive[t]-ChooseDischarge[t]*DISCHARGERATEMAX*TIMESTEP/(NUMBERPASSIVE*STORAGEPASSIVEAREA)+ChargeRate[t]*TIMESTEP/(NUMBERPASSIVE*STORAGEPASSIVEAREA))
###This constraint isn't needed based on the earlier volume balance in NL expression for Active height###
#@NLconstraint(m,[t=1:T-1],HeightActive[t+1]==HeightActive[t]-DischargeRate[t]*TIMESTEP/(STORAGEACTIVEAREA)-ChargeRate[t]*TIMESTEP/(STORAGEACTIVEAREA))

#Water Depth Limits
@constraint(m,[t=1:T],HeightPassive[t]<=HEIGHTPASSIVEMAX)
@constraint(m,[t=1:T],HeightPassive[t]>=HEIGHTPASSIVEMIN)

@constraint(m,[t=1:T],HeightActive[t]<=HEIGHTACTIVEMAX)
@constraint(m,[t=1:T],HeightActive[t]>=HEIGHTACTIVEMIN)

# Charge and Discharge Rate Limits
#@NLconstraint(m,[t=1:T],DischargeRate[t]<=DISCHARGERATEMAX)
@constraint(m,[t=1:T],ChargeRate[t]<=CHARGERATEMAX)
#Energy Balance Constraints
#@NLconstraint(m,[t=1:T],WINDENERGY[t]==SellWindEnergy[t]+SaveWindEnergy[t])

@objective(m,Max,sum(SellWindEnergy[t]*WINDPRICE[t] for t=1:T)+sum(HydroElectricityGen[t]*WINDPRICE[t] for t=1:T))
solve(m)

###################################################################
println("--------------------------------------")
println("Maximized Revenue ",getobjectivevalue(m))
println("Wind Energy We Sell ",getvalue(SellWindEnergy[1:5]))
println("Hydro Electricity Sold ",getvalue(HydroElectricityGen[1:5]))
#println("Grid Energy Purchased ",getvalue(GridEnergy[1:5]))
println("Passive Storage Height ",getvalue(HeightPassive))
println("Active Storage Height ",getvalue(HeightActive[1:5]))
println("Choose to Discharge ", getvalue(ChooseDischarge))
####################################################################
