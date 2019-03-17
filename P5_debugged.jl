#####################################
########### P5 Data Clean ###########
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
ampl_solver="snopt"


using CSV

#Import Wind Data (Energy and Price)
data = CSV.read("P4_test_data_wind_price.csv")

#display(data)
#Wind Data
WINDENERGY=data[:,4]
WINDPRICE=data[:,3]


#Active (lower) and Passive(higher) Parameters
NUMBERPASSIVE=4 #Number of (passive) Windmill Storage Units
HEIGHTDROP=200 #[m] the height from base of passive storage units to lower reservoir.
STORAGEACTIVEAREA=pi*40^2 #[m^2] The cross sectional area of the lower (active) reservoir.
HEIGHTPASSIVEMAX=40 #[m] The maximum height water can fill the upper, passive reservoirs.
HEIGHTPASSIVEMIN=0.5 #[m] The minimum height water can fill the upper, passive reservoirs.
STORAGEPASSIVEMAX=6057 #[m^3] The maximum volume water can fill one upper, passive reservoirs.
STORAGEACTIVEMAX=153943 #[m^3] The maximum volume water can fill the lower, active reservoir.
HEIGHTACTIVEMAX=30.626 #[m] The maximum height water can fill the lower, active reservoir.
HEIGHTACTIVEMIN=0.5 #[m] The minimum height water can fill the lower, active reservoir.
STORAGEPASSIVEAREA=151.425 #[m^2]

#Estimated Based on Oregon State Data
#This can be changed in sensitivity analysis, but we expect around this much. We base this on what would be used for a project this size, exact details are unknown.
OVERALLTURBINEEFFICIENCY=0.89 #This includes the Penstock, Turbine, Nozzle, Generator, and Motor Efficiency.
OVERALLPUMPEFFICIENCY=0.865 #This reflects Penstock, pump, motor, and transformer efficiencies.

#Pretty much never changes
GRAVITY=9.81 #m/s^2
DENSITYWATER=1000 #[kg/m^3]

#Known Penstock Dimensions for Max Bogl Project
PENSTOCKRADIUS=0.8 #[m] can range up to 1.25 m.
PENSTOCKLENGTH=3536 #[m] From Google Earth Estimation of the actual site! Average path traveled by water.
PENSTOCKROUGHNESS=140 #Unitless, Polyethylene roughness

#What we plot over
#T=100
T=length(WINDENERGY)
TIMESTEP=1 #[hr] Everything in terms of hours!

m=Model(solver=AmplNLSolver(joinpath(PATH_TO_SOLVERS,ampl_solver),["outlev=2"]))

#Note that changing expressions use LowerCase while unchanging parameters are all UpperCase

#The water height of the passive reservoirs, chosen as the centerpoint of the flow balances.
@variable(m,HeightPassive[t=1:T]>=0,start=HEIGHTPASSIVEMAX)

#Energy we purchase from the grid [MWh]
@variable(m,GridEnergy[t=1:T]>=0,start=0)

#WindEnergy we use to pump water! [MWh]
@variable(m, SaveWindEnergy[t=1:T] >= 0,start=0)

#Follows discharge rate maximums, but can always be lower or zero. Chosen by solver! [m^3]
@variable(m,DischargeFlow[t=1:T]>=0, start=23900)

@expression(m,HYDRAULICRADIUS,PENSTOCKRADIUS/2) #m This is D/4
#Active Reservoir Height [m]
@NLexpression(m,HeightActive[t=1:T],HEIGHTACTIVEMAX-(NUMBERPASSIVE*HeightPassive[t]*STORAGEPASSIVEAREA)/STORAGEACTIVEAREA)

#Discharge Rate Equation (Hazen-Williams)- the maximum discharge rate for a stored height of water. Takes into account frictional losses.
@NLexpression(m,DischargeRate[t=1:T], 3600*0.849*NUMBERPASSIVE*pi*(PENSTOCKRADIUS^2)*PENSTOCKROUGHNESS*(HYDRAULICRADIUS^0.63)*((HEIGHTDROP+HeightPassive[t])/(PENSTOCKLENGTH))^0.54)

#Charge Rate
@NLexpression(m,PumpEnergy[t=1:T],GridEnergy[t]+SaveWindEnergy[t]) #[Mwh]

@NLexpression(m,ChargeRate[t=1:T],3600*10^6*OVERALLPUMPEFFICIENCY*PumpEnergy[t]/(DENSITYWATER*GRAVITY*(HEIGHTDROP+HeightPassive[t]))) #[m^3/hr]

#Electricity Generation [MWh]
@NLexpression(m,HydroElectricityGen[t=1:T],OVERALLTURBINEEFFICIENCY*DischargeFlow[t]*DENSITYWATER*GRAVITY*(HEIGHTDROP+HeightPassive[t])/(3.6*10^9))

#This fulfills degrees of freedom and the wind energy balance.
@NLexpression(m,SellWindEnergy[t=1:T],WINDENERGY[t]-SaveWindEnergy[t])

#The maximum charge rate is estimated based on a Francis Pump/Turbine pump curve. The relationship is linearized for the region head=200-->240m
@NLexpression(m,CHARGERATEMAX[t=1:T],NUMBERPASSIVE*45*(HEIGHTDROP+HeightPassive[t]))

#Initial Water Depths
@NLconstraint(m,HeightPassive[1]==HEIGHTPASSIVEMAX)
#@NLconstraint(m,HeightActive[1]==HEIGHTACTIVEMIN)

#Water Depth States
@NLconstraint(m,[t=1:T-1],HeightPassive[t+1]==HeightPassive[t]-DischargeFlow[t]/(NUMBERPASSIVE*STORAGEPASSIVEAREA)+ChargeRate[t]*TIMESTEP/(NUMBERPASSIVE*STORAGEPASSIVEAREA))
#@NLconstraint(m,[t=1:T-1],HeightActive[t+1]==HeightActive[t]-DischargeRate[t]*TIMESTEP/(STORAGEACTIVEAREA)-ChargeRate[t]*TIMESTEP/(STORAGEACTIVEAREA))

#Water Depth Limits
@NLconstraint(m,[t=1:T],HeightPassive[t]<=HEIGHTPASSIVEMAX)
@NLconstraint(m,[t=1:T],HeightPassive[t]>=HEIGHTPASSIVEMIN)

@NLconstraint(m,[t=1:T],HeightActive[t]<=HEIGHTACTIVEMAX)
@NLconstraint(m,[t=1:T],HeightActive[t]>=HEIGHTACTIVEMIN)

# Charge and Discharge Rate Limits
#@NLconstraint(m,[t=1:T],DischargeRate[t]<=DISCHARGERATEMAX)
@NLconstraint(m,[t=1:T],DischargeFlow[t]<=DischargeRate[t]*TIMESTEP)
@NLconstraint(m,[t=1:T],ChargeRate[t]<=CHARGERATEMAX[t])

#Objective
@NLobjective(m,Max,sum(((SellWindEnergy[t]+HydroElectricityGen[t]-GridEnergy[t])*WINDPRICE[t]) for t=1:T))

#Solve the model
solve(m)
### Number of results to print
j=10
###################################################################
println("Maximized Revenue ",getobjectivevalue(m))
#println("Wind Energy We Sell ",getvalue(SellWindEnergy[1:j]))
#println("Wind Energy We Save ",getvalue(SaveWindEnergy[1:j]))
#println("Hydro Electricity Sold ",getvalue(HydroElectricityGen[1:T]))
#println("Grid Energy Purchased ",getvalue(GridEnergy[1:j]))
#println("Passive Storage Height ",getvalue(HeightPassive[1:j]))
#println("Active Storage Height ",getvalue(HeightActive[1:j]))
#println("Discharge Trend ",getvalue(DischargeFlow))
#println("Charge Trend ", getvalue(ChargeRate))
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


display(EnergyPlot)
display(HeightPlot)
display(PricePlot)

using DelimitedFiles
check_data = hcat(getvalue(GridEnergy), getvalue(HydroElectricityGen), getvalue(SellWindEnergy), getvalue(HeightActive), getvalue(HeightPassive))
writedlm("checkdata.csv",check_data, ',')
