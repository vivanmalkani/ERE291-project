#####################################
########### P5 Battery Comparison ###########
#####################################



using JuMP
using Clp
using Statistics
using AmplNLWriter
using Pkg
using Plots
using DataFrames
using CSV

#Import Wind Data (Energy and Price)
data = CSV.read("P5_Hatchett.csv")

#display(data)
#Wind Data
WINDENERGY=data[:,4]
WINDPRICE=data[:,3]

#Capacity of hydro turbines adapted to this battery
Capacity = 15.9 #MW

#Mox Bogl provided storage limit
StorageLimit = 70 #MWh

#starting off empty
InitialStorage = StorageLimit
C_cap = 250000#$/MWh
#Lifetime constant, MWh lifetime stored into battery per MWh of capacity
L_c = 3000 #MWh

#Losses from battery storage
W = 0.0645 #MWh/MWh

#What we plot over, using test case of January
T= 744
TIMESTEP=1 #[hr] Everything in terms of hours!


m=Model(solver=ClpSolver())

#Energy Sold from Wind
@variable(m, DirectSale[1:T] >= 0)

#Energy put into storage
@variable(m, StorageIn[1:T] >= 0)

#Energy sold from storage
@variable(m, StorageOut[1:T] >= 0)

#Energy held in storage in time t
@variable(m, InStorage[1:T+1] >= 0)

######################################
######## Objective Functions #########
######################################

# Cost Maximization
@objective(m, Max, sum(DirectSale[t]*data.USD_per_MWh[t] + (1-W)*StorageOut[t]*data.USD_per_MWh[t] for t = 1:T) -
sum(StorageOut[t] for t=1:T)*(C_cap/L_c)*(1-W)/StorageLimit - sum(StorageOut[t]*W*data.USD_per_MWh[t] for t = 1:T))

########### Constraints #############

# Storage initialization constraint
@constraint(m, InStorage[1] == InitialStorage)

# Max Power Output of Battery
@constraint(m, [t=1:T], StorageOut[t] <= Capacity)

# Storage conservation of energy constraint
@constraint(m, [t=1:T], InStorage[t+1] == InStorage[t] + StorageIn[t] - StorageOut[t])

# Storage size constraint
@constraint(m, [t=1:T], InStorage[t] <= StorageLimit)

# Available power constraints
@constraint(m, [t=1:T], DirectSale[t] == (WINDENERGY[t] - StorageIn[t]))

#Solve the model
solve(m)

###################################################################
println("Maximized Revenue ",getobjectivevalue(m))
println("Revenue Without Storage ", sum(WINDPRICE[t]*WINDENERGY[t] for t=1:T))
####################################################################

daily_data = hcat(getvalue(DirectSale),getvalue(StorageOut),WINDENERGY[1:T,:])
daily_data = convert(DataFrame, daily_data)
daily_data = hcat(daily_data, data[1:T,:])
rename!(daily_data, :x1 => :DirectSale)
rename!(daily_data, :x2 => :StorageOut)
rename!(daily_data, :x3 => :WindEnergy)

daily_data = by(daily_data, :date,
DirectSale = :DirectSale => sum,
StorageOut = :StorageOut => sum,
WindEnergy = :WindEnergy => sum)


BatteryPlot = plot(1:length(daily_data.date),
                 [daily_data.WindEnergy, daily_data.StorageOut, daily_data.DirectSale],
                 ylabel = "Energy (MWh)",
                 xlabel = "Day of Month (January)",
                 title="Battery System Comparison",
                 label = ["Wind Energy Available" "Storage Energy Sold" "Direct Wind Sales"],
                 lw = 2)


display(BatteryPlot)
png(BatteryPlot, "BatteryPlot")
