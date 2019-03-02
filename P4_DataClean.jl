#####################################
########### P4 Data Clean ###########
#####################################

using CSV

data = CSV.read("2018Jan_PRC_LMP_DAM.csv")
#test test test
#Let's try this




display(data)


#Storage Parameters
NUMBERWINDMILL=4
HEIGHTDROP=200 #[m]
HEIGHTPASSIVEMAX=40 #[m]
HEIGHTPASSIVEMIN=0.5 #[m]
HEIGHTACTIVEMAX=???
HEIGHTACTIVEMIN=0.5 #[m]
STORAGEPASSIVEMAX=6057 #[m^3]
STORAGEACTIVEMAX=160,000-STORAGEPASSIVEMAX #[m^3]
OVERALLTURBINEEFFICIENCY=0.56 #This includes the Penstock, Maniforld, Nozzle, Runner (a bit low), drive, and generator (also lowball)
OVERALLPUMPEFFICIENCY=0.91 #This reflects waterways losses, pump losses, motor losses, transformer losses.
STORAGEPASSIVEAREA=PASSIVESTORAGEMAX/HEIGHTPASSIVE #[m^2]
HEADLOSSFRACTION=0.05 #Effective height pressure lost from friction in penstock.
LAGTIME=30 #[s]
FLOWLIMIT=2 #[m^3/s]
DENSITYWATER=1000 #[kg/m^3]
DISCHARGECOEFFICIENT=0.6
