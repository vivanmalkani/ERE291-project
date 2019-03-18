
%3D plot

%Import the Data
MaximizedRevenue=csvread("EfficiencySensitivityHatchett.csv");

%same data range used in julia
TurbineEfficiencyRange=linspace(0.5,0.95,15);
PumpEfficiencyRange=linspace(0.5,0.95,15);

%Mesh Plot
figure(1)
OptimalRevenue=mesh(TurbineEfficiencyRange,PumpEfficiencyRange,MaximizedRevenue);
xlabel ('Overall Turbine Efficiency');
ylabel ('Overall Pump Efficiency');
zlabel('Optimized Project Revenue')
title ('Objective Function Sensitivity to Pump and Turbine Efficiency');
print('EfficiencySensitivity3DPlot','-dpng')

