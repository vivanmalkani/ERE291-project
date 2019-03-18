
%3D plot

%Import the Data
MaximizedRevenue=csvread("BaseCaseStartingPoints.csv");

%same data range used in julia
WindEnergyStoreRange=linspace(0,10,10);
GridEnergyRange=linspace(0,10,10);

%Mesh Plot
figure(1)
OptimalRevenue=mesh(WindEnergyStoreRange,GridEnergyRange,MaximizedRevenue);
xlabel ('Initial Value for Storing Wind Energy [MWh]');
ylabel ('Initial Value for Using Grid Energy [MWh]');
zlabel('Optimized Project Revenue')
title ('Objective Function Sensitivity to StartingPoints');
print('StartingPointsSensitivity3DPlot','-dpng')

