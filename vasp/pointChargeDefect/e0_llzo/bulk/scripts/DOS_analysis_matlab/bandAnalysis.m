% clc; clear;
%% This is the matlab file to do band structure analysis
% By Howard on June, 13

% 1. Parameters: 
n_Li = 28;
n_La = 12;
n_Zr = 8;
n_O  = 48;

%% Read in DOS information
% read_dos;

%% Analyze DOS structure
anly_dos

%% Plot
ifg = 0;
% Plot total dos, including spin up and spin down
ifg = ifg + 1;
figure(ifg)
plot(dos_tot(:,2),dos_tot(:,1)-eng_Femi,-dos_tot(:,3),dos_tot(:,1)-eng_Femi)
grid on

% Plot total integrated-dos, including spin up and spin down
ifg = ifg + 1;
figure(ifg)
plot(dos_tot(:,4),dos_tot(:,1),-dos_tot(:,5),dos_tot(:,1))
grid on

% Plot dos of each element
ifg = ifg + 1;
figure(ifg)
plot(sum(dos_tot_Li(:,[2:2:end]),2),dos_tot(:,1)-eng_Femi,'b',sum(dos_tot_La(:,[2:2:end]),2),dos_tot(:,1)-eng_Femi,'g')
hold on
plot(sum(dos_tot_Zr(:,[2:2:end]),2),dos_tot(:,1)-eng_Femi,'r',sum(dos_tot_O(:,[2:2:end]),2),dos_tot(:,1)-eng_Femi,'b')
hold on
plot(-sum(dos_tot_Li(:,[3:2:end]),2),dos_tot(:,1)-eng_Femi,'b',-sum(dos_tot_La(:,[3:2:end]),2),dos_tot(:,1)-eng_Femi,'g')
hold on
plot(-sum(dos_tot_Zr(:,[3:2:end]),2),dos_tot(:,1)-eng_Femi,'r',-sum(dos_tot_O(:,[3:2:end]),2),dos_tot(:,1)-eng_Femi,'b')
hold off
grid on
axis([-40,40,-3,8])
legend('Li','La','Zr','O')