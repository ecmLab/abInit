%% This code is to get the band gap value of all calculations
clc;clear;
% Parameters
LaU = [0.0,4.0,4.5,5.0,5.5,6.5,7.0,7.5,8.0,8.5,9.5,10.0];
ZrU = [0.0,0.5,1.5,2.0,2.5,3.5,4.0];
strLaU = ["0.0",'4.0','4.5','5.0','5.5','6.5','7.0','7.5','8.0','8.5','9.5','10.0'];
strZrU = ["0.0",'0.5','1.5','2.0','2.5','3.5','4.0'];
cutEngUp = 8;       % Up bound for the integrated DOS, in eV
dE       = 0.01;    % energy increacement

% Load DOS data of HSE
dos_hse = load('../results/dosData/HSE_atom.txt');
engy  = dos_hse(:,1);
LaDos = dos_hse(:,4);
ind1  = find(engy>0 & LaDos>0);
ind2  = find(engy>0 & engy<cutEngUp);
bg_hse = engy(ind1(1));            % Bandgap of HSE calculation
cutEngRg = cutEngUp - bg_hse;      % The width of the integrated DOS
area_hse = sum(dos_hse(ind2,2:2:8),1)*dE; % Integrated DOS within the range

% Load DOS data of GGA + U
bg_gga   = zeros(length(LaU),length(ZrU));
area_gga_tot = zeros(length(LaU),length(ZrU));
area_gga_La = zeros(length(LaU),length(ZrU));
area_gga_Zr = zeros(length(LaU),length(ZrU));

for iLa = 1 : length(LaU)
    for iZr = 1 : length(ZrU)
        str = strcat('../results/dosData/','La',strLaU(iLa),'Zr',strZrU(iZr),'_atom.txt');
        try
% Compute the band gap
          dos_gga = load(str);
          engy  = dos_gga(:,1);
          LaDos = dos_gga(:,4);
          ind1  = find(engy>0 & LaDos>0);
          bg_gga(iLa,iZr) = engy(ind1(1));

% Compute the area error
          cutEngUp = bg_gga(iLa,iZr) + cutEngRg;  % The uplimit of GGA energy by keeping the same energy range
          ind2  = find(engy>0 & engy<cutEngUp);
          tmp   = sum(dos_gga(ind2,2:2:8),1);     % Area of each energy within the energy range
          area_gga_La(iLa,iZr) = tmp(2)*dE;
          area_gga_Zr(iLa,iZr) = tmp(3)*dE;
          area_gga_tot(iLa,iZr) = sum(tmp)*dE;
          
        catch ME
          fprintf('WEBREAD without success: %s\n', ME.message);
          continue;
        end
    end
end

%% Plot
LaU_y = repmat(LaU(2:end)',1,length(ZrU));
ZrU_x = repmat(ZrU,length(LaU(2:end)),1);
% Contour plot of Bandgap error, respect to HSE value
figure(1)
% bg_err = (bg_hse - bg_gga(2:end,:)) / bg_hse;
bg_err = (bg_hse - bg_gga(2:end,:));
contourf(ZrU_x,LaU_y,bg_gga(2:end,:))

% Contour plot of total energy area error, respect to HSE value
figure(2)
% area_tot_err = (area_gga_tot(2:end,:) - sum(area_hse)) / sum(area_hse);
area_tot_err = area_gga_tot(2:end,:) - sum(area_hse);
contourf(ZrU_x,LaU_y,abs(area_tot_err))

% Contour plot of La energy area error, respect to HSE value
figure(3)
% area_la_err = (area_gga_La(2:end,:) - area_hse(2)) / area_hse(2);
area_la_err = abs(area_gga_La(2:end,:) - area_hse(2));

contourf(ZrU_x,LaU_y,area_la_err)

% Contour plot of Zr energy area error, respect to HSE value
figure(4)
% area_zr_err = (area_gga_Zr(2:end,:) - area_hse(3)) / area_hse(3);
area_zr_err = abs(area_gga_Zr(2:end,:) - area_hse(3));
contourf(ZrU_x,LaU_y,area_zr_err)

% Contour plot of Zr+La energy area error
aa_la = area_la_err;
aa_la(aa_la>0.2) = 0.4;
figure(5)
contourf(ZrU_x,LaU_y,aa_la)

aa_zr = area_zr_err;
aa_zr(aa_zr>0.2) = 0.4;
figure(6)
contourf(ZrU_x,LaU_y,aa_zr)

