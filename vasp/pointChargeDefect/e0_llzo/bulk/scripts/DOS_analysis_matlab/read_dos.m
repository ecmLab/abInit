%% This code is to get the total dos and dos of each element

% Load original doscar
fileIn = fopen('DOSCAR','r');    % Input file
nln  = 0; 

% 1. Skip the beginning parts
for itmp = 1 : 5
    stmp = fgetl(fileIn);
end

% 2. get system level information
stmp = fgetl(fileIn);
ntmp = str2num(stmp);
eng_max = ntmp(1);            % The max energy
eng_min = ntmp(2);            % The min energy
nSmp    = ntmp(3);            % The number of sampling points
eng_Femi = ntmp(4);           % The Femi energy

% 3. get the total DOS
dos_tot = zeros(nSmp,5);
for itmp = 1 : nSmp
    stmp = fgetl(fileIn);
    dos_tot(itmp,:) = str2num(stmp);
end

% 4. Get the lm-DOS of all possible orbitals (sx2,px6,dx10,fx14) of each Li elements
dos0_Li = zeros(nSmp,33,n_Li);           % For Li element, 
for itmp = 1 : n_Li
    stmp = fgetl(fileIn);               % transition line
    for jtmp = 1 : nSmp
        stmp = fgetl(fileIn);
        dos0_Li(jtmp,:,itmp) = str2num(stmp);
        
        nln  = nln + 1;
        nlen(nln) = length(str2num(stmp));
    end
end

% 5. Get the lm-DOS of all possible orbitals (sx2,px6,dx10,fx14) of each La elements
dos0_La = zeros(nSmp,33,n_La);           % For La element, 
for itmp = 1 : n_La
    stmp = fgetl(fileIn);               % transition line
    for jtmp = 1 : nSmp
        stmp = fgetl(fileIn);
        dos0_La(jtmp,:,itmp) = str2num(stmp);
        
        nln  = nln + 1;
        nlen(nln) = length(str2num(stmp));
    end
end

% 6. Get the lm-DOS of all possible orbitals (sx2,px6,dx10,fx14) of each Zr elements
dos0_Zr = zeros(nSmp,33,n_Zr);           % For Zr element, 
for itmp = 1 : n_Zr
    stmp = fgetl(fileIn);               % transition line
    for jtmp = 1 : nSmp
        stmp = fgetl(fileIn);
        dos0_Zr(jtmp,:,itmp) = str2num(stmp);
        
        nln  = nln + 1;
        nlen(nln) = length(str2num(stmp));
    end
end

% 7. Get the lm-DOS of all possible orbitals (sx2,px6,dx10,fx14) of each O elements
dos0_O = zeros(nSmp,33,n_O);           % For O element, 
for itmp = 1 : n_O
    stmp = fgetl(fileIn);               % transition line
    for jtmp = 1 : nSmp
        stmp = fgetl(fileIn);
        dos0_O(jtmp,:,itmp) = str2num(stmp);
        
        nln  = nln + 1;
        nlen(nln) = length(str2num(stmp));
    end
end

fclose(fileIn);

