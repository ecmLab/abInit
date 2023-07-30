%% This code is to analys the DOS and get band structures

% Extract useful DOS from lm-DOS results
% 1. Get dos of Li elements, Li element only has s and p orbitals, with spin up and spin down
dos_Li = zeros(nSmp,5,n_Li);
dos_tot_Li = zeros(nSmp,5);      % Dos of all elements
for itmp = 1 : n_Li
    dos_Li(:,1:3,itmp) = dos0_Li(:,1:3,itmp);   % Energy range, s-spin up, s-spin down
    dos_Li(:,4,itmp) = sum(dos0_Li(:,[4:2:9],itmp),2); % p-spin up
    dos_Li(:,5,itmp) = sum(dos0_Li(:,[5:2:9],itmp),2); % p-spin down
end
dos_tot_Li(:,1) = dos_tot(:,1);
dos_tot_Li(:,2:end) = sum(dos_Li(:,2:end,:),3);

    
% 2. Get dos of La elements, La element has s, p, d, f orbitals, with spin up and spin down
dos_La = zeros(nSmp,9,n_La);
dos_tot_La = zeros(nSmp,9);
for itmp = 1 : n_La
    dos_La(:,1:3,itmp) = dos0_La(:,1:3,itmp);       % Energy range, s-spin up, s-spin down
    dos_La(:,4,itmp)   = sum(dos0_La(:,[4:2:9],itmp),2);   % p-spin up
    dos_La(:,5,itmp)   = sum(dos0_La(:,[5:2:9],itmp),2);   % p-spin down
    dos_La(:,6,itmp)   = sum(dos0_La(:,[10:2:19],itmp),2); % d-spin up
    dos_La(:,7,itmp)   = sum(dos0_La(:,[11:2:19],itmp),2); % d-spin down
    dos_La(:,8,itmp)   = sum(dos0_La(:,[20:2:33],itmp),2); % f-spin up
    dos_La(:,9,itmp)   = sum(dos0_La(:,[21:2:33],itmp),2); % f-spin down
end
dos_tot_La(:,1) = dos_tot(:,1);
dos_tot_La(:,2:end) = sum(dos_La(:,2:end,:),3);

% 3. Get dos of Zr elements, Zr element has s, p, d orbitals, with spin up and spin down
dos_Zr = zeros(nSmp,7,n_Zr);
dos_tot_Zr = zeros(nSmp,7);
for itmp = 1 : n_Zr
    dos_Zr(:,1:3,itmp) = dos0_Zr(:,1:3,itmp);       % Energy range, s-spin up, s-spin down
    dos_Zr(:,4,itmp)   = sum(dos0_Zr(:,[4:2:9],itmp),2);   % p-spin up
    dos_Zr(:,5,itmp)   = sum(dos0_Zr(:,[5:2:9],itmp),2);   % p-spin down
    dos_Zr(:,6,itmp)   = sum(dos0_Zr(:,[10:2:19],itmp),2); % d-spin up
    dos_Zr(:,7,itmp)   = sum(dos0_Zr(:,[11:2:19],itmp),2); % d-spin down
end
dos_tot_Zr(:,1) = dos_tot(:,1);
dos_tot_Zr(:,2:end) = sum(dos_Zr(:,2:end,:),3);

% 4. Get dos of O elements, O element has s, p orbitals, with spin up and spin down
dos_O = zeros(nSmp,5,n_O);
dos_tot_O = zeros(nSmp,5);
for itmp = 1 : n_O
    dos_O(:,1:3,itmp) = dos0_O(:,1:3,itmp);       % Energy range, s-spin up, s-spin down
    dos_O(:,4,itmp)   = sum(dos0_O(:,[4:2:9],itmp),2);   % p-spin up
    dos_O(:,5,itmp)   = sum(dos0_O(:,[5:2:9],itmp),2);   % p-spin down
end
dos_tot_O(:,1) = dos_tot(:,1);
dos_tot_O(:,2:end) = sum(dos_O(:,2:end,:),3);
