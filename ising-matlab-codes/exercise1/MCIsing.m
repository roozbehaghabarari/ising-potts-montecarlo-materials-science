%  Monte Carlo for 2D Ising Model
%            Coded by                                  
%       Roozbeh Aghabarari       
%  
%           Contact Me
%     ro.aghabarari@gmail.com     
%    www.roozbehaghabarari.com
%
%  input:
%       n = size of cell.  total number of sites = n^2
%       nstep = number of Monte Carlo Steps (MCS)
%       pol = polarization of initial lattice
%       temp = temperature
%       B = magnetic field
%
%  output:
%       aen = average energy
%       amag = average magnetization
%       en = energy at each MCS
%       magn = magnetization at each MCS
%       s = final configuration
%       so = initial configuration
%       de = %error in final energy (a test to make sure it works)
%       pacc = % of trials accepted
%
function[aen,amag,en,magn,s,so,de,pacc]= MCIsing(n,nstep,pol,temp,B)
% create initial lattice of spins
n=20;
nstep=100;
pol=1;
temp=0.25;
B=0;
N=n*n;
s = initIsing(n,pol);
% save initial structure
so = s;
% calculate initial energy and polarization
energy = 0;
mag = 0;
for i=1:n
    for j=1:n
        sig = s(pbc(i+1,n),j) + s(i,pbc(j+1,n)) + s(pbc(i-1,n),j) + s(i,pbc(j-1,n));
        energy = energy - s(i,j)*(sig/2 + B);
        mag = mag + s(i,j);
    end
end
%  initialize for averages and output
en=zeros(nstep,1);
magn = zeros(nstep,1);
MCS = 0;
aen = 0;
amag = 0;
iacc = 0;
% define periodic boundary conditions in pbc
%
% do Monte Carlo
%  number of Monte Carlo steps = nstep
%  number of configurations = N*nstep
nconfig = N*nstep;
for k=1:nconfig
%  pick a site to change
    i = ceil(rand*n);
    j = ceil(rand*n);
% flip the spin and calculate energy change
    sij = s(i,j);
    sig = s(pbc(i+1,n),j) + s(i,pbc(j+1,n)) + s(pbc(i-1,n),j) + s(i,pbc(j-1,n));
    del = 2*sij*(sig + B);
% Metropolis algorithm:  only change state if accepted
% if rejected, state stays the same but is added to list 
% of configurations
    if del < 0 |  exp(-del/temp) > rand
        s(i,j) = -sij;
        energy = energy + del;
        mag = mag - 2*sij;
        iacc = iacc + 1;
    end
%  accumulate for avearging
    aen = aen + energy;
    amag = amag + mag;
%  every MCS add energy and magnetization to list
    if mod(k,N)==0
        MCS = MCS + 1;
        en(MCS) = energy/N;
        magn(MCS) = mag/N;
    end
end
%  calculate average energy and magnetization
aen = aen/nconfig;
amag = amag/nconfig;
% output energy and magnetization per site
aen = aen/N;
amag = amag/N;
% test the final energy to make sure no errors in updating
entest = 0;
for i=1:n
    for j=1:n
        sig = s(pbc(i+1,n),j) + s(i,pbc(j+1,n)) + s(pbc(i-1,n),j) + s(i,pbc(j-1,n));
        entest = entest - s(i,j)*(sig/2 + B);
    end
end
de = (energy-entest)/entest;
% calculate the percent acceptance of the moves
pacc=100*iacc/nconfig;

% Displaying the results

subplot(2,2,1)
imagesc(so);colormap(gray);axis equal;axis off;title('initial','fontsize',18)

subplot(2,2,2)
imagesc(s);colormap(gray);axis equal;axis off;title('final','fontsize',18)

subplot(2,2,3)
plot(en);xlabel('MCS','fontsize',18);ylabel('E','fontsize',18)    

subplot(2,2,4)
plot(magn);xlabel('MCS','fontsize',18);ylabel('M','fontsize',18)