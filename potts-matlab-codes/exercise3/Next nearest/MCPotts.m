%  Monte Carlo for 2D Potts Model
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
%       Q = the range of possible spin values (1 to Q)
%       Es= stored energy
%
%  output:
%       en = energy at each MCS
%       so = initial configuration
%       s1 = configuration at 1/3 of the run
%       s2 = configuration at 2/3 of the run
%       s = final configuration
%       de = %error in final energy (a test to make sure it works)
%
function[en,so,s1,s2,s,de]= MCPotts(n,nstep,q)
n = 30;
nstep= 100;
q = 4;
Es=0;
% create initial lattice of spins
N=n*n;
s = ceil(q.*rand(n));
% save initial structure
so = s;
np = floor(nstep/3);
%inline Kronecker delta
delta = inline('n==0');
% calculate initial energy and polarization
energy = 0;
for i=1:n
    for j=1:n
        sij = s(i,j);
        sig = delta(sij-s(pbc(i+1,n),j)) + delta(sij-s(i,pbc(j+1,n))) + delta(sij-s(pbc(i-1,n),j)) + delta(sij-s(i,pbc(j-1,n))) + delta(sij-s(pbc(i+1,n),pbc(j+1,n))) + delta(sij-s(pbc(i-1,n),pbc(j+1,n))) + delta(sij-s(pbc(i-1,n),pbc(j-1,n))) + delta(sij-s(pbc(i+1,n),pbc(j-1,n)));
        energy = energy + (8-sig)/2 + Es*sig;
    end
end
%  initialize for averages and output
en=zeros(nstep,1);
MCS = 0;
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
    sijo = s(i,j);
    sij = sijo;
    sij = sij + ceil((q-1)*rand);
    sij = sij - q*floor((sij-1)/q);
    sign = delta(sij-s(pbc(i+1,n),j)) + delta(sij-s(i,pbc(j+1,n))) + delta(sij-s(pbc(i-1,n),j)) + delta(sij-s(i,pbc(j-1,n))) + delta(sij-s(pbc(i+1,n),pbc(j+1,n))) + delta(sij-s(pbc(i-1,n),pbc(j+1,n))) + delta(sij-s(pbc(i-1,n),pbc(j-1,n))) + delta(sij-s(pbc(i+1,n),pbc(j-1,n)));
    sigo = delta(sijo-s(pbc(i+1,n),j)) + delta(sijo-s(i,pbc(j+1,n))) + delta(sijo-s(pbc(i-1,n),j)) + delta(sijo-s(i,pbc(j-1,n))) + delta(sijo-s(pbc(i+1,n),pbc(j+1,n))) + delta(sijo-s(pbc(i-1,n),pbc(j+1,n))) + delta(sijo-s(pbc(i-1,n),pbc(j-1,n))) + delta(sijo-s(pbc(i+1,n),pbc(j-1,n)));
    del = -(sign - sigo);
% Metropolis algorithm:  only change state if accepted
% if rejected, state stays the same but is added to list 
% of configurations
    if del <= 0 
        s(i,j) = sij;
        energy = energy + del;
    end
%  every MCS add energy and magnetization to list
    if mod(k,N)==0
        MCS = MCS + 1;
        en(MCS) = energy;
        if MCS == np 
            s1 = s;
        end
        if MCS == 2*np
            s2 = s;
        end
    end
end

% test the final energy to make sure no errors in updating
entest = 0;
for i=1:n
    for j=1:n
        sij = s(i,j);
        sig = delta(sij-s(pbc(i+1,n),j)) + delta(sij-s(i,pbc(j+1,n))) + delta(sij-s(pbc(i-1,n),j)) + delta(sij-s(i,pbc(j-1,n))) + delta(sij-s(pbc(i+1,n),pbc(j+1,n))) + delta(sij-s(pbc(i-1,n),pbc(j+1,n))) + delta(sij-s(pbc(i-1,n),pbc(j-1,n))) + delta(sij-s(pbc(i+1,n),pbc(j-1,n)));
        entest = entest + (8-sig)/2 + Es*sig;
    end
end
de = (energy-entest)/entest;

% Displaying the results

subplot(2,2,1)
imagesc(so);axis equal;axis off;title('initial','fontsize',18);
colormap('default');

subplot(2,2,2)
imagesc(s1),axis equal;axis off;title('1/3 of the run','fontsize',18);

subplot(2,2,3)
imagesc(s);axis equal;axis off;title('final','fontsize',18);  

subplot(2,2,4)
plot(en);xlabel('MCS','fontsize',18);ylabel('E','fontsize',18);