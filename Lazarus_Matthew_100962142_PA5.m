%Matthew Lazarus 100962142

clear all

clearvars
clearvars -GLOBAL
close all

% set(0,'DefaultFigureWindowStyle','docked')
global C
C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665; %metres (32.1740 ft) per s²

nx = 50;
ny = 50;

G = sparse(nx*ny,nx*ny);
for count = 1:nx*ny
    G(count,count)=1;
end
for col = 1:nx
    if(col~=1 &&col~=nx)
        for row = 1:ny            
            if(count~=1 && row ~=1)
                n = row + (col -1)*ny;
                nyBefore = n-1;
                nyAfter = n+1;
                nxBefore = row+(col-2)*ny;
                nxAfter = row+col*ny;
                if(row>10 && row < 20 && col < 20 && col > 10)
                    G(n,n) = -4; 
                else
                    G(n,n) = -4;
                end
                G(n, nyBefore) =1;
                G(n, nyAfter)=1;
                G(n, nxBefore)=1;
                G(n, nxAfter) =1; 
            end        
        end
    end
end 

figure(1);
title("Sparsity of Matrix G");
spy(G);

[E,D] = eigs(G,9,'SM');

for count = 1:length(D)
    diagD(count) = D(count,count);
end
figure
plot(diagD)
title("Eigenvalues");

V=zeros(ny,nx);
rowCount = 1;
colCount = 1;
for count = 1:size(E,2)
    for count2 = 1:size(E,1)        
        if(rowCount>ny)
           rowCount = 1;
           colCount = colCount+1;
        end
        V(colCount, rowCount) = E(rowCount+(colCount-1)*ny,count);
        rowCount= rowCount + 1;
    end
    figure
    surf(V);
    title("Eigenvector");
    rowCount = 1;
    colCount=1;
end


