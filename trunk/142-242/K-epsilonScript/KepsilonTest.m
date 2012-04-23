%SETTINGS

%the dimensions of the K-epsilon matrix
m = 10; %rows
n = 10; %cols

A = 1;    %normalization
dT = 1;   %epsilon time step
hBar = 1; %Planck
M = 1;    %particle mass


%initialize matrix to zero
Ke=complex(zeros(m,n));



for I = 1:m
    for J = 1:n
        %compute i,j element of K-epsilon matrix
        Ke(I,J) = 1/A*exp( (1i)*dT/hBar*Legrangian(I,J) );
    end
end
