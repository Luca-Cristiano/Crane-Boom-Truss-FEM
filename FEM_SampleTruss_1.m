% FEM Sample Truss

% Clear Command Window and Workspace
clear all, clc;

% Knowns - P represents the load, E is the elastic modulus, A is the beam cross-sectional area
P = 1;
E = 1;
A = 1;
L = [12; 6*sqrt(3); 6]; %L array is hand calculated befrore simulation.
theata = [0; 150; 60]; %theata array is hand calculated befrore simulation.

% K coeffients
k = E.*A./L;

% Generate K Matricies
C = cosd(theata);
S = sind(theata);
CS = C.*S; %".*" is a known function in matlab which multiplies matrix C and S element by element and returns the result in CS

% Calculating the stiffness matrix for member 1
k1 = k(1) * [C(1)^2, CS(1), -C(1)^2, -CS(1), 0, 0;
             CS(1), S(1)^2, -CS(1), -S(1)^2, 0, 0;
             -C(1)^2, -CS(1), C(1)^2, CS(1), 0, 0;
             -CS(1), -S(1)^2, CS(1), S(1)^2, 0, 0;
             0, 0, 0, 0, 0, 0;
             0, 0, 0, 0, 0, 0];

% Calculating the stiffness matrix for member 2
k2 = k(2) * [0, 0, 0, 0, 0, 0;
             0, 0, 0, 0, 0, 0;
             0, 0, C(2)^2, CS(2), -C(2)^2, -CS(2);
             0, 0, CS(2), S(2)^2, -CS(2), -S(2)^2;
             0, 0, -C(2)^2, -CS(2), C(2)^2, CS(2);
             0, 0, -CS(2), -S(2)^2, CS(2), S(2)^2];
             
% Calculating the stiffness matrix for member 3
k3 = k(3) * [C(3)^2, CS(3), 0, 0, -C(3)^2, -CS(3);
             CS(3), S(3)^2, 0, 0, -CS(3), -S(3)^2;
             0, 0, 0, 0, 0, 0;
             0, 0, 0, 0, 0, 0;
             -C(3)^2, -CS(3), 0, 0, C(3)^2, CS(3);
             -CS(3), -S(3)^2, 0, 0, CS(3), S(3)^2];

%assembling the global stiffness matrix 
K = k1 + k2 + k3;

% Solve for unknown displacements (u2x, u3x, u3y)
u = [K(3,3), K(3,5:6); K(5,3), K(5,5:6); K(6,3), K(6,5:6)] \  [0; P; 0];

% Create entire displacemnet matrix U (u1x = u1y = u2y = 0)
U = [0; 0; u(1); 0; u(2); u(3)];

% Solve reation forces
F = K * U;

% Solve member forces
fe1 = k(1) * ((U(3)-U(1))*C(1) + (U(4)-U(2))*S(1));
fe2 = k(2) * ((U(5)-U(3))*C(2) + (U(6)-U(4))*S(2));
fe3 = k(3) * ((U(5)-U(2))*C(3) + (U(6)-U(1))*S(3));

Fe = [fe1;fe2;fe3];

% Format and print results
fprintf('Nodal Displacements:\n');
j = 1;
for i = 1:2:length(U)
    fprintf('u%ix = %.3f P/EA\n', j, U(i));
    fprintf('u%iy = %.3f P/EA\n', j, U(i+1));
    j = j+1;
end

fprintf('\nNodal Forces:\n');
j = 1;
for i = 1:2:length(F)
    fprintf('f%ix = %.3f P\n', j, F(i));
    fprintf('f%iy = %.3f P\n', j, F(i+1));
    j = j+1;
end

fprintf('\nMember Forces:\n');
for i = 1:length(Fe)
    fprintf('f(%i) = %.3f P\n', i, Fe(i));
end
