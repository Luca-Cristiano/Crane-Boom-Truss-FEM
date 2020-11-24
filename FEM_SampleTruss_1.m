% FEM Sample Truss

% Clear Command Window and Workspace
clear all, clc;


% Constants for the beam (A is area, E is elastic modulus, L is beam lengths, theta)
A = 1; 
E = 1;
P = 1;

function [] = calculateBridgeValues(numNodes, load, eMod, area, lengths, angles, nodes) 
    numDofs = numNodes*2
    % allAngles = [90, 0, atand(5/3), 180 - atand(5/3)]
    % allLengths = 
    % nodeAngles are always between 0 and 180 degrees inclusive and are
    % from node 1 to node to 2 in nodes matrix
    K = zeros(numDofs)
    k = eMod.*area./lengths;
    mainLoopCount = length(lengths)
    for x = 1:mainLoopCount
        K = generateBeamStiffnessMatrix(nodes(x, 0), nodes(x,1), angles(x), numDofs, K, k(x))
    end
    for x = 1:mainLoopCount
        
    end 
end

function stiffnessMatrix = generateBeamStiffnessMatrix(node1, node2, angle, numDofs, globalK, stiffnessCoefficient)
    c = cosd(angle)
    s = sind(angle)
    cs = c*s
    if node1 < node2
        index1 = node1*2 - 1
        index2 = node1*2
        index3 = node2*2 - 1
        index4 = node2*2
    else 
        index1 = node2*2 - 1
        index2 = node2*2
        index3 = node1*2 - 1
        index4 = node1*2
    end
    indices = [index1 index2 index3 index4]
    specificK = zeros(numDofs)
    for x = 1:4
        for y = 1:4
            if (mod(x,2) == 1) &&(mod(y,2) == 1) && (x == y)
                specificK(indices(x), indices(y)) = c^2
            elseif (mod(x,2) == 1) &&(mod(y,2) == 1)
                specificK(indices(x), indices(y)) = -c^2
            elseif (mod(x,2) == 0) &&(mod(y,2) == 0) && (x == y)
                specificK(indices(x), indices(y)) = s^2
            elseif (mod(x,2) == 0) &&(mod(y,2) == 0)
                specificK(indices(x), indices(y)) = -s^2
            elseif (x == 1 && y == 2) || (x == 2 && y == 1) || (x == 3 && y == 4) || (x == 4 && y == 3)
                specificK(indices(x), indices(y)) = cs
            else 
                specificK(indices(x), indices(y)) = -cs
            end
        end
    end
    stiffnessMatrix = globalK + stiffnessCoefficient*specificK
end

function force = determineForces(node1, node2, angle, stiffnessCoefficient, 


             
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