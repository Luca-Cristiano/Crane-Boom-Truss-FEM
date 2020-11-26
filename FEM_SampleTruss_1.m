% FEM Sample Truss

% Clear Command Window and Workspace
clear all, clc;


% Constants for the beam (A is area, E is elastic modulus, L is beam lengths, theta)
A = 1; 
E = 1;
P = 1;

%CHANGE THESE VALUES
%length array for every member beam
%L = [6; 6*sqrt(3); 12];
L = [5;3;sqrt(34);6;sqrt(34);6;sqrt(34);6;sqrt(34);6;sqrt(34);6;sqrt(34);6;sqrt(34);6;sqrt(34);6;sqrt(34);6;sqrt(34)]
%theta array corresponds to lengths, note it is from node 1 to node 2
%theta = [60; 150; 0]; 
theata = [90;0;atan(5/3);0;atan(-5/3);0;atan(5/3);0;atan(-5/3);0;atan(5/3);0;atan(-5/3);0;atan(5/3);0;atan(-5/3);0;atan(5/3);0;atan(-5/3)]
%node array corresponds to beam, angles are from node #1 to node #2
%nodes = [1 3; 2 3; 1 2]
nodes = [1 2; 1 3; 2 3; 2 4; 3 4; 3 5; 4 5; 4 6;5 6; 5 7; 6 7; 6 8; 7 8; 7 9; 8 9; 8 10; 9 10; 9 11; 10 11; 10 12; 11 12]
%number of total nodes in the system
%numNodes = 3
numNodes = 12
%index where external load P is exerted, note note this index is determined
%by using node*2-1 if the force is in the x direction or by using node*2
%if the force is in the direction
%forceIndex = 5
forceIndex = 24 
%these indexes are using the nodal indexes in terms of directions, for the
%x direction index is equal to node*2 - 1, for y direction index is equal
%to node*2
%zeroDisplacementIndices = [1 2 4]
zeroDisplacementIndices = [1 2 3]
%CHANGE THESE VALUES

calculateBridgeValues(numNodes, P, E, A, L, theata, nodes, forceIndex, zeroDisplacementIndices)

function [] = calculateBridgeValues(numNodes, load, eMod, area, lengths, angles, nodes, pIndex, ignoredIndices) 
    numDofs = numNodes*2
    K = zeros(numDofs)
    k = eMod.*area./lengths;
    mainLoopCount = length(lengths)
    for x = 1:mainLoopCount
        x
        K = generateBeamStiffnessMatrix(nodes(x, 1), nodes(x,2), angles(x), numDofs, K, k(x))
    end
    
    %TESTED UP TO THIS POINT
    externalForces = []
    indexCounter = 0
    for x = 1:numDofs
        if ismember(x, ignoredIndices)
        elseif x == pIndex
            externalForces(end+1) = load
        else
            externalForces(end+1) = 0
        end
    end
    externalForces = externalForces'
    
    kArrayPreSolve = K
    for x = numDofs:-1:1
        if ismember(x, ignoredIndices)
            kArrayPreSolve(x, :) = []
        end
    end
    
    for x = numDofs:-1:1
        if ismember(x, ignoredIndices)
            kArrayPreSolve(:, x) = []
        end
    end
    
    u = kArrayPreSolve \ externalForces
    U = []
    uIndex = 1
    for x = 1:numDofs
        if ismember(x, ignoredIndices)
            U(end+1)= 0
        else
            U(end+1)= u(uIndex)
            uIndex = uIndex + 1
        end   
    end
    
    memberForces = []    
    for x = 1:mainLoopCount
        f = determineForce(nodes(x, 1), nodes(x,2), angles(x), k(x), U)
        memberForces = [memberForces; f]
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
        error('Node indexing is incorrect')
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

function force = determineForce(node1, node2, angle, stiffnessCoefficient, globalDisplacementMatrix)
    if node1 < node2
        x1 = node1*2 - 1
        y1 = node1*2
        x2 = node2*2 - 1
        y2 = node2*2
    else 
        error('Node indexing is incorrect')
    end
    force = stiffnessCoefficient*((globalDisplacementMatrix(x2)-globalDisplacementMatrix(x1))*cosd(angle) + (globalDisplacementMatrix(y2)-globalDisplacementMatrix(y1))*sind(angle))
end


