% FEM Sample Truss

% Clear Command Window and Workspace
clear all, clc;


% Constants for the beam (A is area, E is elastic modulus, L is beam lengths, theta)
A = 1; 
E = 1;
P = 1;

L = [12; 6*sqrt(3); 6]; %L array is hand calculated befrore simulation.
theata = [0; 150; 60]; 
nodes = [1 2; 2 3; 1 3]
numNodes = 3
calculateBridgeValues(numNodes, P, E, A, L, theata, nodes, 3, [1 2 3])

function [] = calculateBridgeValues(numNodes, load, eMod, area, lengths, angles, nodes, pIndex, ignoredIndices) 
    numDofs = numNodes*2
    % allAngles = [90, 0, atand(5/3), 180 - atand(5/3)]
    % allLengths = 
    % nodeAngles are always between 0 and 180 degrees inclusive and are
    % from node 1 to node to 2 in nodes matrix
    K = zeros(numDofs)
    k = eMod.*area./lengths;
    mainLoopCount = length(lengths)
    for x = 1:mainLoopCount
        K = generateBeamStiffnessMatrix(nodes(x, 1), nodes(x,2), angles(x), numDofs, K, k(x))
    end
    
    externalForces = zeros(length(numDofs)-3,1)
    externalForces(pIndex, 1) = load
    
    % Solve for unknown displacements (u2x, u3x, u3y)
    x = numDofs - length(ignoredIndices)
    kArrayPreSolve = zeros(x)
    for x = 1:numDofs - length(ignoredIndices)
        kArrayPreSolve(x)= K(ignoredIndices(end)+x,ignoredIndices(end)+1:end)
    end
    u = kArrayPreSolve \  externalForces
    
    U = zeros(numDofs)
    for x = 1:numDofs-3
        U(x) = u(3+x)
    end

    reactionForces = K * U;
    memberForces = []  
    for x = 1:mainLoopCount
        f = determineForce(nodes(x, 0), nodes(x,1), angles(x), k(x), U)
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

function force = determineForce(node1, node2, angle, stiffnessCoefficient, globalDisplacementMatrix) 
    x1 = node1*2 - 1
    y1 = node1*2
    x2 = node2*2 - 1
    y2 = node2*2
    force = stiffnessCoefficient*((globalDisplacementMatrix(x2)-globalDisplacementMatrix(x1))*cosd(angle) + (globalDisplacementMatrix(y2)-globalDisplacementMatrix(y1))*sind(angle))
end


