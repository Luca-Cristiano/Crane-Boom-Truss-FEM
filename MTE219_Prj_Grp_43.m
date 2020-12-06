%Clear Command Window and Workspace
clear all, clc;

%lengths of every beam
L = [5;2;sqrt(29);6;sqrt(41);6;sqrt(29);6;sqrt(41);6;sqrt(29);6;sqrt(41);6;sqrt(29);6;sqrt(41);6;sqrt(29);6;sqrt(41)]
%angles for every corresponding beam, angles are from node i to node j
theata = [90;0;atand(5/2);0;atand(-5/4);0;atand(5/2);0;atand(-5/4);0;atand(5/2);0;atand(-5/4);0;atand(5/2);0;atand(-5/4);0;atand(5/2);0;atand(-5/3)]
%nodes corresponding to every beam, angles are based on these nodes
nodes = [1 2; 1 3; 2 3; 2 4; 3 4; 3 5; 4 5; 4 6;5 6; 5 7; 6 7; 6 8; 7 8; 7 9; 8 9; 8 10; 9 10; 9 11; 10 11; 10 12; 11 12]
%number of total nodes in the system
numNodes = 12
%node index direction where cantilever load is exerted
forceIndex = 24 
%indices of node directions where there is zero displacement
%corresponds to fixed and roller supports
zeroDisplacementIndices = [1 2 3]

%Constants related to cantilever loading
%P is the current cantilever load (kg)
%W is the force caused by the load
P = 1
g = -9.8
W = g*P
%loop to calculate bridge values and increment load until bridge fails
while 1
    calculateBridgeValues(numNodes, W, L, theata, nodes, forceIndex, zeroDisplacementIndices)
    P = P + 0.1
    W = g*P
end
    
function [] = calculateBridgeValues(numNodes, load, lengths, angles, nodes, pIndex, ignoredIndices) 
    %Geometry values for failure calculations
     t = 0.0015875
     d = 0.0047625
     b = 0.011484
     w = 0.01449
     l = 0
  
    area = t * w
     
    %Normal and Shear Strength and Elastic Modulus (Material Properties)
    bassNormalSTR = 65141501.9
    bassShearSTR = 4200000
    hardNormalSTR = 1139000000
    hardShearSTR = 54000000
    eMod = 10391670000
     
    %Initialization of highest forces and stress variables
    tHigh = 0
    cHigh = 0
    fHgih = 0
    tRupture = 0
    cRupture = 0
    bearing = 0
    tearout = 0
    buckling = 0
    shear = 0
    
    %initialize global K matrix
    numDofs = numNodes*2
    K = zeros(numDofs)
  
    %spring k values for every beam
    k = eMod.*area./lengths;
    
    %loop to generate global K matrix for truss
    for x = 1:length(lengths)
        K = generateBeamStiffnessMatrix(nodes(x, 1), nodes(x,2), angles(x), numDofs, K, k(x))
    end
    
    %external forces vector generated to solve for node displacements later
    %vector doesn't include rows where the displacement is 0
    externalForces = []
    for x = 1:numDofs
        if ismember(x, ignoredIndices)
        elseif x == pIndex
            externalForces(end+1) = load
        else
            externalForces(end+1) = 0
        end
    end
    externalForces = externalForces'
    
    %removes columns that correspond to zero displacements from K matrix
    kArrayPreSolve = K
    for x = numDofs:-1:1
        if ismember(x, ignoredIndices)
            kArrayPreSolve(x, :) = []
        end
    end
    
    %removes rows that correspond to zero displacements from K matrix
    for x = numDofs:-1:1
        if ismember(x, ignoredIndices)
            kArrayPreSolve(:, x) = []
        end
    end
    
    %solves for displacements of each node and creates global displacement 
    %matrix that includes zero displacements
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
    
    %loop to determine the forces of each memeber using displacements,
    %stiffness coefficients and node numbers
    memberForces = []    
    for x = 1:length(lengths)
        f = determineForce(nodes(x, 1), nodes(x,2), angles(x), k(x), U)
        memberForces = [memberForces; f]
    end 
    memberForces (1) = memberForces(1)*-1
    
    %determines resultant force at each node of truss
    nodeResultantX  = zeros(numNodes, 1)
    nodeResultantY = zeros(numNodes, 1)
    for x = 1:length(nodes)
        node1 = nodes(x, 1)
        node2 = nodes(x, 2)
        xComponent = abs(memberForces(x)*cosd(angles(x)))
        yComponent = abs(memberForces(x)*sind(angles(x)))
        nodeResultantX(node1, 1) = nodeResultantX(node1, 1) + xComponent
        nodeResultantX(node2, 1) = nodeResultantX(node2, 1) + xComponent
        nodeResultantY(node1, 1) = nodeResultantY(node1, 1) + yComponent
        nodeResultantY(node2, 1) = nodeResultantY(node2, 1) + yComponent
    end
    nodeResultant = sqrt(nodeResultantX.^2 + nodeResultantY.^2)./2
    maxShearForce = max(nodeResultant)
    
    %Shear stress calculation - max force on pin over cross-sectional area
    shear = maxShearForce/(pi*(d/2)^2)
    
    %Error message displayed when shear force surpasses the shear strength
    %of the hardwood
    if shear > hardShearSTR
        printFinalValues(memberForces, U, load)
        error('Truss failed due to shear')
    end
    
    %Setting up the highest tension and compression forces for future
    %calculations
    tHigh = max(memberForces)
    cHigh = min(memberForces)
    
    %Checking to make sure there are memebers in compression and tension
    if tHigh < 0
        tHigh = 0
        printFinalValues(memberForces, U, load)
        error('No members in tension')
    end
    
    if cHigh > 0
        cHigh = 0
        printFinalValues(memberForces, U, load)
        error('No members in compression')
    end
    
    %Obtaining the magnitude of the compression force
    cHigh = abs(cHigh)
    
    %Calculating the max rupture stress from compression and tension
    tRupture = tHigh/(t*(w-d))
    cRupture = cHigh/(t*w)
    
    %Error message appears if the rupture stress is larger than the normal
    %stress of the basswood
    if cRupture > bassNormalSTR || tRupture > bassNormalSTR
        printFinalValues(memberForces, U, load)
        error ('Truss failed due to rupture')
    end
    
    %Creating a variable for the highest force
    if cHigh >= tHigh
        fHigh = cHigh
    else
        fHigh = tHigh
    end
    
    %Bearing Failure highest magnitude of the force over thickness *
    %diameter for the link and the pin
    bearing = fHigh / (t*d)
    
    %Checking the calculated bearing stress with the basswood normal
    %strength
    if bearing > bassNormalSTR
        printFinalValues(memberForces, U, load)
        error ('Truss failed due to Bearing Failure (Link)')
    end
    
    %Because the holes in the links are the same size as the pins and the
    %same thickness is in contact with the pin and link the pervious
    %bearing stress is compared to the normal strength of the hardwood
    if bearing > hardNormalSTR
        printFinalValues(memberForces, U, load)
        error ('Truss failed due to Bearing Failure (Pin)')
    end
    
    %Finding the tearout for the basswood link which is the highest tension
    %in the system over 2*thickness*the center of the hole to the edge of
    %the link
    tearout = tHigh/(2*b*t)
   
    %Compares the calculated tearout stress and compares it to the shear
    %strength of the Basswood
    if tearout > bassShearSTR
        printFinalValues(memberForces, U, load)
        error ('Truss failed due to tearout')
    end
    
    %Initializes a variable for the amount of compression forces in the
    %system
    compNum = 0
    
    %Calculates the amount of members in compression
    for x = 1: length(memberForces)
        if memberForces(x) < 0
            compNum = compNum + 1  
        end
    end
    
    %Creates a 2 row matrix where the first row will contain the lengths
    %of the members while the second row contains the compression forces
    compLAndF = zeros(2,compNum)
    compCount = 1
    
    %Populates the matrix
    for x = 1: length(memberForces)
        if memberForces(x) < 0
            compLAndF(1, compCount) = lengths (x)
            compLAndF(2, compCount) = abs(memberForces(x))
            compCount = compCount + 1;
        end
    end
    
    %Calculates the second moment of interia
    I = (t^3*w)/12
    
    %Calculates the buckling force for each compression member, pi^2*E*I/l^2. If the
    %buckling force is smaller than the compression force on a member the
    %beam buckles.
    for x = 1: compNum
        l = compLAndF(1, x)/100
        
        buckling = (pi^2 * eMod * I)/(l^2)
        %compLAndF(2, x) = compLAndF(2, x)
        if compLAndF(2, x) > buckling
            printFinalValues(memberForces, U, load)
            error ('Truss failed due to buckling')
        end
    end
end

%function that generates individual stiffness matrix for a beam and adds it
%to the current global stiffness matrix
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

%function used to determine the force of each member
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

%function that outputs final values when the bridge fails
function values = printFinalValues(memberForces, U, load)
    fprintf('Member Forces:\n');
    for i = 1:length(memberForces)
        fprintf('f(%i) = %.3f N\n', i, memberForces(i))
    end
    
    fprintf('\nNode Displacements\n');
    for i = 1:length(U)
        fprintf('u(%i) = %.3f m\n', i, U(i))
    end
    
    fprintf('\nFinal Cantilever Load Supported is %f kg\n\n', -load/9.8 - .1)    
end


