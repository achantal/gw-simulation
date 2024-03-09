# gw-simulation
Simulation of gray whale population abundance - PCFG &amp; ENP

# % Create the original matrix
    p1 = 0.56; % survivorship of calves
    p2 = 0.89; % survivorship of juveniles (estimate)
    p3 = 0.87; % survivorship of adults
    d1 = 1; % duration of calf stage
    d2 = 7; % duration of juvenile stage
    d3 = 60; % duration of adult stage (estimate)

    P1 = ((1-(p1^(d1-1)))/(1-(p1^d1)))*p1; % calculate row 1, col 1 in matrix
    P2 = ((1-(p2^(d2-1)))/(1-(p2^d2)))*p2; % calculate row 2, col 2 in matrix
    P3 = ((1-(p3^(d3-1)))/(1-(p3^d3)))*p3; % calculate row 3, col 3 in matrix

    G1 = ((p1^d1)*(1-p1))/(1-(p1^d1)); % calculate row 2, col 2 in matrix
    G2 = ((p2^d2)*(1-p2))/(1-(p2^d2)); % calculate row 3, col 2 in matrix

    F2 = 0; % f2 is 0 because juveniles do not contribute to fecundity
    F3 = 0.5; % f3 is 0.5 because adults calve every 2 years (avg) and study duration is 1 year

    A1_original = [P1 F2 F3; G1 P2 0; 0 G2 P3] %create original matrix based on above calculations


# % Define parameters (population vector)
    numYears = 50; % Number of years for simulation
    initialPopulations1 = [12468; 3971; 210]; % Initial population size for stages Adults (170), Juveniles (50), Calves (10) - estimates

# % Create population vectors [Adults; Subadults; Calves]
    populationVector_original1 = initialPopulations1;
    populationVector_perturb_adult1 = initialPopulations1;
    populationVector_perturb_subadult1 = initialPopulations1;

# % Create data arrays for tracking the population over time
    adultsData_original1 = zeros(numYears, 1); % Empty array for original matrix
    subadultsData_original1 = zeros(numYears, 1); % Empty array for original matrix
    calvesData_original1 = zeros(numYears, 1); % Empty array for original matrix

# % Simulate population growth
    for year = 1:numYears
        % Calculate population in the next year for all matrices
        populationVector_original1 = A1_original * populationVector_original1;
       
        % Pertub survival rates at specific years
        if year == 5
            p3 = 0.95 % decrease of survival rate of adults at year 10
            p2 = 0.92
            p1 = 0.60
        elseif year == 10
            p3 = 0.89 %increase of survival rate of adults at year 15
            p2 = 0.87
            p1 = 0.54
        elseif year == 15
            p3 = 0.93 % increase of survival rate of adults at year 20
            p2 = 0.89
            p1 = 0.58
        elseif year == 20
            p3 = 0.84 % increase of survival rate of adults at year 20
            p2 = 0.85
            p1 = 0.52
        end

        % Recalculate the adult survival rate (P3)
        P3 = ((1-(p3^(d3-1)))/(1-(p3^d3)))*p3;

        %Update the original matrix (A1) with the new survival rate
        A1_original(3,3) = P3;

        % Recalculate the juvenile survival rate (P2)
        P2 = ((1-(p2^(d2-1)))/(1-(p2^d2)))*p2;

        %Update the original matrix (A1) with the new survival rate
        A1_original(2,2) = P2;

        % Recalculate the calf survival rate (P1)
        P1 = ((1-(p1^(d1-1)))/(1-(p1^d1)))*p1;

        %Update the original matrix (A1) with the new survival rate
        A1_original(1,1) = P1
        
        % Store data for analysis and visualization
        adultsData_original1(year) = populationVector_original1(1);
        subadultsData_original1(year) = populationVector_original1(2);
        calvesData_original1(year) = populationVector_original1(3);

end


# %Plotting
    years = 1:numYears;
    figure;
    plot(years, adultsData_original1, "b", years, subadultsData_original1, "g", years, calvesData_original1, "r");
    xlabel('Years')
    ylabel('Population')
    legend('Original - Adults', 'Original - Juveniles', 'Original - Calves');
    title('PCFG Gray Whale Population Projection (3x3 Lefkovitch Matrix) with perturbed adults');

# %Find dominant eigenvalue
    [V, D] = eigs(A1_original, 1, 'lm')
    lambda = D(1,1) % Calculate the dominant eigenvalue (the eigenvalue that is greater in absolute value than all other eigenvalues)
