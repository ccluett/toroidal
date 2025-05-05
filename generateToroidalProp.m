clear all; clc;  close all; 

%generateToroidalProp.m

% Create base parameters
Params = struct();

Params.D = 800;    % Diameter [mm] 800
Params.L = Params.D/5;    % Total axial span [mm] 160
Params.Z = 3;      % Number of blades

% Section generation parameters
Params.foilType = 'NACA66MOD';  % Available types: 'NACA66MOD', 'NACA16' (not implemented), 'NACA4DIGIT' (not implemented)

% Optional section parameters (if not specified, defaults will be used)
Params.meanlineType  = 'NACA A=0.8';      % Options: 'NACA A=0.8' (default), 'NACA A=0.8 (MODIFIED)', 'PARABOLIC'
Params.thicknessType = 'NACA 66';         % Options: 'NACA 66' (default), 'NACA 66 (MODIFIED)'

% Note: Currently CL values are not used as geometry is driven by offset table.
%       Future development will include proper lift coefficient calculations.

% Export parameters
Params.plotMesh = 1;             % Plot propeller mesh
Params.plotPoints = 1;          % Plot propeller points
Params.plotRollAngleDistribution = 0; % Plot Transition
Params.plotBladeSections = 0;    % Plot blade sections
Params.plotExpandedSections = 0; % Plot expanded blade sections
Params.export = 0;              % Enable export
Params.generatePointCloud = 0;
Params.exportRhino = 1;
Params.generateNURBS = 0;       % Enable NURBS generation
Params.exportIGES = 0;          % Enable IGES export
Params.exportPath = '/Users/ccluett/Desktop/propeller1';  % Optional path for exports
Params.exportName = 'exportedPropeller';  % Base name for exported files

% Optional parameters will default to:
Params.nSpan = 50;                  % Number of spanwise points
Params.nChord = 25;                 % Number of chordwise points per side
Params.sections.frontEnd = 0.4;     % End of front section (l/L)
Params.sections.rearStart = 0.8;    % Start of rear section (l/L)
Params.pitch.front = 1.0;           % No modification to front pitch
Params.pitch.transition = 1.0;      % No modification to transition pitch
Params.pitch.rear = 1.0;            % No modification to rear pitch
Params.pitch.overall = 1.0;         % No overall pitch scaling
Params.chord.scale = 1.0;           % No chord scaling
Params.thickness.scale = 1.0;       % No thickness scaling
Params.camber.scale = 1.0;          % No camber scaling
Params.skew.scale = 1;              % No skew scaling
Params.rake.scale = 1.0;              % No rake scaling
Params.hubRadius = .2;             % Hub radius ratio (Rhub/R = 0.2)

% Generate the propeller
generateToroidalPropeller(Params);

function generateToroidalPropeller(Params)
    % Check inputs
    if nargin < 1
        error('Parameters must be provided');
    end
    
    % Validate required parameters
    required = {'L', 'D', 'Z'};
    for i = 1:length(required)
        if ~isfield(Params, required{i})
            error('Missing required parameter: %s', required{i});
        end
    end
    
    % Set default parameters if not provided
    Params = setDefaultParams(Params);
    
    % Get dimensions
    L = Params.L;     % Total axial span [mm]
    D = Params.D;     % Diameter [mm]
    R = D/2;         % Radius [mm]
    Z = Params.Z;     % Number of blades
    
    % Load base propeller data and scale it according to parameters
    [baseOffsets, scaledOffsets] = loadAndScaleOffsets(Params);
    
    % Store offsets in Params for use by plotting functions
    Params.baseOffsets = baseOffsets;
    Params.scaledOffsets = scaledOffsets;
    
    % % Generate blade points
    % nSpan = Params.nSpan;    % Points along axial span
    % nChord = Params.nChord;  % Points along each side of the airfoil
    % 
    % % Initialize arrays for blade surface coordinates
    % nPoints = nSpan * nChord * 2;  % Doubled for back and face
    % xCoords = zeros(nPoints, 1);
    % yCoords = zeros(nPoints, 1);
    % zCoords = zeros(nPoints, 1);
    % 
    % Generate Blade Surface
    % point = 1;
    % for iSpan = 1:nSpan
    %     % Normalized axial position
    %     lRatio = (iSpan-1)/(nSpan-1);
    % 
    %     % Get geometric parameters at this axial position
    %     geomParams = interpolateParameters(lRatio, scaledOffsets, Params);
    % 
    %     % Generate points for both surfaces
    %     for iChord = 1:nChord
    %         sRatio = (iChord-1)/(nChord-1);
    % 
    %         % Calculate blade back coordinates
    %         [xb, yb, zb] = calculatePoint(lRatio, sRatio, geomParams, L, R, true);
    %         xCoords(point) = xb;
    %         yCoords(point) = yb;
    %         zCoords(point) = zb;
    %         point = point + 1;
    % 
    %         % Calculate blade face coordinates
    %         [xf, yf, zf] = calculatePoint(lRatio, sRatio, geomParams, L, R, false);
    %         xCoords(point) = xf;
    %         yCoords(point) = yf;
    %         zCoords(point) = zf;
    %         point = point + 1;
    %     end
    % end
    
    % Generate blade points using cosine spacing
    [xCoords, yCoords, zCoords] = generateBladePoints(Params, scaledOffsets, L, R);

    % Generate all blades
    allBlades = generateAllBlades(xCoords, yCoords, zCoords, Z);
    
    % Plot the propeller if requested
    if Params.plotPoints
        plotPropellerPoints(allBlades, Params);
    end
    
    if Params.plotMesh
        plotPropellerMesh(allBlades, Params);
    end
    
    if Params.plotRollAngleDistribution
        plotRollAngleDistribution(Params, baseOffsets);
    end
    
    if Params.plotBladeSections
        plotBladeSections(allBlades, Params);
    end
    
    if Params.plotExpandedSections
        plotExpandedSections(allBlades, Params);
    end
    
    % Export
    if isfield(Params, 'export') && Params.export
        % Export point cloud
        if isfield(Params, 'generatePointCloud') && Params.generatePointCloud
            exportPointCloudData(allBlades, Params, Params.exportName);
        end
        
        % Export to Rhino
        if isfield(Params, 'exportRhino') && Params.exportRhino
            rhinoFilename = fullfile(Params.exportPath, [Params.exportName '_rhino.txt']);
            exportRhino(rhinoFilename, allBlades, Params);
        end
        
        % Generate NURBS surfaces if requested
        if isfield(Params, 'generateNURBS') && Params.generateNURBS
            generateNURBSPropeller(allBlades, Params, Params.exportName);
            
            % Export to IGES if requested
            if Params.exportIGES
                exportIGES(allBlades, Params, Params.exportName);
            end
        end
    end
end

function [x, y, z] = calculatePoint(lRatio, sRatio, params, L, R, isBack)
    % Convert degrees to radians
    alpha = deg2rad(params.alpha);
    beta = atan(params.PD/(pi*params.rlR));
    psi = deg2rad(params.psi);
    phi = deg2rad(params.phi);
    theta_s = deg2rad(params.theta_s);

    % Dimensional values
    b = params.bD * 2*R;
    s = sRatio * b; % chordwise coordinate along chord line
    l = lRatio * L;
    r_l = params.rlR * R;

    % Get direct rake (x_l) from interpolation:
    x_l = params.xlD * 2 * R;  % xlD is x_l/D, so multiply by 2*R to get x_l
    x_T = x_l + r_l*theta_s*tan(beta);  % Total rake including skew-induced component

    % Get blade section coordinates using the section generator (using default CLItilde)
    [y_b, y_f] = params.sectionGenerator.generateSection(sRatio, params.tb, params.fb, b); 
    
    % Select y_val depending on back or face surface
    if isBack
        y_val = y_b;
    else
        y_val = y_f;
    end

    % Equation 18 and 20 (leading edge to generatrix formulation):
    % r = r_l - (-(b/2) + s)*sin(alpha) + y_val*sin(psi);
    % 
    % c1 = (b / 2) - ((r * theta_s) / (cos(beta) * cos(alpha)));
    % 
    % theta = phi + (1 / r) * (((-c1 + s)*cos(alpha))*cos(beta) ...
    %                  + (y_val*cos(psi))*sin(beta) );
    % 
    % x = l + x_l + ((-(c1) + s)*cos(alpha))*sin(beta) - (y_val*cos(psi))*cos(beta);

    
    % Equation 16 (direct chord-based formulation):
        % Need to use x_T, total rake
    x = l + x_T + ((-(b/2) + s)*cos(alpha))*sin(beta) - (y_val*cos(psi))*cos(beta);

    r = r_l - (-(b/2) + s)*sin(alpha) + y_val*sin(psi); % tip roll in or out '-' or '+'

    theta = phi + theta_s + (1 / r) * (((-(b/2) + s)*cos(alpha))*cos(beta) ...
                     + (y_val*cos(psi))*sin(beta) );

    % Convert to Cartesian (y,z) (Equation 17):
    y = r*cos(theta);
    z = r*sin(theta);
end

function params = interpolateParameters(lRatio, offsets, Params)
    % Interpolate geometric parameters at given axial position
    params = struct();
    
    % Get table positions
    lTable = offsets(:,1);
    
    % Interpolate each parameter
    params.rlR = interp1(lTable, offsets(:,2), lRatio, 'pchip');
    params.bD = interp1(lTable, offsets(:,3), lRatio, 'pchip');
    params.PD = interp1(lTable, offsets(:,4), lRatio, 'pchip');
    params.tb = interp1(lTable, offsets(:,5), lRatio, 'pchip');
    params.fb = interp1(lTable, offsets(:,6), lRatio, 'pchip');
    params.theta_s = interp1(lTable, offsets(:,7), lRatio, 'pchip');
    params.xlD = interp1(lTable, offsets(:,8), lRatio, 'pchip');
    params.phi = interp1(lTable, offsets(:,9), lRatio, 'pchip');
    params.psi = interp1(lTable, offsets(:,10), lRatio, 'pchip');
    params.alpha = interp1(lTable, offsets(:,11), lRatio, 'pchip');
    
    % Pass through the section generator
    params.sectionGenerator = Params.sectionGenerator;
end

function allBlades = generateAllBlades(x, y, z, Z)
    % Generate all blades by rotating the first blade
    allBlades = struct('x', [], 'y', [], 'z', []);
    
    for i = 1:Z
        angle = 2*pi*(i-1)/Z;
        rotMatrix = [1 0 0;
                    0 cos(angle) -sin(angle);
                    0 sin(angle) cos(angle)];
        
        % Rotate coordinates
        rotated = [x y z] * rotMatrix';
        
        % Store coordinates
        allBlades.x = [allBlades.x; rotated(:,1)];
        allBlades.y = [allBlades.y; rotated(:,2)];
        allBlades.z = [allBlades.z; rotated(:,3)];
    end
end

function Params = setDefaultParams(Params)
    % Discretization defaults
    if ~isfield(Params, 'nSpan')
        Params.nSpan = 50;
    end
    if ~isfield(Params, 'nChord')
        Params.nChord = 15;
    end
    
    % Set default section boundaries
    if ~isfield(Params, 'sections')
        Params.sections.frontEnd = 0.4;    % End of front section (l/L)
        Params.sections.rearStart = 0.8;   % Start of rear section (l/L)
    end
    

    if ~isfield(Params, 'hubRadius')
        Params.hubRadius = 0.2;
    end
    if ~isfield(Params, 'pitch')
        Params.pitch = struct();
    end
    if ~isfield(Params.pitch, 'front')
        Params.pitch.front = 1.0;
    end
    if ~isfield(Params.pitch, 'transition')
        Params.pitch.transition = 1.0;
    end
    if ~isfield(Params.pitch, 'rear')
        Params.pitch.rear = 1.0;
    end
    if ~isfield(Params.pitch, 'overall')
        Params.pitch.overall = 1.0;
    end
    if ~isfield(Params, 'export')
        Params.export = false;
    end
    if ~isfield(Params, 'generateNURBS')
        Params.generateNURBS = false;
    end
    if ~isfield(Params, 'exportPath')
        Params.exportPath = pwd;  % Current directory
    end
    if ~isfield(Params, 'exportName')
        Params.exportName = 'propeller';
    end
    if ~isfield(Params, 'exportIGES')
        Params.exportIGES = false;
    end
    if ~isfield(Params, 'plotBladeSections')
        Params.plotBladeSections = false; 
    end
    if ~isfield(Params, 'plotExpandedSections')
        Params.plotExpandedSections = false;
    end
    if ~isfield(Params, 'plotRollAngleDistribution')
        Params.plotRollAngleDistribution = false;
    end

    % Initialize other geometric parameters with defaults
    paramNames = {'chord', 'thickness', 'camber', 'skew', 'rake'};
    for i = 1:length(paramNames)
        name = paramNames{i};
        if ~isfield(Params, name)
            Params.(name) = struct();
        end
        if ~isfield(Params.(name), 'scale')
            Params.(name).scale = 1.0;
        end
        if ~isfield(Params.(name), 'distribution')
            Params.(name).distribution = 'linear';
        end
        if ~isfield(Params.(name), 'customFcn')
            Params.(name).customFcn = @(x) 1;
        end
    end
    % Set default foil type if not provided
    if ~isfield(Params, 'foilType')
        Params.foilType = 'NACA66mod';
    end
    
    % Create section generator
    Params.sectionGenerator = BladeSection(Params.foilType);
end

function [baseOffsets, scaledOffsets] = loadAndScaleOffsets(Params)
    % Load base propeller definition
    baseOffsets = getBaseOffsets();
    scaledOffsets = baseOffsets;
    
    % Get normalized positions
    lL = baseOffsets(:,1);  % l/L values
    
    % If a new hub radius ratio is specified, adjust the r_l/R values
    % while preserving the distribution shape
    if isfield(Params, 'hubRadius')
        if Params.hubRadius <= 0 || Params.hubRadius >= 1
            error('Hub radius ratio must be between 0 and 1');
        end
        originalHubRatio = baseOffsets(1,2);  % First row, second column is r_l/R
        r_l_R = baseOffsets(:,2);  % All r_l/R values
        
        % Calculate scaling that preserves distribution shape
        % Map values between hub and tip while maintaining ratios
        r_l_R_normalized = (r_l_R - originalHubRatio) / (1 - originalHubRatio);
        r_l_R_scaled = Params.hubRadius + (1 - Params.hubRadius) * r_l_R_normalized;
        
        scaledOffsets(:,2) = r_l_R_scaled;
    end
    
    % Define section ranges and scale other parameters as before
    frontSection = lL <= Params.sections.frontEnd;
    transitionSection = lL > Params.sections.frontEnd & lL < Params.sections.rearStart;
    rearSection = lL >= Params.sections.rearStart;
    
    % Scale pitch (P/D) - column 4
    P_D = baseOffsets(:,4);
    scaledP_D = P_D * Params.pitch.overall;
    scaledP_D(frontSection) = scaledP_D(frontSection) * Params.pitch.front;
    scaledP_D(transitionSection) = scaledP_D(transitionSection) * Params.pitch.transition;
    scaledP_D(rearSection) = scaledP_D(rearSection) * Params.pitch.rear;
    scaledOffsets(:,4) = scaledP_D;
    
    % Scale other parameters
    scaledOffsets(:,3) = scaleParameter(baseOffsets(:,3), lL, Params.chord);        % b/D
    scaledOffsets(:,5) = scaleParameter(baseOffsets(:,5), lL, Params.thickness);    % t/b
    scaledOffsets(:,6) = scaleParameter(baseOffsets(:,6), lL, Params.camber);       % f/b
    scaledOffsets(:,7) = scaleParameter(baseOffsets(:,7), lL, Params.skew);         % θs
    scaledOffsets(:,8) = scaleParameter(baseOffsets(:,8), lL, Params.rake);         % x1/D
end


function scaledValues = scaleParameter(baseValues, lL, scaleParams)
    % Apply scaling based on distribution type
    switch scaleParams.distribution
        case 'linear'
            scaledValues = baseValues * scaleParams.scale;
        case 'custom'
            scaledValues = baseValues .* scaleParams.scale .* scaleParams.customFcn(lL);
        otherwise
            error('Unknown distribution type: %s', scaleParams.distribution);
    end
end

function baseOffsets = getBaseOffsets()
    % Base propeller definition from Table 2
    % Format: [l/L r_l/R b/D P/D t/b f/b θs x_l/D φ ψ α]
    baseOffsets = [
        0       0.200 0.121 1.312 0.301 0.014  0.00 -0.062 -27.56   0.00  0.00;
        0.062   0.300 0.143 1.331 0.243 0.030 -0.89 -0.082 -25.52   2.44  0.12;
        0.120   0.400 0.160 1.336 0.197 0.044 -1.10 -0.093 -23.44   3.38  0.62;
        0.177   0.500 0.172 1.321 0.158 0.055 -0.90 -0.097 -21.26   4.26  1.44;
        0.234   0.600 0.180 1.281 0.125 0.063 -0.36 -0.093 -18.93   5.44  2.61;
        0.294   0.700 0.181 1.198 0.099 0.069  0.43 -0.081 -16.33   8.31  4.37;
        0.359   0.800 0.168 1.078 0.077 0.069  1.46 -0.059 -13.26  14.67  6.89;
        0.438   0.900 0.142 0.932 0.061 0.062  2.83 -0.026  -9.33  28.78  9.39;
        0.489   0.950 0.123 0.863 0.055 0.052  3.78 -0.002  -6.55  42.83 10.45;
        0.525   0.975 0.112 0.830 0.053 0.043  4.45  0.014  -4.56  55.45 10.98;
        0.573   0.995 0.103 0.807 0.052 0.030  5.34  0.035  -1.86  75.91 11.41;
        0.600   1.000 0.101 0.806 0.054 0.021  5.87  0.045  -0.24  90.00 11.51;
        0.644   0.995 0.107 0.853 0.058 0.008  6.69  0.061   2.45 118.06 11.28;
        0.684   0.975 0.119 0.951 0.065 -0.004 7.44  0.073   4.99 139.40 10.52;
        0.715   0.950 0.131 1.040 0.073 -0.011 7.99  0.079   6.95 151.75  8.97;
        0.757   0.900 0.151 1.183 0.087 -0.020 8.75  0.082   9.78 164.02  5.54;
        0.815   0.800 0.176 1.396 0.114 -0.027 9.72  0.072  13.77 173.49  2.33;
        0.858   0.700 0.187 1.543 0.141 -0.028 10.40 0.054  16.83 176.68  1.44;
        0.894   0.600 0.186 1.637 0.169 -0.026 10.94 0.030  19.44 177.66  0.98;
        0.925   0.500 0.177 1.669 0.197 -0.023 11.39 0.002  21.74 178.50  0.63;
        0.952   0.400 0.162 1.663 0.225 -0.019 11.77 -0.027 23.82 179.13  0.38;
        0.977   0.300 0.145 1.652 0.253 -0.015 12.11 -0.058 25.74 179.62  0.17;
        1.000   0.200 0.124 1.631 0.281 -0.010 12.41 -0.091 27.56 180.00  0.00
    ];
end

% Function to generate cosine-spaced points in range [a,b]
function points = getCosineSpacing(a, b, n)
    % Generates n points with cosine spacing between a and b
    % This creates denser spacing near both endpoints
    
    % Generate angles from 0 to pi
    theta = linspace(0, pi, n);
    
    % Use (1 - cos) to get clustering at both ends
    % This gives values from 0 to 2, so we divide by 2 to normalize
    normalized = (1 - cos(theta))/2;
    
    % Scale to desired range [a,b]
    points = a + (b-a)*normalized;
end

% Function to generate blade surface points with cosine spacing
function [xCoords, yCoords, zCoords] = generateBladePoints(Params, scaledOffsets, L, R)
    nSpan = Params.nSpan;    % Points along axial span
    nChord = Params.nChord;  % Points along each side of the airfoil
    
    % Initialize arrays for blade surface coordinates
    nPoints = nSpan * nChord * 2;  % Doubled for back and face
    xCoords = zeros(nPoints, 1);
    yCoords = zeros(nPoints, 1);
    zCoords = zeros(nPoints, 1);
    
    % Generate cosine-spaced points along span (l/L)
    spanPoints = getCosineSpacing(0, 1, nSpan);
    
    % Generate cosine-spaced points along chord
    chordPoints = getCosineSpacing(0, 1, nChord);
    
    point = 1;
    for iSpan = 1:nSpan
        % Use cosine-spaced position along span
        lRatio = spanPoints(iSpan);
        
        % Get geometric parameters at this axial position
        geomParams = interpolateParameters(lRatio, scaledOffsets, Params);
        
        % Generate points for both surfaces
        for iChord = 1:nChord
            % Use cosine-spaced position along chord
            sRatio = chordPoints(iChord);
            
            % Calculate blade back coordinates
            [xb, yb, zb] = calculatePoint(lRatio, sRatio, geomParams, L, R, true);
            xCoords(point) = xb;
            yCoords(point) = yb;
            zCoords(point) = zb;
            point = point + 1;
            
            % Calculate blade face coordinates
            [xf, yf, zf] = calculatePoint(lRatio, sRatio, geomParams, L, R, false);
            xCoords(point) = xf;
            yCoords(point) = yf;
            zCoords(point) = zf;
            point = point + 1;
        end
    end
end

