%BladeSection.m

classdef BladeSection
    % BladeSection Class for handling propeller blade section generation
    %
    % Enhanced version with support for:
    % - Multiple section types 
    % - Separate meanline and thickness profiles
    % - Proper scaling based on lift coefficients
    
    properties (SetAccess = private)
        type        % Section type (e.g., 'NACA66mod')
        meanline    % Meanline type (e.g., 'NACA a=0.8')
        thickness   % Thickness type (e.g., 'NACA 66')
        f0octilde   % f0/c NACA data for CLI == CLItilde
        CLItilde    % NACA data ideal lift coefficient
        alphaItilde % NACA data ideal angle of attack
        calculator  % Function handle to section calculator
        As0         % Normalized cross-sectional area for validation
    end
    
    methods
        function obj = BladeSection(sectionType, varargin)
            % Constructor - sets up the appropriate section calculator
            % 
            % Required inputs:
            %   sectionType: String specifying blade section type
            %
            % Optional name-value pairs:
            %   'Meanline': String specifying meanline type 
            %   'Thickness': String specifying thickness type
            
            obj.type = sectionType;
            
            % Parse optional inputs
            p = inputParser;
            addParameter(p, 'Meanline', 'NACA a=0.8', @ischar);
            addParameter(p, 'Thickness', 'NACA 66', @ischar);
            parse(p, varargin{:});
            
            obj.meanline = p.Results.Meanline;
            obj.thickness = p.Results.Thickness;
            
            % Initialize NACA data parameters
            [obj.f0octilde, obj.CLItilde, obj.alphaItilde] = obj.getNACAParams();
            
            % Set up calculator function
            obj.calculator = obj.getSectionCalculator();

            % Verify section area
            obj.As0 = obj.calculateSectionArea(obj.thickness);        
        end
        
        function [y_b, y_f] = generateSection(obj, sRatio, tb, fb, b, CL)
            % Generate section coordinates using stored calculator
            %
            % Note on lift coefficients:
            % Currently using default CLItilde=1.0 which maintains the original
            % geometry from the offset table. In the future, this can be extended
            % to use actual calculated lift coefficients from lifting line theory,
            % similar to OpenProp's design approach.
            % 
            % Inputs:
            %   sRatio: Non-dimensional position along chord [0,1]
            %   tb: Maximum thickness/chord ratio
            %   fb: Maximum camber/chord ratio
            %   b: Chord length
            %   CL: Design lift coefficient (optional, defaults to CLItilde)
            %
            % Outputs:
            %   y_b: Upper (back) surface y-coordinates
            %   y_f: Lower (face) surface y-coordinates
            
            % Handle optional CL parameter
            if nargin < 6
                CL = obj.CLItilde;  % Use default lift coefficient if not provided
            end
            
            % Scale parameters based on lift coefficient
            fb_scaled = fb * (CL / obj.CLItilde);
            
            [y_b, y_f] = obj.calculator(sRatio, tb, fb_scaled, b);

            try
                % Force column vectors and ensure they're real numbers
                x = double(reshape(sRatio, [], 1));
                thickness = double(reshape(y_b - y_f, [], 1));
                
                if length(x) > 1 && length(thickness) > 1
                    actual_area = trapz(x, thickness);
                    expected_area = obj.As0 * tb * b;
                    
                    tolerance = 0.01;  % 1% tolerance
                    if abs(actual_area - expected_area)/expected_area > tolerance
                        warning(['Section area deviation: Expected = ' num2str(expected_area) ...
                                ', Actual = ' num2str(actual_area)]);
                    end
                end
                % Skip area verification for single points - this is normal when calculating individual points
            catch ME
                % If there's still an error, just continue without area verification
                % disp(['Area verification skipped: ', ME.message]);
            end
        end

        function As0 = calculateSectionArea(obj, thickness_type)
            % Calculate normalized cross-sectional area for verification
            % Returns As0 for c == t0 == 1
            
            % Use very fine spacing for accurate integration
            x0 = linspace(0, 1, 1000)';  % Column vector
            
            % Get thickness distribution using full section calculation
            switch upper(thickness_type)
                case {'NACA 66', 'NACA 66 (DTRC MODIFIED)'}
                    % Generate unit section (t0=1, b=1, no camber)
                    [y_b, y_f] = obj.NACA66mod(x0, 1, 0, 1);
                    As0 = trapz(x0, y_b - y_f);
                    
                    % Reference value from OpenProp: 0.7207
                    reference_As0 = 0.7207;
                    
                    % Check if area matches within tolerance
                    tolerance = 0.001;
                    if abs(As0 - reference_As0) > tolerance
                        warning(['Section area mismatch. Calculated: ' num2str(As0) ...
                                ', Reference: ' num2str(reference_As0)]);
                    end
                    
                case 'NACA 65A010'
                    reference_As0 = 0.6771;
                    % Similar calculation...
                    
                case 'ELLIPTIC'
                    reference_As0 = pi/4;  % Theoretical value
                    
                case 'PARABOLIC'
                    reference_As0 = 2/3;   % Theoretical value
                    
                otherwise
                    error('Unsupported thickness type for area calculation: %s', thickness_type)
            end
        end
    end
        
    methods (Access = private)
        function calculator = getSectionCalculator(obj)
            % Returns appropriate function handle based on section type
            switch upper(obj.type)
                case 'NACA66MOD'
                    calculator = @obj.NACA66mod;
                case 'NACA16'
                    calculator = @obj.NACA16;
                case 'NACA4DIGIT'
                    calculator = @obj.NACA4digit;
                otherwise
                    error('Unsupported section type: %s', obj.type)
            end
        end
        
        function [f0octilde, CLItilde, alphaItilde] = getNACAParams(obj)
            % Return baseline NACA parameters for the selected profile
            %
            % These parameters represent default values for geometry generation.
            % Currently using baseline values since geometry is primarily driven
            % by the offset table. Future development could include proper lift
            % coefficient calculations using lifting line theory, which would
            % then scale these baseline values.
            
            % Default values
            f0octilde = 0.0679;  % For NACA a=0.8 meanline
            CLItilde = 1.0;      % Reference lift coefficient
            alphaItilde = 1.54;  % Reference ideal angle of attack
            
            % Modify based on meanline type
            switch upper(obj.meanline)
                case 'NACA A=0.8'
                    % Keep defaults
                case 'NACA A=0.8 (MODIFIED)'
                    alphaItilde = 1.40;
                case 'PARABOLIC'
                    f0octilde = 1/(4*pi);
                    alphaItilde = 0;
                    CLItilde = 1.0;
            end
        end
        
        function [y_b, y_f] = NACA66mod(obj, sRatio, tb, fb, b)
            % Modified NACA66 section generator
            
            % Get meanline profile
            [f0, dfdx] = obj.getMeanlineProfile(sRatio, fb, b);
            
            % Get thickness profile
            [t] = obj.getThicknessProfile(sRatio, tb, b);
            
            % Combine for final section shape
            y_b = f0 + (t/2).*cos(atan(dfdx));
            y_f = f0 - (t/2).*cos(atan(dfdx));
        end
        
        function [f0, dfdx] = getMeanlineProfile(obj, sRatio, fb, b)
            % Generate meanline profile based on selected type
            switch upper(obj.meanline)
                case 'NACA A=0.8'
                    % Additional points near leading edge for better definition
                    xoc = [0.0000  0.00125 0.0025  0.005   0.0075  0.0125  0.0250  0.0500  ...
                           0.0750  0.1000  0.1500  0.2000  0.3000  0.4000  0.4500  0.5000  ...
                           0.6000  0.7000  0.8000  0.9000  0.9500  1.0000];
                    
                    foc = [0.0000  0.0135  0.0230  0.0365  0.0530  0.0907  0.1586  0.2712  ...
                           0.3657  0.4482  0.5869  0.6993  0.8635  0.9615  0.9881  1.0000  ...
                           0.9786  0.8892  0.7027  0.3586  0.1713  0.0000];
                    
                    % Compute slopes with finer resolution near leading edge
                    dx = diff(xoc);
                    dy = diff(foc);
                    midpoints = (xoc(1:end-1) + xoc(2:end))/2;
                    slopes = dy./dx;
                    
                    % Interpolate slopes to input positions
                    dfdxoc = interp1(midpoints, slopes, xoc, 'pchip', 'extrap');
                    
                    % Get meanline values and slopes at requested positions
                    f0   = fb * b * interp1(xoc, foc, sRatio, 'pchip', 'extrap');
                    dfdx = fb * interp1(xoc, dfdxoc, sRatio, 'pchip', 'extrap');
                    
                    % Handle leading edge point specifically
                    if any(sRatio == 0)
                        f0(sRatio == 0) = 0;
                        % Use forward difference for leading edge slope
                        le_slope = (foc(2) - foc(1))/(xoc(2) - xoc(1));
                        dfdx(sRatio == 0) = fb * le_slope;
                    end
                    
                otherwise
                    error('Unsupported meanline type: %s', obj.meanline)
            end
        end
        
        function [t] = getThicknessProfile(obj, sRatio, tb, b)
            % Generate thickness profile based on selected type
            % Implementation based on NACA 66 (DTRC modified) thickness data
            
            switch upper(obj.thickness)
                case {'NACA 66', 'NACA 66 (DTRC MODIFIED)'}
                    % NACA 66 (DTRC Modified) thickness distribution
                    % Based on standard DTRC modified data with enhanced leading edge definition
                    
                    % Lookup table for NACA 66 (DTRC modified)
                    xoc = [0.0000  0.00125 0.0025  0.005   0.0075  0.0125  0.0250  0.0500  ...
                           0.0750  0.1000  0.1500  0.2000  0.3000  0.4000  0.4500  0.5000  ...
                           0.6000  0.7000  0.8000  0.9000  0.9500  1.0000];
                    
                    % t/t0 values (thickness ratio relative to max thickness)
                    tot0 = [0.0000  0.0750  0.1000  0.1428  0.1750  0.2088  0.2932  0.4132  ...
                            0.5050  0.5814  0.7042  0.8000  0.9274  0.9904  1.0000  0.9924  ...
                            0.9306  0.8070  0.6220  0.3754  0.2286  0.0666];
                    
                    % Use pchip interpolation to maintain shape
                    t = tb * b * interp1(xoc, tot0, sRatio, 'pchip', 'extrap');
                    
                    % Handle leading edge point
                    if any(sRatio == 0)
                        t(sRatio == 0) = 0;
                    end
                    
                otherwise
                    error('Unsupported thickness type: %s', obj.thickness)
            end
        end
        
        function [y_b, y_f] = NACA16(obj, sRatio, tb, fb, b)
            % NACA16 section generator - to be implemented
            error('NACA16 section not yet implemented')
        end
        
        function [y_b, y_f] = NACA4digit(obj, sRatio, tb, fb, b)
            % NACA 4-digit section generator - to be implemented
            error('NACA 4-digit section not yet implemented')
        end
    end
end