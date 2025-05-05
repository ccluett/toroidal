function plotPropellerPoints(blades, Params)
    % Create figure with white background
    fig = figure('Color', 'white', ...
           'Name', 'Toroidal Propeller Visualization', ...
           'NumberTitle', 'off', ...
           'Units', 'normalized', ...
           'Position', [0.1 0.1 0.8 0.8]);
    
    % Disable OpenGL clipping
    opengl('hardware');
    set(gcf, 'renderer', 'opengl');
    
    % Create main axes
    ax = axes('Parent', fig);
    hold(ax, 'on');

    % Calculate points per blade
    pointsPerBlade = size(blades.x, 1)/Params.Z;
    nChord = Params.nChord;
    nSpan = Params.nSpan;
    
    % Create handles to store plot objects
    h1 = gobjects(Params.Z, nSpan);
    h2 = gobjects(Params.Z, nSpan);
    
    % For each blade
    for i = 1:Params.Z
        % Extract single blade data
        startIdx = round(1 + (i-1)*pointsPerBlade);
        endIdx = round(i*pointsPerBlade);
        x = blades.x(startIdx:endIdx);
        y = blades.y(startIdx:endIdx);
        z = blades.z(startIdx:endIdx);
        
        % Separate back and face surfaces
        xBack = x(1:2:end);
        yBack = y(1:2:end);
        zBack = z(1:2:end);
        
        xFace = x(2:2:end);
        yFace = y(2:2:end);
        zFace = z(2:2:end);
        
        % Reshape into grid format
        xBackGrid = reshape(xBack, [nChord, nSpan]);
        yBackGrid = reshape(yBack, [nChord, nSpan]);
        zBackGrid = reshape(zBack, [nChord, nSpan]);
        
        xFaceGrid = reshape(xFace, [nChord, nSpan]);
        yFaceGrid = reshape(yFace, [nChord, nSpan]);
        zFaceGrid = reshape(zFace, [nChord, nSpan]);
        
        % Plot back surface with single scatter call
        backColor = [0.1 0.3 0.6];
        faceColor = [0.8 0.2 0.1];
        
        % Calculate transition region markers
        spanPositions = linspace(0, 1, nSpan);
        transitionStart = 0.5;
        transitionEnd = 0.7;
        
        % Create color arrays for all points
        colors = zeros(nChord * nSpan, 3);
        sizes = ones(nChord * nSpan, 1) * 25;
        
        for j = 1:nSpan
            normPos = spanPositions(j);
            if normPos < transitionStart
                currentColor = backColor;
            elseif normPos > transitionEnd
                currentColor = faceColor;
            else
                t = (normPos - transitionStart)/(transitionEnd - transitionStart);
                currentColor = backColor * (1-t) + faceColor * t;
            end
            
            % Set colors for all points in this span
            idx = ((j-1)*nChord + 1):(j*nChord);
            colors(idx, :) = repmat(currentColor, nChord, 1);
            
            % Set smaller markers in transition region
            if normPos >= transitionStart && normPos <= transitionEnd
                sizes(idx) = 15;
            end
        end
        
        % Plot all points with single scatter calls
        h1(i,:) = scatter3(xBackGrid(:), yBackGrid(:), zBackGrid(:), ...
            sizes, colors, 'filled', 'MarkerFaceAlpha', 0.7);
        h2(i,:) = scatter3(xFaceGrid(:), yFaceGrid(:), zFaceGrid(:), ...
            sizes, colors, 'filled', 'MarkerFaceAlpha', 0.7);
    end
    
    % Add rotation axis
    axisLength = max(blades.x) - min(blades.x);
    line([min(blades.x)-0.1*axisLength max(blades.x)+0.1*axisLength], [0 0], [0 0], ...
        'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
    
    % Set plot properties
    axis equal;
    xlabel('X', 'FontWeight', 'bold');
    ylabel('Y', 'FontWeight', 'bold');
    zlabel('Z', 'FontWeight', 'bold');
    title('Toroidal Propeller', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % Initial view
    view(45, 30);
    
    % Add buttons
    btnWidth = 80;
    btnHeight = 30;
    spacing = 10;
    
    % Front view button (x-axis alignment)
    uicontrol('Style', 'pushbutton', 'String', 'Front View', ...
        'Position', [20 20 btnWidth btnHeight], ...
        'Callback', @(src,event) view(ax, [90 0]));
    
    % Side view button (y-axis alignment)
    uicontrol('Style', 'pushbutton', 'String', 'Side View', ...
        'Position', [20+btnWidth+spacing 20 btnWidth btnHeight], ...
        'Callback', @(src,event) view(ax, [0 0]));
    
    % Top view button (z-axis alignment)
    uicontrol('Style', 'pushbutton', 'String', 'Top View', ...
        'Position', [20+2*(btnWidth+spacing) 20 btnWidth btnHeight], ...
        'Callback', @(src,event) view(ax, [0 90]));
    
    % Toggle zoom mode button (moved to the right)
    zoomState = false;
    hZoomBtn = uicontrol('Style', 'togglebutton', 'String', 'Zoom All', ...
        'Position', [20+3*(btnWidth+spacing) 20 btnWidth btnHeight], ...
        'Callback', @toggleZoomMode);
    
    % Enable 3D rotation
    rotate3d on;
    
    % Function to toggle between zoom modes
    function toggleZoomMode(src, ~)
        zoomState = get(src, 'Value');
        if zoomState
            % Show all points during zoom
            set(src, 'String', 'Zoom Normal');
            set(ax, 'Clipping', 'off');
            xlim('manual');
            ylim('manual');
            zlim('manual');
        else
            % Normal zoom behavior
            set(src, 'String', 'Zoom All');
            set(ax, 'Clipping', 'on');
            xlim('auto');
            ylim('auto');
            zlim('auto');
        end
    end
end

