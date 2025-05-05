function exportRhino(filename, allBlades, Params)
    fid = fopen(filename, 'wt');
    
    [numPoints, ~] = size(allBlades.x);
    pointsPerBlade = numPoints / Params.Z;
    nChord = Params.nChord;
    nSpan = Params.nSpan;
    
    % For first blade only
    x = allBlades.x(1:pointsPerBlade);
    y = allBlades.y(1:pointsPerBlade);
    z = allBlades.z(1:pointsPerBlade);
    
    % Separate back and face
    xBack = x(1:2:end);
    yBack = y(1:2:end);
    zBack = z(1:2:end);
    
    xFace = x(2:2:end);
    yFace = y(2:2:end);
    zFace = z(2:2:end);
    
    % Reshape into grids
    xBackGrid = reshape(xBack, [nChord, nSpan]);
    yBackGrid = reshape(yBack, [nChord, nSpan]);
    zBackGrid = reshape(zBack, [nChord, nSpan]);
    
    xFaceGrid = reshape(xFace, [nChord, nSpan]);
    yFaceGrid = reshape(yFace, [nChord, nSpan]);
    zFaceGrid = reshape(zFace, [nChord, nSpan]);
    
    % Create debug file to inspect point positions
    debugFile = [filename '_debug.txt'];
    debugFid = fopen(debugFile, 'wt');
    
    % Export back surface splines with debug information
    for i = 1:nSpan
        % Write section header to debug file
        fprintf(debugFid, '\nSpan Section %d:\n', i);
        
        fprintf(fid, '_Curve\n');
        fprintf(fid, '_Degree=3\n');
        
        % Write points with extra checks for leading edge region
        for j = 1:nChord
            % Calculate point-to-point distances for debugging
            if j > 1
                dx = xBackGrid(j,i) - xBackGrid(j-1,i);
                dy = yBackGrid(j,i) - yBackGrid(j-1,i);
                dz = zBackGrid(j,i) - zBackGrid(j-1,i);
                dist = sqrt(dx^2 + dy^2 + dz^2);
                
                % Log suspicious points (large changes in distance)
                if j > 2
                    prev_dx = xBackGrid(j-1,i) - xBackGrid(j-2,i);
                    prev_dy = yBackGrid(j-1,i) - yBackGrid(j-2,i);
                    prev_dz = zBackGrid(j-1,i) - zBackGrid(j-2,i);
                    prev_dist = sqrt(prev_dx^2 + prev_dy^2 + prev_dz^2);
                    
                    if dist > prev_dist * 3 || dist < prev_dist / 3
                        fprintf(debugFid, 'Suspicious point at chord=%d, dist=%f (prev=%f)\n', j, dist, prev_dist);
                    end
                end
            end
            
            % Write point coordinates to both files
            fprintf(fid, '%.6f,%.6f,%.6f\n', xBackGrid(j,i), yBackGrid(j,i), zBackGrid(j,i));
            fprintf(debugFid, 'Point %d: (%.6f,%.6f,%.6f)\n', j, xBackGrid(j,i), yBackGrid(j,i), zBackGrid(j,i));
        end
        fprintf(fid, 'Enter\n');
    end
    
    % Export face surface splines
    for i = 1:nSpan
        fprintf(fid, '_Curve\n');
        fprintf(fid, '_Degree=3\n');
        
        for j = 1:nChord
            fprintf(fid, '%.6f,%.6f,%.6f\n', xFaceGrid(j,i), yFaceGrid(j,i), zFaceGrid(j,i));
        end
        fprintf(fid, 'Enter\n');
    end
    
    fclose(fid);
    fclose(debugFid);
end