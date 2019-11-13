function export_data2ros(insertion_depth, coil_currents, fileID)

    Ix = coil_currents(1,:);
    Iy = coil_currents(2,:);
    Iz = coil_currents(3,:);
  
    fprintf(fileID,'omnimagnet:\n');
    fprintf(fileID,'    trajectory:\n');
    
    fprintf(fileID,'        Ix: [');
    for ii = 1:length(Ix)

        if ii == length(Ix)
            fprintf( fileID, ' %2.4f', Ix(ii));
            fprintf( fileID,']\n');
        else
             fprintf( fileID, ' %2.4f,', Ix(ii));
        end
    end
    
    fprintf(fileID,'        Iy: [');
    for ii = 1:length(Iy)

        if ii == length(Iy)
            fprintf( fileID, ' %2.4f', Iy(ii));
            fprintf( fileID,']\n');
        else
             fprintf( fileID, ' %2.4f,', Iy(ii));
        end
    end
    
    fprintf(fileID,'        Iz: [');
    for ii = 1:length(Iz)

        if ii == length(Iz)
            fprintf( fileID, ' %2.4f', Iz(ii));
            fprintf( fileID,']\n');
        else
             fprintf( fileID, ' %2.4f,', Iz(ii));
        end
    end
    
    fprintf(fileID,'        insertion_depth: [');
    for ii = 1:length(insertion_depth)

        if ii == length(insertion_depth)
            fprintf( fileID, ' %2.4f', insertion_depth(ii));
            fprintf( fileID,']\n');
        else
             fprintf( fileID, ' %2.4f,', insertion_depth(ii));
        end
    end
    
    
    fclose(fileID);

end