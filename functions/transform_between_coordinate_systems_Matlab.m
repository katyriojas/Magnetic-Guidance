function transformed_points = transform_between_coordinate_systems_Matlab(input_points,dim,vsize,input_orientation,output_orientation)

if length(vsize)~=length(dim)
    fprintf('ERROR: Vectors with the dimension of the image and the voxel size do not match.\n');
    transformed_points = [];
    return
end

wrong_orientation_change = 0;
for axis = 1:length(dim)
    if input_orientation(axis)=='R' && ~(output_orientation(axis)=='R' || output_orientation(axis)=='L')
        wrong_orientation_change = wrong_orientation_change+1;
    end
    if input_orientation(axis)=='L' && ~(output_orientation(axis)=='R' || output_orientation(axis)=='L')
        wrong_orientation_change = wrong_orientation_change+1;
    end
    if input_orientation(axis)=='P' && ~(output_orientation(axis)=='P' || output_orientation(axis)=='A')
        wrong_orientation_change = wrong_orientation_change+1;
    end
    if input_orientation(axis)=='A' && ~(output_orientation(axis)=='P' || output_orientation(axis)=='A')
        wrong_orientation_change = wrong_orientation_change+1;
    end
    if input_orientation(axis)=='S' && ~(output_orientation(axis)=='S' || output_orientation(axis)=='I')
        wrong_orientation_change = wrong_orientation_change+1;
    end
    if input_orientation(axis)=='I' && ~(output_orientation(axis)=='S' || output_orientation(axis)=='I')
        wrong_orientation_change = wrong_orientation_change+1;
    end
end
if wrong_orientation_change
    fprintf('ERROR: This program cannot transform between the provided orientations.\n');
    transformed_points = [];
    return
end
transformed_points = zeros(size(input_points));
for n_point = 1:size(input_points,2)
    for axis = 1:length(dim)
        if input_orientation(axis)~=output_orientation(axis)
            transformed_points(axis,n_point) = vsize(axis)*(dim(axis)+1)-input_points(axis,n_point);
        else
            transformed_points(axis,n_point) = input_points(axis,n_point);
        end
    end
end