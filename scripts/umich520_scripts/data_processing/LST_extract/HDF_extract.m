% Specify the path to your HDF file
hdfFilePath = 'C:/cygwin64/home/Caleb/MOD11A1.A2023250.h08v05.061.2023252003330.hdf';

% Use the hdfinfo function to get information about the HDF file
fileInfo = hdfinfo(hdfFilePath);

% Display the structure of the HDF file
disp(fileInfo);

% Read data from the HDF file using hdfread function
% Replace 'DatasetName' with the name of the dataset you want to read
data = hdfread(hdfFilePath, 'MYD21A1D');

% Display the size of the data array
disp(size(data));

% Optionally, you can visualize or perform further processing on the data


% Open the HDF file for reading
fileID = H5F.open(hdfFilePath, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

% Get the root group identifier
rootGroupID = H5G.open(fileID, '/');

% Get the number of objects (datasets, groups, etc.) in the root group
numObjects = H5G.get_num_objs(rootGroupID);

% Loop through each object in the root group
for i = 0:numObjects-1
    % Get the name of the i-th object
    objectName = H5G.get_objname_by_idx(rootGroupID, i);
    disp(['Object Name: ', objectName]);
end

% Close the root group and file
H5G.close(rootGroupID);
H5F.close(fileID);
