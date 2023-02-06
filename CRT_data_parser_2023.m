%=========================================================================%
%                       CRT_data_parser.m
%
%   The m-file parses *.dat files collected using the program "crtunnel"
%   in the University of Minnesota closed return wind tunnel.  The data
%   files should be placed in a sub-directory of the directory where this
%   m-file resides.  The directory will be selected by the user once the 
%   m-file is ran. When run, this m-file will generate a matrix where each 
%   row corresponds to the particular data run (angle of attack and 
%   elevator setting).  The matrix has 18 columns where:
%
%   Column  1:      Row number
%   Column  2:      Angle of attack (deg)
%   Column  3:      Elevator deflection (deg)
%   Column  4:      Rudder deflection (rad)
%   Column  5:      Air density (kg/m^3)
%   Column  6:      Air speed (m/s)
%   Column  7:      Normal force (N)
%   Column  8:      Standard deviation of normal force (N)
%   Column  9:      Transverse force (N)
%   Column 10:      Standard deviation of transverse force (N)
%   Column 11:      Axial force (N)
%   Column 12:      Standard deviation of axial force (N)
%   Column 13:      Normal Moment (N-m)
%   Column 14:      Standard deviation of normal moment (N-m)
%   Column 15:      Transverse moment (N-m)
%   Column 16:      Standard deviation of transverse moment (N-m)
%   Column 17:      Axial moment (N-m)
%   Column 18:      Standard deviation of axial moment (N-m)
%
%   Note that column 4 will be zero, as there is no change in rudder
%   deflection. Further, if a user comment "Elevator = __" does not exist,
%   the parser will set every data value to zero and continue along to the
%   next data file. 
%
%   Created:            1/13/2011
%   Written by:         Garrison Hoe, Demoz Gebre-Egziabher & Hamid
%                       Mokhtarzadeh
%   Last Modified:      2/18/2017 M.S. Hemati
%                       -Format 2017
%                       2/19/2018 D. Gebre-Egziabher.
%                       -Format for 2018 (Retained AOA/rudder fields
%                           in comment line).
%                       2/3/2020 C. Sigfrid (TA).
%                       -Format for 2020 (3 new lines added to each dat
%                           file. Adjusted parser accordingly)
%                       2/11/2022 R. D. Halverson (TA)
%                       - Format for 2022 (Help header edits, new index
%                           values for change in crtunnel format, error 
%                           checking for bad comments)
%   Fixes
%                    
%                       2/26/2015  Peter, Sean, Demoz G. 
%                       - Format 2015
%
%                       2/11/2012 (Demoz Gebre-Egziabher)
%                       -Added GUI based file/directory selection
%                       -Increased data matrix size to 18 columns!!!!!
%                        (this will affect files that use data downstream)
%                       -Added angle of attack, elevator and rudder
%                       deflection scanning form *.dat file comment line
%
%=========================================================================%

clear all;
close all;
clc;

%   Define constants

max_elevator = 18;
d2r = pi/180;                           %   Degrees to radians
r2d = 180/pi;                           %   Radians to degrees
lb2N = 4.44822162;                      %   Pounds to Newtons
in2m = 0.0254;                          %   Inches to meters

%   Change working directory to where files are located

file_name_wild_card = '*.dat';
dialog_box_name = 'Select First File in Data Folder';
[tempFileName,file_dir_str,numFiles] = uigetfile(file_name_wild_card,dialog_box_name);

file_dir_listing_str = [file_dir_str file_name_wild_card];

file_list = dir(file_dir_listing_str);
[num_of_data_files temp] = size(file_list);

%   Sort the names so that '-2.dat' comes before '-10.dat' and so on.

temp = zeros(num_of_data_files,1);

for k=1:num_of_data_files
    
    temp1 = file_list(k).name;
    temp2 = strfind(temp1,'-');
    temp3 = strfind(temp1,'.');
    temp4 = temp1(temp2+1:temp3-1);
    temp(k) = str2num(temp4);
    
end

[new_order,old_order_index] = sort(temp,'ascend');

%   Define place holder
    
data_matrix = zeros(num_of_data_files,18);
    
%   Start cycling through each file
    
for k = 1:num_of_data_files
    
    %   Open data files for reading sequentially
    
    data_file_name=[file_dir_str file_list(old_order_index(k)).name];
    fid=fopen(data_file_name,'r');
    
    if(fid == -1)
        error(['Cannot find file named ',file_list(k).name]);
    end

    data_matrix(k,1) = k;
    
    %   Get angle of attack, elevator deflection and rudder deflection
    
    for k2 = 1:3                    %   Read lines 1 through 3
        line_read = fgets(fid);
    end
    
    user_comment_str = line_read;
    equal_sign_loc = findstr(user_comment_str,'=');
    
    if isempty(equal_sign_loc)
        continue
    end

    [a1,b1] = strtok(user_comment_str(equal_sign_loc(1)+1:end),',');
    data_matrix(k,2) = str2num(a1);         %   Angle of attack
%    data_matrix(k,3) = str2num(a1);
    
%    [a1,b1] = strtok(user_comment_str(equal_sign_loc(1)+1:end),',');
     [a1,b1] = strtok(user_comment_str(equal_sign_loc(2)+1:end),',');
     data_matrix(k,3) = str2num(a1);         %   Elevator deflection
%     
%    [a1,b1] = strtok(user_comment_str(equal_sign_loc(2)+1:end),',');
%     [a1,b1] = strtok(user_comment_str(equal_sign_loc(3)+1:end),',');
%     data_matrix(k,4) = str2num(a1);         %   Rudder deflection

    
    
    
    %   Get air density

    for k2 = 4:14                   %   Read lines 4 through 11
        line_read = fgets(fid);
    end

    [a1,b1] = strtok(line_read,'=');
    [a2,b2] = strtok(b1);
    data_matrix(k,5) = str2num(b2);

    %   Get air speed

    for k2 = 15:18                    %   Read lines 12 through 15
        line_read = fgets(fid);
    end

    [a1,b1] = strtok(line_read,'=');
    [a2,b2] = strtok(b1);
    data_matrix(k,6) = str2num(b2);

    %   Normal Force

    for k2 = 18:24                  %   Read lines 16 through 22
        line_read = fgets(fid);
    end

    [a1,b1] = strtok(line_read,'=');
    [a2,b2] = strtok(b1);
    [a3,b3] = strtok(b2,'lb');
    data_matrix(k,7:8) = lb2N*str2num(a3);

    %   Traverse Force

    [a1,b1] = strtok(fgets(fid),'=');
    [a2,b2] = strtok(b1);
    [a3,b3] = strtok(b2,'lb');
    data_matrix(k,9:10) = lb2N*str2num(a3);

    %   Axial Force

    [a1,b1] = strtok(fgets(fid),'=');
    [a2,b2] = strtok(b1);
    [a3,b3] = strtok(b2,'lb');
    data_matrix(k,11:12) = lb2N*str2num(a3);

    %   Normal Moment

    [a1,b1] = strtok(fgets(fid),'=');
    [a2,b2] = strtok(b1);
    [a3,b3] = strtok(b2,'in');
    data_matrix(k,13:14) = (lb2N*in2m)*str2num(a3);

    %   Traverse Moment

    [a1,b1] = strtok(fgets(fid),'=');
    [a2,b2] = strtok(b1);
    [a3,b3] = strtok(b2,'in');
    data_matrix(k,15:16) = (lb2N*in2m)*str2num(a3);

    %   Axial Moment

    [a1,b1] = strtok(fgets(fid),'=');
    [a2,b2] = strtok(b1);
    [a3,b3] = strtok(b2,'in');
	data_matrix(k,17:18) = (lb2N*in2m)*str2num(a3);

    %   Angle of Attack
    [a1,b1] = strtok(fgets(fid),'=');
    [a2,b2] = strtok(b1);
    [a3,b3] = strtok(b2,'deg');
    data_matrix(k,2) = str2num(a3);
    

    fclose(fid);
  %  disp(k)
end
    
%   Save the data file

temp = clock;
file_suffix = [num2str(temp(1)),'_', num2str(temp(2)),'_',num2str(temp(3)),...
                '_',num2str(temp(4)),'_',num2str(temp(5))];
            
save_file_name = ['CRT_data_',file_suffix,'.mat'];
eval(['save ',save_file_name,' data_matrix;']);