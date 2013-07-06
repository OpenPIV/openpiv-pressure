function [final,iterations] = piv_poisson(dirname,vector_grid_spacing_pixels,mm_per_pixel)
%
% PIV_POISSON - calculates velocity flactuations, mean velocity and
%		pressure field for the free flow outside the channel tube !!! (only)
%
%
%	flag = 1 - calculate Pressure field
%	flag = 0 - calculate ONLY velocities, averages, fluctuations and Reynold stresses.
%	Last modified at:10 November, 1998
%
% 
% number of layers in the final 3D matrix is,
%		5n+7
% where n is the number of realizations.
%
% The values are complex to provide the x componenet
% as real and y component as imaginary
% 
% 7
% ^^^
% 1	Location
% 2	Mean
% 3	Raynolds Stress
% 4	rms in x
% 5	rms in y
% 6	Mean Pressure
% 7	rms Pressure
%  
%  5
%  ^^^
%  1	Velocity
%  2	Fluctuation
% 3	Fluctuation Square
% 4	Pressure
%  5	Pressure Fluctuation
%
%
% Copyright (c) 1998, Alex Liberzon & Roi Gurka, Dep. of Agr.Eng, Technion
% A special thanks to Asaf Arnon, Dep. of Agr.Eng., Technion
%
% Last modified at August, 19 by Alex:
%
% Last mofifications:
% 1. reading all txt files in the given directory.
% 2. Very, very nice vectorization of reading all data and
%    arrange in the FINAL matrix, see tmp1, minx, pos, vel...
%
%	mm_per_pixel = 0.01/120;		% 0.01 m = 120 pixels = scaling
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: All data matrices have to be of the same size!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Some initializations

% global right dudx dvdy dudy dvdx

% "ro" is a density of the air in ??? units

% ro = 1.293 e-9;


% Check the number of inputs/outputs
if nargout ~= 3 || nargin ~= 4
    disp('Usage:    [data,iterations] = piv_poisson(dirname,vector_grid_spacing_pixels,scale,flag)  ');
    return;
end



ro = 1; % density, kg/m^3, % ro = 1.293 e-9;
iter = 0;

data  = struct('x',[],'y',[],'u',[],'v',[],'dudx',[],'dudy',[],'dvdx',...
    [],'dvdy','pressure',[]);


% mm_per_pixelrogation size, total size and step sizes initialization:

dx = vector_grid_spacing_pixels*mm_per_pixel;		% vector_grid_spacing_pixels = spacing grid
dy = vector_grid_spacing_pixels*mm_per_pixel;		% dx,dy = scaled grid in X,Y directions


% read the list of TXT files
d = get_list_txt_files(dirname);



file_num = length(d);
data = repmat(data,file_num,1);



% Load the first file in order to get the rectangular matrix shape:

fid = fopen(fullfile(dirname,d(1).name),'r');
[tmp1,count] = fscanf(fid,'%f %f %f %f %f');
fclose(fid);
tmp1  = reshape(tmp1, 5, count/5)';
tmp1 = sortrows(tmp1,1);	% Alex, 23/08/98

minx = min(data(1).x);
maxx = max(data(1).x);
miny = min(data(1).y);
maxy = max(data(1).y);

winx = (maxx-minx)/vector_grid_spacing_pixels + 1;
winy = (maxy-miny)/vector_grid_spacing_pixels + 1;


% Assign the first data matrix

[data(1).x,data(1).y,data(1).u,data(1).v]  = ...
    parse_tmp_matrix(tmp1,winx,winy);




% Report the number of files, or exit with error, if no file has been found

if file_num
    disp(['Found ',num2str(file_num),' file(s)']);
else
    error('No file has been found');
end



for ind = 1:file_num
    fid = fopen(fullfile(dirname,d(ind).name),'r');
    [tmp,count] = fscanf(fid,'%f %f %f %f %f');
    fclose(fid);
    tmp  = reshape(tmp, 5, count/5)';
    tmp = sortrows(tmp,1);	% Alex, 23/08/98
    [data(ind).x,data(ind).y,data(ind).u,data(ind).v]  = ...
    parse_tmp_matrix(tmp1,winx,winy);

end


% pos = tmp1(:,1,1) + 1j*tmp1(:,2,1);
% pos = reshape(pos, winy, winx);
% right = pos;
% 
% 
% vel = tmp1(:,3,:) + 1j*tmp1(:,4,:);
% vel = reshape(vel,winy,winx,size(tmp1,3));
% 
% final = cat(3,pos,vel);


% Ensemble average velocities:
data.mean_u = 0*data(1).u;
data.mean_v = data.mean_u;
for i = 1:file_num
    data.mean_u = data.mean_u + data(i).u;
    data.mean_v = data.mean_v + data(i).v;
end
data.mean_u = data.mean_u/file_num;
data.mean_v = data.mean_v/file_num;

% % Velocity Fluctuations matrices:
% 
% len = size(final,3);
% % Fluctuations of velcocity:
% for ind = 2:file_num+1
%     fluct = final(:,:,ind) - final(:,:,len);
%     final = cat(3,final,fluct);
%     par = cat(2, par,'uf,vf');
% end



% Reynolds stress, 24/08/98

% [row,col,len] = size(final);
% 
% ii = len - file_num+1;
% 
% if file_num > 1
%     reynolds = mean(real(final(:,:,ii:len)).*imag(final(:,:,ii:len)),3);
% else
%     reynolds = zeros(row,col);
% end
% 
% final = cat(3,final,reynolds);
% par = cat(2,par,'Rs');
% 
% %%%% 10-Nov-98, RMS velocity %%%%%%%%%%%
% 
% for ind = file_num+3:2*file_num+2
%     turb_int_u = real(final(:,:,ind)).^2;
%     turb_int_v=imag(final(:,:,ind)).^2;
%     final = cat(3,final,turb_int_u + 1j*turb_int_v);
%     par = cat(2,par,'uf/U,vf/V');
% end
% 
% if file_num > 1
%     turb_int_avg_u=sqrt(mean(real(final(:,:,2*file_num+4:end)),3));
%     turb_int_avg_v=sqrt(mean(imag(final(:,:,2*file_num+4:end)),3));
% else
%     turb_int_avg_u = zeros(row,col);
%     turb_int_avg_v = zeros(row,col);
%     
% end
% final = cat(3,final,turb_int_avg_u + 1j*turb_int_avg_v);
% % final = cat(3,final,turb_int_avg_v);
% par   = cat(2,par,'mean. turb.int');

   
    
    for k = 1:file_num	
        
 
        [data(k).dudx,data(k).dudy] = gradient(data(k).u);
        [data(k).dvdx,data(k).dvdy] = gradient(data(k).v);
        
        
        %
        
        rhsv = data(k).dudx.^2 + data(k).dvdy.^2 + ...
            2*data(k).dudy.*data(k).dvdx;
        
        right = cat(3,right,rhsv);
        
        % Pressure matrix initialization with boundary conditions:
        
        P = zeros(row,col);
        %	P(1,:)=zeros(1,col);
        %	P(:,1)=zeros(row,1);
        %	P(:,col)=P(:,1);
        PP = P;
        lamda = 1.8; % weighting factor
        
        % Poisson equation solution by Liebmann's (iterative) method
        
        tol = 1e-3;		% error is 1%
        maxerr = inf;		% initial error
        iter = 0;
        while maxerr > tol
            iter = iter+1;
            
            disp(['Iteration no. ',num2str(iter)]);
            for c = 2:col-1
                for r = 2:row-1
                    P(r,c)=(P(r+1,c) + P(r-1,c) + P(r,c+1) + P(r,c-1) + ro*rhsv(r,c)*dx^2)/4;
                    P(r,c) = lamda*P(r,c) + (1-lamda)*PP(r,c);
                end
                    r = row;		% Neuman boundary condition:
                	P(r,c)=(2*P(r-1,c)+P(r,c+1)+P(r,c-1)+ro*rhsv(r,c)*dx^2)/4;
                	P(r,c)=lamda*P(r,c)+(1-lamda)*PP(r,c);
            end
            
            maxerr = max(max(abs((P-PP)./P)));
            disp(['Maximum error is  ',num2str(maxerr)]);
            PP = P;
        end			% while the error larger than tolerance
        
        final = cat(3,final,P);
        par = cat(2,par,'P');
        
        
    end			% for all velocity matrices
    
    
    % Pressure mean value, 24/08/98
    
    if file_num > 1
        final = cat(3,final,mean(final(:,:,end-file_num+1:end),3));
    else
        final = cat(3, final,final(:,:,end));		% in testing we put only one file
    end
    
    par = cat(2,par,'Mean P');
    
    % Pressure fluctuations, 24/08/98
    
    [~,~,len] = size(final);
    
    for ind = len-file_num:len-1
        fluct = final(:,:,ind) - final(:,:,len);
        final = cat(3,final,fluct);
        par = cat(2,par,'pf');
    end
    
    press_intens = sqrt(mean(final(:,:,end-file_num+1:end).^2,3));
    final = cat(3,final,press_intens);
    par = cat(2,par,'P int.');
    

outname = strtok(d(1).name,'.');

save(outname,'final','iter','par');


% cd(wd);


