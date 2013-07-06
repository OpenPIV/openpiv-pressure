function output = piv_poisson(dirname,mm_per_pixel)
%
% OUTPUT = PIV_POISSON(DIRNAME, MM_PER_PIXEL)
% 
% estimates pressure in the flow field of the air impinging jet analysed 
% using OpenPIV software ONLY (!). The Neumann boundary conditions on the 
% lowest horizontal surface and Dirichlet boundary conditions are HARDCODED
% and should NOT be used without proper modifications for other flow cases.
% 
% Inputs:
%       dirname = directory that contains TXT files analysed by OpenPIV 
%       mm_per_pixel = scaling of millimeters to pixels
%
% Outputs:
%      output = structure that contains the following fields:
%           x,y = in milimmeters
%           data(n_files).{u,v} = velocities 
%           data(n_files).pressure = instantaneous pressure, result of the
%           Poisson equation
%           mean_pressure, mean_u, mean_v
%           inputs: directory, spacing, scaling (for preservation)
%           
% Example:
%       dirname = '/Users/alex/Desktop/imp_pep_3d';
%       mm_per_pixel = 0.01/120;
%       out = piv_poisson(dirname,mm_per_pixel)
%       surf(out.x(1,:),out.y(:,1),out.mean_pressure)
%       xlabel('x [mm]'); ylabel('y, [mm]'); 
%       zlabel('Pressure [Pa]');
%       
%       
% Assumptions: 
%       \rho = air density is taken as 1 for simplicity
%       \delta x = \delta y, i.e spacing in x and y directions is equal
%       flow is an impinging jet coming from top and hitting the flat 
%       horizontal plate at the bottom, see the setup in the reference.
%
% Reference:
%  Gurka R., Liberzon A., Hefetz D., Rubinstein D. and Shavit U.,
% "Computation of Pressure Distribution Using PIV Velocity Data", 
%  3rd International Workshop on Particle Image Velocimetry, Santa Barbara, 
%  California, September 16-18, 1999.

%
% Copyright (c) 1998, Alex Liberzon and Roi Gurka
%
%
% Last modified at July 6, 2013 by Alex:

% Check the number of inputs/outputs
if nargout ~= 1 || nargin ~= 2
    help piv_poisson
    return
end



ro = 1; % density, kg/m^3, % ro = 1.293 e-9;

data  = struct('u',[],'v',[],'pressure',[]);


% mm_per_pixelrogation size, total size and step sizes initialization:



% read the list of TXT files
d = get_list_txt_files(dirname);



n_files = length(d);
data = repmat(data,n_files,1);



% Load the first file in order to get the rectangular matrix shape:

fid = fopen(fullfile(dirname,d(1).name),'r');
[tmp1,count] = fscanf(fid,'%f %f %f %f %f');
fclose(fid);
tmp1  = reshape(tmp1, 5, count/5)';
tmp1 = sortrows(tmp1,1);	% Alex, 23/08/98

minx = min(tmp1(:,1));
maxx = max(tmp1(:,1));
miny = min(tmp1(:,2));
maxy = max(tmp1(:,2));

if max(diff(tmp1(:,1))) ~= max(diff(tmp1(:,2)))
    error('The solution is valid for the equally spaced matrix in x and y directions only');
else
    spc = max(diff(tmp1(:,1)));
end

dx = spc*mm_per_pixel;	
n_rows = (maxx-minx)/spc + 1;
n_cols = (maxy-miny)/spc + 1;


% Assign the first data matrix

[x,y,data(1).u,data(1).v]  = ...
    parse_tmp_matrix(tmp1,n_rows,n_cols);


% Report the number of files, or exit with error, if no file has been found
if n_files
    disp(['Found ',num2str(n_files),' file(s)']);
else
    error('No file has been found');
end



for ind = 1:n_files
    fid = fopen(fullfile(dirname,d(ind).name),'r');
    [tmp,count] = fscanf(fid,'%f %f %f %f %f');
    fclose(fid);
    tmp  = reshape(tmp, 5, count/5)';
    tmp = sortrows(tmp,1);	% Alex, 23/08/98
    [~,~,data(ind).u,data(ind).v]  = ...
        parse_tmp_matrix(tmp,n_rows,n_cols);
    
end




for k = 1:n_files
    
    
    [dudx,dudy] = gradient(data(k).u);
    [dvdx,dvdy] = gradient(data(k).v);
    
    rhsv = dudx.^2 + dvdy.^2 + 2*dudy.*dvdx;
    
    % Pressure matrix initialization with boundary conditions:
    P = zeros(n_rows,n_cols);
    PP = P;
    % Optimal lambda:
    % Ref: http://www.enm.bris.ac.uk/admin/courses/ANA&PDEs/pde-lect8.pdf
    tmpc = cos(pi/(31+1)) + cos(pi/(30+1));
    lambda = 4/(2+sqrt(4-tmpc^2));
    % lambda = 1.8; % weighting factor
    
    % Poisson equation solution by Liebmann's (iterative) method
    
    tol = 1e-3;		 % error is 1%
    maxerr = inf;	 % initial error
    iter = 0;
    while maxerr > tol
        iter = iter + 1;
        
        disp(['Iteration no. ',num2str(iter)]);
        for c = 2:n_cols-1
            for r = 2:n_rows-1
                P(r,c)=(P(r+1,c) + P(r-1,c) + P(r,c+1) + P(r,c-1) + ro*rhsv(r,c)*dx^2)/4;
                P(r,c) = lambda*P(r,c) + (1-lambda)*PP(r,c);
            end
            % Neuman boundary condition along the bottom wall for
            % the impinging jet problem.
            % TODO: make it generic
            r = n_rows;
            P(r,c)=(2*P(r-1,c)+P(r,c+1)+P(r,c-1)+ro*rhsv(r,c)*dx^2)/4;
            P(r,c)=lambda*P(r,c)+(1-lambda)*PP(r,c);
        end
        
        maxerr = max(max(abs((P-PP)./P)));
        disp(['Maximum error is  ',num2str(maxerr)]);
        PP = P;
    end			% as long the error larger than tolerance, continue
    
    data(k).pressure = P;
end		% for n_files


% Ensemble averages:

mean_u = data(1).u;
mean_v = data(1).v;
mean_pressure = data(1).pressure;

for i = 2:n_files
    mean_pressure = mean_pressure + data(i).pressure;
    mean_u = mean_u + data(i).u;
    mean_v = mean_v + data(i).v;
end
mean_pressure = mean_pressure/n_files;
mean_u = mean_u/n_files;
mean_v = mean_v/n_files;

% preserve the inputs with the output


output = struct('data',data,'mean_pressure',mean_pressure,'U',mean_u,...
    'V',mean_v,'directory', dirname,...
    'scaling',mm_per_pixel,'x',x,'y',y);


