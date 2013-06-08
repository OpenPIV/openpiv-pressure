function [final,iter] = piv_poisson(dirname,spc,inter,flag)
%
% PIV_POISSON - calculates velocity flactuations, mean velocity and
%		pressure field for the free flow outside the channel tube !!! (only)
%
%
%	flag = 1 - calculate Pressure field
%	flag = 0 - calculate ONLY velocities, averages, fluctuations and Reynold stresses.
%	Last modified at:10 November, 1998
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
%	inter = 0.01/120;		% 0.01 m = 120 pixels = scaling
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: All data matrices have to be of the same size!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Some initializations

% global right dudx dvdy dudy dvdx

% "ro" is a density of the air in ??? units

% ro = 1.293 e-9;

ro = 1;
iter = 0;

% warning off

% Check the number of inputs/outputs


if nargout ~= 2 || nargin ~= 4
    disp('Usage:    [final,iter] = piv_poisson(dirname,spc,scale,flag)  ');
    return;
end

% Keep current directory, change to dirname, read all txt files

% wd = cd;
% cd(dirname);
% fullfile removes the need to change directories [AL, 8.6.13]
d1 = dir(fullfile(dirname,'*.txt'));
d2 = dir(fullfile(dirname,'*_noflt.txt'));
d3 = dir(fullfile(dirname,'*_flt.txt'));

% the following loops remove _flt and _noflt.txt from the list
% of names. 

d = d1;
ind = false(1,length(d1));

for i = 1:length(d1)
    for j = 1:length(d2)
        if strcmp(d1(i).name,d2(j).name), ind(i) = true; end
    end
    for j = 1:length(d3)
        if strcmp(d1(i).name,d3(j).name), ind(i) = true; end
    end
end
d(ind) = [];



file_num = length(d);

% Report the number of files, or exit with error, if no file has been found

if file_num
    disp(['Found ',num2str(file_num),' file(s)']);
else
    error('No file has been found');
end

% Load all files and find a mininum/maximum of X and Y values

fid = fopen(fullfile(dirname,d(1).name),'r');
[tmp1,count] = fscanf(fid,'%f %f %f %f %f');
fclose(fid);
tmp1  = reshape(tmp1, 5, count/5)';
tmp1 = sortrows(tmp1,1);	% Alex, 23/08/98


minx = min(tmp1(:,1));
maxx = max(tmp1(:,1));
miny = min(tmp1(:,2));
maxy = max(tmp1(:,2));


for ind = 2:file_num
    fid = fopen(fullfile(dirname,d(ind).name),'r');
    [tmp,count] = fscanf(fid,'%f %f %f %f %f');
    fclose(fid);
    tmp  = reshape(tmp, 5, count/5)';
    tmp = sortrows(tmp,1);	% Alex, 23/08/98
    tmp1 = cat(3,tmp1,tmp);
end


% Interrogation size, total size and step sizes initialization:

dx = spc*inter;		% spc = spacing grid
dy = spc*inter;		% dx,dy = scaled grid in X,Y directions

winx = (maxx-minx)/spc + 1;
winy = (maxy-miny)/spc + 1;

pos = tmp1(:,1,1) + 1j*tmp1(:,2,1);
pos = reshape(pos, winy, winx);
right = pos;


vel = tmp1(:,3,:) + 1j*tmp1(:,4,:);
vel = reshape(vel,winy,winx,size(tmp1,3));

final = cat(3,pos,vel);

% Concatenate matrix of average velocities:

% Removed some unnecessary check, [AL, 8.6.13]
% if file_num > 1
    final = cat(3,final,mean(final(:,:,2:end),3));
% else
    % final = cat(3, final, final(:,:,end));		% in testing we put only one file
% end

% Velocity Fluctuations matrices:

[~,~,len] = size(final);

for ind = 2:file_num+1
    fluct = final(:,:,ind) - final(:,:,len);
    final = cat(3,final,fluct);
end


% Reynolds stress, 24/08/98

[row,col,len] = size(final);

ii = len - file_num+1;

if file_num > 1
    reynolds = mean(real(final(:,:,ii:len)).*imag(final(:,:,ii:len)),3);
else
    reynolds = zeros(row,col);
end

final = cat(3,final,reynolds);

%%%% 10-Nov-98, RMS velocity %%%%%%%%%%%

for ind = file_num+3:2*file_num+2
    turb_int_u = real(final(:,:,ind)).^2;
    turb_int_v=imag(final(:,:,ind)).^2;
    final = cat(3,final,turb_int_u + 1j*turb_int_v);
end

if file_num > 1
    turb_int_avg_u=sqrt(mean(real(final(:,:,2*file_num+4:end)),3));
    turb_int_avg_v=sqrt(mean(imag(final(:,:,2*file_num+4:end)),3));
else
    turb_int_avg_u = zeros(row,col);
    turb_int_avg_v = zeros(row,col);
    
end
final = cat(3,final,turb_int_avg_u);
final = cat(3,final,turb_int_avg_v);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% IF FLAG == 0, DON'T CONTINUE WITH PRESSURE CALCULATIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag
    
    % du/dx, dv/dy, du/dy, dv/dx calculations:
    
    dudx=zeros(row,col);
    dudy=zeros(row,col);
    dvdx=zeros(row,col);
    dvdy=zeros(row,col);
    
    for k = 2:file_num + 1	% for all velocity matrices
        
        u=real(final(:,:,k));
        v=imag(final(:,:,k));
        
        % boundary columns:
        
        dudx(:,1)=(u(:,2)-u(:,1))/dx;	% forward differences
        dvdx(:,1)=(v(:,2)-v(:,1))/dx;
        dudx(:,col)=(u(:,col)-u(:,col-1))/dx;	% backward --//--
        dvdx(:,col)=(v(:,col)-v(:,col-1))/dx;
        
        % the rest of the columns:
        
        for c = 2:col-1
            dudx(:,c) = (u(:,c+1)-u(:,c-1))/2/dx;	% central differences
            dvdx(:,c) = (v(:,c+1)-v(:,c-1))/2/dx;	% central differences
        end
        
        % boundary rows:
        
        dudy(1,:)=(u(2,:)-u(1,:))/dy;	% forward differnces
        dvdy(1,:)=(v(2,:)-v(1,:))/dy;
        dudy(row,:)=(u(row,:)-u(row-1,:))/dy;	% backward ...
        dvdy(row,:)=(v(row,:)-v(row-1,:))/dy;
        
        % the rest of the rows
        
        for r = 2:row-1
            dudy(r,:)=(u(r+1,:)-u(r-1,:))/2/dy;	% central ...
            dvdy(r,:)=(v(r+1,:)-v(r-1,:))/2/dy;
        end
        
        %
        
        rhsv = dudx.^2+dvdy.^2+2*dudy.*dvdx;
        
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
                % 	    r = row;		% Neuman boundary condition:
                % 	    P(r,c)=(2*P(r-1,c)+P(r,c+1)+P(r,c-1)+ro*rhsv(r,c)*dx^2)/4;
                % 	    P(r,c)=lamda*P(r,c)+(1-lamda)*PP(r,c);
            end
            
            maxerr = max(max(abs((P-PP)./P)));
            disp(['Maximum error is  ',num2str(maxerr)]);
            PP = P;
        end			% while the error larger than tolerance
        
        final = cat(3,final,P);
        
    end			% for all velocity matrices
    
    
    % Pressure mean value, 24/08/98
    
    if file_num > 1
        final = cat(3,final,mean(final(:,:,end-file_num+1:end),3));
    else
        final = cat(3, final,final(:,:,end));		% in testing we put only one file
    end
    
    % Pressure fluctuations, 24/08/98
    
    [~,~,len] = size(final);
    
    for ind = len-file_num:len-1
        fluct = final(:,:,ind) - final(:,:,len);
        final = cat(3,final,fluct);
    end
    
    press_intens = sqrt(mean(final(:,:,end-file_num+1:end).^2,3));
    final = cat(3,final,press_intens);
    
end			% if flag == 1

outname = strtok(d(1).name,'.');

save(outname,'final','iter');


% cd(wd);
