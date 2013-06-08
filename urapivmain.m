function [varargout] = URAPIVMAIN(dirname,itt,spc,s2n_method,s2nlimit,sclt,mean10,crop)
% URAPIVmain- The first trial to implement LOCAL MEDIAN BASED
% PIV DATA VALIDATION. (See Westerweel's lecture for reference)
% This file is based on PIV_MAIN function, see help below.
%
%
% function [res] = URAPIVmain(dirname,itt,spc,s2n_method,s2nlimit,sclt,mean10)
%
% dirname = directory name
% itt = interrogation area in pixels,
% spc = grid spacing in pixels,
% s2n_method = 1 - for peak detectability (two highest peaks division)
%					2 - for peak-to-mean ratio
% s2limit = threshold for Signal-To-Noise ratio (1 leave S2N unused)
% sclt = scaling from pixel to m/sec 
% mean10 = threshold of GLOBAL mean filter
% crop = vector of lenght 4x1, with the following order:
%	[left right up bottom] crop values as number of lines*itt:
% For example, if you want to crop 2 lines from left (2 x itt pixels), 2 from
% right, 3 from up and 3 from bottom: crop should be [2 2 3 3];
% Default [0 0 0 0];

% Written by: Alex Liberzon and Roi Gurka,
%             
% at: 30-Jul-98
%
% Last modified at: 01-Nov-98.
% at f:\matlab\user\piv_software\ on Alex's PC.


%%%%%%%% START OF THE FUNCTION %%%%%%%%%%%%

% OPEN THE DIRECTORY AND READ THE FILE NAMES %
%	wd = cd;
%	cd(dirname);

if nargin ~=8		% only right number of inputs - no mistakes allowed
	error('Usage: urapivmain(dirname,itt,spc,s2n_method,s2nlimit,sclt,mean10,crop)');
end

if s2n_method ~=1 & s2n_method ~= 2
	error('Signal-To-Noise method is 1 for peak detectability'); 
	disp('or 2 for peak/mean ratio');
end



	direc=dir([dirname,filesep,'*.bmp']); filenames={};
	[filenames{1:length(direc),1}] = deal(direc.name);
	filenames=sortrows(char(filenames{:}));
	amount = max(str2num(filenames(:,end-8:end-6)));
	filebase = filenames(1,1:end-10);

%%%%%%% START THE MAIN LOOP - FOR ALL FILES IN THE DIRECTORY  %%%%%%%%%%%

 for fileind = 1:2:2*amount-1	% main loop, for each pair of files

%----- Read Base & Cross images------ %

 a = double(imread([dirname,filesep,filenames(fileind,:)],'bmp'))/255;
 b = double(imread([dirname,filesep,filenames(fileind+1,:)],'bmp'))/255;
 [sx,sy]=size(a);
 
% Just for a case that B is not the same size as A:

 [sxb,syb]=size(b);
 sx = min(sx,sxb); sy = min(sy,syb);
 


a = a(1+crop(1)*itt:spc*floor(sx/spc)-crop(2)*itt,1+crop(3)*itt:spc*floor(sy/spc)-crop(4)*itt);
b = b(1+crop(1)*itt:spc*floor(sx/spc)-crop(2)*itt,1+crop(3)*itt:spc*floor(sy/spc)-crop(4)*itt); 
 [sx,sy]= size(a);
 
 
% figure,imshow(a),figure,imshow(b); break

% AL & RG changed the following lines for the URI's satisfaction
% 13/09/98 1309_full
%  a = a(1:spc*floor(sx/spc)-9*(itt/2),1:spc*floor(sy/spc));
%  b = b(1:spc*floor(sx/spc)-9*(itt/2),1:spc*floor(sy/spc));

 
% MODIFIED BY A.L. at 27-Sep-98:  
 reslenx = (sx-itt)/spc+1;
 resleny = (sy-itt)/spc+1;
% 

 res = zeros(reslenx*resleny,5);
 resind = 0; 
 a2 = zeros(itt);
 b2 = zeros(itt);
 Nfft = 2*itt;



%%%%%% Start the loop for each interrogation block %%%%%%%

    for k=1:spc:sx-itt+1
      for m=1:spc:sy-itt+1

 disp([k,m])
        	
			a2 = a(k:k+itt-1,m:m+itt-1);
         b2 = b(k:k+itt-1,m:m+itt-1);
         a2 = a2 - mean2(a2); 
			b2 = b2 - mean2(b2);
 			b2 = b2(itt:-1:1,itt:-1:1); 
        	ffta=fft2(a2,Nfft,Nfft);
			fftb=fft2(b2,Nfft,Nfft);
      	c = real(ifft2(ffta.*fftb));

% Find your majour peak - displacement:

peak1 = max(c(:));
% The following IF.. statement means that interrogation area
% contains only black values and therefore we have to throw the
% result out. Maybe we should think more about this option?
if ~peak1
   peak1 = eps;
   pixi = 1;        % It doesn't matter, we can put here pixi = Nfft;
   pixj = 1;        % or other "flag" of bad data!
   disp('Black, black, black !!!'); 
else
   [pixi,pixj]=find(c==peak1);
end


% Modified at 25-Sep-98, by Alex:
%
% The decision is not to use PEAK DETECTABILITY FUNCTION
% that is actually the ration between the maximum peak and the 
% second highest peak, but use the ratio between the highest peak
% and the mean value of the rest ??? of the matrix. A couple of questions exist
% 1. Should we use the mean value of the whole matrix? of the matrix without the
%    maximum value? 
% 2. Should we use the ABS? or only POSITIVE values?
%
% 

      tmp = c;
		tmp(pixi,pixj) = 0;

 if pixi==1 | pixj==1 | pixi==Nfft | pixj==Nfft		% if the peak on the border

       peak2 = peak1;

 else 			% look for second peak

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if s2n_method == 1		% PEAK DETECTABILITY METHOD NEEDS SECOND PEAK SEARCH:
        	tmp(pixi-1:pixi+1,pixj-1:pixj+1) = NaN;
        	
        	peak2 = max(tmp(:));
        	[x2,y2] = find(tmp==peak2);
        	tmp(x2,y2) = NaN;
        	
    if x2 > 1 & y2 > 1 & x2 < Nfft & y2 < Nfft


        	while peak2 < max(max(c(x2-1:x2+1,y2-1:y2+1)))
        	
      		peak2 = max(tmp(:));
        		[x2,y2] = find(tmp==peak2);
        		
				if x2 == 1 | y2==1 | x2 == Nfft | y2 == Nfft 
					peak2 = peak1;	% will throw this one out
					break;
				end
			
				tmp(x2,y2) = NaN;
			end		% end of while
      else			% second peak on the border means "second peak doesn't exist"
				peak2 = peak1;
  end    % if x2 >1 ......end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif s2n_method == 2	% PEAK-TO-MEAN VALUE RATIO:

	% If the highest peak is not on the border, then:
	%-----------------------------------------------%
	% Mean value of the matrix without maximum peak %
	%-----------------------------------------------%	
		peak2 = mean2(abs(tmp));		

end		% end of second peak search, both methods.
		

end				% end of if highest peak on the border

 

% Now let's find Signal-To-Noise Ratio:       	
        
         if ~peak2
        	     s2n = Inf;		% Just to protect from zero.
         else 
        	     s2n = peak1/peak2; 
         end

% !!!!!!!!!!!! NOTICE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% We have to find the appropriate S2Nlimit value!

         if s2n < s2nlimit
         		peakx = itt;		% Alex, 19-10-98
         		peaky = itt;
%			disp(' BAD Signal-to-Noise Ratio'); pause(1);  %debugging option
         else

%%%%%%%% Gaussian Interpolation %%%%%%%%%%%%%

%%%%%%%% ELIMINATE THE PROBLEM WITH NEGATIVE VALUES %%%%%%%%
%		if any(any(c<0))				% if there is at least one negative value
%			   c  = c - min(c(:)) + 1;	% 18-10-1998, by Alex Liberzon
% 		end

% disp([pixi pixj c(pixi,pixj)])
% keyboard
 	      	f0 = log(c(pixi,pixj));
	      	f1 = log(c(pixi-1,pixj));
	      	f2 = log(c(pixi+1,pixj));
	      	peakx = pixi+ (f1-f2)/(2*f1-4*f0+2*f2);  
%
  	      	f0 = log(c(pixi,pixj));
  	      	f1 = log(c(pixi,pixj-1));
  	      	f2 = log(c(pixi,pixj+1));
 	      	peaky = pixj+ (f1-f2)/(2*f1-4*f0+2*f2); 

	      	if ~isreal(peakx) | ~isreal(peaky)
					peakx = itt;						% Alex , 19-10-98
					peaky = itt;
			disp('NEGATIVE value');
	      	end

         end 		

%%%%%%% Sub-pixel displacement to velocity scaling %%%%%%%%%%

	      u = (itt-peaky)*sclt;
	      v = (itt-peakx)*sclt;	% Note the directions !
	      x = m+itt/2-1;
	      y = k+itt/2-1;
			resind = resind + 1;
	      res(resind,:) = [x y u v s2n];
%
     end				
  end	

% NO_FILT_RES will be written in another TXT file, containing the result of 
% cross-correlation and first signal-to-noise filtering, but before the
% local filtering:

 no_filt_res = res;

	
%%%%%%%%%% END OF THE LOOP FOR EACH INTERROGATION BLOCK %%%%%%%%%

%%%%%%%%%%% THROW OUT OUTLAYERS BY GLOBAL MEAN FILTER %%%%%%%%%%%%%%

% Reshape U and V matrices in two-dimensional grid and produce 
% velocity vector in U + i*V form (real and imaginary parts):

  u = reshape(res(:,3),resleny,reslenx);
  v = reshape(res(:,4), resleny,reslenx);
  vector = u + sqrt(-1)*v;

 
  % Modified by A.L. at 27-Sep-98 to use NANMEAN function.
  %	meanvel = sqrt(res(:,3).^2+res(:,4).^2);
 	
	vector(abs(vector)>mean(abs(vector(find(vector))))*mean10) = 0;

   	u = real(vector);
	v = imag(vector);


%%%%%% END OF THE GLOBAL MEAN FILTER %%%%%%%%%%%%%%%%

%%%%%% LOCAL FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%

% The heart of the local filter is its kernel:

kernel = [-1 -1 -1; -1 8 -1; -1 -1 -1];

% The filtering is by founding the outlayers:

tmpv = abs(conv2(v,kernel,'same'));
tmpu = abs(conv2(u,kernel,'same'));

% WE HAVE TO DECIDE WHICH LIMIT TO USE:
% 1. Mean + 3*STD for each one separately OR
% 2. For velocity vector length (and angle)
% 3. OR OTHER.

lmtv = mean(tmpv(find(tmpv))) + 3*std(tmpv(find(tmpv)));
lmtu = mean(tmpu(find(tmpu))) + 3*std(tmpu(find(tmpu)));

u_out = find(tmpu>lmtu);
v_out = find(tmpv>lmtv);

% Here we can add the union of the indeces:
%
% outlayers = union(u_out,v_out);

% Let's throw the outlayers out:

u(u_out) = 0; u(v_out) = 0;
v(v_out) = 0; v(u_out) = 0;

vector = u + sqrt(-1)*v;


%%%%%%% END OF THE LOCAL FILTER %%%%%%%%%%%%%%

  res(:,3) = reshape(real(vector),resleny*reslenx,1);
  res(:,4) = reshape(imag(vector),resleny*reslenx,1);

filt_res = res;
  
%%%%%%% FILL EMPTY PLACES WITH AVERAGE OF NEIGHBORS %%%%%%%

  [indx,indy] = find(~vector);

 while ~isempty(indx)

  for z=1:length(indx)
  		k = [max(3,indx(z))-2:min(resleny-2,indx(z))+2];
  		m = [max(3,indy(z))-2:min(reslenx-2,indy(z))+2];
%%%
      	tmpvec = vector(k,m);
			tmpvec = tmpvec(find(tmpvec));
if ~isempty(tmpvec)
  vector(indx(z),indy(z)) = mean(real(tmpvec))+ sqrt(-1)*mean(imag(tmpvec));
else
	vector(indx(z),indy(z)) = 0;
end

  end	% of for

  [indx,indy] = find(~vector);
 end % of while


  res(:,3) = reshape(real(vector),resleny*reslenx,1);
  res(:,4) = reshape(imag(vector),resleny*reslenx,1);
  
  
%%%%%%%% WRITE TEXT FILE OF RESULTS %%%%%%%%%%%%%%%5
  
  fid = fopen([dirname,filesep,filenames(fileind,1:end-6),'.txt'],'w');
  fprintf(fid,'%3d %3d %7.4f %7.4f %7.4f\n',res');
  fclose(fid);


% Comment next 3 lines, if you don't need the unfiltered results of
% the cross-correlation  
	fid = fopen([dirname,filesep,filenames(fileind,1:end-6),'_noflt.txt'],'w');
	fprintf(fid,'%3d %3d %7.4f %7.4f %7.4f\n',no_filt_res');
	fclose(fid);
% Same thing if you don't want filtered result, but without filling of
% mean values
 	fid = fopen([dirname,filesep,filenames(fileind,1:end-6),'_flt.txt'],'w');
	fprintf(fid,'%3d %3d %7.4f %7.4f %7.4f\n',filt_res');
	fclose(fid);
 	 
end

%%%%%%%%% END OF THE MAIN LOOP - ALL BMP FILES ARE PROCEEDED %%%%%%%%%%
 
%%%%%%%% END OF THE FUNCTION %%%%%%%%%%%%%%%
