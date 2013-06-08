function [] = remove_white(dirname)
% REMOVE_WHITE removes defected pixels in the base and the cross image data from
% every image file in 'dirname' directory.
% Base and cross images have to be stored in 'dirbase' directory.
%
% Written by AL and RG at 01-Nov-98
% Modified at:
% on Alex's PC

   direc=dir([dirname,filesep,'*.bmp']); filenames={};
	[filenames{1:length(direc),1}] = deal(direc.name);
	filenames=sortrows(char(filenames{:}));
	amount = max(str2num(filenames(:,end-8:end-6)));
	filebase = filenames(1,1:end-10);

  	base = imread([dirname,filesep,'defected',filesep,'1.bmp'],'bmp');
   cross = imread([dirname,filesep,'defected',filesep,'1c.bmp'],'bmp');

%%%%%%% START THE MAIN LOOP - FOR ALL FILES IN THE DIRECTORY  %%%%%%%%%%%

 for fileind = 1:2:2*amount-1	% main loop, for each pair of files

%----- Read Base & Cross images------ %

 a = imread([dirname,filesep,filenames(fileind,:)],'bmp');
 [b,map] = imread([dirname,filesep,filenames(fileind+1,:)],'bmp');

 % figure(1), imshow(b);

 a = uint8(abs(double(a) - double(base)));
 b = uint8(abs(double(b) - double(cross)));

% figure(2),imshow(b);
% figure(3), imshow(cross);

% pause


 imwrite(a,map,[dirname,filesep,filenames(fileind,:)],'bmp');
 imwrite(b,map,[dirname,filesep,filenames(fileind+1,:)],'bmp');

end