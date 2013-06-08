'test_comp' includes 'urapivmain.m' to analyze the pictures 
and a folder '10_pic' includes 10 pairs of images (Q=30, ratio=0.2)


Open Matlab 5.1

call directory "test_comp"

write:
	urapivmain('10_pic',64,32,2,17,4.44,4,[0 0 0 0]);

		10 txt files hane been created in the directory '10_pic'

		and also another 20 not filtered txt files in the same directory
