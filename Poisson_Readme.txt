General
^^^^^^^

number of layers in the final 3D matrix is,
		5n+7
where n is the number of realizations.

The values are complex to provide the x componenet
as real and y component as imaginary

 7
^^^
1	Location
2	Mean
3	Raynolds Stress
4	rms in x
5	rms in y
6	Mean Pressure
7	rms Pressure

 5
^^^
1	Velocity
2	Fluctuation
3	Fluctuation Square
4	Pressure
5	Pressure Fluctuation

=============

To compute the mean values
^^^^^^^^^^^^^^^^^^^^^^^^^^
for Dana I created a directory 6-5-99-txt containing
the 39 filtered interpolated txt files. The dir location
was parallel to the 6-5-99 dir.

I was above the dirs and used the command:
[final,iter]=piv_poisson('6-5-99-txt',32,32,3.0000e-005,0);

To see the results I used:
>> a=final;
>> mesh(imag(-a(:,:,41)));
>> mesh(real(-a(:,:,41)));

>> quiver(real(a(:,:,1)),imag(a(:,:,1)),real(a(:,:,41)),imag(a(:,:,41)));
>> zoom on

To save results
^^^^^^^^^^^^^^^
>> x=real(a(:,:,1));
>> y=imag(a(:,:,1));
>> u=real(a(:,:,41));
>> v=imag(a(:,:,41));
>> pwd

ans =

D:\Users\Dana

>> cd 6-5-99-txt
>> save x.txt x
>> save x.txt x -ascii
>> save y.txt y -ascii
>> save u.txt u -ascii
>> save v.txt v -ascii
>> save xyuv.txt x y u v -ascii

To Read results later
^^^^^^^^^^^^^^^^^^^^^
>> x=load('x.txt');
>> y=load('y.txt');
>> u=load('u.txt');
>> v=load('v.txt');
>> quiver(x,y,u,v);



