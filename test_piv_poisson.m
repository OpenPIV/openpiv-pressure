        dirname = './imp_3';
        mm_per_pixel = 0.01/120;
        out = piv_poisson(dirname,mm_per_pixel)
        surf(out.x(1,:),out.y(:,1),out.mean_pressure)
        xlabel('x [mm]'); ylabel('y, [mm]'); 
        zlabel('Pressure [Pa]');