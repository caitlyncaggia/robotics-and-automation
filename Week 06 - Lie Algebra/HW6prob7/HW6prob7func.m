function gdot = HW6prob7func(t, g)

xi = [1; -1; 1/4];
theta = g(3);
Rot = [ cos(theta) , -sin(theta) , 0 ; 
         sin(theta), cos(theta), 0; 
         0, 0, 1 ];
gdot = Rot * xi;

end