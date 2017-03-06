load determinant.mat
h = 1;
g = 1;
count = 1;
C_GN_airy = zeros(100,1);
x_axis = zeros(100,1);
for j = 0 : 0.1 : 10
    T = 2*pi/sqrt( g*j*tanh(j*h) );
    sig = 2*pi/T;
    y  = solve(detA-0,sigma);
    cgn  = y(4) ./ k ;
    airy = sqrt( g/j*tanh(j*h) );
    C_GN_airy(count,1) = subs(cgn,j)/airy ;
    x_axis(count,1) = sqrt((h/g)*g*j*tanh(j*h))/2/pi ;
    count = count + 1;
end