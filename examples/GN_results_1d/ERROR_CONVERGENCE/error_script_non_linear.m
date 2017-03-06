load('error_new.mat')
load('non_l_1.mat')
load('non_l_2.mat')
p1 = eta_L2err_1(1,:);
p2 = eta_L2err(2,:);
p1 = p1';
p2 = p2' ;
for i = 1 : 3
    p1_err(i,1) = abs(p1(4,1) - p1(i,1));
    p2_err(i,1) = abs(p2(4,1) - p2(i,1));
end
slope_1 = (log(p1_err(3,1)) - log(p1_err(2,1)))/(log(h(3,1)) - log(h(2,1)));
slope_2 = (log(p2_err(3,1)) - log(p2_err(2,1)))/(log(h(3,1)) - log(h(2,1)));
  
  