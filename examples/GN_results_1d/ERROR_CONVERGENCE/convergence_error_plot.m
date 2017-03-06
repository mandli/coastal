load('error_new.mat')
x1 = h;
y1 = p1;
marker1 = 'o';
lstyle1 = '--';
colorrgb1 = [1 0 0];
y2 = p2;
marker2 = '+';
lstyle2 = '-';
colorrgb2 = [ 0 1 0];
y3 = p3;
marker3 = 's';
lstyle3 = ':';
colorrgb3 = [0 0 1];
error_plot(x1,y1,y2,y3,marker1,marker2,marker3,lstyle1,lstyle2,lstyle3,colorrgb1,colorrgb2,colorrgb3);