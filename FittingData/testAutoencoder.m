% establish axis mapping:
% log plot
x_label = [0.2,0.3,0.4,0.5,1,1.5,2];
log_x_label = log(x_label);
pixel_x = [199,319,403,469,673,792,877];
fit_x = fit(pixel_x',log_x_label','poly1');
y_label = [60,100,150,200,300,400,500,1000,1500,2000,3000];
log_y_label = log(y_label);
pixel_y = [1190,1039,919,835,715,631,565,361,241,157,37];
fit_y = fit(pixel_y',log_y_label','poly1');
pixel = @(px,py)[exp(fit_x(px)),exp(fit_y(py))];

FanSpeed = [790,870,1010,1350,1630,2890];
v = FanSpeed'/1000;
Efficiency = [70,75,78,81,78,75,70]/100;
u = [-sqrt(81-70), -sqrt(81-75),-sqrt(81-78),sqrt(81-81),sqrt(81-78),sqrt(81-75),sqrt(81-70)]';
Pixel = [242,889;305,902;370,923;447,976;492,1046;512,1099;521,1140;...
         270,833;333,845;398,866;476,920;521,989;540,1042;549,1082;...
         314,744;377,757;442,778;520,832;565,901;584,954;593,994;...
         400,574;463,586;528,607;605,660;650,731;670,783;679,823;...
         455,463;518,475;583,496;661,550;706,619;725,672;734,712;...
         625,123;688,136;753,157;831,210;876,279;895,334;904,372];
Axis_QP=pixel(Pixel(:,1),Pixel(:,2));
Q = reshape(Axis_QP(:,1),length(Efficiency),length(FanSpeed));
P = reshape(Axis_QP(:,2),length(Efficiency),length(FanSpeed));
[Efficiency_grid,FanSpeed_grid] = meshgrid(Efficiency,FanSpeed);

LnQ = log(Q);
LnP = log(P);

[U,V] = meshgrid(u,v);

ft = fittype( 'p00 + p10*x + p01*y + p20*x^2 + p02*y^2', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.383306318625529 0.617279232316449 0.575494859702814 0.530051704765016 0.275069755821935];
[xData, yData, zData] = prepareSurfaceData( u, v, LnQ );
[LnQ_fit, LnQ_gof] = fit( [xData, yData], zData, ft, opts );
[xData, yData, zData] = prepareSurfaceData( u, v, LnP );
[LnP_fit, LnP_gof] = fit( [xData, yData], zData, ft, opts );
LnQ_Coeff = coeffvalues(LnQ_fit);
LnP_Coeff = coeffvalues(LnP_fit);
Eff_max = 0.81;



Q = 0.8;
v = 1200/1000;
delta = LnQ_Coeff(4)^2-4*LnQ_Coeff(5)*(LnQ_Coeff(1)+LnQ_Coeff(2)*v+LnQ_Coeff(3)*v^2-log(Q));
u = (-LnQ_Coeff(4)+sign(delta)*sqrt(abs(delta)))/2/LnQ_Coeff(5);
P = exp(LnP_Coeff*[1;v;v^2;u;u^2])
Eff = Eff_max-u^2/100