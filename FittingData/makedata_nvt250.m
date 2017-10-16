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
v = FanSpeed';
Efficiency = [70,75,78,81,78,75,70]'/100;
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

[V,U] = meshgrid(v,u);

LnQ_V = mean(log(Q./V),2);
LnP_VV = mean(log(P./V.^2),2);

figure(3)
plot(LnQ_V,LnP_VV)

plot((2*LnQ_V-LnP_VV),u)

% ft = fittype( 'p00 + p10*x + p01*y + p20*x^2 + p02*y^2', 'independent', {'x', 'y'}, 'dependent', 'z' );
ft = fittype( 'k/(x+b)+h', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% opts.StartPoint = [0.383306318625529 0.617279232316449 0.575494859702814 0.530051704765016 0.275069755821935];
opts.StartPoint = [1,1,1];
[Model_fit, Model_gof] = fit( LnQ_V, LnP_VV, ft, opts );
Model_Coeff = coeffvalues(Model_fit);

Eff_max = 0.81;
Eff_min = min(Efficiency);

save('nvt250.mat','Model_Coeff','LnQ_V','Efficiency');
%%
% load('nvt250.mat');
q = linspace(-1,1,100000);
v = 500:100:2000;
[V,Q] = meshgrid(v,q);
fan_1 = @(q,v,s) v.^2*exp(s(2)); %q<=0
fan_2 = @(q,v,s) v.^2.*exp(s(3)./(log(q./v)+s(1))+s(2)); %0<q<v*exp(-s(1))
% p = 0 for q>=v*exp(-s(1))
fan = @(q,v,s)full(sparse(find(q<=0),1,fan_1(q(q<=0),v(q<=0),s),length(q),1)+sparse(find(q>0 & q<v*exp(-s(1))),1,fan_2(q(q>0 & q<v*exp(-s(1))),v(q>0 & q<v*exp(-s(1))),s),length(q),1));
P = fan(Q(:),V(:),Model_Coeff);
P = reshape(P,size(Q));
plot(Q,P);
% Q = 0.8;
% v = 1200/1000;
% delta = LnQ_Coeff(4)^2-4*LnQ_Coeff(5)*(LnQ_Coeff(1)+LnQ_Coeff(2)*v+LnQ_Coeff(3)*v^2-log(Q));
% u = (-LnQ_Coeff(4)+sign(delta)*sqrt(abs(delta)))/2/LnQ_Coeff(5);
% P = exp(LnP_Coeff*[1;v;v^2;u;u^2])
% Eff = Eff_max-u^2/100