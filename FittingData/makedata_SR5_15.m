FileName='SR5_15.txt';

Cb2 = zeros(11,11);
QbQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AbAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
Cb2 = Get2DTable(FileName, 6,16)';
Cb = Cb2.*(QbQc'.^2*AbAc.^-2);

save('SR5_15.mat','Cb','Cb2','QbQc','AbAc');