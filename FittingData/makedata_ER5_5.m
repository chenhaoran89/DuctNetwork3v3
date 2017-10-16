FileName='ER5_5.txt';

Cb2 = zeros(11,4);
QbQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AbAc = [0.25 0.5 0.75 1];
Cb2 = Get2DTable(FileName, 6,9)';
Cb = Cb2.*(QbQc'.^2*AbAc.^-2);

save('ER5_5.mat','Cb','Cb2','QbQc','AbAc');