FileName='SD5_9.txt';

Cb2 = zeros(11,11);
QbQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AbAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
Cb2 = Get2DTable(FileName, 6,16)';
Cb = Cb2.*(QbQc'.^2*AbAc.^-2);

Cs2 = zeros(11,11);
QsQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AsAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
Cs2 = Get2DTable(FileName, 23,33)';
Cs = Cs2.*(QsQc'.^2*AsAc.^-2);

save('SD5_9.mat','Cb','Cb2','QbQc','AbAc','Cs','Cs2','QsQc','AsAc');