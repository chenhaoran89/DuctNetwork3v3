FileName='ER5_1.txt';

Cb = zeros(11,3,3);
Cb2 = zeros(11,3,3);
Cs = zeros(11,3,3);
Cs2 = zeros(11,3,3);

QbQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
QsQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AbAc = [0.25      0.5      1.0];
AsAc = [0.5      0.75      1.0];

for ii=1:3
    Cb2(:,:,ii) = Get3DTable(FileName, (ii-1)*3+6,(ii-1)*3+8)';
    Cb(:,:,ii) = Cb2(:,:,ii).*(QbQc'.^2*AbAc.^-2);
    Cs2(:,:,ii) = Get3DTable(FileName, (ii-1)*3+21,(ii-1)*3+23)';
    Cs(:,:,ii) = Cs2(:,:,ii).*(QsQc'.^2*(AsAc(ii)*ones(1,3)).^-2);
end

save('ER5_1.mat','Cb','Cb2','QbQc','AbAc','Cs','Cs2','QsQc','AsAc');