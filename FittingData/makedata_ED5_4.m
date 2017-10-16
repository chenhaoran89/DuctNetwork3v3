FileName='ED5_4.txt';

Cb1 = zeros(11,11,11);
Cb1_o = zeros(11,11,11);
Cb2 = zeros(11,11,11);
Cb2_o = zeros(11,11,11);

QbQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
QsQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AbAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AsAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];

for ii=1:11
    Cb1_o(:,:,ii) = Get3DTable(FileName, (ii-1)*11+6,(ii-1)*11+16)';
    Cb1(:,:,ii) = Cb1_o(:,:,ii).*(QbQc'.^2*AbAc.^-2);
    Cb2_o(:,:,ii) = Get3DTable(FileName, (ii-1)*11+133,(ii-1)*11+143)';
    Cb2(:,:,ii) = Cb2_o(:,:,ii).*(QsQc'.^2*(AsAc(ii)*ones(1,11)).^-2);
end

save('ED5_4.mat','Cb1','Cb1_o','QbQc','AbAc','Cb2','Cb2_o','QsQc','AsAc');