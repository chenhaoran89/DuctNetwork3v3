FileName='ED5_3.txt';

Cb_part1 = zeros(11,11,11);
Cs_part1 = zeros(11,11,11);
Cb_part2 = zeros(11,11,11);
Cs_part2 = zeros(11,11,11);
Cb2_part1 = zeros(11,11,11);
Cs2_part1 = zeros(11,11,11);
Cb2_part2 = zeros(11,11,11);
Cs2_part2 = zeros(11,11,11);

QbQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
QsQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AbAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AsAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];

for ii=1:11
    Cb2_part1(:,:,ii) = Get3DTable(FileName, (ii-1)*11+6,(ii-1)*11+16)';
    Cb_part1(:,:,ii) = Cb2_part1(:,:,ii).*(QbQc'.^2*AbAc.^-2);
    Cs2_part1(:,:,ii) = Get3DTable(FileName, (ii-1)*11+133,(ii-1)*11+143)';
    Cs_part1(:,:,ii) = Cs2_part1(:,:,ii).*(QsQc'.^2*(AsAc(ii)*ones(1,11)).^-2);
    Cb2_part2(:,:,ii) = Get3DTable(FileName, (ii-1)*11+260,(ii-1)*11+270)';
    Cb_part2(:,:,ii) = Cb2_part2(:,:,ii).*(QbQc'.^2*AbAc.^-2);
    Cs2_part2(:,:,ii) = Get3DTable(FileName, (ii-1)*11+387,(ii-1)*11+397)';
    Cs_part2(:,:,ii) = Cs2_part2(:,:,ii).*(QsQc'.^2*(AsAc(ii)*ones(1,11)).^-2);
end

save('ED5_3.mat','Cb_part1','Cb2_part1','Cs_part1','Cs2_part1','Cb_part2','Cb2_part2','Cs_part2','Cs2_part2','QbQc','AbAc','QsQc','AsAc');