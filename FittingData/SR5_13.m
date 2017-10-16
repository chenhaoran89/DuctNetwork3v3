function [dP, dPdQ, dPdS]=SR5_13(q, s, Selection)
% q = [Qs,Qb,Qc];
% s = [Hs,Ws,Hb,Wb,Hc,Wc,rho]
ModifiedTable = 0;
if (ModifiedTable == 1)
    gExp = 0.5*s(7)*(q(3)/(s(5)*s(6)))^2;
    dgdq = [0,0,s(7)*q(3)/(s(5)*s(6))^2];
    dgds = [0,0,0,0,-2/s(5),-2/s(6),1/s(7)]*gExp;
    switch Selection
        case 'b'
            GridVec = {DuctNetwork.Table_SR5_13.QbQc,DuctNetwork.Table_SR5_13.AbAc};
            ZExp = [q(2)/q(3);(s(3)*s(4))/(s(5)*s(6))];
            dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0];
            dZds = [0,0,0,0,0,0,0;0,0,1/s(3),1/s(4),-1/s(5),-1/s(6),0]*(s(3)*s(4))/(s(5)*s(6));
            [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_SR5_13.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
        case 's'
            GridVec = {DuctNetwork.Table_SR5_13.QsQc,DuctNetwork.Table_SR5_13.AsAc};
            ZExp =[q(1)/q(3);(s(1)*s(2))/(s(5)*s(6))];
            dZdq =[1/q(3),0,-q(1)/q(3)^2;0,0,0];
            dZds =[0,0,0,0,0,0,0;1/s(1),1/s(2),0,0,-1/s(5),-1/s(6),0]*(s(1)*s(2))/(s(5)*s(6));
            [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_SR5_13.Cs, ZExp, dZdq, dZds, gExp, dgdq, dgds);
        otherwise
            dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0,0,0];
    end
else
    switch Selection
        case 'b'
            gExp = 0.5*s(7)*(q(2)/(s(3)*s(4)))^2;
            dgdq = [0,0,s(7)*q(2)/(s(3)*s(4))^2];
            dgds = [0,0,-2/s(3),-2/s(4),0,0,1/s(7)]*gExp;
            GridVec = {DuctNetwork.Table_SR5_13.QbQc,DuctNetwork.Table_SR5_13.AbAc};
            ZExp = [q(2)/q(3);(s(3)*s(4))/(s(5)*s(6))];
            dZdq = [0,1/q(3),-q(2)/q(3)^2;0,0,0];
            dZds = [0,0,0,0,0,0,0;0,0,1/s(3),1/s(4),-1/s(5),-1/s(6),0]*(s(3)*s(4))/(s(5)*s(6));
            [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_SR5_13.Cb2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
        case 's'
            gExp = 0.5*s(7)*(q(1)/(s(1)*s(2)))^2;
            dgdq = [0,0,s(7)*q(1)/(s(1)*s(2))^2];
            dgds = [-2/s(1),-2/s(2),0,0,0,0,1/s(7)]*gExp;
            GridVec = {DuctNetwork.Table_SR5_13.QsQc,DuctNetwork.Table_SR5_13.AsAc};
            ZExp =[q(1)/q(3);(s(1)*s(2))/(s(5)*s(6))];
            dZdq =[1/q(3),0,-q(1)/q(3)^2;0,0,0];
            dZds =[0,0,0,0,0,0,0;1/s(1),1/s(2),0,0,-1/s(5),-1/s(6),0]*(s(1)*s(2))/(s(5)*s(6));
            [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_SR5_13.Cs2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
        otherwise
            dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0,0,0];
    end
end
end
