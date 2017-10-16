function [dP, dPdQ, dPdS]=ER5_3(q, s, Selection)
% q = [Qs,Qb,Qc];
% s = [Hs,Ws,Hb,Wb,rho];
ModifiedTable = 0;
if (ModifiedTable == 1)
    gExp =  0.5*s(5)*(q(3)/(s(1)*s(2)))^2;
    dgdq =  [0,0,s(5)*q(3)/(s(1)*s(2))^2];
    dgds =  [-2/s(1),-2/s(2),0,0,1/s(5)]*gExp;
    switch Selection
        case 'b'
            GridVec = {DuctNetwork.Table_ER5_3.QbQc};
            ZExp = [q(2)/q(3)];
            dZdq = [0,1/q(3),-q(2)/q(3)^2];
            dZds = [0,0,0,0,0];
            [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ER5_3.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds);
        case 's'
            GridVec = {DuctNetwork.Table_ER5_3.QsQc};
            ZExp = [q(1)/q(3)];
            dZdq = [1/q(3),0,-q(1)/q(3)^2];
            dZds = [0,0,0,0,0];
            [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ER5_3.Cs, ZExp, dZdq, dZds, gExp, dgdq, dgds);
        otherwise
            dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
    end
else
    switch Selection
        case 'b'
            gExp =  0.5*s(5)*(q(2)/(s(3)*s(4)))^2;
            dgdq =  [0,0,s(5)*q(2)/(s(3)*s(4))^2];
            dgds =  [0,0,-2/s(3),-2/s(4),1/s(5)]*gExp;
            GridVec = {DuctNetwork.Table_ER5_3.QbQc};
            ZExp = [q(2)/q(3)];
            dZdq = [0,1/q(3),-q(2)/q(3)^2];
            dZds = [0,0,0,0,0];
            [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ER5_3.Cb2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
        case 's'
            gExp =  0.5*s(5)*(q(1)/(s(1)*s(2)))^2;
            dgdq =  [0,0,s(5)*q(3)/(s(1)*s(2))^2];
            dgds =  [-2/s(1),-2/s(2),0,0,1/s(5)]*gExp;
            GridVec = {DuctNetwork.Table_ER5_3.QsQc};
            ZExp = [q(1)/q(3)];
            dZdq = [1/q(3),0,-q(1)/q(3)^2];
            dZds = [0,0,0,0,0];
            [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, DuctNetwork.Table_ER5_3.Cs2, ZExp, dZdq, dZds, gExp, dgdq, dgds);
        otherwise
            dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0,0];
    end
end
end
