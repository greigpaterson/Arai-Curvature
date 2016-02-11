function Par = LMA(XY,ParIni)

%--------------------------------------------------------------------------
%  
%     Geometric circle fit (minimizing orthogonal distances)
%     based on the Levenberg-Marquardt scheme in the 
%     "algebraic parameters" A,B,C,D  with constraint B*B+C*C-4*A*D=1
%        N. Chernov and C. Lesort, "Least squares fitting of circles",
%        J. Math. Imag. Vision, Vol. 23, 239-251 (2005)
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%             ParIni = [a b R] is the initial guess (supplied by user)
%
%     Output: Par = [a b R] is the fitting circle:
%                           center (a,b) and radius R
%
%--------------------------------------------------------------------------

factorUp=10;  factorDown=0.04;
lambda0=0.01;  epsilon=0.000001;  
IterMAX = 50;  AdjustMax = 20;
Xshift=0;  Yshift=0;  dX=1;  dY=0;

n = size(XY,1);      % number of data points

%     starting with the given initial guess

anew = ParIni(1) + Xshift;
bnew = ParIni(2) + Yshift;

Anew = 1/(2*ParIni(3));
aabb = anew*anew + bnew*bnew;
Fnew = (aabb - ParIni(3)*ParIni(3))*Anew;
Tnew = acos(-anew/sqrt(aabb));
if (bnew > 0) 
    Tnew = 2*pi - Tnew;
end
%        if 1+4*Anew*Fnew < epsilon
%            fprintf(1,' +++ violation:  %f\n',1+4*Anew*Fnew);
%        end

VarNew = VarCircle(XY,ParIni);

%     initializing lambda and iter

lambda = lambda0;  finish = 0;

for iter=1:IterMAX

    Aold = Anew;  Fold = Fnew;  Told = Tnew;  VarOld = VarNew;

    H = sqrt(1+4*Aold*Fold);
    aold = -H*cos(Told)/(Aold+Aold) - Xshift;
    bold = -H*sin(Told)/(Aold+Aold) - Yshift;
    Rold = 1/abs(Aold+Aold);
%    fprintf(1,'%2d  (%f, %f)  %f  %.8f\n',iter,aold,bold,Rold,sqrt(VarOld));

%           computing moments

    DD = 1 + 4*Aold*Fold;
    D = sqrt(DD);
    CT = cos(Told);
    ST = sin(Told);

    H11=0; H12=0; H13=0; H22=0; H23=0; H33=0; F1=0; F2=0; F3=0;
    
    for i=1:n
        Xi = XY(i,1) + Xshift;
        Yi = XY(i,2) + Yshift;
        Zi = Xi*Xi + Yi*Yi;
        Ui = Xi*CT + Yi*ST;
        Vi =-Xi*ST + Yi*CT;

        ADF = Aold*Zi + D*Ui + Fold;
        SQ = sqrt(4*Aold*ADF + 1);
        DEN = SQ + 1;
        Gi = 2*ADF/DEN;
        FACT = 2/DEN*(1 - Aold*Gi/SQ);
        DGDAi = FACT*(Zi + 2*Fold*Ui/D) - Gi*Gi/SQ;
        DGDFi = FACT*(2*Aold*Ui/D + 1);
        DGDTi = FACT*D*Vi;

        H11 = H11 + DGDAi*DGDAi;
        H12 = H12 + DGDAi*DGDFi;
        H13 = H13 + DGDAi*DGDTi;
        H22 = H22 + DGDFi*DGDFi;
        H23 = H23 + DGDFi*DGDTi;
        H33 = H33 + DGDTi*DGDTi;

        F1 = F1 + Gi*DGDAi;
        F2 = F2 + Gi*DGDFi;
        F3 = F3 + Gi*DGDTi;
    end

    for adjust=1:AdjustMax

%             Cholesly decomposition

        G11 = sqrt(H11 + lambda);
        G12 = H12/G11;
        G13 = H13/G11;
        G22 = sqrt(H22 + lambda - G12*G12);
        G23 = (H23 - G12*G13)/G22;
        G33 = sqrt(H33 + lambda - G13*G13 - G23*G23);

        D1 = F1/G11;
        D2 = (F2 - G12*D1)/G22;
        D3 = (F3 - G13*D1 - G23*D2)/G33;

        dT = D3/G33;
        dF = (D2 - G23*dT)/G22;
        dA = (D1 - G12*dF - G13*dT)/G11;

%            updating the parameters

        Anew = Aold - dA;
        Fnew = Fold - dF;
        Tnew = Told - dT;
%        fprintf(1,'%2d   %.8f\n',iter,lambda);    

        if (1+4*Anew*Fnew < epsilon && lambda>1)
%            fprintf(1,'     violation:  %f\n',1+4*Anew*Fnew);
            Xshift = Xshift + dX;
            Yshift = Yshift + dY;

            H = sqrt(1+4*Aold*Fold);
            aTemp = -H*cos(Told)/(Aold+Aold) + dX;
            bTemp = -H*sin(Told)/(Aold+Aold) + dY;
            rTemp = 1/abs(Aold+Aold);

            Anew = 1/(rTemp + rTemp);
            aabb = aTemp*aTemp + bTemp*bTemp;
            Fnew = (aabb - rTemp*rTemp)*Anew;
            Tnew = acos(-aTemp/sqrt(aabb));
            if bTemp > 0 
               Tnew = 2*pi - Tnew;
            end
            VarNew = VarOld;
            break;
        end

        if 1+4*Anew*Fnew < epsilon
           lambda = lambda * factorUp;
           continue;
        end

        DD = 1 + 4*Anew*Fnew;
        D = sqrt(DD);
        CT = cos(Tnew);
        ST = sin(Tnew);

        GG = 0;

        for i=1:n
            Xi = XY(i,1) + Xshift;
            Yi = XY(i,2) + Yshift;
            Zi = Xi*Xi + Yi*Yi;
            Ui = Xi*CT + Yi*ST;

            ADF = Anew*Zi + D*Ui + Fnew;
            SQ = sqrt(4*Anew*ADF + 1);
            DEN = SQ + 1;
            Gi = 2*ADF/DEN;
            GG = GG + Gi*Gi;
        end
    
        VarNew = GG/(n-3);

        H = sqrt(1+4*Anew*Fnew);
        anew = -H*cos(Tnew)/(Anew+Anew) - Xshift;
        bnew = -H*sin(Tnew)/(Anew+Anew) - Yshift;
        Rnew = 1/abs(Anew+Anew);
       
%             checking if improvement is gained

        if VarNew <= VarOld      %   yes, improvement
           progress = (abs(anew-aold) + abs(bnew-bold) + abs(Rnew-Rold))/(Rnew+Rold);
           if progress < epsilon
              Aold = Anew;
              Fold = Fnew;
              Told = Tnew;
              VarOld = VarNew;
              finish = 1;
              break;
           end
           lambda = lambda * factorDown;
           break;
        else                     %   no improvement
           lambda = lambda * factorUp;
           continue;
        end
    end
    if finish == 1
       break;
    end
end

H = sqrt(1+4*Aold*Fold);
Par(1) = -H*cos(Told)/(Aold+Aold) - Xshift;
Par(2) = -H*sin(Told)/(Aold+Aold) - Yshift;
Par(3) = 1/abs(Aold+Aold);

end  % LMA

function Var = VarCircle(XY,Par)

%--------------------------------------------------------------------------
%  
%             computing the sample variance of distances 
%             from data points (XY) to the circle Par = [a b R]
%
%--------------------------------------------------------------------------
 

n = size(XY,1);      % number of data points
Dx = XY(:,1) - Par(1);  Dy = XY(:,2) - Par(2);
D = sqrt(Dx.*Dx + Dy.*Dy) - Par(3);

Var = D'*D/(n-3);

end  %  VarCircle
