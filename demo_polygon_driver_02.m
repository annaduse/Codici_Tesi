function demo_polygon_driver_02 %versione per sigma rumore solo gaussiano

%--------------------------------------------------------------------------
% OBJECT.
% This function runs the experiments for hyperinterpolation over polygons.
% Several hyperinterpolants are use on possibly noisy functions.
%
% It makes the following numerical experiments:
% a) Fixed noise and several well-chosen lambda parameters (where useful).
% b) Fixed lambda and variable noise.
%
% At the end of the code
% 1. it produces the pertinent Latex files, making the relavant table.
% 2. it plots graphics comparing sparsity vs L2 errors on some hyp.variants;
%    all the latter data are saved on a file.
%--------------------------------------------------------------------------
% In "Hybrid hyperinterpolation over general regions", we consider:
%
% 1. f(x,y)=(1-x.^2-y.^2).*exp(x.*cos(y)).
% 2. Maximum hyperinterpolant degree: 15 (rule: ade=30).
% 3. Reference cubature rule degree of precision: 40.
% 4. Experiment 1. Noise: a=0, sigma=0.075, lambda=[10 20 40 80].
% 5. Experiment 2. Noise: a=0, sigma=[0.05 0.1 0.4 0.8], lambda=10.
% 6. Experiment 3. Noise: a=0, sigma=0.075, average Lambda: 4.404e-03.
% 7. Cubature points for hyperinterpolation: 1197; note: no rule compress.
% 8. Cubature points for testing L2 errors: 2065; note: no rule compress.
% 9. Total amount of hyperinterpolation coefficients: 136.
%
% a) This routine runs the experiments cited in the paper above.
% Polygon: non convex polygon with vertices:
%   vertices=(1/4)*[1 2; 1 0; 3 2; 3 0; 4 2; 3 3; 3 0.85*4; 2 4; 0 3; 1 2];
% b) The process took 207 seconds on a MacBook Pro, with Apple M1 processor
% and 16 GB of RAM (stage 3 may be time consuming).
%--------------------------------------------------------------------------
% Dates:
% Written on January 1, 2023: A. Sommariva.
% Modified on April 24, 2023: A. Sommariva.
% Joint work with Cong-Pei An, Jia-Shu Ran.
%--------------------------------------------------------------------------
% COPYRIGHT
%--------------------------------------------------------------------------
% Copyright (C) 2023-
%
% Authors:
% Alvise Sommariva
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------


clear;

% ................................SETTINGS ................................


% a) PART 1: fixed noise
LV=15; % Hyperinterpolation degree.
dimPS=(LV+1)*(LV+2)/2;
if dimPS < 10
    kV=[1 dimPS];
else
    kV=[];
    kk=1;
    while kk*10 <= dimPS
        kV(end+1)=kk*10;
        kk=kk*2;
    end
end

a1=0.0; sigma1=0.5;

% b) PART 2: fixed lambda
kS=10; % lambda is choosen equal to the "k"-th hyp. coeff. in magnitude.
a2=0.0; sigma2V=[0.05 0.1 0.4 0.8]; % noise parameters

% c) PART 3: sparsity plots
kP=2:2:dimPS;
a3=0.0; sigma3=0.5;

% ......................FIRST STAGE OF EXPERIMENTS   ......................

fprintf(2,'\n \n \n \t * <strong>Stage 1.</strong> \n \n \n')

AE2V=[]; betaV=[]; lambdaHIST=[];

XYW=[]; XYWR=[];

% note: cubature rule is computed only once.
for k=kV
    [AEinf,AE2,beta,lambdaV,XYW,XYWR]=demo_polygon(k,a1,sigma1,XYW,XYWR,...
        LV);
    AE2V=[AE2V AE2]; betaV=[betaV beta]; lambdaHIST=[lambdaHIST lambdaV];
end

AE2T=AE2V(1,:); BETAT=betaV(1,:); % tikhonov data storage
AE2F=AE2V(2,:); BETAF=betaV(2,:); % filtered hyp. data storage
AE2L=AE2V(3,:); BETAL=betaV(3,:); % lasso hyp. data storage
AE2H=AE2V(4,:); BETAH=betaV(4,:); % Hybrid hyp. data storage
AE2HA=AE2V(5,:); BETAHA=betaV(5,:); % Hard hyp. data storage
AE20=AE2V(6,:); BETA0=betaV(6,:); % Hyp. data storage


% Making tables
AE2V1=AE2V; betaV1=betaV; BETAH1=BETAH;
HypType=categorical({'tikhonov'; 'filtered'; 'lasso'; 'hybrid'; ...
    'hard'; 'hyperint.'; 'hyb.spars.'});
tablemat=[AE2V1; BETAH1];
T = table(HypType,tablemat); disp(T);

fprintf('\n \t * Lambda interval (for hyp. variants): [%1.2e,%1.2e] \n',...
    min(min(lambdaHIST(kV,:))),max(max(lambdaHIST(kV,:))));

for k=1:length(kV)
    kVL=kV(k);
    fprintf('\n \t * Lambda %1.0f column (average): %1.4e',...
        k,mean(lambdaHIST(kVL,:)));
end



% ...................... SECOND STAGE OF EXPERIMENTS  .....................

fprintf(2,'\n \n \n \t * <strong>Stage 2.</strong> \n ')
fprintf(2,'\n \n \t -> Average lambda: %1.3e \n \n',...
    mean(lambdaV(kS,:)));

% 2. Making experiments
AE2V=[]; betaV=[];

for sigma2=sigma2V
    % fprintf('\n \t sigma: %3.5e',sigma)
    [AEinf,AE2,beta,lambdaV,XYW,XYWR]=demo_polygon(kS,a2,sigma2,XYW,XYWR);
    AE2V=[AE2V AE2]; betaV=[betaV beta];
end

% 3. Making tables
AE2V2=AE2V; betaV2=betaV; BETAH2=betaV(4,:);
HypType=categorical({'tikhonov'; 'filtered'; 'lasso'; 'hybrid'; 'hard'; ...
    'hyperint.'; 'hyb.spars.'});
tablemat=[AE2V2; BETAH2];
T = table(HypType,tablemat); disp(T)



% ....................... PRODUCE LATEX TABLE  ............................

fprintf(2,'\n \n \t -> LaTeX table \n \n');

fileID=fopen('results_latex_polygon_stage1.txt','w');

strC={'Tikhonov'; 'Filtered'; 'Lasso'; 'Hybrid'; 'Hard'; 'Hyperint.'; ...
    '$\|{\mathbf{\beta}}\|_{0}$'};
for k=1:7
    strL=strC{k};

    strNUM='';

    for s=1:size(AE2V1,2)

        if k< 7
            strNUML=['$',num2str(AE2V1(k,s),'%.4f'),'$ &'];
        else
            strNUML=[ '$',num2str(BETAH1(s),'%.1f'),'$ &'];
        end

        strNUM=[strNUM strNUML];

    end


    for s=1:size(AE2V2,2)

        if k< 7
            strNUML=['$',num2str(AE2V2(k,s),'%.4f'),'$ &'];
        else
            strNUML=[ '$',num2str(BETAH2(s),'%.1f'),'$ &'];
        end

        strNUM=[strNUM strNUML];

    end

    if k ==7
        fprintf('\n \t'); disp('\hline');
        fprintf(fileID,'\n \t \\hline');
    end

    strNUM=[extractBefore(strNUM,length(strNUM)),'\\'];

    strdisp=['{\em{',strL,'}} &',strNUM];
    fprintf('\n \t '); disp(strdisp);

    strdispfile=replace(strdisp,'\','\\');
    fprintf(fileID,strcat('\n \t',' ',strdispfile));
end

fclose(fileID);





% ....................... THIRD STAGE OF EXPERIMENTS  .....................

fprintf(2,'\n \n \n \t * <strong>Stage 3.</strong> \n')

AE2P=[]; betaP=[]; JzHIST=[]; HzHIST=[];

% Making experiments.
for k=kP
    % fprintf('\n \t k: %3.0f',k)
    [AEinf,AE2,beta,lambdaV,XYW,XYWR,JzL,HzL,Wg_norm2]=demo_polygon(k,...
        a3,sigma3,XYW,XYWR);
    AE2P=[AE2P AE2]; betaP=[betaP beta];
    JzHIST=[JzHIST JzL]; HzHIST=[HzHIST HzL];
end


% .................... PLOT FIGURE SPARSITY vs L2 ERROR  ..................

fprintf(2,'\n \n \t -> Plotting sparsity/L2 errors figure.')

fig=figure(1);
plot(betaP(3,:),JzHIST(3,:),'<','linewidth',1,'markersize',6,'color','r')
hold on
plot(betaP(3,:),HzHIST(3,:),'-<','linewidth',1,'markersize',6,'MarkerFaceColor','r','color','r')
plot(betaP(4,:),JzHIST(4,:),'s','linewidth',1,'markersize',6,'color','k')
plot(betaP(4,:),HzHIST(4,:),'-s','linewidth',1,'markersize',6,'MarkerFaceColor','k','color','k')
%plot(betaP(5,:),JzHIST(5,:),'d','linewidth',1,'markersize',6,'color','b')
%plot(betaP(5,:),HzHIST(5,:),'-d','linewidth',1,'markersize',6,'MarkerFaceColor','b','color','b')
legend('Lasso (J)','Lasso (H)','Hybrid (J)','Hybrid (H)');
xlabel('Sparsity (a)');
ylabel('Values');
grid on;
yscale_symlog;
hold off

dtm=datetime; 
str=datestr(dtm); str=replace(str,' ','_'); str=replace(str,':','');
str=replace(str,'-','');

figname=['figure1',str];
print(fig,figname,'-depsc');
savefig([figname,'.fig']);



fig=figure(2);
semilogy(betaP(3,:),AE2P(3,:),'r<','linewidth',2,'markersize',6)
hold on
semilogy(betaP(4,:),AE2P(4,:),'ks','linewidth',2,'markersize',6)
%semilogy(betaP(5,:),AE2P(5,:),'bd','linewidth',2,'markersize',6)
legend('Lasso', 'Hybrid');
xlabel('Sparsity (a)');
ylabel('L_2 errors');
grid on;
hold off

figname=['figure2',str];
print(fig,figname,'-depsc');
savefig([figname,'.fig']);















