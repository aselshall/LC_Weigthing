function [Err,Err1,OF,mw]=a2Case04(x,AD,ZOSN,ZOSS,NR,KB,nm,row,rangemin,rangemax,PrMsMin,PrMsMax,SOF)  %Comment for testing only
% 27/March/2019
%clc; clear *; close all;    %Uncomment for postprocessing 
%This function calculates the objective function (OF)]
% [5] Process decision variables 
% [6] Process information per node to link decision variables to SUTRA withdrawal nodes 
% [7] Wrting withdrawal and recharge as input to SUTRA (rech.trans file)
% [8] Run SUTRA
% [9] Read Concentration observation package 'ModelName.obs'
% [10] Match SUTRA observations with pumping well nodes and spring nodes
% [11] PW(Drawdown and salinity) and Spring(Discharge and Salinity) of optimal solution 


%% [5] Processing decision variables 

%(5.2) Define prior (uniform paramter distrubutions) 
OptRange=[rangemin;rangemax];                      %Normal distribution range
NDV=length(x);                         %Number of decision variable (#PW or #PumpingClusters)
ParaRange=repmat([PrMsMin;PrMsMax],1,NDV)';        %Parameter range

%(5.3)Scaling CMA-ES parameter space to SUTRA parameter space
UN=size(x,1);   %Number of unknowns
xi=NaN(UN,1);   %Model parameters
x=real(x);

for m=1:length(x)
    xi(m,1) = interp1(OptRange,ParaRange(m,:)',x(m,1),'linear','extrap');
end

%(5.6) Boundary handling (optimization mode only)
BH=1;   %(1) Reject out-of-bound, (2) Rebound to zero
%Reject values outside the lower and upper endpoints of the uniform distribution
if BH==1
    if min(x(:))<rangemin || max(x(:))>rangemax     %Outside  parameter range
        OF=1e6*max(abs(x));%Huge penalty
        Err=100;
        Err1=100;
        mw=ones(1,UN);
        return
    end

elseif BH==2
    if min(x(:))<rangemin %Outside  parameter range lower bound
        xi(xi<0)=0;
    end
end  %End BH


%% [8] Run model

%(8.1) Weight mean (ZOS) 
NR=NR';
mw=(xi/sum(xi));        %Model weight per single-model ensemble
mwm=repelem(mw./NR,NR); %Model weight per ensemble member
%mwm=repmat(1/20,20,1);  %Equal weight for all realizations for testing purpose only
WMN=sum(ZOSN.*permute(repmat(mwm,[1,size(ZOSN,2),size(ZOSN,1)]),[3 2 1]),3);  %Weighted mean ZOS
WMS=sum(ZOSS.*permute(repmat(mwm,[1,size(ZOSS,2),size(ZOSS,1)]),[3 2 1]),3);  %Weighted mean ZOS

%(8.2) Mean segement(ZOS)
DMN=mean(WMN,2,'omitnan');
DMS=mean(WMS,2,'omitnan');

%(8.3) Delta ZOS_mean_Segement(North, South)
DM=DMN-DMS;

%(8.4) LC per period
MD=max(reshape(DM,nm,[]),[],1,'omitnan')';

%% [9] Predictors
%           1       2   3   4   5       6       7       8       9       10          11          12          13      14      15      16             
%           e_size	KB	LCN	LCS	LCN_NB	LCN_B	LCS_NB	LCS_B	Err_KB	Match_LCN	Match_LCS	Match_Tot	Err_LCN	Err_LCS	Err_Tot	RMSE
% obs	    1       15	32	12	 17	    15      12      0       0           32          12          44      0       0       0       0
% 3210	    41      15	3	41	 2      1       27      14      0.933       2           11          13      0.938	0.083	0.705	5.13
% 321X	    33      15	34	10	 22     12      7       3       0.2         25          3           28      0.219	0.75	0.364	3.71
% 32XX	    28      15	23	21	 17     6       12      9       0.6         17          6           23      0.469	0.5     0.477	3.92
% 3XXX	    20      15	35	9	 22     13      7       2       0.133       28          5           33      0.125	0.583	0.25	3.68
% XXX0	    8       15	0	44	 0      0       29      15      1           0           12          12      1       0       0.727	13.52

%(9.0) General
AD=AD';
KB=KB';
e_size=size(ZOSN,3);    %Table(1)
KBB=sum(KB);            %Table(2)

%(9.2) LC-N / LC-S count 
LCN(:,1)=(AD>=0);
LCS(:,1)=(AD<0);
LCN(:,2)=(MD>=0);
LCS(:,2)=(MD<0);
LCN(:,3)=(AD>0 & MD>=0);
LCS(:,3)=(AD<0 & MD<0);
LCNM=sum(LCN(:,2));     %Table(3)
LCSM=sum(LCS(:,2));     %Table(4) 
LCNA=sum(LCN(:,1));     %Table(3)_data
LCSA=sum(LCS(:,1));     %Table(4)_data

%(9.3) LC-N / LC-S count with bloom  
LCN_NB=sum((MD>=0 & KB==0));	%Table(5)
LCN_B=sum((MD>=0 & KB==1));     %Table(6)
LCS_NB=sum((MD<0 & KB==0));     %Table(7)
LCS_B=sum((MD<0 & KB==1));      %Table(8)
LCS_B_A=sum((AD<0 & KB==1));    %Table(8)_data

%(9.4) KB error 
Err_KB=LCS_B/KBB;               %Table(9)

%(9.5) LC match 
Match_LCN=sum(LCN(:,3));        %Table(10)
Match_LCS=sum(LCS(:,3));        %Table(11)
Match_Tot=Match_LCN+Match_LCS;  %Table(12)

%(9.6) LC Error 
Err_LCN=((LCNA-Match_LCN)/LCNA);    %Table(13)
Err_LCS=((LCSA-Match_LCS)/LCSA);    %Table(14)
Err_Tot= (length(AD) - sum((AD<0)==(MD<0))) /length(AD);  %Table(15)

%(9.7) RMSE
%resm.loc[member,'RMSE']=np.round(np.sqrt(np.mean(np.square(LC-LCO)))*1e2,decimals=2)
RMSE=sqrt(mean(((MD-AD).^2)))*100;            %RMSE (cm) Table(16)
rms=mean(((MD-AD).^2))*100;

%(9.8) Ratio constrain
LCSA_Ratio=sum(LCS(:,1))/length(AD);
LCSM_Ratio=sum(LCS(:,2))/length(MD);
LCS_Ratio_Delta=abs(LCSA_Ratio-LCSM_Ratio);



%% [10] Objective function

%(10.1) LC Ratio constrain
if LCSM_Ratio<0.2 || LCSM_Ratio>0.35 
    OF=1e15;
    Err=100;
    Err1=100;
    return
end


%(10.2) Ojbective function
if length(NR)>4
    Err_TotT=0.2;
else
    Err_TotT=0.23;
end

if SOF==1
    OF=(LCS_B^5)*Err_LCS;
elseif SOF==2
    OF=(LCS_B^5)*(Err_Tot+1)*(rms+1);
elseif SOF==3
    OF=(LCS_B^5)*((LCS_Ratio_Delta*10)+1)*(rms+1);
elseif SOF==4
    OF=(LCS_B^5)*((LCS_Ratio_Delta*10)+1)*(rms+1)*(Err_Tot+1);
elseif SOF==5
    OF=((LCS_B^5))*Err_LCS*((LCS_Ratio_Delta*10)+1)*(rms+1);
elseif SOF==6
    OF=(LCS_B^5)*((Err_LCS+1)^4)*((LCS_Ratio_Delta*10)+1)*(rms+1)*((Err_Tot+1)^2);
elseif SOF==7
    OF=(LCS_B+1)*(Err_LCS+1)*(LCS_Ratio_Delta+1)*(RMSE+1)*(Err_Tot+1);
elseif SOF==8
    if LCS_B>2  || Err_LCS>0.5  || Err_Tot>Err_TotT
        PEN=10;
    else
        PEN=1;
    end
    OF=(LCS_B+1)*((Err_LCS+1)^10)*(LCS_Ratio_Delta+1)*(Err_Tot+1)*RMSE*PEN;   
end

Err=Err_LCS;
Err1=LCS_B;



%(10.3) Display and save good results
if LCS_B>2  || Err_LCS>0.5  || Err_Tot>=Err_TotT
    return
end

%File for saving
Run=row.Ensemble(~isspace(row.Ensemble));
sfile=['opt_res_R' num2str(Run) '.mat'];
if ~ isfile(sfile)
    count=0;
    SO_Info=[];SO_MW=[]; SO_MD=[]; SO_AD=AD; RMSE_MAX=100;
    save(sfile,'count','SO_Info','SO_MW','SO_MD','SO_AD','RMSE_MAX')
end

%Display and save results
load(sfile,'count','SO_Info','SO_MW','SO_MD','RMSE_MAX')
if RMSE<RMSE_MAX
    count=count+1;
    SO_Info(count,:)=[OF e_size KBB LCNM LCSM LCN_NB LCN_B LCS_NB LCS_B Err_KB Match_LCN Match_LCS Match_Tot Err_LCN Err_LCS Err_Tot RMSE];
    SO_MW(count,:)=mw';
    SO_MD(count,:)=MD';
    
    %Display results
    disp([num2str(count) ': ' ...
        'LCN(Data,Model,Match,Err): ' num2str(LCNA) ' ' num2str(LCNM) ' ' num2str(Match_LCN) ' ' num2str(Err_LCN,3) ...
        ' LCS(Data,Model,Match,Err): ' num2str(LCSA) ' ' num2str(LCSM) ' ' num2str(Match_LCS) ' ' num2str(Err_LCS,3)...
        ' LCS_B(Data,Model): ' num2str(LCS_B_A) ' ' num2str(LCS_B) ' Err_Tot: '  num2str(Err_Tot) ' RMSE : ' num2str(RMSE)])
    % %Display model weights
    % disp([num2str(count) ': ' num2str(mw')])
    RMSE_MAX=RMSE;
    save(sfile,'count','SO_Info','SO_MW','SO_MD','RMSE_MAX')
end
close all






end

