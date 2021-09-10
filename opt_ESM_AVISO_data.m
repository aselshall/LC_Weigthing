function [LCA,KB,row,ZOSN,ZOSS,NR]=opt_ESM_AVISO_data(file,folder)

%[0] Model data
file=[folder '/' file];
load([file '.mat'],'ZOSN','ZOSS','LCA','Member','KB','row')
%fprintf('KB: %s | Ensemble: %s| Sections: %s-%s | model: %s | ecKBM: %i | Err_LCS: %f | Err_Tot: %f | Err_AM_Ratio: %f\n',...
%    row.KB,row.Ensemble,row.N,row.S,row.model,row.ecKBM,row.Err_LCS, row.Err_Tot, row.Err_Ratio)

%[1] Select ensemble
%3XXX: NR=[6 3 4 7];
%3210: NR=1 1 1 3 1 5 3 3 6 3 1 2 1 3 3 4];
[NR, ~] = hist(Member(:,2),unique(Member(:,2))); %#ok<HIST>


%[2] ZOS all members
ZOSN=permute(ZOSN, [2 3 1]);             %Convert numpy to MATLAB
ZOSS=permute(ZOSS, [2 3 1]);             %Convert numpy to MATLAB

end