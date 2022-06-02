% function P_output (Population,time,Algorithm,Problem,M,Run,Replacesort,cosine,Theta,W,K1,K2)
function P_output (Population,FunctionValue,time,Algorithm,Problem,M,Run,MaxGen,K)
	eval(['save Data/',Algorithm,'/',Problem,'/',num2str(MaxGen),'_',num2str(K),'_',Algorithm,'_',Problem,'_',num2str(M),'_',num2str(Run),' Population FunctionValue'])
end