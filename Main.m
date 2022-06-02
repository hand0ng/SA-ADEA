function Main()
clear;format short;addpath public;addpath dace;addpath Data;addpath draw;
RunNum = 30;
MaxGen = 20%:5:40;
K = 5%:2:11; %K-Medoide
Problems = {'DTLZ1','DTLZ2','DTLZ3','DTLZ4'};
Objectives = [3,5,8,10]; % number of objectives

for i = 1:length(MaxGen)
    for j = 1:length(K)
        disp(['Para:',num2str(MaxGen(i)),'-',num2str(K(j))]);
        for Prob = 1:length(Problems)
            for Obj = 1:length(Objectives)
                Algorithm = {'SAADEA'};
                Problem   = Problems{Prob};
                M         = Objectives(Obj);
                for Run = 1:RunNum
                    Start(Algorithm,Problem,M,Run,MaxGen(i),K(j));
                    disp([Problem,'_',num2str(M),'_',num2str(Run),' ',datestr(now)]);
                end
            end
        end
    end
end
end