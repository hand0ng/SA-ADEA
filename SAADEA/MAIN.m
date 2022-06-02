% Decomposition Main File
function MAIN(Algorithm,Problem,M,Run,MaxGen,K)
format compact;tic;
global Theta W Zmin Zmax Replace Nadir Realmin
%%
% basic settings
rand('seed', sum(100 * clock));
[~,N,p1,p2] = P_settings(Algorithm,Problem,M);
Generations = MaxGen;
delta   = 1; % the probability of choosing parents locally
Tm      = ceil(N/10); % the mating neighborhood size
Theta   = 5*ones(N,1);
Replace = zeros(N,1);

% weight vector initialization
[N,Ws] = F_weight(p1,p2,M);
for i = 1:N
    Ws(i,:) = Ws(i,:)./norm(Ws(i,:));
end
W = Ws;

% calculat neighboring angle for angle normalization
cosineW = W*W';
[scosineW, ~] = sort(cosineW, 2, 'descend');
acosW = acos(scosineW(:,2));
refW = rad2deg(acosW);

% perform the mating neighborhood
NB     = pdist2(W,W);
[~,NB] = sort(NB,2);
NB     = NB(:,1:Tm);

%% build models
[~,Boundary,Coding] = P_objective('init',Problem,M,1);
D     = size(Boundary,2);
NI    = 11*D-1;
A.Pop = repmat(Boundary(1,:)-Boundary(2,:),NI,1).*lhsamp(NI,D)+repmat(Boundary(2,:),NI,1);
[~,index] = unique(A.Pop,'rows');
A.Pop     = A.Pop(index,:);
A.FV  = P_objective('value',Problem,M,A.Pop);
B     = A;
THETA = 5.*ones(M,D);
Model = cell(1,M);
for i = 1 : M
    % The parameter 'regpoly1' refers to one-order polynomial
    % function, and 'regpoly0' refers to constant function. The
    % former function has better fitting performance but lower
    % efficiency than the latter one
    dmodel     = dacefit(A.Pop,A.FV(:,i),'regpoly1','corrgauss',THETA(i,:),1e-5.*ones(1,D),100.*ones(1,D));
    Model{i}   = dmodel;
    THETA(i,:) = dmodel.theta;
end
Realmin = min(A.FV,[],1);
%%
% population initialization
if (length(A.Pop)>N)
    duplicated = randi(length(A.Pop),N,1);
    Population = A.Pop(duplicated,:);
else
    [PopulationTmp,~,~] = P_objective('init',Problem,M,N - length(A.Pop));
    Population = [A.Pop;PopulationTmp];
end
for i = 1: N
    for j = 1 : M
        [FunctionValue(i,j),~,MSE(i,j)] = predictor(Population(i,:),Model{j});
    end
end

% extream point
Nadir = max(FunctionValue,[],1);
Zmin = min(FunctionValue,[],1);
NonDominated    = P_sort(FunctionValue,'first')==1;
BFunctionValue  = FunctionValue(NonDominated,:);
BPopulation     = Population(NonDominated,:);
[Zmax, ZmaxInd] = max(BFunctionValue,[],1);
ZmaxFV  = BFunctionValue(ZmaxInd,:);
ZmaxPop = BPopulation(ZmaxInd,:);

%% run SAADEA
MaxFEs = 300;
FEs    = NI;
while FEs < MaxFEs
%     FEs
    %% run ADEA
    for Gene = 1 : Generations
        uFunctionValue = (FunctionValue - repmat(Zmin, [size(FunctionValue,1) 1]));
        uFunctionValue = uFunctionValue./repmat(sqrt(sum(uFunctionValue.^2,2)), [1 M]);
        cosine         = sum(uFunctionValue.*W,2);
        Theta = 0.06*M*(refW+rad2deg(acos(cosine)));
        for sub = 1 : N
            % choose the candidate matingpool
            if rand < delta
                P = NB(sub,:);
            else
                P = 1:N;
            end
            % binary tounament selection for matingpool
            [MatingPool] = F_mating(sub, P, Population, FunctionValue, Coding, W);
            % generate an offspring
            PopulationOF    = P_generator(MatingPool,Boundary,Coding,N);
            for j = 1 : M
                [FunctionValueOF(j),~,MSEOF(j)] = predictor(PopulationOF,Model{j});
            end
            % update the ideal point
            Zmin = min(Zmin,FunctionValueOF);
            % get the nadir point
            Nadir = max(Nadir,FunctionValueOF);
            % the translation
            uFunctionValueOF = (FunctionValueOF - repmat(Zmin, [size(FunctionValueOF,1) 1]));
            % the association
            uFunctionValueOF = uFunctionValueOF./repmat(sqrt(sum(uFunctionValueOF.^2,2)), [1 M]);
            cosineOFF    = uFunctionValueOF*W'; % calculate the cosine values between each solution and each vector
            [maxc, maxcidx] = max(cosineOFF, [], 2);
            R = maxcidx;
            % update the population
            for j = 1 : length(R)
                if F_scalar(FunctionValue(R(j),:),R(j)) > F_scalar(FunctionValueOF,R(j))
                    Population(R(j),:)    = PopulationOF;
                    FunctionValue(R(j),:) = FunctionValueOF;
                    MSE(R(j),:) = MSEOF;
                end
            end
        end
        % update the extream point
        CombineFV      = [FunctionValue;ZmaxFV];%otherOffspringFV];
        CombinePop     = [Population;ZmaxPop];%otherOffspring];
        NonDominated   = P_sort(CombineFV,'first')==1;
        BFunctionValue = CombineFV(NonDominated,:);
        BPopulation    = CombinePop(NonDominated,:);
        [Zmax, ZmaxInd] = max(BFunctionValue,[],1);
        ZmaxFV  = BFunctionValue(ZmaxInd,:);
        ZmaxPop = BPopulation(ZmaxInd,:);
        % update the reference vectors
        if(mod(Gene, ceil(Generations*0.1)) == 0) && Gene~=Generations
            W = Ws;
            W = W.*repmat((Zmax - Zmin)*1.0,N,1);
            for i = 1:N
                W(i,:) = W(i,:)./norm(W(i,:));
            end
            NB     = pdist2(W,W);
            [~,NB] = sort(NB,2);
            NB     = NB(:,1:Tm);
            cosineW = W*W';
            [scosineW, ~] = sort(cosineW, 2, 'descend');
            acosW = acos(scosineW(:,2));
            refW = rad2deg(acosW);
            ZmaxFV  = [];
            ZmaxPop = [];
        end
    %     clc; 
    %     fprintf('Progress %4s%%\n',num2str(round(Gene/Generations*100,-1)));
    end
    %% update model
    pop.Pop = Population;
    pop.FV = FunctionValue;
    pop.S = sqrt(MSE);
    [A,B] = F_updatearchive(pop,A,B,Problem,M,K);
    FEs = FEs + K;
    for i = 1 : M
        % The parameter 'regpoly1' refers to one-order polynomial
        % function, and 'regpoly0' refers to constant function. The
        % former function has better fitting performance but lower
        % efficiency than the latter one
        dmodel     = dacefit(A.Pop,A.FV(:,i),'regpoly1','corrgauss',THETA(i,:),1e-5.*ones(1,D),100.*ones(1,D));
        Model{i}   = dmodel;
        THETA(i,:) = dmodel.theta;
    end
    % update function value
    for i = 1: N
        for j = 1 : M
            [FunctionValue(i,j),~,MSE(i,j)] = predictor(Population(i,:),Model{j});
        end
    end
    % update extream point
    Nadir = max(FunctionValue,[],1);
    Zmin = min(FunctionValue,[],1);
    NonDominated    = P_sort(FunctionValue,'first')==1;
    BFunctionValue  = FunctionValue(NonDominated,:);
    BPopulation     = Population(NonDominated,:);
    [Zmax, ZmaxInd] = max(BFunctionValue,[],1);
    ZmaxFV  = BFunctionValue(ZmaxInd,:);
    ZmaxPop = BPopulation(ZmaxInd,:);
end

NonDominated  = P_sort(B.FV,'first')==1;
Population    = B.Pop(NonDominated,:);
FunctionValue = B.FV(NonDominated,:);
%%
P_output(Population,FunctionValue,toc,Algorithm,Problem,M,Run,MaxGen,K);
end


