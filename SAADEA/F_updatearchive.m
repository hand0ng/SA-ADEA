function [A, B] = F_updatearchive_new(pop,A,B,Problem,M,K)
% function F_updatearchive()
	global W LCBmin Realmin
    % K-Medoids
    pop.uFV = pop.FV./repmat(sqrt(sum(pop.FV.^2,2)),1,size(pop.FV,2));
    pop.LCB = pop.FV - 2 * pop.S;
    LCBmin = min(pop.LCB,[],1);
    [idx,~,~,~,midx] = kmedoids(pop.uFV,K);
    for i = 1 : K
        Pop = pop.Pop(idx==i,:);
        LCB = pop.LCB(idx==i,:);
        [~, minidx] = min(F_scalar_LCB(LCB,midx(i)));
        Poptemp(i,:) = Pop(minidx,:);
    end
    FVtemp = P_objective('value',Problem,M,Poptemp);
    Realmin = min([Realmin;FVtemp],[],1);
    % update A
    Class = F_associate(A.Pop, A.FV, W);
    DelIndex = zeros(1,K);
    for i = 1:K
         uFVtemp = FVtemp(i,:);
         uFVtemp = (uFVtemp - repmat(Realmin, [size(uFVtemp,1) 1]));
         uFVtemp = uFVtemp./repmat(sqrt(sum(uFVtemp.^2,2)), [1 M]);
         cosineOFF    = uFVtemp*W'; % calculate the cosine values between each solution and each vector
         [~, classfocus] = max(cosineOFF, [], 2);
         if ~isempty(Class(classfocus).c)
             Distance = sum((repmat(FVtemp(i,:),Class(classfocus).num,1) - Class(classfocus).indFV).^2,2);
             [~, maxidx] = max(Distance);
             tmp = Class(classfocus).c(maxidx);
             DelIndex(i) = tmp;
             Class(classfocus).c(maxidx) = [];
             Class(classfocus).indFV(maxidx,:) = [];
             Class(classfocus).num = Class(classfocus).num - 1;
         else
             num_set = cat(1,Class().num);
             num =  max(num_set);
             id = find(num_set==num);
             classfocus = id(unidrnd(length(id)));
             index = unidrnd(num);
             tmp = Class(classfocus).c(index);
             DelIndex(i) = tmp;
             Class(classfocus).c(index) = [];
             Class(classfocus).indFV(index,:) = [];
             Class(classfocus).num = Class(classfocus).num - 1;
         end
    end
    A.Pop(DelIndex,:) = [];
    A.FV(DelIndex,:)  = [];
    A.Pop = [A.Pop;Poptemp];
    A.FV  = [A.FV;FVtemp];
    
    % update B
    B.Pop = [B.Pop;Poptemp];
    B.FV  = [B.FV;FVtemp];
end