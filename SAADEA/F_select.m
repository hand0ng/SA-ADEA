% The selection
function [Population, FunctionValue] = F_select(Class, Population, FunctionValue, cosine)
    global W Replace
    VN = size(W, 1);
    for sub = 1:VN
        if(~isempty(Class(sub).c))
            c = Class(sub).c;
            ind = Class(sub).ind;
            indFunctionValue = Class(sub).indFV;
            maxc = Class(sub).maxc;
            [~, selind] = max(maxc);
            indd1 = F_d1(indFunctionValue(selind,:),sub);
            subd1 = F_d1(FunctionValue(sub,:),sub);         

            if cosine(sub) < maxc(selind) && indd1 < subd1
                Population(sub,:)    = ind(selind,:);
                FunctionValue(sub,:) = indFunctionValue(selind,:);
                Replace(sub)=Replace(sub) + 1;
            end
                
%             %PBI calculation
%             subScalarValue = F_scalar(subFunctionValue, sub);
%             indScalarValue = F_scalar(indFunctionValue, sub);
%             [~, selind] = min(indScalarValue);
% 
%             if indScalarValue(selind) < subScalarValue
%                 Population(sub,:)    = ind(selind,:);
%                 FunctionValue(sub,:) = indFunctionValue(selind,:);
%                 Replace(sub)=Replace(sub) + 1;
%             end

%             candi = [ind; Population(sub,:)];
%             candiFunctionValue = [indFunctionValue;subFunctionValue];
%             NonDominated  = P_sort(candiFunctionValue,'first')==1;
%             candi = candi(NonDominated,:);
%             candiFunctionValue = candiFunctionValue(NonDominated,:);
%             
%             candiScalarValue = F_scalar(candiFunctionValue, sub);
%             [~, selind] = min(candiScalarValue);
%             Population(sub,:)    = candi(selind,:);
%             FunctionValue(sub,:) = candiFunctionValue(selind,:);
        end;
    end;
end