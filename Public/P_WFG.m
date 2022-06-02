function [Output,Boundary,Coding] = P_WFG(Operation,Problem,M,Input)
    k = find(~isstrprop(Problem,'digit'),1,'last');
    switch str2double(Problem(k+1))
        case 1
            [Output,Boundary,Coding] = WFG1(Operation,M,Input);
        case 2
            [Output,Boundary,Coding] = WFG2(Operation,M,Input);
        case 3
            [Output,Boundary,Coding] = WFG3(Operation,M,Input);
        case 4
            [Output,Boundary,Coding] = WFG4(Operation,M,Input);
        case 5
            [Output,Boundary,Coding] = WFG5(Operation,M,Input);
        case 6
            [Output,Boundary,Coding] = WFG6(Operation,M,Input);
        case 7
            [Output,Boundary,Coding] = WFG7(Operation,M,Input);
        case 8
            [Output,Boundary,Coding] = WFG8(Operation,M,Input);
        case 9
            [Output,Boundary,Coding] = WFG9(Operation,M,Input);
        otherwise
            error([Problem,'Not Exist']);
    end
end