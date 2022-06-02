function [Generations,N,varargout] = P_settings(Algorithm,Problem,M)
    [Generations,N] = set_problem(Problem,M);
    Parameter = set_algorithm(Algorithm,Problem,M);
    varargout = num2cell(Parameter);
end

function [Generations,N] = set_problem(Problem,M)
    k = find(~isstrprop(Problem,'digit'),1,'last');
    D = str2double(Problem(k+1:end));
    Problem = Problem(1:k);
    if M==15
        switch Problem
            case {'DTLZ'}
                Generations = [1000 500 1000 500 500 1000 500];
                Generations = Generations(D);
            case {'WFG'}
                Generations = [1000 1000 1000 1000 1000 1000 1000 1000 1000];
                Generations = Generations(D);
        end
        N = 135;
    else
    switch Problem
        case {'DTLZ'}
            if M < 2 || M > 10
                error('Objective Number Not Supported !');
            end
            if D < 1 || D > 7
                error([Problem,'Not Exist']);
            end
            Generations = [1000 500 1000 500 500 1000 500];
            Generations = Generations(D);
        case {'SDTLZ'}
            if M < 2 || M > 10
                error('Objective Number Not Supported !');
            end
            if D < 1 || D > 7
                error([Problem,'Not Exist']);
            end
            Generations = [1000 500 1000 500];
            Generations = Generations(D);
        case {'WFG'}
            if M < 2 || M > 10
                error('Objective Number Not Supported !');
            end
            if D < 1 || D > 9
                error([Problem,'Not Exist']);
            end
            Generations = [1000 1000 1000 1000 1000 1000 1000 1000 1000];
            Generations = Generations(D);
            Generations = 20;
        case {'UF'}
            if D < 1 || D > 9
                error([Problem,'Not Exist']);
            end
            Generations = [1000 1000 1000 1000 1000 1000 1000 1000 1000];
            Generations = Generations(D);
        case {'Refine'}
            Generations = 20;
        case {'RefineHCR'}
            Generations = 20;
        case {'RefineHCR2FEED'}
            Generations = 20;
        otherwise
            error([Problem,'Not Exist']);
    end
    N = [50 105 120 126 132 112 156 90 275];
    N = N(M-1);
    end
end

function Parameter = set_algorithm(Algorithm,Problem,M)
    Parameter = NaN;
    k = find(~isstrprop(Problem,'digit'),1,'last');
    D = str2double(Problem(k+1:end));
    Problem = Problem(1:k);
    switch Algorithm
        case {'RVEA'}
            p1 = [49 13  7  5  4  3  3  2  3]; 
            p2 = [ 0  0  0  0  1  2  2  2  2];
            if M==15
                p1 = 2;
                p2 = 1;
            else
                p1 = p1(M-1);
                p2 = p2(M-1);
            end
            Parameter(1) = p1; Parameter(2) = p2; %H1��H2
        case {'Decomposition'}
            p1 = [49 13  7  5  4  3  3  2  3]; 
            p2 = [ 0  0  0  0  1  2  2  2  2];
            if M==15
                p1 = 2;
                p2 = 1;
            else
                p1 = p1(M-1);
                p2 = p2(M-1);
            end
            Parameter(1) = p1; Parameter(2) = p2; %H1��H2
        case 'null'
        otherwise
            error([Algorithm,'Not Exist']);
    end
end

