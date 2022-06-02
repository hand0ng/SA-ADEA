function [Output,Boundary,Coding] = WFG5(Operation,M,Input)
    Boundary = NaN; Coding = NaN;
    K = M-1;
    switch Operation
        case 'init'
            D        = M + 9;
            MinValue = zeros(1,D);
            MaxValue = 2 : 2 : 2*D;
            Population = rand(Input,D);
            Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
            Output   = Population;
            Boundary = [MaxValue;MinValue];
            Coding   = 'DE';
        case 'value'
            PopDec = Input;
            [N,D]  = size(PopDec);
            
            L = D - K;
            D = 1;
            S = 2 : 2 : 2*M;
            A = ones(1,M-1);
            
            z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
            
            t1 = zeros(N,K+L);
            t1 = s_decept(z01,0.35,0.001,0.05);

            t2 = zeros(N,M);
            for i = 1 : M-1
                t2(:,i) = r_sum(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
            end
            t2(:,M) = r_sum(t1(:,K+1:K+L),ones(1,L));

            x = zeros(N,M);
            for i = 1 : M-1
                x(:,i) = max(t2(:,M),A(:,i)).*(t2(:,i)-0.5)+0.5;
            end
            x(:,M) = t2(:,M);

            h = concave(x);
            PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
            
            Output = PopObj;
        case 'true'
            h = UniformPoint(Input,M);
            h = h./repmat(sqrt(sum(h.^2,2)),1,M);
            h = repmat(2:2:2*M,size(h,1),1).*h;
            Output = h;
    end
end

function Output = s_decept(y,A,B,C)
    Output = 1+(abs(y-A)-B).*(floor(y-A+B)*(1-C+(A-B)/B)/(A-B)+floor(A+B-y)*(1-C+(1-A-B)/B)/(1-A-B)+1/B);
end

function Output = r_sum(y,w)
    Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
end

function Output = concave(x)
    Output = fliplr(cumprod([ones(size(x,1),1),sin(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),cos(x(:,end-1:-1:1)*pi/2)];
end

function W = UniformPoint(N,M)
    H1 = 1;
    while nchoosek(H1+M,M-1) <= N
        H1 = H1 + 1;
    end
    W = nchoosek(1:H1+M-1,M-1) - repmat(0:M-2,nchoosek(H1+M-1,M-1),1) - 1;
    W = ([W,zeros(size(W,1),1)+H1]-[zeros(size(W,1),1),W])/H1;
    if H1 < M
        H2 = 1;
        while nchoosek(H1+M-1,M-1)+nchoosek(H2+M,M-1) <= N
            H2 = H2 + 1;
        end
        W2 = nchoosek(1:H2+M-1,M-1) - repmat(0:M-2,nchoosek(H2+M-1,M-1),1) - 1;
        W  = [W;([W2,zeros(size(W2,1),1)+H2]-[zeros(size(W2,1),1),W2])/H2];
    end
end