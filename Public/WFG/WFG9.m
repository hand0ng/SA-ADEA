function [Output,Boundary,Coding] = WFG9(Operation,M,Input)
    Boundary = NaN; Coding = NaN;
    K = M-1;
    switch Operation
        case 'init'
            D        = M + 9;
            MinValue    = zeros(1,D);
            MaxValue    = 2 : 2 : 2*D;
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
            % Same as <t1(:,i)=b_param(z01(:,i),r_sum(z01(:,i+1:end),ones(1,K+L-i)),0.98/49.98,0.02,50)>
            Y = (fliplr(cumsum(fliplr(z01),2))-z01)./repmat(K+L-1:-1:0,N,1);
            t1(:,1:K+L-1) = z01(:,1:K+L-1).^(0.02+(50-0.02)*(0.98/49.98-(1-2*Y(:,1:K+L-1)).*abs(floor(0.5-Y(:,1:K+L-1))+0.98/49.98)));
            % ------------------------------------------------------------------------------------------
            t1(:,end)     = z01(:,end);

            t2 = zeros(N,K+L);
            t2(:,1:K)     = s_decept(t1(:,1:K),0.35,0.001,0.05);
            t2(:,K+1:end) = s_multi(t1(:,K+1:end),30,95,0.35);

            t3 = zeros(N,M);
            for i = 1 : M-1
                t3(:,i) = r_nonsep(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),K/(M-1));
            end
            % Same as <t3(:,M)=r_nonsep(t2(:,K+1:end),L)>
            SUM = zeros(N,1);
            for i = K+1 : K+L-1
                for j = i+1 : K+L
                    SUM = SUM + abs(t2(:,i)-t2(:,j));
                end
            end
            t3(:,M) = (sum(t2(:,K+1:end),2)+SUM*2)/ceil(L/2)/(1+2*L-2*ceil(L/2));
            % -------------------------------------------

            x = zeros(N,M);
            for i = 1 : M-1
                x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
            end
            x(:,M) = t3(:,M);

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

function Output = b_param(y,Y,A,B,C)
    Output = y.^(B+(C-B)*(A-(1-2*Y).*abs(floor(0.5-Y)+A)));
end

function Output = r_sum(y,w)
    Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
end

function Output = s_decept(y,A,B,C)
    Output = 1+(abs(y-A)-B).*(floor(y-A+B)*(1-C+(A-B)/B)/(A-B)+floor(A+B-y)*(1-C+(1-A-B)/B)/(1-A-B)+1/B);
end

function Output = s_multi(y,A,B,C)
    Output = (1+cos((4*A+2)*pi*(0.5-abs(y-C)/2./(floor(C-y)+C)))+4*B*(abs(y-C)/2./(floor(C-y)+C)).^2)/(B+2);
end

function Output = r_nonsep(y,A)
    Output = zeros(size(y,1),1);
    for j = 1 : size(y,2)
        Temp = zeros(size(y,1),1);
        for k = 0 : A-2
            Temp = Temp+abs(y(:,j)-y(:,1+mod(j+k,size(y,2))));
        end
        Output = Output+y(:,j)+Temp;
    end
    Output = Output./(size(y,2)/A)/ceil(A/2)/(1+2*A-2*ceil(A/2));
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