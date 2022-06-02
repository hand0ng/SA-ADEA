function [ScalarValue] = F_scalar_LCB(FunctionValue, sub)
    global W LCBmin
	[N, ~] = size(FunctionValue);
    normW   = sqrt(sum(W(sub,:).^2,2));
	normP   = sqrt(sum((FunctionValue-repmat(LCBmin,N,1)).^2,2));
    CosineP = sum((FunctionValue-repmat(LCBmin,N,1)).*repmat(W(sub,:),N,1),2)./normW./normP;
    ScalarValue = normP.*CosineP;
end