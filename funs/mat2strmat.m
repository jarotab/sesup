function S = mat2strmat(M,n,varargin)

if nargin < 2
    n = 3;
end

Ms = mat2str(M,n);
idx = strfind(Ms,';');
idX = [[1 idx+1]; [idx, numel(Ms)]];

S = repmat("",size(M,1),1);
for i = 1:size(M,1)
    S(i,:) = string(Ms(idX(1,i):idX(2,i)));    
end

end