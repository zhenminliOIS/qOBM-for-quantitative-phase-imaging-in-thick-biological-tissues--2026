function [A] = cspad(A,x)
%Patrick Ledwig 4/5/2019 
%Add zero padding after circle shift operation
% A - matrix
% x - amount shifted
x=round(x);
Ndims=ndims(A);
dimorder=1:Ndims;
for n=1:Ndims
    A=permute(A,circshift(dimorder,n-1));
    if x(n)>0
        A(1:x(n),:)=0;
    elseif x(n)<0
        A(end+x(n)+1:end,:)=0;
    end
    A=permute(A,circshift(dimorder,1-n));
end
end

