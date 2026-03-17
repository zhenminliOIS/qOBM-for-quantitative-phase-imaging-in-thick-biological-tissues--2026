function S=brainDiffuse(DIMS,filePath,label)
% DIMS=[2048,2048];
    load([filePath,'\Source_',label,'.mat'],'fitobject');
    u0 = linspace(-1,1,DIMS(1));
    du0 = u0(2)-u0(1);
    u0 = u0(1:end)+(du0/2);
    v0 = linspace(-1,1,DIMS(2));
    dv0 = v0(2)-v0(1);
    v0 = v0(1:end)+(dv0/2);
    [U0,V0] = ndgrid(u0,v0);
    S = fitobject(U0,V0);
    S = S'./(sum(S(:)));
    S(S<0) = 0; S(isnan(S)) = 0;
end