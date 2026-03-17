function out = normalizeI(in)
out = (in-min(in(:)))./(max(in(:))-min(in(:)));
end