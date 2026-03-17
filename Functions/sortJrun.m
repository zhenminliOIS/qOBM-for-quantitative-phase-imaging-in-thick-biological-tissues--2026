function Jorder = sortJrun(Jrun)

Jorder = cell(4,0);

counter = 1;
for m = 1:length(Jrun)   
    lengths = sum(~cellfun(@isempty,Jorder),2);
    Jorder{counter,lengths(counter)+1} = Jrun(m);
    counter = counter + 1;
    if counter > 4; counter = 1; end
end

end