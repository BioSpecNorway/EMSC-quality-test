function ZY = split2levels(NameRange,saisir)

NumLevels = size(NameRange,1);
ZYExt = cell(1,NumLevels);
ZY = cell(1,NumLevels);

for i = 1:NumLevels
Range = NameRange{i};
ZYExt{i} = IndicatorVariables(saisir,Range);
    ZY{i} = ZYExt{i}.i;
end

end
