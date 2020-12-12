

%==[1] Load the file.
parmFile = 'pikparms.mat';
pp = load(parmFile);

disp('If any columms go from higher to lower, then they need to be flipped.');
pp.alphaLims

%==[2] Check that limits are from low to high, and that orientation agrees.
for ii=1:size(pp.alphaOrient,2)
  if (pp.alphaLims(1,ii) > pp.alphaLims(3,ii))
    pp.alphaLims([1 3],ii) = pp.alphaLims([3, 1], ii);
    if (pp.alphaOrient(ii) == 1)
      pp.alphaOrient(ii) = -1;
    end
  end
end
disp('They should have been flipped now.');
pp.alphaLims


save(parmFile,'-struct','pp');
