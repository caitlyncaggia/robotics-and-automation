gAO = [sqrt(2)/2 -sqrt(2)/2 7;sqrt(2)/2 sqrt(2)/2 4; 0 0 1];
gBO = [0 -1 2; 1 0 7; 0 0 1];

gOA = inv(gAO);
gOB = inv(gBO);


gCA = [cosd(-15) -sind(-15) 1; sind(-15) cosd(-15) 4; 0 0 1];

gCO = gCA * gAO;
gCB = gCO * gOB;