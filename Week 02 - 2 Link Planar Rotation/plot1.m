

figure(1);
axis([-5 5 -5 5]);
for ii=1:size(x,1)
  drawHilare(x(ii,:));
  pause(0.1);
end
