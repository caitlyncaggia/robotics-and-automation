function SCARA_EX1(l1, l2)

alpha1 = (-90:10:90) * pi()/180;
alpha2 = (-180:10:180) * pi()/180;

alpha1Cnt = numel(alpha1);
alpha2Cnt = numel(alpha2);
total = alpha1Cnt*alpha2Cnt;
x     = zeros(1,total);
y     = zeros(1,total);
theta = zeros(1,total);

ll = 1;
for ii=1:alpha1Cnt
  for jj=1:alpha2Cnt
    x(ll) = l1*cos(alpha1(ii))+l2*cos(alpha1(ii)+alpha2(jj));
    y(ll) = l1*sin(alpha1(ii))+l2*sin(alpha1(ii)+alpha2(jj));
    theta(ll) = alpha1(ii)+alpha2(jj);
    ll = ll+1;
  end
end

    
method = 'plot';
switch method
  case 'mesh'
    x     = reshape(x, [alpha2Cnt, alpha1Cnt]);
    y     = reshape(y, [alpha2Cnt, alpha1Cnt]);
    theta = reshape(theta, [alpha2Cnt, alpha1Cnt]);
    figure(1);
    surf(x, y, theta);
    hold on;
    patch(x,y,'g');
    hold off;
    axis equal;
    ha = gca;
    zticklabels = get(ha, 'ZTickLabel');
    zticklabels = num2str(round((180/pi)*str2num(zticklabels)),'%3d');
    set(ha, 'ZTickLabel',zticklabels);
  case 'plot'
    figure(1);
    plot3(x, y, theta,'r*');
    hold on;
    plot(x, y, 'g*');
    axis equal
    hold off;
    %surf(x, y, theta);
    
    figure(2);
    plot(x,y, 'g*');
    axis equal 
end

end
