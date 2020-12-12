function planar_r3(l1, l2, l3, method)

if (nargin < 4)
    method = 'plot';
end

alpha1 = 2*(-90:10:90) * pi()/180;
alpha2 = 2*(-90:10:90) * pi()/180;
alpha3 = 2*(-90:10:90) * pi()/180;
alpha1Cnt = numel(alpha1);
alpha2Cnt = numel(alpha2);
alpha3Cnt = numel(alpha3);
total = alpha1Cnt*alpha2Cnt*alpha3Cnt;

x = zeros(1,total);
y = zeros(1,total);
theta = zeros(1,total);
ll = 1;

for ii=1:alpha1Cnt
    for jj=1:alpha2Cnt
        for kk=1:alpha3Cnt
            x(ll) = l1*cos(alpha1(ii))+l2*cos(alpha1(ii)+alpha2(jj)) ...
                +l3*cos(alpha1(ii)+alpha2(jj)+alpha3(kk));
            y(ll) = l1*sin(alpha1(ii))+l2*sin(alpha1(ii)+alpha2(jj)) ...
                +l3*sin(alpha1(ii)+alpha2(jj)+alpha3(kk));
            theta(ll) = alpha1(ii)+alpha2(jj)+alpha3(kk);
            ll = ll+1;
        end
    end
end

switch method
    case 'mesh'
        x = reshape(x, [alpha2Cnt*alpha3Cnt, alpha1Cnt]);
        y = reshape(y, [alpha2Cnt*alpha3Cnt, alpha1Cnt]);
        theta = reshape(theta, [alpha2Cnt*alpha3Cnt, alpha1Cnt]);
        figure(1);
        patch(x, y,'g*');
    case 'plot'
        plot3(x, y, theta,'r*');
        hold on;
        plot(x, y, 'g*');
        axis equal
        hold off;
        %surf(x, y, theta);
end

end