function sp1_draw(g)

wasHeld = ishold;

rho = 10/100;
w   = 20/100;
rad = w/sqrt(2);

x = g(1);
y = g(2);
th = g(3);

x1 = x + rad*cos(th+pi/4);
y1 = y + rad*sin(th+pi/4);

pts = [ [x1, y1] ];

x2 = x + rad*cos(th-pi/4);
y2 = y + rad*sin(th-pi/4);

pts = [ pts ; [x2, y2] ];

x2 = x + rad*cos(th-3*pi/4);
y2 = y + rad*sin(th-3*pi/4);

pts = [ pts ; [x2, y2] ];

x2 = x + rad*cos(th+3*pi/4);
y2 = y + rad*sin(th+3*pi/4);

pts = [ pts ; [x2, y2] ];
pts = [ pts ; [x1, y1] ];

hold on;
plot (pts(:,1), pts(:,2),'k');

x1 = x;
y1 = y;
x2 = x + (w/4)*cos(th);
y2 = y + (w/4)*sin(th);

pts = [ [x1, y1] ; [x2, y2] ];

plot (pts(:,1), pts(:,2), 'r');

x1 = x + (w/8)*cos(th+pi/2);
y1 = y + (w/8)*sin(th+pi/2);
x2 = x + (w/8)*cos(th-pi/2);
y2 = y + (w/8)*sin(th-pi/2);

pts = [ [x1, y1] ; [x2, y2] ];

plot (pts(:,1), pts(:,2), 'r');

if (~wasHeld)
  hold off;
end
