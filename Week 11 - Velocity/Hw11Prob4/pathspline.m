function coef = pathspline(alpha_i, alpha_f, T)

[x, y] = size(alpha_i);

coef = zeros(x,4);
coef(:,1) = alpha_i;
coef(:,2) = zeros(x,1);
coef(:,3) = (3/(T^2))*(alpha_f-alpha_i);
coef(:,4) = (2/(T^3))*(alpha_i-alpha_f);

end