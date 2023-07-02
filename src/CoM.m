function c = CoM(mass, theta, shapeSize, oil_density, origin)

% Compute the center of mass of a tank
% mass: the mass of oil in the tank
% theta: the pitch angle of the flight
% shapeSize: three-demension vector, the size of the tank
% oil_density: the density of oil in the tank
% origin: the original point of the tank's coordinate frame w.r.t flight
% coordinate frame

assert(length(shapeSize)==3);
a = shapeSize(1);
b = shapeSize(2);
c = shapeSize(3);
V = a*b*c;
M = oil_density * V;

sign = 0;
if theta < 0
    sign = -1;
elseif theta >= 0
    sign = 1;
end


% situation 1: A(x_A, 0, 0.5c_i), B(0.5a_i*sign(theta), 0, z_B)
% the CoM of (M-m) is r((x_A+a_i)/3, 0, (z_B+c_i)/3)
syms xA zB
eqns = [zB - 0.5*c == -tan(theta) * (0.5*a*sign - xA)];
eqns = [eqns, (M - mass) / M == 0.5 * (0.5*a-sign*xA)*(0.5*c-zB)/(a*c)];
eqns = [eqns, xA >= -0.5*a, xA <= 0.5*a, zB >= -0.5 * c, zB <= 0.5 * c];
S = solve(eqns, [xA zB], 'Real', true, 'ReturnConditions', true);

if length(S.xA) == 1 && length(S.zB) == 1
    c = -(M-mass) / mass * [(eval(S.xA) + a * sign) / 3, 0, (eval(S.zB) + c) / 3]';
    c(1) = c(1) + origin(1);
    c(2) = c(2) + origin(2);
    c(3) = c(3) + origin(3);
    return;
end


% situation 2:  A(x_A, 0, 0.5c_i), B(x_B, 0, -0.5c_i)
syms xB
eqns = c == -tan(theta) *  (xA - xB);
eqns = [eqns, mass / M == 0.5 + 0.5 / a * (xA + xB) * sign];
eqns = [eqns, xA >= -0.5*a, xA <= 0.5*a, xB >= -0.5*a, xB <= 0.5*a];
S = solve(eqns, [xA xB], 'Real', true);

if length(S.xA) == 1 && length(S.xB) == 1
    x_A = eval(S.xA);
    x_B = eval(S.xB);
    alpha = (a + 2*x_A*sign) / (x_B - x_A) * sign;
    r1 = [(x_B-x_A) / 3 - 0.5*a*sign, 0, -c/6]';
    r2 = [(x_B - 0.5*a*sign)/2, 0, 0]';
    c = 1/(1+alpha) * r1 + alpha/(1+alpha) * r2;
    c(1) = c(1) + origin(1);
    c(2) = c(2) + origin(2);
    c(3) = c(3) + origin(3);
    return;
end


% situation 3: A(-0.5a_i, 0, z_A), B(0.5a_i, 0, z_B)
syms zA
eqns = zA == zB + a * tan(theta);
eqns = [eqns, mass / M == 0.5 + (zA + zB) / (2 * c)];
eqns = [eqns, zA >= -0.5*c, zA <= 0.5*c, zB >= -0.5*c, zB <= 0.5*c];
S = solve(eqns, [zA, zB], 'Real', true);

if length(S.zA) == 1 && length(S.zB) == 1
    z_A = eval(S.zA);
    z_B = eval(S.zB);
    alpha = (z_A + 0.5*c) / (z_B + 0.5*c);
    r1 = [-0.5*a, 0, z_A+z_B-0.5*c]' / 3;
    r2 = [0.5*a, 0, z_B-c]' / 3;
    c = alpha/(1+alpha) * r1 + 1/(1+alpha) * r2;
    c(1) = c(1) + origin(1);
    c(2) = c(2) + origin(2);
    c(3) = c(3) + origin(3);
    return;
end


% situation 4:  A(x_A, 0, -0.5c_i), B(-0.5a_i*sign(theta), 0, z_B)
eqns = tan(theta) * (xA + 0.5*a*sign) == (0.5*c + zB);
eqns = [eqns, mass / M == 0.5 * (0.5*a + xA*sign)*(zB + 0.5*c) / (a*c)];
eqns = [eqns, xA >= -0.5*a, xA <= 0.5*a, zB >= -0.5*c, zB <= 0.5*c];
S = solve(eqns, [xA zB], 'Real', true);

if length(S.xA) == 1 && length(S.zB) == 1
    x_A = eval(S.xA);
    z_B = eval(S.zB);
    c = [x_A - a*sign, 0, z_B-c]' / 3;
    c(1) = c(1) + origin(1);
    c(2) = c(2) + origin(2);
    c(3) = c(3) + origin(3);
    return;
end

% no solution found
error("No solution found when compute CoM")

end

