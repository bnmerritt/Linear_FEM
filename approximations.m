clear;clc;

% Homework question 53-55

% Solve with everything in terms of xsi and across -1 to 1
syms x y x1
degree = 2; %user defined
F_Leg = zeros(1,degree+1);
M_Leg = zeros(degree+1);
F_Ber = zeros(1,degree+1);
M_Ber = zeros(degree+1);
F_Lag = zeros(1,degree+1);
M_Lag = zeros(degree+1);

x = (1/2)*(y + 1);
% f_x = x^3 - (8*x^2)/5 + (3*x)/5;
f_x = sin(pi*x);
nodes = linspace(-1,1,degree+1);

N_Leg(1) = x^0;
N_Leg(2) = 2*x - 1;
N_Leg(3) = (1/2)*(3*(2*x - 1)^2 - 1);
N_Leg(4) = (1/6)*(15*(2*x-1)^3 - 9*(2*x-1));
for basis_idx = 0:degree
    N_Ber(basis_idx+1) = nchoosek(degree,basis_idx) * (x)^basis_idx * (1-(x))^(degree-basis_idx);
    val = x^0;
    for j = 0:degree
        if basis_idx == j
            N_Lag(basis_idx+1) = val;
        else
            val = val*((2*x-1) - nodes(j+1)) / (nodes(basis_idx+1) - nodes(j+1));
            N_Lag(basis_idx+1) = val;
        end
    end
end
for i = 1:(degree+1)
    for j = 1:(degree+1)
        M_Leg(i,j) = double(int(N_Leg(i)*N_Leg(j),-1,1)/2);
        M_Ber(i,j) = double(int(N_Ber(i)*N_Ber(j),-1,1)/2);
        M_Lag(i,j) = double(int(N_Lag(i)*N_Lag(j),-1,1)/2);
    end
    F_Leg(i) = double(int(f_x*N_Leg(i),-1,1)/2);
    F_Ber(i) = double(int(f_x*N_Ber(i),-1,1)/2);
    F_Lag(i) = double(int(f_x*N_Lag(i),-1,1)/2);
end
d_Leg = inv(M_Leg)*F_Leg';
d_Ber = inv(M_Ber)*F_Ber';
d_Lag = inv(M_Lag)*F_Lag';

% interpolating polynomial Lagrange
P = x*0;
for j = 1:length(nodes)
    f = subs(f_x,y,nodes(j));
    fun_basis = N_Lag(j);
    P = P + fun_basis*f;
end
figure(1);
fplot(f_x);
hold on
fplot(P);
xlim([-1,1]);

%Error Lagrange
uh_Leg = x*0;
for i = 1:length(d_Leg)
    uh_Leg = uh_Leg + d_Leg(i)*N_Leg(i);
end
error_uh = double(int(abs(f_x - uh_Leg), 0,1));
error_interpolating_polynomial = double(int(abs(f_x - P), 0,1));