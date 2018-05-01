%{
This function runs the Cos-FFT method for option pricing, according to 
Fang, F., & Oosterlee, C. W. (2008). A novel pricing method for European 
options based on Fourier-cosine series expansions. SIAM Journal on 
Scientific Computing, 31(2), 826-848.

Inputs needed are:
N, the precision; a scalar
S0, the price of the underlying asset; a scalar
K, the strike price; a scalar
r, the interest rate; a scalar
lambda, eta, mu, u, u0, and rho, parameters needed to calculate the characteristic function
    phi_hes and the boundaries; scalars
T, the time of maturity; a scalar
t0, the start date; a scalar
put, an indicator for whether the option is a put (put==1) or a call
    (put==0)

Outputs produced are:
v, the price of the option; a scalar
phi_hes, the characteristic function; a (1xN) vector
omega, needed for plotting; a (1xN) vector
Uk, needed for further calculations; a (1xN) vector
%}

function [v,phi_hes, omega, Uk]=optionpriceCOSFFT(N,S0,K,r,lambda,eta,u,u0,rho,T,t0,mu,put)

%% Calculated inputs

dt=T-t0;

%formulas from table 11, Heston model

c1 = mu*T + (1-exp(-lambda*T)) * (u -u0)/(2*lambda) - 0.5*u*T;
c2= 1/(8*lambda^3)...
*(eta*T*lambda*exp(-lambda*T)*(u0-u)*(8*lambda*rho-4*eta)...
+ lambda*rho*eta*(1-exp(-lambda*T))*(16*u-8*u0)...
+ 2*u*lambda*T*(-4*lambda*rho*eta+eta^2+4*lambda^2)...
+eta^2*((u-2*u0)*exp(-2*lambda*T)+u*(6*exp(-lambda*T)-7)+2*u0)...
+8*lambda^2*(u0-u)*(1-exp(-lambda*T)));
w=0;

x=log(S0/K); %calculate x as in section 3.1

%boundaries, from section 5.3
a=c1 - 12*sqrt(abs(c2)); %lower boundary
b=c1 + 12*sqrt(abs(c2));%upper boundary


%% Calculations

% 1) calculate phi_hes, which we only have to do once
phi_hes=zeros(1,N); %set up phi, we'll have 1 value for every k
Uk= zeros(1,N); %set up Uk, we'll have 1 value for every k
omega=zeros(1,N); %set up omega, 1 value for every k

for k=0:N-1
    omega(k+1) = k*pi/(b-a); %fill the vector omega
    
    %call or put option? for calculation of chi and psi
    if put==1 %put option
        c=a;
        d=0;
    else %call option
        c=0;
        d=b;
    end
    
    %calculation of chi(c,d), formula 22
    chi_k = 1/(1+ omega(k+1)^2) * ...
        (cos(k*pi*(d-a)/(b-a))*exp(d)...
        -cos(k*pi*(c-a)/(b-a))*exp(c)...
        + omega(k+1)*sin(k*pi*(d-a)/(b-a))*exp(d)...
        - omega(k+1)*sin(k*pi*(c-a)/(b-a))*exp(c));
    
    %calculation of psi(c,d), formula 23
    if k==0
    psi_k = d-c;  
    else
    psi_k = (sin(k*pi*(d-a)/(b-a))...
        -sin(k*pi*(c-a)/(b-a)))*(b-a)/(k*pi);
    end
    
    %calculation of Uk, formula 29
    if put==1
        Uk(k+1)=2/(b-a)*(-chi_k + psi_k); %put option
    else
        Uk(k+1)=2/(b-a)*(chi_k-psi_k); %call option
    end
    
    %calculate inputs for phi_hes
    D = sqrt ( (lambda - 1i*rho*eta*omega(k+1))^2 ...
    + (omega(k+1)^2 + 1i*omega(k+1))*eta^2); % formula 34, part 2
    
    G = (lambda - 1i*rho*eta*omega(k+1) - D) / (lambda - 1i*rho*eta*omega(k+1) + D); %formula 34, part 2
    
    %calculate phi_hes, the characteristic function
    phi_hes(k+1)= exp(1i*omega(k+1)*mu*dt ...
        + u0/(eta^2)*((1-exp(-D*dt))/(1-G*exp(-D*dt))) ... 
        * (lambda-1i*rho*eta*omega(k+1)-D)) ...
        * exp(lambda*u/(eta^2) ...
        * (dt * (lambda - 1i*rho*eta*omega(k+1) - D) ...
        - 2*log( (1-G*exp(-D*dt))/(1-G)))); %formula 34, part 2
    
end
phi_hes(1)=phi_hes(1)/2; %apply the formula for k=0 here


% 2) calculate the function for all x in the vector range, according to
% formula 34

range=a:0.01:b; %setting up the output
allk=0:N-1;
longSum = phi_hes.*Uk.*exp((x-a)/(b-a)*1i*pi*allk); %an intermediate step in the formula
v = K*exp(-r*dt)*real(sum(longSum));
