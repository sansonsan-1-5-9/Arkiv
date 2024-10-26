clear all
clc
tic

cpuT=cputime;

N = 1001;       %grid boxes
a = 1319;       %wave propagation 1286 expremintal, 1319 mandair sim
Hr = 32;      %reservoir pressure head for BC
D = .0221;       %pipe diameter
A = .25*pi*D^2;       %crossectional area 
g = 9.8066502;       %gravitational const 9.8066502
L = 37.23;       %Pipelength
f = .0345;       %Darcy-Weisbach friction factor
NS = N+1;        %number of grid points
tc = 0.009;        %closing time
tend = .5;     
rho0 = 999.057;     %watrer density at 15.5 Celsius
dx = L/N;       
dt = dx/a; %at the Courant boundary, could be smaller
nu = 1.18181818181818e-06;

timestep=1+floor((tend/dt));
time=linspace(0,tend,timestep);
Pmax=size(time);
Pmax=Pmax(2);

%valve orifice
tau = 1;        %valve opening
t = 0;      %time in the beginning
Em = 1.5;       %taken from example in Streeter


%Parameters for Friction computation Zielke
R = f*dx/(2*g*D*A^2);
B = a/(g*A);
Rc = 16*nu*dx/(g*A*D^2);
A_z = .5/sqrt(pi);
B_z = 210.08;       %laminar flow Vardy
Q_Z = zeros(round(tend/dt),NS);
h_u = zeros(1,NS);

%Weighting Function
tauZ=4*nu*time/((D)^2);
tauZ=tauZ+tauZ(2);
W=A_z*exp(-B_z*tauZ)./(sqrt(tauZ));



Nvector = 1:NS;
x = 0:dx:L;

%initial values and steady state
V0 = .1;
Q0(1:NS) = V0*A;
H0(1:NS) = Hr - f*dx*(Nvector-1).*Q0.*abs(Q0)/(2*g*A^2*D);
H00 = g*rho0*H0;        %Pressure head in Pa

figure(1)
subplot(2,1,1)
plot(x,H0);
title('initial pressure through pipe')
axis([0 L Hr*.999 Hr]);
subplot(2,1,2)
plot(x,Q0);
axis([0 L 0 Q0(1)*1.1]);
title('initial volume flow through pipe')

Q = Q0;
H = H0;

Qnew = zeros(1,NS);
Hnew = zeros(1,NS);

P =1;
P985=H00(985)
Massflow_init=Q0(1)*rho0

%---------------Variables for Heat Equation-----------------------%
% Z_dVdx2 = zeros(round(tend/dt),N+1);
% Z_dPdx = zeros(round(tend/dt),N+1);
% %steadystate condition
% ph=rho0*g*H;
% Vh=Q/A;
% 
% Z_dVdx2(1,1) = (Vh(1)-2*Vh(2)+Vh(3))/(dx^2);
% Z_dVdx2(1,2:N) = (Vh(3:N+1)-2*Vh(2:N)+Vh(1:N-1))/(dx^2);
% Z_dVdx2(1,1) = (Vh(N+1)-2*Vh(N)+Vh(N-1))/(dx^2);
% Z_dPdx(1,1) = (ph(2)-ph(1))/dx;
% Z_dPdx(1,2:N) = (ph(3:N+1)-ph(1:N-1))/(2*dx);
% Z_dPdx(1,N+1) = (ph(N+1)-ph(N))/dx;


%figure(2)

for i = 0:dt:tend
    %Parameters for Zielke
    Q_Z(P,:)=Q;
    for j=1:(P-1)
        h_u = h_u+(Q_Z(j+1,:)-Q_Z(j,:))*W(P-j);
    end
    h_u = Rc*h_u;
    
    %All Cp, Cv
    Cp = H(1:NS-1) + B*Q(1:NS-1) - R*Q(1:NS-1).*abs(Q(1:NS-1)) - h_u(1:NS-1);
    Cm = H(2:NS) - B*Q(2:NS) + R*Q(2:NS).*abs(Q(2:NS)) + +h_u(2:NS);
    
    %Interior points
    Hnew(2:NS-1) = .5*(Cp(1:NS-2)+Cm(2:NS-1));
    Qnew(2:NS-1) = (Hnew(2:NS-1)-Cm(2:NS-1))/B;
    
    %upper boundary
    Hnew(1) = Hr;
    Qnew(1) = (Hnew(1)-Cm(1))/B;
    
    if t <= tc   
        tau=(1-(1/tc)*t)^Em; 
    else
        tau=0;
    end
    
    
    %lower boundary
    Cv = (Q0(1)*tau)^2/(2*H0(NS));
    
    if B*Cv+2*Cp(NS-1) < 0
        Qnew(NS) = -B*Cv - sqrt(abs((B*Cv)^2+2*Cv*Cp(NS-1)));
    else
        Qnew(NS) = -B*Cv + sqrt((B*Cv)^2+2*Cv*Cp(NS-1));
    end
    
    Hnew(NS) = Cp(NS-1) - B*Qnew(NS);
    
    
    disp(i);
    %Values for plot
    Hend(P) = Hnew(NS);
    Qend(P)=Qnew(NS);
    Hmid(P) = Hnew(NS/2);
    hfmid(P)=h_u(NS/2);
    P985(P)=rho0*g*Hnew(985);
    M985(P)=rho0*Qnew(985);
    Vend(P) = Qnew(NS)/A;
    Qin(P) = Qnew(1);
    
    P = P + 1;
    Q = Qnew;
    H = Hnew;
    h_u = zeros(1,NS);
    
    
%-------------------Gradients fro heatEQ-------------------------%
%     ph=rho0*g*H;
%     Vh=Q/A;
%     Z_dVdx2(P,1) = (Vh(1)-2*Vh(2)+Vh(3))/(dx^2);
%     Z_dVdx2(P,2:N) = (Vh(3:N+1)-2*Vh(2:N)+Vh(1:N-1))/(dx^2);
%     Z_dVdx2(P,1) = (Vh(N+1)-2*Vh(N)+Vh(N-1))/(dx^2);
%     Z_dPdx(P,1) = (ph(2)-ph(1))/dx;
%     Z_dPdx(P,2:N) = (ph(3:N+1)-ph(1:N-1))/(2*dx);
%     Z_dPdx(P,N+1) = (ph(N+1)-ph(N))/dx;
    
    
    t=t+dt;
    
    
%     plot(x,H);
%     grid on
%     title(['Transient head t= ', num2str(t,3)]);
%     axis([0 L .5*Hr Hr*1.5]);
%     pause(.00001);
end

Zielke_trans_values(:,1) = time;
Zielke_trans_values(:,2) = Hend;

e_cpuT=cputime-cpuT;
disp(e_cpuT)
toc

figure(3)
plot(time,Hend)
grid on;
title(['Transient head just upstream valve a=', num2str(a,4)])
figure(4)
plot(time,Qend)
title('volume flow as a function of time')
axis([0 0.2*tend -Qend(1) Qend(1)]);
hold on

figure(5)
hold on
plot(time,hfmid/max(hfmid))
plot(time,W/max(W))
xlim([0 .1])
title(['Unsteady Headloss at midpoint a=', num2str(a,4)])







