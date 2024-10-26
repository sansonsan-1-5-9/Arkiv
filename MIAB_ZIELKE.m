close all
clear all
clc
tic


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
tend = 1.5;     
rho0 = 999.057;     %watrer density at 15.5 Celsius
dx = L/N;       
dt = dx/a;
nu = 1.18181818181818e-06;

NZ = 3;     %Number of nodes computed by Zielke


timestep=1+floor((tend/dt));
time=linspace(0,tend,timestep);

%valve orifice
tau = 1;        %valve opening
t = 0;      %time in the beginning
Em = 1.5;       %taken from example in Streeter


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
Cp = zeros(1,N);
Cm = zeros(1,N);

%%% Parrameters for MIAB Friction model
kx = .03;
kt = .04;

dt_ = zeros(1,N+1);
R = f*dx/(2*g*D*A^2);


%Parameters for Friction computation Zielke
1.18181818181818e-06;
R = f*dx/(2*g*D*A^2);
B = a/(g*A);
Rc = 16*nu*dx/(g*A*D^2);
A_z = .5/sqrt(pi);
B_z = 210.08;       %laminar flow Vardy
Q_Z = zeros(round(tend/dt),NZ+1);
h_u = zeros(1,NZ+1);


%Weighting Function
tauZ=4*nu*time/((D)^2);
tauZ=tauZ+tauZ(2);
W=A_z*exp(-B_z*tauZ)./(sqrt(tauZ));



%Head at end for plot
%time = 1:dt:tend;
% Hend = zeros(1,numel(time));
% Qend = zeros(1,numel(time));
P =1;



% figure(2)


for i = 0:dt:tend
    
    phi_p =sign(Q(2:NS).*(Q(2:NS)-(Q(1:NS-1)))/dx);
    phi_m =sign(Q(1:NS-1).*(Q(2:NS)-(Q(1:NS-1)))/dx);
    for j = 1:N
        if phi_p(j) == 0
            phi_p(j) = 1;
        end
        if phi_m(j) == 0
            phi_m(j) = 1;
        end
    end
    Bp = 2*a/(g*A)*(1+kt)./(-phi_p*kx + sqrt(kx^2 + 4*(1+kt)));
    Bm = 2*a/(g*A)*(1+kt)./(phi_m*kx + sqrt(kx^2 + 4*(1+kt)));
    
    alpha_p = -phi_p*kx + sqrt(kx^2 + 4*(1+kt));
    alpha_m = phi_m*kx + sqrt(kx^2 + 4*(1+kt));
    
    %%%% Finding dx,dt,.... for interpolation
    dxp = 2*dx*alpha_m./(alpha_p+alpha_m);
    dxm = 2*dx - dxp;
    dt_(1) = -dx*alpha_m(1)/(-2*a);
    dt_(2:N+1) = (dxp.*alpha_p)/(2*a);
    dx2p = dt*dxp./dt_(2:N+1);
    dx2m = dt*dxm./dt_(1:N);
    
    %%%% Finding H-,H+....
    %Upstream boundary
    if dt_(1) < dt
        Hm(1) = H(2) + (H(3)-H(2))*(dx2m(1)-dx)/dx;
        Qm(1) = Q(2) + (Q(3)-Q(2))*(dx2m(1)-dx)/dx;
    elseif dt_(1) > dt
        Hm(1) = H(1) + (H(2)-H(1))*dx2p(1)/dx;
        Qm(1) = Q(1) + (Q(2)-Q(1))*dx2p(1)/dx;
    else
        Hm(1) = H(2);
        Qm(1) = Q(2);
    end
        %%% M� finne Qnew,Hnew(1) for � finne Hp,Qp(1)
    Cm(1) = Hm(1)-Bm(1)*Qm(1)+R*Qm(1)*abs(Qm(1));
    Hnew(1) = Hr;
    Qnew(1)=(Hnew(1)-Cm(1))/Bm(1);
    if dx2p(1) > dx
        Hp(1) = Hr;
        Qp(1) = Q(1)+(Qnew(1)-Q(1))*(dt-dt_(2))/dt;
    elseif dx2p(1)<dx
        Hp(1) = H(1) + (H(2)-H(1))*(dx-dx2p(1))/dx;
        Qp(1) = Q(1) + (Q(2)-Q(1))*(dx-dx2p(1))/dx;
    else
        Hp(1)=Hr;
        Qp(1)=Q(1);
    end
    %Interior Points H+H-Q+-
    for k = 2:NS-NZ
        if dx2p(k)>dx
            Hp(k)=H(k-1)+(H(k)-H(k-1))*(2*dx-dx2p(k))/dx;
            Qp(k)=Q(k-1)+(Q(k)-Q(k-1))*(2*dx-dx2p(k))/dx;
        elseif dx2p(k)<dx
            Hp(k)=H(k)+(H(k+1)-H(k))*(dx-dx2p(k))/dx;
            Qp(k)=Q(k)+(Q(k+1)-Q(k))*(dx-dx2p(k))/dx;
        else
            Hp(k) = H(k);
            Qp(k) = Q(k);
        end
        if dx2m(k)<dx
            Hm(k)=H(k)+(H(k+1)-H(k))*dx2m(k)/dx;
            Qm(k)=Q(k)+(Q(k+1)-Q(k))*dx2m(k)/dx;
        elseif dx2m(k)>dx
            Hm(k)=H(k+1)+(H(k+2)-H(k+1))*(dx2m(k)-dx)/dx;
            Qm(k)=Q(k+1)+(Q(k+2)-Q(k+1))*(dx2m(k)-dx)/dx;
        else
            Hm(k)=H(k+1);
            Qm(k)=Q(k+1);
        end
    end
    
    %%%% Computing Hnew, Qnew
    Cp(1:NS-NZ) = Hp + Bp(1:NS-NZ).*Qp - R*Qp.*abs(Qp);
    Cm(1:NS-NZ) = Hm - Bm(1:NS-NZ).*Qm + R*Qm.*abs(Qm);
    %interior points 2:N
    Qnew(2:NS-NZ)=(Cp(1:NS-NZ-1)-Cm(2:NS-NZ))./(Bm(2:NS-NZ)+Bp(1:NS-NZ-1));
    Hnew(2:NS-NZ)=Cp(1:NS-NZ-1)-Bp(1:NS-NZ-1).*Qnew(2:NS-NZ);
    %Boundaries already computed!!!
    
    
    %%%Computing last NZ points and downstream boundary using Zielke's
    %%%method
    
    Q_Z(P,:)=Q(NS-NZ:NS);
    for j=1:(P-1)
        h_u = h_u+(Q_Z(j+1,:)-Q_Z(j,:))*W(P-j);
    end
    
    h_u = Rc*h_u;
    
    %Cp, Cv for computing nodes NS-NZ+1 --> NS
    Cp(NS-NZ:N) = H(NS-NZ:N) + B*Q(NS-NZ:N) - R*Q(NS-NZ:N).*abs(Q(NS-NZ:N)) - h_u(1:NZ);
    Cm(NS-NZ:N) = H(NS-NZ+1:NS) - B*Q(NS-NZ+1:NS) + R*Q(NS-NZ+1:NS).*abs(Q(NS-NZ+1:NS))+h_u(2:NZ+1);
    
    
    Hnew(NS-NZ+1:N) = .5*(Cp(NS-NZ:N-1)+Cm(NS-NZ+1:N));
    Qnew(NS-NZ+1:NS-1) = (Hnew(NS-NZ+1:N)-Cm(NS-NZ+1:N))/B;
    
    %%&& Lower boundary
    if t <= tc   
        tau=1-(1/tc)*t; 
    else
        tau=0;
    end 
    
    Cv = (Q0(1)*tau)^2/(2*H0(NS));
    
    if B*Cv+2*Cp(NS-1) < 0
        Qnew(NS) = -B*Cv - sqrt(abs((B*Cv)^2+2*Cv*Cp(NS-1)));
    else
        Qnew(NS) = -B*Cv + sqrt((B*Cv)^2+2*Cv*Cp(NS-1));
    end
    
    Hnew(NS) = Cp(NS-1) - B*Qnew(NS);
    
    
    
    %Values for plot
    Hend(P) = Hnew(NS);
    Hmid(P) = Hnew(NS/2);
    Qend(P) = Qnew(NS);
    Vend(P) = Qnew(NS)/A;
    Qin(P) = Qnew(1);
    Hin(P) = Hnew(2);
    
    disp(t)
    Q = Qnew;
    H = Hnew;
    t=t+dt;
    P = P+1;
    h_u = zeros(1,NZ+1);
    
%     grid on
%     subplot(2,1,1)
%     plot(x,H);
%     title(['Pressure through pipe at t= ', num2str(t,3)]);
%     axis([0 L .5*Hr Hr*1.5]);
%     subplot(2,1,2)
%     plot(x,Q);
%     axis([0 L -Q0(1) Q0(1)*1.1]);
%     title('Volume flow through pipe')
%     pause(.00001);
end


toc

MIAB_Zielke_trans(:,1) = time;
MIAB_Zielke_trans(:,2) = Hend;

figure(3)
plot(time,Hend)
grid on;
title(['Transient head just upstream valve a=', num2str(a,4)])
figure(4)
plot(time,Qend)
title('volume flow as a function of time')
axis([0 0.2*tend -Qend(1) Qend(1)]);
% hold on
