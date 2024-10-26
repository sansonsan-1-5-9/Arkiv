close all
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
tend = 0.101;     
rho0 = 999.057;     %watrer density at 15.5 Celsius
%grid
dx = L/N;       
dt = dx/a; 


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

%%% Parrameters for MIAB Friction model
kx = .03;
kt = .04;

dt_ = zeros(1,N+1);
R = f*dx/(2*g*D*A^2);



%---------------Variables for Heat Equation-----------------------%
M_dVdx2 = zeros(round(tend/dt),N+1);
M_dPdx = zeros(round(tend/dt),N+1);
%steadystate condition
ph=rho0*g*H;
Vh=Q/A;

M_dVdx2(1,1) = (Vh(1)-2*Vh(2)+Vh(3))/(dx^2);
M_dVdx2(1,2:N) = (Vh(3:N+1)-2*Vh(2:N)+Vh(1:N-1))/(dx^2);
M_dVdx2(1,1) = (Vh(N+1)-2*Vh(N)+Vh(N-1))/(dx^2);
M_dPdx(1,1) = (ph(2)-ph(1))/dx;
M_dPdx(1,2:N) = (ph(3:N+1)-ph(1:N-1))/(2*dx);
M_dPdx(1,N+1) = (ph(N+1)-ph(N))/dx;



%-----------------Time Variables-----------------------%
P =1;
timestep=1+floor((tend/dt));
time=linspace(0,tend,timestep);


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
    for k = 2:N-1
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
    %downstream
    if dt_(N+1)<dt
        Hp(N)=H(N-1)+(H(N)-H(N-1))*(2*dx-dx2p(N))/dx;
        Qp(N)=Q(N-1)+(Q(N)-Q(N-1))*(2*dx-dx2p(N))/dx;
    elseif dt_(N+1)>dt
        Hp(N)=H(N)+(H(N+1)-H(N))*(dx-dx2p(N))/dx;
        Qp(N)=Q(N)+(Q(N+1)-Q(N))*(dx-dx2p(N))/dx;
    else
        Hp(N)=H(N);
        Qp(N)=Q(N);
    end
        %%% M� finne Qnew Hnew(N+1) for � finne Hm Qm (N)
    Cp(N)=Hp(N)+Bp(N)*Qp(N)-R*Qp(N)*abs(Qp(N));
    if t <= tc   
        tau=(1-(1/tc)*t)^Em; 
    else
        tau=0;
    end
    Cv = (Q0(1)*tau)^2/(2*H0(NS));
    if (Bp(N)*Cv)^2+2*Cv*Cp(N)>=0
        Qnew(N+1)=-Bp(N)*Cv+sqrt(abs((Bp(N)*Cv)^2+2*Cv*Cp(N)));
    else
        Qnew(N+1)=-Bp(N)*Cv-sqrt(abs((Bp(N)*Cv)^2+2*Cv*Cp(N)));
    end
    Hnew(N+1)=Cp(N)-Bp(N)*Qnew(N+1);
    % Hm(N) downstream
    if dx2m(N)<dx
        Hm(N)=H(N)+(H(N+1)-H(N))*dx2m(N)/dx;
        Qm(N)=Q(N)+(Q(N+1)-Q(N))*dx2m(N)/dx;
    elseif dx2m(N)>dx
        Hm(N)=H(N+1)+(Hnew(N+1)-H(N+1))*(dt-dt_(N+1))/dt;
        Qm(N)=Q(N+1)+(Qnew(N+1)-Q(N+1))*(dt-dt_(N+1))/dt;
    else
        Hm(N)=H(N+1);
        Qm(N)=Q(N+1);
    end
    
    %%%% Computing Hnew, Qnew
    Cp = Hp + Bp.*Qp - R*Qp.*abs(Qp);
    Cm = Hm - Bm.*Qm + R*Qm.*abs(Qm);
    %interior points 2:N
    Qnew(2:N)=(Cp(1:N-1)-Cm(2:N))./(Bm(2:N)+Bp(1:N-1));
    Hnew(2:N)=Cp(1:N-1)-Bp(1:N-1).*Qnew(2:N);
    %Boundaries already computed!!!
    
    
    disp(t);
    %Values for plot
    Hend(P) = Hnew(NS);
    Qend(P)=Qnew(NS);
    H200(P)=Hnew(200);
    Hmid(P) = Hnew((NS)/2);
    P985(P)=rho0*g*Hnew(985);
    M985(P)=rho0*Qnew(985);
    Pend(P)=rho0*g*Hnew(NS);
    Mend(P)=rho0*Qnew(NS);
    Vend(P) = Qnew(NS)/A;
    Qin(P) = Qnew(1);
    P = P + 1;
    
    Q = Qnew;
    H = Hnew;
    
    ph=rho0*g*H;
    Vh=Q/A;
    
    M_dVdx2(P,1) = (Vh(1)-2*Vh(2)+Vh(3))/(dx^2);
    M_dVdx2(P,2:N) = (Vh(3:N+1)-2*Vh(2:N)+Vh(1:N-1))/(dx^2);
    M_dVdx2(P,1) = (Vh(N+1)-2*Vh(N)+Vh(N-1))/(dx^2);
    M_dPdx(P,1) = (ph(2)-ph(1))/dx;
    M_dPdx(P,2:N) = (ph(3:N+1)-ph(1:N-1))/(2*dx);
    M_dPdx(P,N+1) = (ph(N+1)-ph(N))/dx;
    
    disp(i)
    t=t+dt;
    
    
%     plot(x,H);
%     grid on
%     title(['Pressure through pipe at t= ', num2str(t,3)]);
%     axis([0 L .5*Hr Hr*1.5]);
%     pause(.00001);
end

e_cpuT=cputime-cpuT;
disp(e_cpuT)
toc

%---------------Variables for 3D simulation--------------------%
MIAB_M985(:,1) = time;
MIAB_M985(:,2) = M985;
MIAB_Mend(:,1) = time;
MIAB_Mend(:,2) = rho0*Qend;
MIAB_P985(:,1) = time;
MIAB_P985(:,2) = P985;
MIAB_Pend(:,1) = time;
MIAB_Pend(:,2) = Pend;


MOC_trans_values(:,1) = time;
MOC_trans_values(:,2) = Hend;

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
plot(time,H200)
set(gca, 'FontName', 'Times New Roman','FontSize',13)
ylabel('Pressure head [m]')
xlabel('Time [s]')
xlim([0 tend])
box off

