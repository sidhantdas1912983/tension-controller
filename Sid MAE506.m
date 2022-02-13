clear all

Kf= 0.1;
KL = 0.1;
Rl=0.5;
m = 0.511;
J = 0.03435;
r = 1.5;
KA = 27.2;
KB = 27.2;

%openloop
A = [0 0 1 0;
    0 0 0 1;
    -(KA+KB)/m -(KA-KB)/m 0 0;
    -(r^2)*(KA-KB)/J -(r^2)*(KA+KB)/J 0 0];
B = [0 ; 0 ;(KL/(Rl*m))+ KA/Kf*m; KA/Kf*J ];
C = [KB -r*KB 0 0].*1e-1;
D = [0];
sys_ol = ss(A,B,C,D);
x0 = [-1; -1; 0; 0];
t = 0:0.01:15;
u = 1*ones(length(t),1);
VA=1*ones(length(t),1);
y_ol=lsim(sys_ol,u,t,x0);
%openloop_information
Ob=obsv(sys_ol);
rank(Ob);
Ct=ctrb(sys_ol);
rank(Ct);
[wn_ol,zeta_ol,p_ol] = damp(sys_ol);
tf(sys_ol);

% closed loop
%K= place(A,B,[-1,-2,-3,-4.0])
K= place(A,B,[-1+j,-1-j,(-7.522/(99.52*2))-61.5872j,(-7.522/(99.52*2))+61.5872j]);
G= -4/(C*inv(A-B*K)*B);
sys_cl= ss(A-B*K , B*G , C-D*G , D*G);
t2 = 0:0.01:200;
u =1*ones(length(t2),1);
y_cl=lsim(sys_cl,u,t2,x0);
%closed loop information
s_cl=stepinfo(sys_cl);
[wn_cl,zeta_cl,p_cl] = damp(sys_cl);
tf(sys_cl);

%LQR_Cheap
Q= [100 0 0 0 ; 0 100 0 0; 0 0 100 0; 0 0 0 100];
R = 1;
K_lqr = lqr(A,B,Q,R);
G_lqr = -4/(C*inv(A-B*K_lqr)*B);
sys_lqr= ss(A-B*K_lqr,B*G_lqr,C-D*K_lqr,0);
t=0:0.01:15;
%ref= zeros(length(t),1);
ref= 1*ones(length(t),1);
y_LQR_cheap = lsim(sys_lqr,ref,t,x0);
%LQR information
s_lqr=stepinfo(sys_lqr);
[wn_lqr,zeta_lqr,p_lqr] = damp(sys_lqr);
tf(sys_lqr);

%LQR_Expensive
Q= [1 0 0 0 ; 0 1 0 0; 0 0 1 0; 0 0 0 1];
R = 100;
K_lqr = lqr(A,B,Q,R);
G_lqr = -4/(C*inv(A-B*K_lqr)*B);
sys_lqr= ss(A-B*K_lqr,B*G_lqr,C-D*K_lqr,0);
t=0:0.01:15;
%ref= zeros(length(t),1);
ref= 1*ones(length(t),1);
y_LQR_exp = lsim(sys_lqr,ref,t,x0);
%LQR information
s_lqr=stepinfo(sys_lqr);
[wn_lqr,zeta_lqr,p_lqr] = damp(sys_lqr);
tf(sys_lqr);

figure
subplot(4,1,1);
plot(t,y_LQR_cheap,'b','LineWidth',1)
xlabel('Time (s)');
ylabel('Tension(gf)');
grid on
subplot(4,1,2);
plot(t,y_LQR_exp,'b','LineWidth',1)
xlabel('Time (s)');
ylabel('Tension(gf)');
grid on
subplot(4,1,3);
plot(t2,y_cl,'gr','LineWidth',1)
xlabel('Time (s)');
ylabel('Tension(gf)');
grid on
subplot(4,1,4);
plot(t,y_ol,'r','LineWidth',1)
xlabel('Time (s)');
ylabel('Tension(gf)');
grid on