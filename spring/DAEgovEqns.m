function [A,F] = DAEgovEqns()
syms x1(t) y1(t) th1(t) x2(t) y2(t) th2(t) x3(t) y3(t) th3(t)...
     dx1(t) dy1(t) dth1(t) dx2(t) dy2(t) dth2(t) dx3(t) dy3(t) dth3(t)...   
     ddx1(t) ddy1(t) ddth1(t) ddx2(t) ddy2(t) ddth2(t) ddx3(t) ddy3(t) ddth3(t)...
     Rax Ray Rcx Rcy Rex Rey Rfx Rfy...
     Lae beta1 AC0x AC0y B0x B0y...
     Am w0 t...
     Ls1 Ls2 k1 k2 Fs2 dap y1_0 th1_0...
     m1 m2 m3 g...
     L1 L2 L3 
% variable corresponding to trianle link 3
L32 = L3/(2*cos(beta1));
a = (1/3)*L32*sin(beta1);
b = (1/2)*L3;
c = norm([a b]);
d = L32*sin(beta1);
psi = atan(a/b);
phi = th3+pi()-beta1;
sigma = beta1-psi;

% Need convert all vectors to formula version
% vectors for amb3
    r_d_a = [L1*cos(th1) L1*sin(th1) 0];            r_d_a = formula(r_d_a);
    r_c1_a = [L1/2*cos(th1) L1/2*sin(th1) 0];       r_c1_a = formula(r_c1_a);
    r_e_a = [Lae*cos(th1) Lae*sin(th1) 0 ];         r_e_a = formula(r_e_a);
% vectors of amb6
    r_c2_c = [L2/2*cos(th2) L2/2*sin(th2) 0];       r_c2_c = formula(r_c2_c);
    r_f_c = [L2*cos(th2) L2*sin(th2) 0];            r_f_c = formula(r_f_c);
%vectors for amb9
    r_g_e = [-L3*cos(th3) -L3*sin(th3) 0];          r_g_e = formula(r_g_e);
    r_c3_e = [c*cos(phi+sigma) c*sin(phi+sigma) 0]; r_c3_e = formula(r_c3_e);
    r_f_e = [L32*cos(phi) L32*sin(phi) 0];          r_f_e = formula(r_f_e);
% other vectors 
    r_f_c3 = [d*2/3*cos(th3+pi()/2) d*2/3*sin(th3+pi()/2) 0];   r_f_c3 = formula(r_f_c3);
    r_g_c3 = [c*cos(th3+pi()+psi) c*sin(th3+pi()+psi) 0];   r_g_c3 = formula(r_g_c3);



    
% Spring number 1 forces
h = (3/2)*c;
eta = th3+pi()+psi;
r_g_b = [x3+2*h/3*cos(eta)-B0x, y3+2*h/3*sin(eta)-B0y, 0];
r_g_b = formula(r_g_b);
r_g_b_mag = sqrt(r_g_b(1)^2+r_g_b(2)^2+r_g_b(3)^2);
Fs1 = (r_g_b_mag-Ls1)*k1*(r_g_b/r_g_b_mag);
Fs1 = formula(Fs1);
% Spring number 2 forces
% Dealt within state vector function ODE
% o_yprime = y1_0+L1/2*sin(th1_0)-Ls2+Am*sin(w0*t);
% r_d_op = [0 o_yprime-y1 0];
% r_d_op = formula(r_d_op);
% %Fs2 = [0 (r_d_op(2)-Ls2)*k2 0];
% Fs2 = [0 (r_d_op(2)-Ls2)*k2 0];
% Fs2 = formula(Fs2);


%damper forces
v3 = [dx3 dy3 0]; v3 = formula(v3);
F_d = dap*(v3 + r_g_c3*dth3);
F_d = formula(F_d);
% Moments of inertias
I1 = 1/3*m1*L1^2;
I2 = 1/3*m2*L2^2;
I3 = 1/4*L3*h^3;
%define LMB and AMB equations (1-9)
eqns = [];
eqns = [eqns; m1*ddx1 == Rax+Rex];
eqns= [eqns; m1*ddy1 == Ray+Rey+Fs2-m1*g];
    
amb3 = [0 0 I1*ddth1] == cross(r_d_a,[0 Fs2 0]) + cross(r_c1_a,[0 -m1*g 0]) + cross(r_e_a, [Rex Rey 0]);
amb3f = formula(amb3);
eqns = [eqns; amb3f(3)];
eqns = [eqns; m2*ddx2 == Rcx+Rfx];
eqns = [eqns; m2*ddy2 == -m2*g+Rcy+Rfy];
    
amb6 = [0 0 I2*ddth2] == cross(r_c2_c,[0 -m2*g 0])+cross(r_f_c,[Rfx Rfy 0]);
amb6f = formula(amb6);
eqns = [eqns; amb6f(3)];
eqns = [eqns; m3*ddx3 == Fs1(1)+Rfx+Rex+F_d(1)];
eqns = [eqns; m3*ddy3 == Fs1(2)+Rfy+Rey-m3*g+F_d(2)];
   
amb9 = [0 0 I3*ddth3] == cross(r_g_e, Fs1)+cross(r_c3_e,[0 -m3*g 0])+cross(r_f_e,[Rfx Rfy 0])+cross(r_g_c3,F_d);
amb9f = formula(amb9);
eqns = [eqns; amb9f(3)];
%define constraint equations (10-17)
%link 1
% r_0 = r_0_c1 - r_c1_a
l1i = 0 == x1-r_c1_a(1);
l1j = 0 == y1-r_c1_a(2);
%link 2
% r_c2_0 = r_e_a + r_f_e + r_c2_F
l2i = x2 == r_e_a(1)+r_f_e(1)-L2/2*cos(th2);
l2j = y2 == r_e_a(2)+r_f_e(2)-L2/2*sin(th2);
% l2i = x2 == AC0x + L2/2*cos(th2);
% l2j = y2 == AC0y + L2/2*sin(th2);
% %link 3
%r_c3_0 = r_e_a + r_c3_e

l3i = x3 == r_e_a(1)+r_c3_e(1);
l3j = y3 == r_e_a(2)+r_c3_e(2);
%total link constraint
% ltoti = AC0x == Lae*cos(th1)+L32*cos(phi)-L2*cos(th2);
% ltotj = AC0y == Lae*sin(th1)+L32*sin(phi)-L2*sin(th2);

% r_A_C = r_e_a + r_c3_e + r_f_c3 + r_c_f
ltoti = AC0x == r_e_a(1)+r_c3_e(1)+r_f_c3(1)-r_f_c(1);
ltotj = AC0y == r_e_a(2)+r_c3_e(2)+r_f_c3(2)-r_f_c(2);
%take derivatives

eqns = [eqns; diff(diff(l1i,t),t)];
eqns = [eqns; diff(diff(l1j,t),t)];

eqns = [eqns; diff(diff(l2i,t),t)];
eqns = [eqns; diff(diff(l2j,t),t)];

eqns = [eqns; diff(diff(l3i,t),t)];
eqns = [eqns; diff(diff(l3j,t),t)];

eqns = [eqns; diff(diff(ltoti,t),t)];
eqns = [eqns; diff(diff(ltotj,t),t)];

%replace 1st and 2nd time derivatives with symbolic variable
M = formula(eqns);
Msub = subs(M,diff(x1,t,t),ddx1);
Msub = subs(Msub,diff(x1,t),dx1);
Msub = subs(Msub,diff(y1,t,t),ddy1);
Msub = subs(Msub,diff(y1,t),dy1);
Msub = subs(Msub,diff(th1,t,t),ddth1);
Msub = subs(Msub,diff(th1,t),dth1);
Msub = subs(Msub,diff(x2,t,t),ddx2);
Msub = subs(Msub,diff(x2,t),dx2);
Msub = subs(Msub,diff(y2,t,t),ddy2);
Msub = subs(Msub,diff(y2,t),dy2);
Msub = subs(Msub,diff(th2,t,t),ddth2);
Msub = subs(Msub,diff(th2,t),dth2);
Msub = subs(Msub,diff(x3,t,t),ddx3);
Msub = subs(Msub,diff(x3,t),dx3);
Msub = subs(Msub,diff(y3,t,t),ddy3);
Msub = subs(Msub,diff(y3,t),dy3);
Msub = subs(Msub,diff(th3,t,t),ddth3);
Msub = subs(Msub,diff(th3,t),dth3);
eqns = Msub;

vars = [ddx1(t); ddy1(t); ddth1(t); ddx2(t); ddy2(t); ddth2(t); ddx3(t); ddy3(t); ddth3(t);...
    Rax; Ray; Rcx; Rcy; Rex; Rey; Rfx; Rfy];

syms ddx1_ ddy1_ ddth1_ ddx2_ ddy2_ ddth2_ ddx3_ ddy3_ ddth3_...
    Rax_ Ray_ Rcx_ Rcy_ Rex_ Rey_ Rfx_ Rfy_

vars_ = [ddx1_; ddy1_; ddth1_; ddx2_; ddy2_; ddth2_; ddx3_; ddy3_; ddth3_;...
    Rax_; Ray_; Rcx_; Rcy_; Rex_; Rey_; Rfx_; Rfy_]; 

[A,F] = equationsToMatrix(subs(eqns,vars,vars_),vars_);

syms x1_ y1_ th1_ x2_ y2_ th2_ x3_ y3_ th3_...
    dx1_ dy1_ dth1_ dx2_ dy2_ dth2_ dx3_ dy3_ dth3_

vs2 = [x1(t); y1(t); th1(t); x2(t); y2(t); th2(t); x3(t); y3(t); th3(t);...
    dx1(t); dy1(t); dth1(t); dx2(t); dy2(t); dth2(t); dx3(t); dy3(t); dth3(t)];

vs2_ = [x1_; y1_; th1_; x2_; y2_; th2_; x3_; y3_; th3_;...
    dx1_; dy1_; dth1_; dx2_; dy2_; dth2_; dx3_; dy3_; dth3_];

A = subs(A,vs2,vs2_);
F = subs(F,vs2,vs2_);
end