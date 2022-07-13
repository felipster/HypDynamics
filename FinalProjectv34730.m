% Work done by Felipe Hanuch!

clear 
close all
clc
%% parameters
%initial angles
th1_0 = -pi/8;
th2_0 = pi/4;
th3_0 = -pi/3;
% linkage masses
m1 = 20;
m2 = 2;
m3 = 3;
g = 9.81;
% linkage lengths
L1 = 4;
Lae = 2;
L2 = 2;
L3 = 1;

%link 3 stuff

beta1 = pi/4;
gamma = pi-2*beta1;
L32 = L3/(2*cos(beta1));
a = (1/3)*L32*sin(beta1);
b = (1/2)*L3;
c = norm([a b]);
d = L32*sin(beta1);
psi = atan(a/b);
phi = th3_0+pi()-beta1;
sigma = beta1-psi;
h = (3/2)*c;

%alpha = th3_0-psi;
eta = th3_0+pi()+psi;

% track profile constants
Am = 0.2;
w0 = 4;
t0 = 0;

% spring 1 and 2 constants
Ls1 = 1;
Ls2 = 1;
k1 = 400;
k2_r = 1000;
% damper constant
dap = 20;

% position of point A
Ap = [0,0];

% position of point B

%initial link points vectors
%   find fixed point of C with initial theta values and lengths of linkages
r_e_a = [Lae*cos(th1_0),Lae*sin(th1_0)];
r_f_e = [L32*cos(phi), L32*sin(phi)];
r_c_f = [-L2*cos(th2_0),-L2*sin(th2_0)];
AC0x= r_e_a(1)+r_f_e(1)+r_c_f(1);
AC0y = r_e_a(2)+r_f_e(2)+r_c_f(2);
AC0 = [AC0x, AC0y];

Cp = Ap+AC0;

%   find initial pose for link 1
r_c1_a = [L1/2*cos(th1_0) L1/2*sin(th1_0) 0];
r01 = [Ap th1_0]+r_c1_a;
x1_0=r01(1);
y1_0 = r01(2);
%   find initial pose for link 2
r_c2_c = [L2/2*cos(th2_0) L2/2*sin(th2_0) 0];
r02 = [Cp th2_0]+r_c2_c;
x2_0=r02(1);
y2_0 = r02(2);
%   find initial pose for link 3
r_c3_e = [c*cos(phi+sigma) c*sin(phi+sigma) 0];
r_e_a = [Lae*cos(th1_0) Lae*sin(th1_0) 0 ];
r03 = [Ap th3_0]+r_e_a+r_c3_e;
x3_0=r03(1);
y3_0 = r03(2);
%   all other relevent vectors:
%r_d_op = [0 Ls2+L1/2*sin(th1_0)-Am*sin(w0*t) 0];
r_d_a = [L1*cos(th1_0) L1*sin(th1_0) 0];
r_f_c = [L2*cos(th2_0) L2*sin(th2_0) 0];
r_g_e = [-L3*cos(th3_0) -L3*sin(th3_0) 0];
r_f_c3 = [d*2/3*cos(th3_0+pi()/2) d*2/3*sin(th3_0+pi()/2) 0];
r_g_c3 = [c*cos(th3_0+pi()+psi) c*sin(th3_0+pi()+psi) 0];
% Make initial spring direction parallel to link2
Bp = Ap+r_e_a(1:2)+r_g_e(1:2)+[-Ls1*cos(th2_0),-Ls1*sin(th2_0)];
B0x = Bp(1);
B0y = Bp(2);
r_g_b = [x3_0+2*h/3*cos(eta)-B0x y3_0+2*h/3*sin(eta)-B0y 0];
% build initial state vector
v0 = [0 0 0 0 0 0 0 0 0]';
X0 = [r01';r02';r03';v0];


%% Plot initial position
figure('Units','centimeters',...
'Position',[12 2 15 15])
hold on
axis equal
% xlim([-1 1])
% ylim([-.1 1])

plot(0,0,'r+')
xlabel('X-axis (m)')
ylabel('Y-Axis (m)')
title('Motorcylce wheel dynamics initial positions')
box on
%% plot swingarm
AD = [Ap 0;r_d_a];
h2 = plot(AD(:,1),AD(:,2),'r','LineWidth',5);
h1 = plot(X0(1),X0(2),'b.','MarkerSize',40);

%% plot the triangle
%cog3

%position vectors
vertsL3 = [];
vertsL3(1,:) = Ap+r_e_a(1:2);
vertsL3(2,:) = vertsL3(1,:)+r_f_e(1:2);
vertsL3(3,:) = vertsL3(1,:)+r_g_e(1:2);
h5 = fill(vertsL3(:,1),vertsL3(:,2),'y');
h4 = plot(X0(7),X0(8),'k.','MarkerSize',20);
R_c3_e = [Ap+r_e_a(1:2); Ap+r_e_a(1:2)+r_c3_e(1:2)];
R_f_c3 = [X0(7:8)'; X0(7:8)'+r_f_c3(1:2)];
R_g_c3 = [X0(7:8)'; X0(7:8)'+r_g_c3(1:2)];
h8 = plot(R_c3_e(:,1),R_c3_e(:,2));
h9 = plot(R_f_c3(:,1),R_f_c3(:,2));
h10 = plot(R_g_c3(:,1),R_g_c3(:,2));
%% plot link 2 CF
%cog2

%plot link CF
CF = [Cp 0];
CF = [CF; [Cp 0]+r_f_c];
h7 = plot(CF(:,1),CF(:,2),'m','Linewidth',3);
h6 = plot(X0(4),X0(5),'c.','MarkerSize',10);
%delete([h1,h2,h4,h5,h6,h7])
hold off
%% find solution
%calculate DAE gov. equations matrix
[A,F] = DAEgovEqns();
S = v2struct(y1_0, L1, th1_0, Ls2, Am, w0,...
      k2_r);
A = subs(A);
F = subs(F);

% Timespan
tfinal = 5;
deltat = 0.001;
t = 0:deltat:tfinal;

% Create global variable to store Forces
global Forces
Forces = [];
fdynamic = @(t,X) fourBarLinkSusp(t,X,S,A,F);

%options = odeset('absTol', 1e-12,'relTol', 1e-12);
options = odeset('RelTol',1e-6);
tic
[t,X] = ode45(fdynamic,t,X0,options);
disp(toc)
x1 = X(:,1);...+p.OAx;
y1 = X(:,2);...+p.OAy;
th1 = X(:,3);
x2 = X(:,4);
y2 = X(:,5);
th2 = X(:,6);
x3 = X(:,7);
y3= X(:,8);
th3 = X(:,9);
xtrack = x1+L1/2*cos(th1);
ytrack = X0(2)+L1/2*sin(th1_0)-Ls2+Am*sin(w0*t);
%% Plot Reaction Forces
Ra_mag = sqrt(Forces(:,1).^2+Forces(:,2).^2);
Rc_mag = sqrt(Forces(:,3).^2+Forces(:,4).^2);
Re_mag = sqrt(Forces(:,5).^2+Forces(:,6).^2);
Rf_mag = sqrt(Forces(:,7).^2+Forces(:,8).^2);
figure4 = figure;
hold on
box on
xlabel('time [s]')
ylabel('Force [N]')
title('Magnitude of Reaction Forces')
t_f = linspace(0,tfinal,size(Forces,1));
plot(t_f,Ra_mag)
plot(t_f,Rc_mag)
plot(t_f,Re_mag)
plot(t_f,Rf_mag)
legend('R_a','R_c','R_e','R_f')
hold off
%% Plot Spring Forces
FFs2 = Forces(:,9);
disp = Forces(:,10);
figure5 = figure;
hold on
box on
xlabel('time [s]')
ylabel('Force [N]')
title('Spring 2 Force')
t_f = linspace(0,tfinal,size(Forces,1));
plot(t_f,FFs2)
legend('Fs2')
hold off

figure6 = figure;
hold on
box on
xlabel('time [s]')
ylabel('displacement [m]')
title('Spring 2 Displacement')
plot(t_f,disp)
legend('r_d_op')
hold off
%% Plot details
figure2 = figure('Units','centimeters',...
'Position',[12 2 15 15]);
hold on
axis equal
xlim([-1 5])
ylim([-3 0.5])
plot(0,0,'r+')
xlabel('X-axis (m)')
ylabel('Y-Axis (m)')
title('Motorcylce wheel dynamics')
box on
idx = 1;

for j = 1:(length(t)-1)/20+1
    i = (j-1)*20+1;
    % time dependent parameters and vectors
    phi = th3(i)+pi()-beta1;
    eta = th3(i)+pi()+psi;
    r_d_a = [L1*cos(th1(i)) L1*sin(th1(i)) 0];
    r_e_a = [Lae*cos(th1(i)) Lae*sin(th1(i)) 0 ];
    r_f_e = [L32*cos(phi), L32*sin(phi)];
    r_g_e = [-L3*cos(th3(i)) -L3*sin(th3(i)) 0];
    r_f_c = [L2*cos(th2(i)) L2*sin(th2(i)) 0];
    r_g_c3 = [c*cos(th3(i)+pi()+psi) c*sin(th3(i)+pi()+psi) 0];
    %plot cog1 swingarm
    
    AD = [Ap 0;r_d_a];
    h2 = plot(AD(:,1),AD(:,2),'r','LineWidth',5);
    h1 = plot(x1(i),y1(i),'b.','MarkerSize',40);
    
    %plot cog3 triangle

    %position vectors
    vertsL3 = [];
    vertsL3(1,:) = Ap+r_e_a(1:2);
    vertsL3(2,:) = vertsL3(1,:)+r_f_e(1:2);
    vertsL3(3,:) = vertsL3(1,:)+r_g_e(1:2);
    h5 = fill(vertsL3(:,1),vertsL3(:,2),'y');
    h4 = plot(x3(i),y3(i),'k.','MarkerSize',40);
    
    %plot cog2 link 2

    CF = [Cp 0];
    CF = [CF; [Cp 0]+r_f_c];
    h7 = plot(CF(:,1),CF(:,2),'m','Linewidth',3);
    h6 = plot(x2(i),y2(i),'c.','MarkerSize',40);
    %plot path
    h8 = plot(x1(1:i),y1(1:i),'m-');
    h9 = plot(x3(1:i),y3(1:i),'c-');
    h10 = plot(x2(1:i),y2(1:i),'b');
    %plot track height
    
    h11 = plot(xtrack(1:i),ytrack(1:i),'g.','MarkerSize',3);
    h3 = plot([-5 5],[ytrack(i) ytrack(i)],'g','LineWidth',3);
    % plot springs
    ro = 0.1;
    [xs1,ys1] = spring(xtrack(i),ytrack(i),Ap(1)+r_d_a(1),Ap(2)+r_d_a(2),5,Ls2,ro); 
    [xs2,ys2] = spring(B0x,B0y,x3(i)+r_g_c3(1),y3(i)+r_g_c3(2),5,Ls2,ro); 
    h12 = plot(xs1,ys1,'b','LineWidth',2);
    h13 = plot(xs2,ys2,'b','LineWidth',2);
    % plot damper portion
    h14 = plot([B0x,x3(i)+r_g_c3(1)],[B0y,y3(i)+r_g_c3(2)],'k','LineWidth',4);
    
    % plot wheel
    h15 = circle(Ap(1)+r_d_a(1),Ap(2)+r_d_a(2),Ls2);%abs(Ap(2)+r_d_a(2)-ytrack(i)));
    
    
    pause(0.1)
    
    frame = getframe(figure2);
    im{idx} = frame2im(frame);
    if i>length(t)-22
        break
    end
    delete([h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14,h15]...,h3,h4])
        )
    idx = idx+1;
end

filename = 'testAnimated.gif'; % Specify the output file name
nImages = idx;
for idx = 1:nImages
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
end

legend('origin','', 'C.O.M 1','','C.O.M 3','','C.O.M 2','path1','path3',...
    'path2','wheel contact path','track','spring1','spring2','damper','wheel')
vertsAframex = [0, B0x,Cp(1),-0.5];
vertsAframey = [0, B0y,Cp(2),Cp(2)];
h16 = fill(vertsAframex,vertsAframey,'k');
function h = circle(x,y,r)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
h = plot(x+xp,y+yp,'k','LineWidth',4);
end