function Xdot = fourBarLinkSusp(t,X,S,A,F)
     [y1_0, L1, th1_0, Ls2, Am, w0,...
      k2_r] = v2struct(S);
    x1_ = X(1); y1_ = X(2); th1_ = X(3); dx1_ = X(10); dy1_ = X(11); dth1_ = X(12);
    x2_ = X(4); y2_ = X(5); th2_ = X(6); dx2_ = X(13); dy2_ = X(14); dth2_ = X(15);
    x3_ = X(7); y3_ = X(8); th3_ = X(9); dx3_ = X(16); dy3_ = X(17); dth3_ = X(18);
    
    o_yprime = y1_0+L1/2*sin(th1_0)-Ls2+Am*sin(w0*t);
    r_d_op = [0 o_yprime-(y1_+L1/2*sin(th1_)) 0];
    if abs(r_d_op(2))-Ls2 >= 0
        k2 = 0;
    else 
        k2 = k2_r;
    end
    Fs2 = -(abs(r_d_op(2))-Ls2)*k2;
    
    A = subs(A);
    A = eval(A);
    F = subs(F);
    F = eval(F);
    u = A\F;
    vdot = u(1:9);
    rdot = X(10:18);
   
    Xdot = [rdot;vdot];
    global Forces
    Forces = [Forces;u(10:end)',Fs2,r_d_op(2)];
end