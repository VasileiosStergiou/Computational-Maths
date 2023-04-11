#Initial Values of constants and initial conditions 
global AM=(4039+4145+4300)/3;
global m = 425000;
global m_z = 357000000;
global m_a = 113000;

global psi0 = 0;
global x0 = AM/1000;
global y0 = 0;

global D_x = 11835;
global D_y = 11835 - AM;
global D_psi = 11835 + AM;
global h = 0.1;

global dy0 = 0;
global dx0 = 0;
global dpsi0 = 0;

global t = [0:0.1:600];
global n = t(length(t))/h;

global f_x = AM;
global f_y = 0;
global n_z = 0;

global dx = zeros(1,n);
global dy = zeros(1,n);
global dpsi = zeros(1,n);
global x = zeros(1,n);
global y = zeros(1,n);
global psi = zeros(1,n);

global Kp_x=60000+5*AM;
global Kd_x=5000000;
global Kp_y=60000;
global Kd_y=5000000-100*AM;
global Kp_psi=50000;
global Kd_psi=7000000;
global x_des=AM/100;
global y_des=-AM/100;
global psi_des=AM/10000;