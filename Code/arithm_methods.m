clc;
#{
-Solves the system system in Ex1a using the Euler Method
-returns solved x y and psi
#}
function setSystem()
  global x y psi dx dy dpsi dx0 dy0 dpsi0 x0 y0 psi0 n
  #reset all vectors
  dx = zeros(1,n);
  dy = zeros(1,n);
  dpsi = zeros(1,n);
  
  x = zeros(1,n);
  y = zeros(1,n);
  psi = zeros(1,n);
  
  #starting values
  dx(1) = dx0;
  dy(1) = dy0;
  dpsi(1) = dpsi0;
  x(1) = x0;
  y(1) = y0;
  psi(1) = psi0;
endfunction

#Calculate x, y, psi using Euler
function [x,y,psi] = System1Euler()
  global h n_z m_z D_psi m m_a f_x D_x f_y D_y dx dy dpsi x y psi n
  setSystem();
  for i=1:1:n;
    #derivatives at step i
    dx(i+1) =  dx(i) + h*((sin(psi(i))*dpsi(i)*dx(i) - cos(psi(i))*dpsi(i)*dy(i)) + ((f_x - D_x* abs(dx(i))*dx(i))/(m+3*m_a)))/cos(psi(i));
    dy(i+1) = dy(i) + h*(cos(psi(i))*dpsi(i)*dx(i) + sin(psi(i))*dpsi(i)*dy(i) + ((f_y - D_y* abs(dy(i))*dy(i))/(m+3*m_a)))/cos(psi(i));
    dpsi(i+1) = dpsi(i) + h*(n_z/m_z - (D_psi*abs(dpsi(i))*dpsi(i))/m_z);
    #functions at step i
    x(i+1) = x(i) + h*dx(i);
    y(i+1) = y(i) + h*dy(i);
    psi(i+1) = psi(i) + h*dpsi(i);
  endfor;
endfunction;

#{
-Solves the system system in Ex1c using the Euler Method
-returns solved x y and psi
#}
function [x,y,psi] = System2Euler()
  global h m_z D_psi m m_a D_x D_y dx dy dpsi x y psi Kp_x Kp_y Kp_psi Kd_x Kd_y Kd_psi x_des y_des psi_des n
  
  setSystem();  
  for i=1:1:n;
    #part2 forces
    f_x=Kp_x*(x_des-x(i))-Kd_x*(dx(i));
    f_y=Kp_y*(y_des-y(i))-Kd_y*(dy(i));
    n_z=Kp_psi*(psi_des-psi(i))-Kd_psi*(dpsi(i));
    #derivatives at step i
    dx(i+1) =  dx(i) + h*((sin(psi(i))*dpsi(i)*dx(i) - cos(psi(i))*dpsi(i)*dy(i)) + ((f_x - D_x* abs(dx(i))*dx(i))/(m+3*m_a)))/cos(psi(i));
    dy(i+1) = dy(i) + h*(cos(psi(i))*dpsi(i)*dx(i) + sin(psi(i))*dpsi(i)*dy(i) + ((f_y - D_y* abs(dy(i))*dy(i))/(m+3*m_a)))/cos(psi(i));
    dpsi(i+1) = dpsi(i) + h*(n_z/m_z - (D_psi*abs(dpsi(i))*dpsi(i))/m_z);
    #functions at step i
    x(i+1) = x(i) + h*dx(i);
    y(i+1) = y(i) + h*dy(i);
    psi(i+1) = psi(i) + h*dpsi(i);
  endfor;
endfunction;

#{
-Solves the system system in Ex1a using the Modified Euler Method
-returns solved x y and psi
#}
function [x,y,psi] = System1ModEuler()
  global h n_z m_z D_psi m m_a f_x D_x f_y D_y psi0 dx dy dpsi x y psi n
  setSystem();
  for i=1:n
    #second derivatives at step i
    ddy = ((cos(psi(i))*dpsi(i)*dx(i) + sin(psi(i))*dpsi(i)*dy(i)) + ((f_y - D_y* abs(dy(i))*dy(i))/(m+3*m_a)))/cos(psi(i));
    ddx = ((sin(psi(i))*dpsi(i)*dx(i) - cos(psi(i))*dpsi(i)*dy(i)) + ((f_x - D_x* abs(dx(i))*dx(i))/(m+3*m_a)))/cos(psi(i));
    #second derivatives at step i
    temp_psi = dpsi(i) + (h/(2*m_z))*(n_z - D_psi*abs(dpsi(i))*dpsi(i));
    temp_dy = dy(i)+(h/2)*ddy;
    temp_dx = dx(i)+(h/2)*ddx;
    #derivatives at step i
    dpsi(i+1) = dpsi(i) + (h/m_z)*(n_z -(D_psi* abs(temp_psi))*(temp_psi));
    dy(i+1) =dy(i) + h*(cos(psi(i))*dpsi(i)*dx(i) + sin(psi(i))*dpsi(i)*temp_dy + ((f_y - D_y* abs(temp_dy)*temp_dy)/(m+3*m_a)))/cos(psi(i));
    dx(i+1) =  dx(i) + h*((sin(psi(i))*dpsi(i)*temp_dx - cos(psi(i))*dpsi(i)*dy(i)) + ((f_x - D_x* abs(temp_dx)*temp_dx)/(m+3*m_a)))/cos(psi(i));
    #functions at step i
    psi(i+1) = psi(i)+ h*dpsi(i);
    y(i+1) = y(i)+ h*dy(i);
    x(i+1) = x(i) + h*dx(i);
  endfor;
endfunction

#{
-Solves the system system in Ex1c using the Modified Euler Method
-returns solved x y and psi
#}
function [x,y,psi] = System2ModEuler()
  global h m_z D_psi m m_a D_x D_y dx dy dpsi x y psi Kp_x Kp_y Kp_psi Kd_x Kd_y Kd_psi x_des y_des psi_des n
  setSystem();
  for i=1:n
    #new forces at step i
    f_x=Kp_x*(x_des-x(i)+0.5*dx(i))-Kd_x*(dx(i));
    f_y=Kp_y*(y_des-y(i)+0.5*dy(i))-Kd_y*(dy(i));
    n_z=Kp_psi*(psi_des-psi(i)+0.5*dpsi(i))-Kd_psi*(dpsi(i));
    #second derivatives at step i
    ddy = ((cos(psi(i))*dpsi(i)*dx(i) + sin(psi(i))*dpsi(i)*dy(i)) + ((f_y - D_y* abs(dy(i))*dy(i))/(m+3*m_a)))/cos(psi(i));
    ddx = ((sin(psi(i))*dpsi(i)*dx(i) - cos(psi(i))*dpsi(i)*dy(i)) + ((f_x - D_x* abs(dx(i))*dx(i))/(m+3*m_a)))/cos(psi(i));
    #replacements at step i
    temp_psi = dpsi(i)/m_z + (h/(2*m_z))*(n_z - D_psi*abs(dpsi(i))*dpsi(i));
    temp_dy = dy(i)+(h/2)*ddy;
    temp_dx = dx(i)+(h/2)*ddx;
    #derivatives at step i
    dpsi(i+1) = dpsi(i) + (h/m_z)*(n_z -(D_psi* abs(temp_psi))*(temp_psi));    
    dy(i+1) =dy(i) + h*(cos(psi(i))*dpsi(i)*dx(i) + sin(psi(i))*dpsi(i)*temp_dy + ((f_y - D_y* abs(temp_dy)*temp_dy)/(m+3*m_a)))/cos(psi(i));
    dx(i+1) =  dx(i) + h*((sin(psi(i))*dpsi(i)*temp_dx - cos(psi(i))*dpsi(i)*dy(i)) + ((f_x - D_x* abs(temp_dx)*temp_dx)/(m+3*m_a)))/cos(psi(i));
    #functions at step i
    psi(i+1) = psi(i)+ h*dpsi(i);
    y(i+1) = y(i)+ h*dy(i);
    x(i+1) = x(i) + h*dx(i);
  endfor;
endfunction

function Problem2Euler()
  global Kp_psi psi_des psi n_z Kd_psi dpsi h m_z n D_psi
  setSystem();
  for i=1:n
    n_z= Kp_psi*(psi_des-psi(i)) - Kd_psi*(dpsi(i));
    dpsi(i+1) = dpsi(i) + h*(n_z/m_z - (D_psi*dpsi(i))/m_z);
    psi(i+1) = psi(i) + h*dpsi(i);
  endfor;
endfunction;

function Problem2ModEuler()
  global Kp_psi psi_des psi n_z Kd_psi dpsi h m_z n D_psi
  setSystem();
  for i=1:n
    n_z=Kp_psi*(psi_des-psi(i))-Kd_psi*(dpsi(i));
    temp_psi = dpsi(i) + (h/(2*m_z))*(n_z - D_psi*dpsi(i));
    dpsi(i+1) = dpsi(i) + (h/m_z)*(n_z -(D_psi*temp_psi));
    psi(i+1) = psi(i)+ h*dpsi(i);
  endfor;
endfunction;

#{
-arguments(1) the name of the window (optional)
-Plots the current x y and psi
#}
function plotfuncs(name='')
  global t x y psi
  a = figure('name', name);
  #plot psi
  subplot(3,1,1);
  plot(t,psi);
  title('\Psi');
  #plot y
  subplot(3,1,2);
  plot(t,y);
  title('y');
  #plot x
  subplot(3,1,3);
  plot(t,x);
  title('x');
endfunction

#{
-writes current x y psi to file
#}
function writetofile(filename = 'testtext.txt')
  global t x y psi
  P=[t' x' y' psi'];
  save("-ascii", filename, "P");
endfunction