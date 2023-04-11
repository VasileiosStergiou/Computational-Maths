clear all;
InitialValues;
arithm_methods;
Psi_s_denom = [m_z (Kd_psi+D_psi) Kp_psi 0];

poles = roots(Psi_s_denom);

a = figure('name', 'Poles');
stem(real(poles),imag(poles),'-.');
ylabel('Imaginary Part');
xlabel('Real Part');

l = real(poles(1));
m = imag(poles(1));

t0=0;
c = AM/10000;
c1 = -c;
c2 = c*l/m;

psi_new = c1*exp(l*t).*cos(m*t) + c2*exp(l*t).*sin(m*t) + c;

global D_psi=11835;

b = figure('name', 'new psi');
plot(t,psi_new);
ylabel('\Psi');
xlabel('t');

Problem2Euler();
b = figure('name', "\psi");
subplot(2,1,1);
plot(t,psi);
ylabel('\Psi');
xlabel('t');

Problem2ModEuler();
subplot(2,1,2)
plot(t,psi);
ylabel('\Psi');
xlabel('t');

rs = [0];
imgs = [0];
test = [10 100 1000 10000 100222 11857777];
for k1=1:1:length(test)
  for k2=1:1:length(test)
    r1= roots([m_z (test(k1)+D_psi) test(k2) 0]);
    rs(length(rs)+1) = real(r1(1));
    rs(length(rs)+1) = real(r1(2));
    rs(length(rs)+1) = real(r1(3));
    
    imgs(length(imgs)+1) = imag(r1(1));
    imgs(length(imgs)+1) = imag(r1(2));
    imgs(length(imgs)+1) = imag(r1(3));
  endfor;
endfor;
figure('name', 'Poles');
stem(rs,imgs,'-.');
ylabel('Imaginary Part');
xlabel('Real Part');