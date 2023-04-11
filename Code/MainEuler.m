clear all;
InitialValues;
arithm_methods;
#Ex1 Euler
#1st calculation
[x,y,psi] = System1Euler(n);
#1st plot
plotfuncs('Part 1a Euler');
#2nd initial conditions
f_x = 0;
f_y = -AM;
n_z = 0;
#2st calculation
[x,y,psi] = System1Euler(n);
#2st plot
plotfuncs('Part 1b Euler');
#3rd initial conditions
f_x = 0;
f_y = 0;
n_z = -AM;
#3st calculation
[x,y,psi] = System1Euler(n);
#3rd plot
plotfuncs('Part 1c Euler');
#Ex1 c Euler
#change some initial values
x0=0;
y0=-AM/1000;
psi0=0;
D_x=11835+AM;
D_y=11835;
D_psi=11835;
#make the calculation with variable forces
[x,y,psi] = System2Euler();
#new system plot
plotfuncs('Part 1d Euler');