//SOLVE Schrodinger equation for H-atom in Screened Anharmonic Oscillator potential (using Finite Difference Method)

clc; clear; clf;    //Clearing console, variables and fig

//Describing Constants : h, k, m, b1, b2, b3
h=197.3; k=100; m=940; b1=0; b2=10; b3=30;

rmin=0.01; rmax=10; n=1000;
r=linspace(rmin,rmax,n);    //linspace = linearly spaced vector 
d=r(2)-r(1);    //Incremental Step Size

//Defining Potential Energy Matrix
V1=zeros(n,n);
for i=1:n
    V1(i,i)=(k*(r(i)^2))/2 + (b1*(r(i)^3))/3;
end
V2=zeros(n,n);
for i=1:n
    V2(i,i)=(k*(r(i)^2))/2 + (b2*(r(i)^3))/3;
end
V3=zeros(n,n);
for i=1:n
    V3(i,i)=(k*(r(i)^2))/2 + (b3*(r(i)^3))/3;
end

//Defining Kinetic Energy using given formula
K=eye(n,n)*(-2);
for i=1:(n-1)
    K(i,i+1)=1;
    K(i+1,i)=1;
end

//Defining Hamiltonian Matrix using given formula
H1=(-(h^2)/(2*m*d^2))*K+V1;
H2=(-(h^2)/(2*m*d^2))*K+V2;
H3=(-(h^2)/(2*m*d^2))*K+V3;

//Evaluating Eigenvalues & Eigenvectors of H matrix using "spec" function
[U1,EV1]=spec(H1);
[U2,EV2]=spec(H2);
[U3,EV3]=spec(H3);
E1=diag(EV1);
E2=diag(EV2);
E3=diag(EV3);
format(6)   //changing number format

disp("Grounded State Energy (in eV) for b=0, 10 and 30",[E1(1) E2(1) E3(1)],"1st Excited State Energy (in eV)for b=0, 10 and 30",[E1(2) E2(2) E3(2)])

//Plotting Wavefunctions at GS State & 1st Excited State
subplot(3,1,1)
plot(r',[abs(U1(:,1))**2,abs(U1(:,2))**2],"linewidth",3)
legend("Ground State","1st Excited State",1)
xlabel("r","fontsize",2);
ylabel("Probability Density","fontsize",2);
title("b=0","position",[5 0.08])

subplot(3,1,2)
plot(r',[abs(U2(:,1))**2,abs(U2(:,2))**2],"linewidth",3)
legend("Ground State","1st Excited State",1)
xlabel("r","fontsize",2);
ylabel("Probability Density","fontsize",2);
title("b=10","position",[5 0.08])

subplot(3,1,3)
plot(r',[abs(U3(:,1))**2,abs(U3(:,2))**2],"linewidth",3)
legend("Ground State","1st Excited State",1)
xlabel("r","fontsize",2);
ylabel("Probability Density","fontsize",2);
title("b=30","position",[5 0.08])
