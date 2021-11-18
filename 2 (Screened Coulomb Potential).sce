//SOLVE Schrodinger equation for H-atom in Screened Coulomb potential (using Finite Difference Method)

clc; clear; clf;    //Clearing console, variables and fig

//Describing Constants : h, e, m, a1, a2, a3
h=1973; e=3.795; m=0.511e6; a1=3; a2=5; a3=7;

rmin=0.01; rmax=10; n=1000;
r=linspace(rmin,rmax,n);    //linspace = linearly spaced vector 
d=r(2)-r(1);    //Incremental Step Size

//Defining Potential Energy Matrix
V1=zeros(n,n);
for i=1:n
    V1(i,i)=-(e^2)*exp(-r(i)/a1)/r(i);
end
V2=zeros(n,n);
for i=1:n
    V2(i,i)=-(e^2)*exp(-r(i)/a2)/r(i);
end
V3=zeros(n,n);
for i=1:n
    V3(i,i)=-(e^2)*exp(-r(i)/a3)/r(i);
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

disp("Grounded State Energy (in eV) for a=3A, 5A and 7A",[E1(1) E2(1) E3(1)],"1st Excited State Energy (in eV)for a=3A, 5A and 7A",[E1(2) E2(2) E3(2)])

//Plotting Probability Densities at GS State & 1st Excited State
subplot(3,1,1)
plot(r',[U1(:,1)**2,U1(:,2)**2],"linewidth",3)
legend("Ground State","1st Excited State",1)
xlabel("r","fontsize",2);
ylabel("Probability Density","fontsize",2);
title("a=3A")

subplot(3,1,2)
plot(r',[U2(:,1)**2,U2(:,2)**2],"linewidth",3)
legend("Ground State","1st Excited State",1)
xlabel("r","fontsize",2);
ylabel("Probability Density","fontsize",2);
title("a=5A")

subplot(3,1,3)
plot(r',[U3(:,1)**2,U3(:,2)**2],"linewidth",3)
legend("Ground State","1st Excited State",1)
xlabel("r","fontsize",2);
ylabel("Probability Density","fontsize",2);
title("a=7A")

//OUTPUT : 
//"Grounded State Energy (in eV) for a=3A, 5A and 7A"
//
//-9.386  -10.95  -11.67
//
//"1st Excited State Energy (in eV)for a=3A, 5A and 7A"
//
//-0.483  -1.272  -1.747
