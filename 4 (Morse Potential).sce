//SOLVE Schrodinger equation for H-atom in Morse potential (using Finite Difference Method)

clc; clear; clf;    //Clearing console, variables and fig

//Describing Constants
h=1973; D=0.755501; m=940*10^6; a=1.44; ro=0.131349;

rmin=0.01; rmax=10; n=1000;
r=linspace(rmin,rmax,n);    //linspace = linearly spaced vector 
d=r(2)-r(1);    //Incremental Step Size

//Defining Morse Potential using given formula
V=zeros(n,n);
for i=1:n
    rp=(r(i)-ro)/r(i);
    V(i,i)=D*(exp(-2*a*rp)-exp(-a*rp));
end

//Defining Kinetic Energy using given formula
K=eye(n,n)*(-2);
for i=1:(n-1)
    K(i,i+1)=1;
    K(i+1,i)=1;
end

//Defining Hamiltonian Matrix using given formula
H=(-(h^2)/(2*m*d^2))*K+V;

//Evaluating Eigenvalues & Eigenvectors of H matrix using "spec" function
[U,EV]=spec(H); //U=Eigenvectors(used to plot Wavefunction) & EV=Eigenvalues(used to find Energies)
E=diag(EV); //Extracting diagonal elements of EV matrix using "diag" function
format(6);   //changing number format

disp("Grounded State Energy (in eV)",E(1),"1st Excited State Energy (in eV)",E(2));

//Plotting Wavefunctions at GS State & 1st Excited State
plot(r',[abs(U(:,1))**2,abs(U(:,2))**2],"linewidth",3);
legend("Ground State","1st Excited State",1);
xlabel("r","fontsize",5);
ylabel("Probability Density","fontsize",5);
xgrid(5);

//Finding Bohr Radius (using SI unit)
E= -13.6; //Energy required to separate electron and proton, eV
e= 1.6*(10^(-19)); //charge of an electron, C
E= E*e; //converting to J
Po= 8.85*(10^(-12)); //Permittivity of free space, F/m
r= e^2/(8*(%pi)*Po*E); //radius, m
r=-r;
disp("Bohr Radius(in m) = ",r);

//NOTE : U here is not the radial function. Instead U/r is the radial function

//SECOND PART OF PROGRAM to find Dissociation Energy & plot potential
//Defining Morse Potential using given formula
rmin=0.001; rmax=10; n=1000;
r=linspace(rmin,rmax,n);    //linspace = linearly spaced vector 
V=zeros(n,1);
for i=1:n
    V(i)=D*( (1-exp(-a*(r(i)-ro))) **(2));
end
scf(2);
plot(r',V,"linewidth",3);
legend("Morse Potential Graph",2);
xlabel("Internuclear Distance (r)","fontsize",5);
ylabel("Morse Potential (V(r))","fontsize",5);
xgrid(5);

//D is the well depth (defined relative to the dissociated atoms)
Dissociation_Energy = V(n)-V(1)
disp("Dissociation Energy = ",Dissociation_Energy)

//OUTPUT
//
//  "Grounded State Energy (in eV)"
//
//  -0.155
//
//  "1st Excited State Energy (in eV)"
//
//  -0.143
//
//  "Bohr Radius(in m) = "
//
//   5.D-11
//
//  "Dissociation Energy = "
//
//   0.723

