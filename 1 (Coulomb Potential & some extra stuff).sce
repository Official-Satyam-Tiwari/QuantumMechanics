//SOLVE Schrodinger equation for H-atom in Coulomb's potential (using Finite Difference Method)

clc; clear; clf;    //Clearing console, variables and fig

//Describing Constants : h, e & m
h=1973; e=3.795; m=0.511e6;

rmin=1e-10; rmax=10; n=1000;
r=linspace(rmin,rmax,n);    //linspace = linearly spaced vector 
d=r(2)-r(1);    //Incremental Step Size

//Defining Coulomb Potential using given formula
V=zeros(n,n);
for i=1:n
    V(i,i)=(-(e^2)/r(i));
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
E=diag(EV) //Extracting diagonal elements of EV matrix using "diag" function
format(6)   //changing number format

//Since 1st eigenvalue is absurd, hence we took n=2 as GS and n=3 as 1st Excited state
disp("Grounded State Energy (in eV)",E(2),"1st Excited State Energy (in eV)",E(3))

//Plotting Probability Densities at GS State & 1st Excited State
plot(r',[U(:,2)**2,U(:,3)**2],"linewidth",3)
legend("Ground State","1st Excited State",1)
xset("font size",4);
xlabel("r","fontsize",5);
ylabel("Probability Density","fontsize",5);
xgrid(5);

b= h^2/(e^2*m)
disp("The Bohr radius is (in angstrom) :",b)

scf(2)
for i=1:5:n
        subplot(2,2,1)
        
        plot(r(i),U(i,2),"r*","linewidth",2,"markersize",2)
        plot(r(i),U(i,3),"b-","linewidth",2,"markersize",2)
        xgrid()
        legend("Ground State","First Excited State")
        ylabel("Wave Function u(r)")
        xlabel("r")
        title('Numerical')
        
        subplot(2,2,2)
        
        plot(r(i),U(i,2)**2,"r*","linewidth",2,"markersize",2)
        plot(r(i),U(i,3)**2,"b-","linewidth",2,"markersize",2)
        xgrid()
        legend("Ground State","First Excited State")
        ylabel("Probability Density u(r)*u(r)")
        xlabel("r")
        title('Numerical')
        
        subplot(2,2,3)
        
        x = [0:0.1:10]
        plot((2/b^(3/2))*x.*exp(-x/b),'r') 
        plot( 1/b^(3/2)*(x.*exp(-x/(2*b)) - 1/(2*b)*x^2.*exp(-x/(2*b)))) 
        xgrid()
        legend("Ground State","First Excited State")
        ylabel("Wavefunction u(r)")
        xlabel("r")
        title('Theoretical')
        
        subplot(2,2,4)
         
        plot(((2/b^(3/2))*x.*exp(-x/b))^2,'r')
        plot((1/b^(3/2)*(x.*exp(-x/(2*b)) - 1/(2*b)*x^2.*exp(-x/(2*b))))^2)
        xgrid()
        legend("Ground State","First Excited State")
        ylabel("Probability Density u(r)*u(r)")
        xlabel("r")
        title('Theoretical')
end
