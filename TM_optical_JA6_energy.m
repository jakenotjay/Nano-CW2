% TM_Optical_v4
% 20/2/2019
% J. Allam
% Generalised OPTICAL transmission through multiple layers
% with peak detection and calculation of eigenfunctions
% based on Levi Exercise 4.2
% this version scans energy rather than wavelength

clear
clf

%change the following file locations:
input_file='TM_O_layers.dat';
output_transmission='TM_O_trans.txt';
output_simple='TM_O_simple.csv';
output_peaks='TM_O_peaks.txt';
output_waves='TM_O_waves.txt';
output_DiffPhase='TM_O_DiffPhase.txt';

% generalised input from file for arbitrary number of layers
load (input_file);

layer = TM_O_layers(:,1)
L = TM_O_layers(:,2)*1e-9  %distance array (nm)
mu = TM_O_layers(:,3)       %refractive index

mumax=max(mu);
mumin=min(mu);
N=max(layer);			%number of layers

E_min=1.0 ;               %minimum enery (eV)
E_max=2.0;               %maximum energy (eV)		
npoints=300;                 %number of points in energy scan

% constants and parameters
dE=(E_max-E_min)/npoints;		%energy increment (nm)
eye=complex(0.,1.);		%square root of -1
hbar=1.0545715e-34;		%Planck's constant (Js)
c0=2.998e8;             %speed of light (m/s)
echarge=1.6021764e-19;	%electron charge (C)

%write input data to figure as a table
f = figure(1);
set (f,'Position',[200 200 400 300]);
cnames = {'Layer','L(nm)','Index'};
t = uitable('Parent',f,'Data',TM_O_layers,'ColumnName',cnames,'Position',[20 20 360 250]);

%main calculation of transmission parameters against energy
   for j=1:npoints
      E(j)=dE*j+E_min+1e-8;   %add (pi*1.0e-5) to avoid divide by zero
        if j<npoints
          Eminusone(j)=E(j);    %to make an array with the same dimensions as the derivative
        end
     	bigP=[1,0;0,1];         %default value of matrix bigP 
       for i=1:N
         k(i)=E(j)*echarge*mu(i)/hbar/c0;

      end
      for n=1:(N-1)
         fac = k(n+1)/k(n);     
         p(1,1)=0.5*(1+fac)*exp(-eye*k(n)*L(n));
         p(1,2)=0.5*(1-fac)*exp(-eye*k(n)*L(n));
         p(2,1)=0.5*(1-fac)*exp(eye*k(n)*L(n));
         p(2,2)=0.5*(1+fac)*exp(eye*k(n)*L(n));
         bigP=bigP*p;
      end
      Trans(j)=(abs(1/bigP(1,1)))^2;                    %transmission probability
      TransPhase(j)=angle(1/bigP(1,1))*180/pi;          %transmitted phase
      Ref(j)=(abs(bigP(2,1)/bigP(1,1)))^2;              %reflection pro ability
      RefPhase(j)=angle(bigP(2,1)/bigP(1,1))*180/pi;    %reflected phase
   end
Tmin=min(Trans);    % min transmission for graph plotting
DPhase=diff(TransPhase);
DiffPhase=(DPhase/dE)*E(j)/2/pi/c0; %convert to angular frequency CHECK THIS
MaxDiffPhase=max(DiffPhase);
MinDiffPhase=min(abs(DiffPhase));
   
%plot potential, transmission and reflection coefficients
figure(2); 				
% generalised generation of potential-distance for 
% arbitrary number of input layers
mup=[mu';mu'];mup=mup(:);
dx=1e-12;				%small distance increment used in potential plot
Lx(1)=L(1);
x1=L;
x2=L;
x1(1)=0;
x2(1)=Lx(1)-dx;
for i=2:N
   for j=2:i
      Lx(i)=L(j)+Lx(j-1);				%distance, x
   end
   x1(j)= Lx(j-1)
   x2(j) = Lx(j)- dx
end
x3=[x1';x2'];
xp=x3(:)*1e9;
maxL=Lx(N)*1e9;

subplot(3,3,1),plot(xp,mup),axis([0,Lx(N)*1e9,0.9*mumin,1.2*mumax]);
xlabel('Position, x (nm)'),ylabel('Refractive index');

subplot(3,3,2),plot(Trans,E),axis([0,2,E_min,E_max]);
xlabel('Transmission coefficient'),ylabel('energy (eV)');

subplot(3,3,3),plot(log10(Trans),E),axis([log10(Tmin),0,E_min,E_max]);
%subplot(3,3,3),plot((Trans),lamda),axis([1e-5,1,700,900]);
xlabel('log10(trans. coeff.)'),ylabel('energy (eV)');

subplot(3,3,4),plot(TransPhase,E),axis([-180,180,E_min,E_max]);
xlabel('TransPhase'),ylabel('energy (eV)');


subplot(3,3,6),plot(Ref,E),axis([0,1,E_min,E_max]);
xlabel('ref. coeff.'),ylabel('energy (eV)');

subplot(3,3,7),plot(RefPhase,E),axis([-180,180,E_min,E_max]);
xlabel('RefPhase'),ylabel('energy (eV)');

figure (4);
plot(E,Ref),axis([E_min,E_max,-0.1,1.2]),xlabel('energy (eV)'),ylabel('ref. coeff.');


%output trans to file 
fileid=fopen(output_transmission,'w');
fprintf(fileid, '%s\t %s\t %s\t %s\t %s\n', 'energy','Trans','TransPhase','Ref','RefPhase');
fprintf(fileid, '%f %e %f %f %f \r\n', [E; Trans; TransPhase; Ref; RefPhase]);
fclose(fileid);


%output simple to file 
fileid=fopen(output_simple,'w');
fprintf(fileid, '%s\t, %s\n', 'energy (eV)','Reflectance');
fprintf(fileid, '%f, %f \r\n', [E; Ref]);
fclose(fileid);


%output differentiated phase to file 
fileid=fopen(output_DiffPhase,'w');
fprintf(fileid, '%s\t %s\n', 'energy','DiffPhase');
fprintf(fileid, '%f %e \r\n', [Eminusone; DiffPhase]);
fclose(fileid);

%find peaks IN REFLECTANCE NOT TRANSMISSION
npeak=0;
for j=2:npoints-1
    Ref(j),E(j)
    if Ref(j)>Ref(j-1)
        if Ref(j)>Ref(j+1)
            npeak=npeak+1;
            Epeak(npeak)=E(j);
            Rpeak(npeak)=Ref(j);
        end
    end
end

%output peaks to file 
fileid=fopen(output_peaks,'w');
fprintf(fileid, '%s\t %s\n', 'Epeak','Rpeak');
fprintf(fileid, '%f %f \r\n', [Epeak; Rpeak]);
fclose(fileid);

%find wavefunctions at peak energies VVV
%prepared for subplots depending on number of peaks detected;
figure(3);
peakplotrows=1;
if npeak>4
        peakplotrows=2;
end
if npeak>6
        peakplotrows=3;
end
if npeak>9
        peakplotrows=4;
end
peakplotcols=round(npeak/peakplotrows+0.5);

%first find complex transmission at peaks for initial waveform on RHS

%output waves to file 
fileid=fopen(output_waves,'w');
fprintf(fileid, '%s\t %s\n', 'position (nm)','wave');

for j=1:npeak
    otherE(j)=Epeak(j);
     	bigP=[1,0;0,1];	%default value of matrix bigP
       for i=1:N
         k(i)=2*pi*mu(i)/E(j)/1e-9;
       end
      for n=1:(N-1)
         fac = k(n+1)/k(n);
         p(1,1)=0.5*(1+fac)*exp(-eye*k(n)*L(n));
         p(1,2)=0.5*(1-fac)*exp(-eye*k(n)*L(n));
         p(2,1)=0.5*(1-fac)*exp(eye*k(n)*L(n));
         p(2,2)=0.5*(1+fac)*exp(eye*k(n)*L(n));
         bigP=bigP*p;
      end

   % initial waveform on RHS layer N;
        A(N)=bigP(1,1); 
        B(N)=0;
        % subsequent waveforms from sucessive matrices;
       for nn=1:(N-1)
           n=N-nn;
           fac = k(n+1)/k(n);
           p(1,1)=0.5*(1+fac)*exp(-eye*k(n)*L(n));
           p(1,2)=0.5*(1-fac)*exp(-eye*k(n)*L(n));
           p(2,1)=0.5*(1-fac)*exp(eye*k(n)*L(n));
           p(2,2)=0.5*(1+fac)*exp(eye*k(n)*L(n));
           A(n)= p(1,1)*A(n+1)+p(1,2)*B(n+1);
           B(n)= p(2,1)*A(n+1)+p(2,2)*B(n+1);
       end
      
layer=1;
xw(1)=0;
wf(1)=0;
xx=0;
for i = 2:(N*100)
    xw(i) = xw(i-1)+maxL/100/N;
    xx=xx+maxL/N/100;
    wf(i) = real(A(layer)*exp(eye*k(layer)*xx*1e-9)+B(layer)*exp(-eye*k(layer)*xx*1e-9));
    if xx > L(layer)*1e9;
        layer=layer+1;
        xx=0;
    end
end


%output successive waves to same file
fprintf(fileid, '%f %f \r\n', [xw; wf]);

%plot each wave as subplot
wfmin=min(wf);
wfmax=max(wf);
subplot(peakplotrows,peakplotcols,j),plot(xw,wf),axis([0,Lx(N)*1e9,wfmin,wfmax]);
xlabel('position (nm)'),ylabel('wavefunction');
ttl = sprintf('peak %3.0f of %3.0f, otherE=%3.1f nm, Ref=%1.5f',j,npeak,Epeak(j),Rpeak(j));
title (ttl);
end  

figure (4);
plot(E,Ref),axis([E_min,E_max,-0.1,1.2]),xlabel('energy (eV)'),ylabel('ref. coeff.');

fclose(fileid);
  

  