function [z] = supercoiling(F,L0,dLk,offset,length_correction)       % define force, length (nm) and the array of linking number (f.e from -50 to 50), z-offset, hidden length due to off-center bead attachment



kbT = 4.1;                                              % kBT - Boltzmann constant in room temperature (pN*nm);                                                          
C = 100;                                                % DNA twist persistence length (nm)
Cp = 24;                                                % twist persistence length of plectonemic DNA (nm)
At = 50;                                                % DNA stretch persistence length (nm)
Ct = C * (1 - C ./ (4 .* At) .* sqrt(kbT ./ (At.*F)))   % twist persistence length as a function of force with constant C = 100 nm


g = (F - sqrt(F*kbT/At));                               % free energy per nm of torsionally unconstrained DNA;
tau = Ct .* (2 .* pi .* dLk ./ L0) .* kbT;              % DNA restoring torque 
tau_buck = sqrt(2 .* kbT .* Cp .* g ./ (1 - Cp/Ct))     % DNA buckling torque
tau_melt = -11;                                         % DNA melting torque
Lk_buck = tau_buck .* L0 /(2 .*pi .* kbT .* Ct)         % maximal linking number that can be absorbed by DNA before buckling
Lk_melt = tau_melt .* L0 /(2 .*pi .* kbT .* Ct)         % maximal linking number that can be absorbed by DNA before melting

%% DNA extension as a function of linking number

% twisted DNA
z_t =  (1 - 0.5 .* sqrt(kbT/(At .* F))) - (C.^2 ./ 16) .* (2 .* pi .* dLk ./ L0).^2 .* ( kbT ./ (At*F)) .^ (3/2);   % extension of twisted DNA (z/L)

% DNA buckling
z_slope = (2 .* pi .* (1 - 0.5 .* sqrt(kbT/(At .* F)) - C.^2/(16 .* Ct.^2) * (kbT ./(At.*F)).^(3/2) .* (tau_buck./kbT).^2))/((tau_buck./kbT).*(1./Cp - 1./Ct))  % slope that described DNA buckling (nm/turn)

G_t = - g + (0.5 * kbT .* Ct) .* (2 .* pi .* dLk ./ L0).^2;                 % Gt / L0
dG_DNA = 2 .* pi .* tau_buck;                                               % free energy of DNA gained per each turn;


%% conditions for DNA buckling at positive twist


for i = 1:length(dLk)
    if dLk(i) >= Lk_buck && i>1
        z_t(i) = z_t(i-1);                                                  % once buckling torque is reached, the change in extension is attributed only to the plectoneme formation
        z(i) = (z_t(i).* L0  - (dLk(i)-Lk_buck) .* z_slope) ./ L0;          % keep in mind that the extension of twisted DNA is expressed in z/L
    
    else  
        z(i) = z_t(i);
    end
    
end

%% DNA buckling and melting at negative twist

  
k = 1;
start = dLk(k);

while start ~= 0                                                                
    k = k+1;
    start = dLk(k);
    k;                                                                      % the position of point in the turns array that starts to go up from negative to positive twist
end
                                                                            

for i = 1:(k-2)
    
  if -tau_buck > tau_melt                                                   % DNA buckles when melting torque is not reached
      
      if dLk(k-i) <= - Lk_buck 
            z(k-i-1) = z(k-i);   
            z(k-i) = (z(k-i).* L0  - (abs(dLk(k-i))-Lk_buck) .* z_slope) ./ L0; 
            
      else
            z(k-i) = z(k-i);
                
      end
      
  else
      
      if tau(k-i) <= tau_melt  
            z(k-i) = z(k-i+1);
      else
            z(k-i) = z(k-i);
      end
      
  end
  
end
 z(1) = z(2);
                                                                 
%% erasing negative extensions

for i = 1:length(dLk)                                                       
    if z(i) < length_correction+0.06 ./L0     %0.06 to correct for the fact that the border conditions do not include noise in the measurement
        z(i) = length_correction+0.06 ./L0;
    else
        z(i) = z(i);
    end
end

z = z + offset;
%% plottting 

%figure(1)

sigma = dLk/(L0/0.34/10.4);
%plot(sigma,z,'o')
xlim ([-0.1 0.1]) 
%hold on;


