% Potts Monte Carlo method,
% see, 

ii = 0;
cc =colormap(jet);
 
Q = [1:1000]; % possible grain orientations

NX = 100;
NY = 100;

[X,Y] = meshgrid(NX,NY);

for i=1:NX
for j=1:NY
    
    Z(i,j) = Q(randi(length(Q)));
    
end
end

for t=1:10000000
    
   % pick a random lattice site
   iRand = randi(NX-2)+1;
   jRand = randi(NY-2)+1;
    
   % analyse the neighbours
   q = [Z(iRand+1,jRand), ...
        Z(iRand-1,jRand), ...
        Z(iRand,jRand+1), ...
        Z(iRand,jRand-1), ...
        Z(iRand+1,jRand+1), ...
        Z(iRand-1,jRand-1), ...
        Z(iRand+1,jRand-1), ...
        Z(iRand-1,jRand+1)];
    qq = q;
    
   % calculate its current energy
   E = 0;
   for i=1:8
      if (Z(iRand,jRand) == q(i))
            % do nothing
      else
            % grain boundary costs energy
            E = E + 1;
      end
   end
    
   % make it unique
   q = unique(q);
   if (size(q,2) == 1)
       continue;
   end
   
   % remove the like indices
   q(find(q==Z(i,j)))=[];
   
   % randomly assign one of the unlike neighbours
   newZ = q(randi(length(q)));
    
   % calculate the energy change
   newE = 0;
   for i=1:8
      if (newZ == qq(i))
            % do nothing
      else
            % grain boundary costs energy
            newE = newE + 1;
      end
   end   
   dE = newE - E;
   
   % calculate the probability of acquiring this state
   P = 0.5*(1-tanh(dE));
   pp = zeros(1,100)+1;
   pp(ceil(P*100):end) = 0;
   changeTakesPlace = pp(randi(length(pp)));
   
   if (changeTakesPlace == 1)
       Z(iRand,jRand) = newZ;

    % apply periodic bc
    Z(1,:) = Z(NX-1,:);
    Z(NX,:) = Z(2,:);
    Z(:,1) = Z(:,NY-1);
    Z(:,NY) = Z(:,2);
   else
       % do nothing
   end
   
   if (mod(t,10000)==0)
    surf(Z),axis equal, view(2); shading flat
    pause(0.01)
   end
end
