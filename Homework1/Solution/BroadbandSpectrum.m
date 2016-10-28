function [SPL13,fcentre,SPLO,fcentreO,OASPL]=BroadbandSpectrum(SPL,f)

SPL13=zeros(33,1);

%% One-thrid octave band sound spectrum
 fcentre  = (10.0^3) * ((2.0) .^ ([-19:13]./3));
 fd = (2^(1/6));
 fupper = fcentre .* fd;
 flower = fcentre ./ fd;
 
 for i=1:size(fcentre,2)
   sum=0;
   for j=1:size(f)
      if f(j)>=flower(i) && f(j)<=fupper(i)
         sum=sum+10^(SPL(j)/10);
      end
   end
     SPL13(i)=10*log10(sum);
 end

%% Octave band sound spectrum  
SPLO=zeros(11,1);

fcentreO = (10.0^3) * ((2.0) .^ ([-6:4]));
fdO=(2^(1/2));
fupperO = fcentreO .* fdO;
flowerO = fcentreO ./ fdO;

for i=1:size(fcentreO,2)
   sum=0;
   for j=1:size(fcentre,2)
      if fcentre(j)>=flowerO(i) && fcentre(j)<=fupperO(i)
         sum=sum+10^(SPL13(j)/10);
      end
   end
     SPLO(i)=10*log10(sum);
end


%% Overall sound spectrum
sum=0;
 for j=1:size(fcentreO,2)
     sum=sum+10^(SPLO(j)/10);
 end
     OASPL=10*log10(sum);
end   
   


 

 