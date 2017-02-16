function r = gaussManualRegimeRule(par)
   if (par(2) > par(1)) 
      r=1;
   else
      r=2;
   end
end