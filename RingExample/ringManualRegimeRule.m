function r = ringManualRegimeRule(par)
   if (par(1) < 0 && par(2) < 0) 
      r=1;
   elseif (par(1) < 0 && par(2) >= 0)
      r=2;
   elseif (par(1) >= 0 && par(2) >= 0)
      r=3;
   else
      r=4;
   end
end