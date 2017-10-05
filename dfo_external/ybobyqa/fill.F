#include "fintrf.h"
C ======================================================================
#if 0
C     
C     fill.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C     fill.f
C     This is the subroutine sincall calls to fill the 
C     mxArray with data.  You can write your own version 
C     to load your data.
C
C     This is a MEX-file for MATLAB.
C     Copyright 1984-2011 The MathWorks, Inc.
C     
C ======================================================================
C     Computational routine - subroutine for filling up data

      subroutine fill(pr, m, n, max)

C     Declarations
      implicit none

      real*8  pr(*)
      mwSize m, n, max, i
     
      m=max/2
      n=1
      do 10 i=1,m
 10      pr(i)=i*(4*3.1415926/max)

      return
      end
