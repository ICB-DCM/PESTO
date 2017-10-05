#include "fintrf.h"
C ======================================================================
C     
C     sincall.F is an example for illustrating how to use mexCallMATLAB.
C     
C     It creates an mxArray and passes its associated pointers to 
C     subfunction fill() to get data filled up. It then calls 
C     mexCallMATLAB to calculate sin function and plot the result.
C
C     NOTE: The subfunction fill() is in the file called fill.F.
C     To build this MEX-file you need to also compile fill.F, 
C     for example:
C
C        mex sincall.F fill.F
C
C     Copyright 1984-2009 The MathWorks, Inc.
C     
C ======================================================================
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

C     Function declarations:
      integer mexCallMATLAB
      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix

C     Pointers to input/output mxArrays:
      mwPointer input(1), output(1)

C     Array information:
      mwSize m, n, max
      
      integer status
      mwPointer NULL
      parameter (NULL = 0)
C-----------------------------------------------------------------------
C     Initialize variables.
      m=1
      n=1
      max=1000

C     Create the matrix to be passed into MATLAB.
      input(1) = mxCreateDoubleMatrix(max, n, 0)
 
C     Pass the pointer and variable and let fill() fill up data.
      call fill(%VAL(mxGetPr(input(1))), m, n, max)
      call mxSetM(input(1), m)
      call mxSetN(input(1), n)

      status = mexCallMATLAB(1, output, 1, input, 'sin')
      status = mexCallMATLAB(0, NULL, 1, output, 'plot')

C     Cleanup the un-freed memory after calling mexCallMATLAB.
      call mxDestroyArray(input(1))
      call mxDestroyArray(output(1))

      return
      end
