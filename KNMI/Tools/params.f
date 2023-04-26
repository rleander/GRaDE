! GLOBAL PARAMETERS 
	integer WWIN,NKMAX,NDIM,MAXYR,YR0

        !       PARAMETERS WHICH CAN BE ALTERED TO MODIFY THE SIMULATION
        !------------------------------------------------------------------
        parameter (WWIN=60) ! half-window width (days)        
C        parameter (WWIN=182) ! half-window width (days)        
C                                ! window is chosen symmetrically from        
C                                ! today-winwidth : today+winwidth        

        !       PARAMETERS WHICH SHOULD ONLY BE CHANGED FOR SPECIFIC REASONS 
        !------------------------------------------------------------------
        parameter (NDIM=4)      ! maximum number of dimensions in the feature 
                                ! vector, change ndim in sub 'SORTNNB' 
                                ! accordingly!     
                                ! Note that the actual choice of the feature 
        parameter (MAXYR=100)   ! dimension of the feature-vector array in years
        parameter (YR0=1900)    ! offset of historical years
        !---------------------------------------------------------------------------
        parameter (NKMAX=(2*WWIN+1)*MAXYR)   ! max number of nearest neighbours 


