    integer NITER_PKW,NITER_AEC_AQ
    integer NITER_POA,NITER_AEC_DRY    
    double precision ::  ALFHP! percentage of H+ allowed to c/e(0.1)
    double precision ::  EPSER
    double precision ::  TINYM,TINYN,MTSBL
    double precision ::  DTMAX,DTAEROMIN
    !double precision ::  Diam_MIN,Diam_MAX
    parameter(EPSER = 1.D-2) ! Relative error for time step adjustment
    parameter(ALFHP = 0.1D0)
    parameter(TINYM = 0.)
    parameter(TINYN = 1.d-6)
    parameter(DTMAX = 500.D0)
    parameter(NITER_AEC_AQ = 1)
    parameter(NITER_AEC_DRY = 1)
    parameter(NITER_PKW = 5)
    parameter(NITER_POA = 10)
    parameter(MTSBL = 1.0D0)
