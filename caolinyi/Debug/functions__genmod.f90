        !COMPILER-GENERATED INTERFACE MODULE: Tue Feb 09 23:07:03 2016
        MODULE FUNCTIONS__genmod
          INTERFACE 
            SUBROUTINE FUNCTIONS(IROW,XC,GRID,PP,QX,QDX,N,A,ALPHA,BETA, &
     &MU,B,KAPPA,DELTA,XCR,FX,FDX)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IROW
              REAL(KIND=8) :: XC(N)
              REAL(KIND=8) :: GRID(N)
              REAL(KIND=8) :: PP(N,N)
              REAL(KIND=8) :: QX(N)
              REAL(KIND=8) :: QDX(N)
              REAL(KIND=8) :: A
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: BETA
              REAL(KIND=8) :: MU
              REAL(KIND=8) :: B
              REAL(KIND=8) :: KAPPA
              REAL(KIND=8) :: DELTA
              REAL(KIND=8) :: XCR
              REAL(KIND=8) :: FX(N)
              REAL(KIND=8) :: FDX(N)
            END SUBROUTINE FUNCTIONS
          END INTERFACE 
        END MODULE FUNCTIONS__genmod
