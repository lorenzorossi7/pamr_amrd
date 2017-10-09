!============================================================================= 
!  operator types for fortran files, mirroring those in ../include/pamr.h
!============================================================================= 

!============================================================================= 
! fine-to-coarse operators
! VC & CC 
!============================================================================= 
        integer PAMR_NO_INJECT 
        parameter ( PAMR_NO_INJECT=0 )
!============================================================================= 
! VC 
!============================================================================= 
        integer PAMR_STRAIGHT_INJECT,PAMR_HW_RESTR,PAMR_FW_RESTR
        parameter ( PAMR_STRAIGHT_INJECT=1 )
        parameter ( PAMR_HW_RESTR=2 )
        parameter ( PAMR_FW_RESTR=3 )
!============================================================================= 
! CC 
!============================================================================= 
        integer PAMR_NN_AVERAGE,PAMR_NN_ADD
        parameter ( PAMR_NN_AVERAGE=4 )
	parameter  ( PAMR_NN_ADD=5 )

!============================================================================= 
! coarse-to-fine operators
! VC & CC 
!============================================================================= 
        integer PAMR_NO_INTERP
        parameter ( PAMR_NO_INTERP=0 )
!============================================================================= 
! VC & CC
!============================================================================= 
        integer PAMR_FIRST_ORDER,PAMR_SECOND_ORDER,PAMR_FOURTH_ORDER
        parameter ( PAMR_SECOND_ORDER=2 )
        parameter ( PAMR_FOURTH_ORDER=4 )
!============================================================================= 
! CC
!============================================================================= 
        integer PAMR_SECOND_ORDER_ENO,PAMR_FOURTH_ORDER_ENO
        integer PAMR_FIRST_ORDER_CONS, PAMR_MC
	integer PAMR_FIRST_ORDER_EXTENSIVE
        parameter ( PAMR_SECOND_ORDER_ENO=1 )
        parameter ( PAMR_FOURTH_ORDER_ENO=3 )
        parameter ( PAMR_FIRST_ORDER_CONS=5 )
        parameter ( PAMR_MC=6 )
	parameter ( PAMR_FIRST_ORDER_EXTENSIVE=7 )

!============================================================================= 
! post-regrid transfer bits; for values > 1 uses interpolation types */
! VC & CC 
!============================================================================= 
        integer PAMR_NO_TRANSFER
        parameter ( PAMR_NO_TRANSFER=0 )

!============================================================================= 
! physical boundary types 
! VC & CC 
!============================================================================= 
        integer PAMR_UNKNOWN,PAMR_EVEN,PAMR_ODD
        parameter ( PAMR_UNKNOWN=0 )
        parameter ( PAMR_EVEN=1 )
        parameter ( PAMR_ODD=2 )

!============================================================================= 
! V->C transfer
! VC 
!============================================================================= 
        integer PAMR_V_TO_C_NO_TRANSFER,PAMR_V_TO_C_SECOND_ORDER
        integer PAMR_V_TO_C_FOURTH_ORDER
        parameter ( PAMR_V_TO_C_NO_TRANSFER=0 )
        parameter ( PAMR_V_TO_C_SECOND_ORDER=1 )
        parameter ( PAMR_V_TO_C_FOURTH_ORDER=2 )

!============================================================================= 
! C->V transfer
! CC 
!============================================================================= 
        integer PAMR_C_TO_V_NO_TRANSFER,PAMR_C_TO_V_SECOND_ORDER
        integer PAMR_C_TO_V_SECOND_ORDER_ENO
        integer PAMR_C_TO_V_FOURTH_ORDER_ENO
        parameter ( PAMR_C_TO_V_NO_TRANSFER=0 )
        parameter ( PAMR_C_TO_V_SECOND_ORDER=1 )
        parameter ( PAMR_C_TO_V_SECOND_ORDER_ENO=3 )
        parameter ( PAMR_C_TO_V_FOURTH_ORDER_ENO=4 )
