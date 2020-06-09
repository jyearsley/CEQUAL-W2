!*
!*                                                                    **
!*                         CE-QUAL-W2-INC                             **
!*                                                                    **
!*                 The INCLUDE file for CE-QUAL-W2                    **
!*                                                                    **
!*                          Version 2.0                               **
!*                         October, 1995                              **
!*                                                                    **
!*                  Developed by: Thomas M. Cole                      **
!*                  Water Quality Modeling Group                      **
!*                  U.S. Army Corps of Engineers                      **
!*                  Waterways Experiment Station                      **
!*                  Vicksburg, Mississippi 39180                      **
!*                  Phone number: (601) 634-3283                      **
!*                                                                    **
!*

!*                                                                    **
!*   This file contains the PARAMETER specifications for CE-QUAL-W2   **
!*   and the EQUIVALENCEs between the constituent array and the var-  **
!*   iables used in the kinetics subroutines.                         **
!*                                                                    **
!*   The PARAMETERS are defined as follows:                           **
!*                                                                    **
!*     CONFN - Control filename                                       **
!*     IMP   - Total number of segments                               **
!*     IMC   - Total number of segments                               **
!*     KMP   - Total number of layers                                 **
!*     KMC   - Total number of layers                                 **
!*     NBP   - Total number of branches                               **
!*     NSP   - Total number of outlet structures                      **
!*     NTP   - Total number of tributaries                            **
!*     NWP   - Total number of withdrawals                            **
!*     NCP   - Total number of constituents                           **
!*     NDP   - Maximum number of dates for I/O parameters             **
!*             (see manual, page __)                                  **
!*                                                                    **
!*   Notes:                                                           **
!*                                                                    **
!*   1. The variables IMC and KMC may be set to 1 in order to save    **
!*      memory if constituents are not being modeled.  If consti-     **
!*      tuents are modeled, then IMC and KMC should be set equal      **
!*      to IMP and KMP.                                               **
!*   2. Do not change the value of NCP unless more water quality      **
!*      state variables have been added to the model.                 **
!*   3. No parameters may be set to zero.  They must be set to at     **
!*      least 1.                                                      **
!*                                                                    **
!*

!**** PARAMETER declarations

      CHARACTER  CONFN*72
      PARAMETER (CONFN='w2_con.npt')
      PARAMETER (IMP=28,KMP=41,IMC=1,KMC=1,NBP=1,NSP=2,NTP=1,NWP=1,      &
                 NDP=50)
      PARAMETER (NCP=21)
!
!**** Type declarations

      REAL  LPOM,   LDOM,   NH4,   NO3
      REAL  LPOMSS, LDOMSS, NH4SS, NO3SS

!**** DIMENSION declarations

      DIMENSION    CSSK(KMC,IMC,NCP), C1(KMC,IMC,NCP), C2(KMC,IMC,NCP)
      DIMENSION    SS(KMC,IMC),       COL(KMC,IMC),    TDS(KMC,IMC),  &
                   LDOM(KMC,IMC),     RDOM(KMC,IMC),   ALGAE(KMC,IMC),  &
                   LPOM(KMC,IMC),     PO4(KMC,IMC),    NH4(KMC,IMC),  &
                   NO3(KMC,IMC),      DO(KMC,IMC),     SED(KMC,IMC),  &
                   TIC(KMC,IMC),      ALKAL(KMC,IMC),  PH(KMC,IMC),  &
                   CO2(KMC,IMC),      HCO3(KMC,IMC),   CO3(KMC,IMC),  &
                   FE(KMC,IMC),       CBOD(KMC,IMC)
      DIMENSION    SSSS(KMC,IMC),     COLSS(KMC,IMC),  LDOMSS(KMC,IMC),  &
                   RDOMSS(KMC,IMC),   ALGSS(KMC,IMC),  LPOMSS(KMC,IMC),  &
                   PO4SS(KMC,IMC),    NH4SS(KMC,IMC),  NO3SS(KMC,IMC),  &
                   DOSS(KMC,IMC),     TICSS(KMC,IMC),  FESS(KMC,IMC),  &
                   CBODSS(KMC,IMC)

!**** EQUIVALENCE declarations

      EQUIVALENCE (SS(1,1),C2(1,1,2)),       (COL(1,1),C2(1,1,3)),  &
                  (TDS(1,1),C2(1,1,4)),      (LDOM(1,1),C2(1,1,5)),  &
                  (RDOM(1,1),C2(1,1,6)),     (ALGAE(1,1),C2(1,1,7)),  &
                  (LPOM(1,1),C2(1,1,8)),     (PO4(1,1),C2(1,1,9)),  &
                  (NH4(1,1),C2(1,1,10)),     (NO3(1,1),C2(1,1,11)),  &
                  (DO(1,1),C2(1,1,12)),      (SED(1,1),C1(1,1,13)),  &
                  (TIC(1,1),C2(1,1,14)),     (ALKAL(1,1),C2(1,1,15)),  &
                  (PH(1,1),C1(1,1,16)),      (CO2(1,1),C1(1,1,17)),  &
                  (HCO3(1,1),C1(1,1,18)),    (CO3(1,1),C1(1,1,19)),  &
                  (FE(1,1),C2(1,1,20)),      (CBOD(1,1),C2(1,1,21))
      EQUIVALENCE (SSSS(1,1),CSSK(1,1,2)),   (COLSS(1,1),CSSK(1,1,3)),  &
                  (LDOMSS(1,1),CSSK(1,1,5)), (RDOMSS(1,1),CSSK(1,1,6)),  &
                  (ALGSS(1,1),CSSK(1,1,7)),  (LPOMSS(1,1),CSSK(1,1,8)),  &
                  (PO4SS(1,1),CSSK(1,1,9)),  (NH4SS(1,1),CSSK(1,1,10)),  &
                  (NO3SS(1,1),CSSK(1,1,11)), (DOSS(1,1),CSSK(1,1,12)),  &
                  (TICSS(1,1),CSSK(1,1,14)), (FESS(1,1),CSSK(1,1,20)),  &
                  (CBODSS(1,1),CSSK(1,1,21))

!**** COMMON declarations

      COMMON /GCNSTC/ CSSK, C1, C2








