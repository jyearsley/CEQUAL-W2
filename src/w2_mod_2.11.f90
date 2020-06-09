!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                                                                          **
!*                                 CE-QUAL-W2                               **
!*                                                                          **
!*                    A Two-Dimensional, Laterally Averaged,                **
!*                     Hydrodynamic and Water Quality Model                 **
!*                                    for                                   **
!*                   Rivers, Lakes, Reservoirs, and Estuaries               **
!*                                                                          **
!*                               Version 2.11                               **
!*                              January, 1998                               **
!*                                                                          **
!*                               Revisions by                               **
!*                                                                          **
!*                              Thomas M. Cole                              **
!*                         Water Quality Modeling Group                     **
!*                         U.S. Army Corps of Engineers                     **
!*                         Waterways Experiment Station                     **
!*                         Vicksburg, Mississippi 39180                     **
!*                         Phone number: (601) 634-3283                     **
!*                       e-mail tcole@lasher.wes.army.mil                   **
!*                                                                          **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                                                                          **
!*        The long arm of the lawyers has found its way into the water      **
!*        quality modeling arena, so:                                       **
!*                                                                          **
!*        This program is furnished by the U.S. Government and is accepted  **
!*        and used by the recipient with the express understanding that the **
!*        U.S. Government gives no warranties, expressed or implied,        **
!*        concerning the accuracy, reliability, usability, or suitability   **
!*        for any particular purpose of the information and data contained  **
!*        in this program or furnished in connection therewith.  The United **
!*        States of America shall be under no liability whatsoever to any   **
!*        person by reason of any use made thereof.  This program belongs   **
!*        to the U.S. Government, therefore, the recipient further agrees   **
!*        not to assert any proprietary rights therein, or to represent     **
!*        this program to anyone as other than a U.S. Government program.   **
!*        Distribution of this program is restricted by the Export          **
!*        Administration Act of 1979, 50 App. USC 2401-2420, as amended,    **
!*        and other applicable laws on regulations.                         **
!*                                                                          **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

 
!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                            Task 1: Initialization                        **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                             Task 1.1: Block Data                         **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      BLOCK DATA COMMON_DATA
        PARAMETER (NCP=21)
        REAL       NONZERO
        INTEGER    SNP
        CHARACTER  CUNIT1*6,    CNAME2*24,   HNAME*24
        DIMENSION  CUNIT1(NCP), CNAME2(NCP), HNAME(4)
        COMMON /NAMESC/ HNAME,  CNAME2,      CUNIT1
        COMMON /SNPUC/  SNP
        COMMON /NONZC/  NONZERO
        DATA HNAME  /'Horizontal_velocity, m/s',  &
                     ' Vertical_velocity, mm/s',  &
                     '   Temperature, øC   ',  &
                     '   Limiting timestep    '/
        DATA CNAME2 /'    Tracer, g/m^3      ',      &                     !1
                     'Suspended solids, g/m^3',      &                     !2
                     '    Coliform, g/m^3    ',      &                     !3
                     'Dissolved solids, g/m^3',      &                     !4
                     '   Labile DOM, g/m^3   ',      &                     !5
                     ' Refractory DOM, g/m^3 ',      &                     !6
                     '      Algae, g/m^3     ',      &                     !7
                     '   Labile POM, g/m^3   ',      &                     !8
                     '   Phosphorus, mg/m^3  ',      &                  !9
                     '    Ammonium, mg/m^3   ',  &                         !10
                     'Nitrate-Nitrite, mg/m^3',  &                         !11
                     'Dissolved oxygen, g/m^3',    &                       !12
                     '    Sediment, g/m^3    ',    &                       !13
                     'Inorganic carbon, g/m^3',    &                       !14
                     '   Alkalinity, g/m^3   ',    &                       !15
                     '           pH          ',    &                       !16
                     ' Carbon dioxide, g/m^3 ',    &                       !17
                     '   Bicarbonate, g/m^3  ',    &                       !18
                     '    Carbonate, g/m^3   ',    &                       !19
                     '      Iron, g/m^3      ',    &                       !20
                     '      CBOD, g/m^3      '/                         !21
        DATA CUNIT1 /8*'g/m^3', 3*'mg/m^3', 4*'g/m^3', ' ',  &
                     3*'g/m^3',   'mg/m^3',   'g/m^3'/ 
        DATA SNP   /20/, NONZERO /1.0E-20/
      END

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                         Task 1.2: Variable Declarations                  **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      PROGRAM  CE_QUAL_W2
      INCLUDE  'w2.inc'

!**** Type declaration

      REAL    MINDLT, NONZERO
      REAL    LAT,    LONG
      REAL    JDAY,   JDMIN
      REAL    JDAYTS, JDAYPR,  JDAYSP                                   !FORTRAN
      REAL    NH4T1,  NH4T2,   NO3T1,  NO3T2
      REAL    NH4K1,  NH4K2,   NO3K1,  NO3K2
      REAL    LPOMD,  LDOMDK,  LPOMDK, LRDDK
      REAL    NH4DK,  NO3DK
      REAL    NH4D,   NO3D
      REAL    NH4RM,  NO3RM
      REAL    NH4R,   KBOD
      REAL    NXTVD,  NXTMSN,  NXTMTS, NXTMPR, NXTMCP, NXTMVP, NXTMRS,   &
              NXTMSC, NXTMSP,  NXELV1, NXELV2
      REAL    ICETH,  ICESW,   ICETHU, ICETH1, ICETH2, ICEMIN, ICET2,    &
              ICETHI, ICE_TOL
      INTEGER CON,    BTH,     RSI,    RSO,    TSR,    PRF,    VPL,      &
              CPL,    VPR,     SNP,    ERR,    WRN,    SPR,    QWS
      INTEGER TSRDP,  SNPDP,   VPLDP,  CPLDP,  PRFDP,  RSODP,  WSCDP,    &
              DLTDP,  SCRDP,   SPRDP
      INTEGER UHS,    DHS,     US,     DS,     CUS
      INTEGER CN,     CN0,     TRCN,   DTCN,   PRCN,   UHCN,   DHCN
      INTEGER FREQUK, YEAR,    GDAY,   SKTI
      LOGICAL SNAPSHOT,        PROFILE,         TIME_SERIES,             &
              VECTOR,          CONTOUR,         RESTART_IN,              &
              RESTART_OUT,     SPREADSHEET,     SCREEN_OUTPUT,           &
              ONE_LAYER
      LOGICAL TRIBUTARIES,     DIST_TRIBS,      WITHDRAWALS,             &
              EVAPORATION,     PRECIPITATION
      LOGICAL UP_FLOW,         DN_FLOW,         UP_HEAD,                 &
              DN_HEAD,         UH_INTERNAL,     DH_INTERNAL,             & 
              UH_EXTERNAL,     DH_EXTERNAL,     HEAD_BOUNDARY
      LOGICAL ISO_TEMP,        VERT_TEMP,       LONG_TEMP,               &
              ISO_CONC,        VERT_CONC,       LONG_CONC
      LOGICAL ADD_LAYER,       SUB_LAYER
      LOGICAL OPEN_FILES,      OPEN_VPR,        OPEN_LPR
      LOGICAL FRESH_WATER,     SALT_WATER,      SUSP_SOLIDS,             &
              CONSTITUENTS,    LIMITING_FACTOR
      LOGICAL NO_WIND,         NO_INFLOW,       NO_OUTFLOW,              &
              NO_HEAT
      LOGICAL WINTER,          ICE_IN,          ALLOW_ICE,               &
              ICE,             ICE_CALC,        DETAILED_ICE,            &
              TERM_BY_TERM
      LOGICAL INTERP_INFLOW,   INTERP_OUTFLOW,  INTERP_MET,              &
              INTERP_DTRIBS,   INTERP_TRIBS,    INTERP_HEAD,             &
              INTERP_WITHDRWL, INTERPOLATE
      LOGICAL VOLUME_BALANCE,  ENERGY_BALANCE,  MASS_BALANCE,            &
              WATER_BALANCE
      LOGICAL PLACE_QIN
      LOGICAL PLACE_QTR,       SPECIFY_QTR
      LOGICAL WARNING_OPEN,    VOLUME_WARNING
      LOGICAL SEL_WITHDRAWAL,  POINT_SINK
      LOGICAL TRANSPORT,       OXYGEN_DEMAND
      LOGICAL UPDATE_KINETICS, UPDATE_RATES
      LOGICAL END_RUN,         BRANCH_FOUND,    UPWIND,                  &
              SHIFT_DEMAND,    LEAP_YEAR,       LIMITING_DLT       
      LOGICAL LASERJET_II,     LASERJET_III,    LASERJET_IV
      CHARACTER*1  EXT1,   SEG1,  ESC
      CHARACTER*2  EXT2,   SEG2,  DEG,   CUNIT3
      CHARACTER*3  CPRC,   HPRC,  GDCH,  EXT3,  SCRC
      CHARACTER*3  ACC,    INACC, TRACC, DTACC, PRACC
      CHARACTER*3  SNPC,   VPLC,  TSRC,  PRFC,  CPLC,  RSOC, RSIC,       &
                   SPRC
      CHARACTER*3  INFIC,  OUTIC, TRIC,  DTRIC, WDIC,  HDIC, METIC
      CHARACTER*3  CCC,    PQC,   VBC,   EVC,   ICEC,  PRC,  LIMC,       &
                   DTRC,   QINC,  QOUTC, WINDC, HEATC, SWC,  SDC,        &
                   MBC,    EBC,   WBC
      CHARACTER*3  MON,    DUM
      CHARACTER*4  LFAC,   EXT4
      CHARACTER*5  LJPC,   WTYPE, EXT5
      CHARACTER*6  CUNIT1, CUNIT2
      CHARACTER*6  SEGMENT
      CHARACTER*7  BLANK,  CONV,  LFPR
      CHARACTER*8  SLTR,   SINK,  SLICE, SLHEAT, TRC
      CHARACTER*8  CTIME   
      CHARACTER*10  CDATE                                                !FORTRAN
      CHARACTER*9  MONTH
      CHARACTER*9  CONV1,  BLANK1
      CHARACTER*12 RSOFN
      CHARACTER*17 CNAME1, DUM1
      CHARACTER*24 HNAME,  CNAME2
      CHARACTER*72 RSIFN,  VPRFN, TSRFN, PRFFN, VPLFN, CPLFN, METFN,     &
                   QINFN,  TINFN, CINFN, QOTFN, QWDFN, QTRFN, TTRFN,     &
                   CTRFN,  QDTFN, TDTFN, CDTFN, PREFN, TPRFN, CPRFN,     &
                   EUHFN,  TUHFN, CUHFN, EDHFN, TDHFN, CDHFN, SNPFN,     &
                   LPRFN,  BTHFN, SPRFN, ELOFN, TITLE
      DOUBLE PRECISION VOLSBR, VOLTBR, VOLEV, VOLPR,  VOLTR,  VOLDT,     &
                       VOLWD,  VOLUH,  VOLDH, VOLIN,  VOLOUT, DLVOL
      DOUBLE PRECISION TSSEV,  TSSPR,  TSSTR, TSSDT,  TSSWD,  TSSUH,     &
                       TSSDH,  TSSIN,  TSSOUT,TSSS,   TSSB,   TSSICE,    &
                       ESBR,   ETBR,   EIBR
      DOUBLE PRECISION CMBRS,  CMBRT
      DOUBLE PRECISION BHRHO,  GMA,    BTA,   A,      C,      D,         &
                       V,      F
      DOUBLE PRECISION Z,      SZ,     ZMIN
      DOUBLE PRECISION BTAT,   GMAT,   AT,    CT,     DT,     VT
      DOUBLE PRECISION TADL,   TADV,   CADL,  CADV
      DOUBLE PRECISION AD1L,   AD2L,   AD3L,  AD1V,   AD2V,   AD3V
      DOUBLE PRECISION DX1,    DX2,    DX3,   ALFA
      DOUBLE PRECISION T1L,    T2L,    T3L,   C1L,    C2L,    C3L

!**** Dimension statements

      DIMENSION ADMX(KMP,IMP), AZ(KMP,IMP),  ADMZ(KMP,IMP), DZ(KMP,IMP), &
                DX(KMP,IMP),   DM(KMP,IMP),  RHO(KMP,IMP),  P(KMP,IMP),  &
                HPG(KMP,IMP),  HDG(KMP,IMP), SB(KMP,IMP),   ST(KMP,IMP), &
                DZQ(KMP,IMP)
      DIMENSION U(KMP,IMP),    SU(KMP,IMP),  W(KMP,IMP),    SW(KMP,IMP), &
                SAZ(KMP,IMP)
      DIMENSION T1(KMP,IMP),   T2(KMP,IMP),  TSS(KMP,IMP),  QSS(KMP,IMP)
      DIMENSION B(KMP,IMP),    BB(KMP,IMP),  BR(KMP,IMP),   BH(KMP,IMP), &
                BHR(KMP,IMP),  HSEG(KMP,IMP)
      DIMENSION CONV(KMP,IMP), CONV1(KMP,IMP)
      DIMENSION NDLT(KMP,IMP)
      DIMENSION FPSS(KMC,IMC), FPFE(KMC,IMC)
      DIMENSION HKT1(IMP),     HKT2(IMP),    AVHKT(IMP),    BKT(IMP),    &
                BHKT1(IMP),    BHKT2(IMP),   BHRKT1(IMP),   BHRKT2(IMP), &
                DLX(IMP),      DLXR(IMP),    ELWS(IMP),     ELWS2(IMP)
      DIMENSION KTI(IMP),      SKTI(IMP),    KB(IMP),       KBMIN(IMP)
      DIMENSION BHRHO(IMP),    A(IMP),       C(IMP),        D(IMP),      &
                F(IMP),        V(IMP),       BTA(IMP),      GMA(IMP),    &
                Z(IMP),        SZ(IMP),      DLXRHO(IMP)
      DIMENSION EV(IMP),       QDT(IMP),     QPR(IMP)
      DIMENSION Q(IMP),        QC(IMP),      QSSUM(IMP)
      DIMENSION ICETH(IMP),    ICE(IMP),     ICESW(IMP)
      DIMENSION SOD(IMP),      SODS(IMP)
      DIMENSION RN(IMP),       PHI0(IMP)
      DIMENSION IPRF(IMP),     ISPR(IMP),    SEGMENT(IMP)
      DIMENSION H(KMP),        AVH(KMP),     EL(KMP)
      DIMENSION DEPTHB(KMP),   DEPTHM(KMP)
      DIMENSION TVP(KMP)
      DIMENSION QWD(NWP),      TWD(NWP),     IWD(NWP),      KWD(NWP),    &
                JBWD(NWP)
      DIMENSION QTR(NTP),      TTR(NTP),     JBTR(NTP),     ITR(NTP),    &
                KTTR(NTP),     KBTR(NTP),    ELTRT(NTP),    ELTRB(NTP),  &
                TRC(NTP)
      DIMENSION QTRFN(NTP),    TTRFN(NTP),   CTRFN(NTP)
      DIMENSION C2I(NCP),      AVCOUT(NCP),  CPRC(NCP)
      DIMENSION CNAME1(NCP),   CNAME2(NCP),  CUNIT1(NCP),   CUNIT2(NCP), &
                CUNIT3(NCP)
      DIMENSION CN(NCP),       ACC(NCP),     INACC(NCP),    TRACC(NCP),  &
                DTACC(NCP),    PRACC(NCP),   CN0(NCP)
      DIMENSION INCN(NCP),     TRCN(NCP),    DTCN(NCP),     PRCN(NCP),   &
                UHCN(NCP),     DHCN(NCP)
      DIMENSION DTRC(NBP),     SWC(NBP)
      DIMENSION QSUM(NBP),     NOUT(NBP)
      DIMENSION KTQIN(NBP),    KBQIN(NBP)
      DIMENSION ELUH(NBP),     ELDH(NBP),    JBUH(NBP),     JBDH(NBP)
      DIMENSION US(NBP),       CUS(NBP),     DS(NBP),       UHS(NBP),    &
                DHS(NBP),      NL(NBP)
      DIMENSION QIN(NBP),      PR(NBP),      QPRBR(NBP),    QDTR(NBP),   &
                EVBR(NBP)
      DIMENSION TIN(NBP),      TPR(NBP),     TDTR(NBP)
      DIMENSION VOLEV(NBP),    VOLPR(NBP),   VOLTR(NBP),    VOLDT(NBP),  &
                VOLUH(NBP),    VOLDH(NBP),   VOLIN(NBP),    VOLOUT(NBP), &
                VOLWD(NBP),    VOLSBR(NBP),  VOLTBR(NBP),   DLVOL(NBP)
      DIMENSION QINFN(NBP),    TINFN(NBP),   CINFN(NBP),    QDTFN(NBP),  &
               TDTFN(NBP),    CDTFN(NBP),   PREFN(NBP),    TPRFN(NBP),   &
               CPRFN(NBP),    EUHFN(NBP),   TUHFN(NBP),    CUHFN(NBP),   &
                EDHFN(NBP),    TDHFN(NBP),   CDHFN(NBP),    QOTFN(NBP)
      DIMENSION TSSEV(NBP),    TSSPR(NBP),   TSSTR(NBP),    TSSDT(NBP),  &
                TSSWD(NBP),    TSSUH(NBP),   TSSDH(NBP),    TSSIN(NBP),  &
                TSSOUT(NBP),   TSSS(NBP),    TSSB(NBP),     TSSICE(NBP), &
                ESBR(NBP),     ETBR(NBP),    EIBR(NBP)
      DIMENSION UP_FLOW(NBP),      DN_FLOW(NBP),     UP_HEAD(NBP),       &
                DN_HEAD(NBP),      UH_EXTERNAL(NBP), DH_EXTERNAL(NBP),   &
                UH_INTERNAL(NBP),  DH_INTERNAL(NBP), DIST_TRIBS(NBP)
      DIMENSION NSTR(NBP),         KBSW(NSP,NBP)
      DIMENSION QSTR(NSP,NBP),     ESTR(NSP,NBP),    WSTR(NSP,NBP),      &
                SINK(NSP,NBP)
      DIMENSION CIN(NCP,NBP),      CDTR(NCP,NBP),    CPR(NCP,NBP)
      DIMENSION CTR(NCP,NTP),      CWD(NCP,NWP)
      DIMENSION CUH(KMC,NCP,NBP),  CDH(KMC,NCP,NBP)
      DIMENSION CMBRS(NCP,NBP),    CMBRT(NCP,NBP)
      DIMENSION A1(KMC,IMC),       A2(KMC,IMC),      A3(KMC,IMC)
      DIMENSION AGR(KMC,IMC),      ARR(KMC,IMC),     AMR(KMC,IMC),       &
                AER(KMC,IMC)
      DIMENSION LPOMD(KMC,IMC),    SEDD(KMC,IMC),    DOMD(KMC,IMC),      &
                SODD(KMC,IMC),     NH4D(KMC,IMC),    NO3D(KMC,IMC),      &
                CBODD(KMC,IMC) 
      DIMENSION OMRM(KMC,IMC),     NH4RM(KMC,IMC),   NO3RM(KMC,IMC)
      DIMENSION ARMR(KMC,IMC),     ARMF(KMC,IMC)
      DIMENSION SETIN(KMC,IMC),    SETOUT(KMC,IMC)
      DIMENSION LFAC(KMC,IMC),     LFPR(KMC,IMC)
      DIMENSION QUH1(KMP,NBP),     QUH2(KMP,NBP),    QDH1(KMP,NBP),      &
                QDH2(KMP,NBP),     TUH(KMP,NBP),     TDH(KMP,NBP)
      DIMENSION TSSUH1(KMP,NBP),   TSSUH2(KMP,NBP),  TSSDH1(KMP,NBP),    &
                TSSDH2(KMP,NBP)
      DIMENSION QINF(KMP,NBP),     QTRF(KMP,NTP)
      DIMENSION AKBR(KMP,NBP)
      DIMENSION QOUT(KMP,NBP),     KOUT(KMP,NBP)
      DIMENSION FETCHU(IMP,NBP),   FETCHD(IMP,NBP)
      DIMENSION ISO_CONC(NCP),     VERT_CONC(NCP),   LONG_CONC(NCP)
      DIMENSION C1S(KMC,IMC,NCP),  CSSB(KMC,IMC,NCP)
      DIMENSION CVP(KMC,NCP)
      DIMENSION SF1L(KMP,IMP),     SF2L(KMP,IMP,2),  SF3L(KMP,IMP,2),    &
                SF4L(KMP,IMP,2),   SF5L(KMP,IMP,2),  SF6L(KMP,IMP,2),    &
                SF7L(KMP,IMP,2),   SF8L(KMP,IMP,2),  SF9L(KMP,IMP,2),    &
                SF10L(KMP,IMP,2),  SF11L(KMP,IMP,2), SF12L(KMP,IMP,2),   &
                SF13L(KMP,IMP,2)
      DIMENSION SF1V(KMP),         SF2V(KMP,2),      SF3V(KMP,2),        &
                SF4V(KMP,2),       SF5V(KMP,2),      SF6V(KMP,2),        &
                SF7V(KMP,2),       SF8V(KMP,2),      SF9V(KMP,2),        &
                SF10V(KMP,2)
      DIMENSION DX1(KMP,IMP),      DX2(KMP,IMP),     DX3(KMP,IMP),       &
                AD1L(KMP,IMP),     AD2L(KMP,IMP),    AD3L(KMP,IMP),      &
                AD1V(KMP,IMP),     AD2V(KMP,IMP),    AD3V(KMP,IMP),      &
                TADL(KMP,IMP),     TADV(KMP,IMP)
      DIMENSION CT(KMP,IMP),       AT(KMP,IMP),      BTAT(KMP,IMP)
      DIMENSION VT(KMP),           DT(KMP),          GMAT(KMP)
      DIMENSION CADL(KMC,IMC,NCP),    CADV(KMC,IMC,NCP)
      DIMENSION CSSUH1(KMP,NCP,NBP),  CSSUH2(KMP,NCP,NBP),               &
                CSSDH1(KMP,NCP,NBP),  CSSDH2(KMP,NCP,NBP)
      DIMENSION SNPD(NDP), VPLD(NDP), PRFD(NDP), CPLD(NDP), RSOD(NDP),   &
                TSRD(NDP), DLTD(NDP), WSCD(NDP), SPRD(NDP), SCRD(NDP)
      DIMENSION SNPF(NDP), VPLF(NDP), PRFF(NDP), CPLF(NDP), RSOF(NDP),   &
                TSRF(NDP), DLTF(NDP), SCRF(NDP), SPRF(NDP)
      DIMENSION WSC(NDP),  DLTMAX(NDP)
      DIMENSION IPR(IMP),  IPRI(IMP), TITLE(7),  HNAME(4),  HPRC(4)
      DIMENSION SEL_WITHDRAWAL(NBP),  POINT_SINK(NSP,NBP),               &
                ICE_IN(NBP),          ALLOW_ICE(IMP),                    &
                TRANSPORT(NCP),       ONE_LAYER(IMP)
      DIMENSION PLACE_QTR(NTP),       SPECIFY_QTR(NTP)

!**** Common declaration

      COMMON /GLOBLC/ JB,    JC,      IU,     ID,     KT,    ELKT,       &
                      DLT,   KB,      KTI
      COMMON /GEOMHC/ EL,    ELWS,    ELWS2,  H,      HKT1,  HKT2,       &
                      DLX
      COMMON /GEOMBC/ B,     BKT,     BH,     BHKT1,  BHKT2, BHRKT1
      COMMON /DEPTHC/ DEPTHM
      COMMON /NAMESC/ HNAME, CNAME2,  CUNIT1
      COMMON /SNPUC/  SNP
      COMMON /TEMPC/  T1,    T2
      COMMON /ICEC/   ICE,   ICETH,   ICE_CALC
      COMMON /HYDRC1/ U,     W,       AZ,     RHO,    NDLT
      COMMON /HYDRC2/ Z
      COMMON /PRNTC1/ IPR,   IEPR,    KBPR
      COMMON /PRNTC2/ TITLE, CPRC,    HPRC,   CONV
      COMMON /PRNTC3/ CUS,   DS
      COMMON /GRDLGC/ LIMITING_FACTOR,OXYGEN_DEMAND
      COMMON /DNSPHC/ FRESH_WATER,    SALT_WATER,     SUSP_SOLIDS
      COMMON /GRTVDC/ CONSTITUENTS,   CN,             NAC
      COMMON /INTERC/ INTERP_INFLOW,  INTERP_OUTFLOW, INTERP_MET,        &
                      INTERP_DTRIBS,  INTERP_TRIBS,   INTERP_HEAD,       &
                      INTERP_WITHDRWL
      COMMON /TVDDSC/ NO_WIND,        NO_INFLOW,      NO_OUTFLOW,        &
                      NO_HEAT
      COMMON /TVDLC1/ PRECIPITATION,  WITHDRAWALS,    TRIBUTARIES,       &
                      DIST_TRIBS
      COMMON /TVDLC2/ UP_FLOW,        DN_FLOW,        UH_INTERNAL,       &
                      UH_EXTERNAL,    DH_INTERNAL,    DH_EXTERNAL
      COMMON /TVDLC3/ OPEN_FILES,     TERM_BY_TERM,   WATER_BALANCE
      COMMON /TVDMTC/ TAIR,   TDEW,   CLOUD,  PHI,    ET,     CSHE,      &
                      SRO,    SRON,   LAT,    LONG,   WINDH,  RHOWCP
      COMMON /TVDFNC/ METFN,  QWDFN,  QOTFN,  QINFN,  TINFN,  CINFN,     &
                      QTRFN,  TTRFN,  CTRFN,  QDTFN,  TDTFN,  CDTFN,     &
                      PREFN,  TPRFN,  CPRFN,  EUHFN,  TUHFN,  CUHFN,     &
                      EDHFN,  TDHFN,  CDHFN,  ELOFN
      COMMON /TVDQC/  QIN,    QTR,    QDTR,   PR,     ELUH,   ELDH,      &
                      QOUT,   QWD,    QSUM
      COMMON /TVDTC/  TIN,    TTR,    TDTR,   TPR,    TUH,    TDH
      COMMON /TVDCC1/ CIN,    CTR,    CDTR,   CPR,    CUH,    CDH
      COMMON /TVDCC2/ INCN,   TRCN,   DTCN,   PRCN,   UHCN,   DHCN
      COMMON /TVDCC3/ NACIN,  NACTR,  NACPR,  NACDT
      COMMON /RTMLTC/ OMT1,   OMT2,   NH4T1,  NH4T2,  NO3T1,  NO3T2,     &
                      AT1,    AT2,    AT3,    AT4,    OMK1,   OMK2,      &
                      NH4K1,  NH4K2,  NO3K1,  NO3K2,  AK1,    AK2,       &
                      AK3,    AK4
      COMMON /DKORGC/ LPOMD,  DOMD
      COMMON /DKSEDC/ SEDD,   SODD,   SOD
      COMMON /DKNITC/ NH4D,   NO3D
      COMMON /DKBODC/ CBODD
      COMMON /GBLRTC/ OMRM,   NH4RM,  NO3RM,  ARMR,   ARMF
      COMMON /GLBLCC/ PALT,   APOM,   O2LIM,  WIND,   WSCDP,  WSC
      COMMON /DKMLTC/ A1,     A2,     A3
      COMMON /CLFRMC/ COLQ10, COLDK
      COMMON /ORGDKC/ SDK,    LPOMDK, LDOMDK, RDOMDK, LRDDK
      COMMON /SETLC1/ SETIN,  SETOUT
      COMMON /SETLC2/ SSS,    POMS,   AS,     FES
      COMMON /PHYTGC/ AGR,    ARR,    AMR,    AER
      COMMON /PHYTC1/ AE,     AM,     AG,     AR,     ASAT,   AHSN,      &
                      AHSP
      COMMON /PHYTC2/ BETA,   EXH2O,  EXSS,   EXOM
      COMMON /PHOSPC/ PO4R,   BIOP,   PARTP
      COMMON /NITROC/ BION,   NH4DK,  NH4R,   NO3DK
      COMMON /PO4C1/  FPSS,   FPFE
      COMMON /OXYGNC/ O2OM,   O2AG,   O2NH4,  O2AR
      COMMON /CARBNC/ CO2R,   BIOC
      COMMON /IRONC/  FER
      COMMON /CBODC/  KBOD,   TBOD,   RBOD
      COMMON /LFACC/  LFAC,   LFPR
      COMMON /SELWC/  NSTR,   QSTR,   ESTR,   WSTR,   KBSW,   NOUT,      &
                      KOUT
      COMMON /GDAYC1/ GDAY,   DAYM,   JDAYG,  LEAP_YEAR, M
      COMMON /GDAYC2/ MONTH
      COMMON /SCRNC1/ JDAY,   DLTS,   ILOC,   KLOC,   MINDLT, JDMIN,     &
                      IMIN,   KMIN,   DLTAV,  NIT,    NV,     YEAR,      &
                      ELTMJD
      COMMON /SCRNC2/ NWD,    NTR
      COMMON /SCRNC3/ ZMIN
      COMMON /TVDSWC/ SEL_WITHDRAWAL, POINT_SINK
      COMMON /NONZC/  NONZERO
      COMMON /ELBALC/ NXELV1, NXELV2, ELWSO1, ELWSO2

!**** Data statements

      DATA CNAME1 /'Tracer',          'Suspended_solids',                &
                   'Coliform',        'Dissolved_solids',                &
                   'Labile_DOM',      'Refractory_DOM',                  &
                   'Algae',           'Labile POM',                      &
                   'Phosphorus',      'Ammonium',                        &
                   'Nitrate-nitrite', 'Dissolved_oxygen',                &
                   'Sediment',        'Inorganic_carbon',                &
                   'Alkalinity',      'pH',                              &
                   'Carbon_dioxide',  'Bicarbonate',                     &
                   'Carbonate',       'Iron',                            &
                   'CBOD'/
      DATA CUNIT2 /21*'g/m 3'/
      DATA CUNIT3 /21*'g '/
      DATA RHOA   /1.25/,   RHOW  /1000.0/,   RHOI    /916.0/,           &
           RK1    /2.12/,   RL1   /333507.0/, RIMT    /0.0/
      DATA G      /9.8/,    VTOL  /1.0E3/,    CP      /4186.0/,          &
           FRAZDZ /0.14/,   AZMIN /1.4E-5/,   DZMIN   /1.4E-7/,          &
           AZMAX  /1.0E-4/, DZMAX /1.0/,      ICE_TOL /0.005/
!*kevin's mod
      DATA BLANK  /'       '/
      DATA BLANK1 /'         '/
      DATA CON /10/, BTH /11/, RSI /12/, VPR /13/, LPR /14/
      DATA RSO /21/, TSR /22/, PRF /23/, VPL /24/, CPL /25/, WRN /26/,   &
           ERR /27/, SPR /28/, QWS /29/
      NB = NBP

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                               Task 1.3: Inputs                           **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

!**** Open control file

      OPEN (CON,FILE=CONFN,STATUS='OLD')

!**** Title cards

      READ (CON,*)
      READ (CON,1000) (TITLE(J),J=1,6)

!**** Time control cards

      READ (CON,1010)  TMSTRT, TMEND, YEAR
      READ (CON,1020)  NDT,    DLTMIN
      READ (CON,1030) (DLTD(J),  J=1,NDT)
      READ (CON,1030) (DLTMAX(J),J=1,NDT)
      READ (CON,1030) (DLTF(J),  J=1,NDT)

!**** Grid definition cards

      READ (CON,1040) (US(JB), DS(JB), UHS(JB), DHS(JB), NL(JB),  &
                       JB=1,NBP)
      READ (CON,1030)  LAT,    LONG,   ELBOT
      
!**** Initial condition cards

      READ (CON,1050)  T2I,    ICETHI, WTYPE
      READ (CON,1060)  VBC,    EBC,    MBC,    WBC,    PQC,   EVC,  &
                       PRC
      READ (CON,1060)  INFIC,  TRIC,   DTRIC,  HDIC,   OUTIC, WDIC,  &
                       METIC
      READ (CON,1060)  WINDC,  QINC,   QOUTC,  HEATC
      READ (CON,1070)  ICEC,   SLICE,  SLHEAT, ALBEDO, HWI,   BETAI,  &
                       GAMMAI, ICEMIN, ICET2
      READ (CON,1080)  SLTR,   THETAI
      READ (CON,1020)  NWSC,   WINDH
      READ (CON,1030) (WSCD(J),J=1,NWSC)
      READ (CON,1030) (WSC(J), J=1,NWSC)
      READ (CON,1030)  AX, DXI, CHEZY, CBHE, TSED

!**** Inflow-outflow cards

      READ (CON,1060) (SWC(JB), JB=1,NBP)
      READ (CON,1090) (NSTR(JB),JB=1,NBP)
      READ (CON,1090)   (KBSW(JS,1), JS=1,NSTR(1))
      DO JB=2,NB
        READ (CON,1130) (KBSW(JS,JB),JS=1,NSTR(JB))
      END DO
      READ (CON,1100)   (SINK(JS,1), JS=1,NSTR(1))
      DO JB=2,NB
        READ (CON,1110) (SINK(JS,JB),JS=1,NSTR(JB))
      END DO
      READ (CON,1030)   (ESTR(JS,1), JS=1,NSTR(1))
      DO JB=2,NB
        READ (CON,1120) (ESTR(JS,JB),JS=1,NSTR(JB))
      END DO
      READ (CON,1030)   (WSTR(JS,1), JS=1,NSTR(1))
      DO JB=2,NB
        READ (CON,1120) (WSTR(JS,JB),JS=1,NSTR(JB))
      END DO
      READ (CON,1090)   (NOUT(JB),  JB=1,NBP)
      READ (CON,1090)   (KOUT(JO,1),JO=1,NOUT(1))
      DO JB=2,NB
        READ (CON,1130) (KOUT(JO,JB),JO=1,NOUT(JB))
      END DO
      READ (CON,1090)  NWD
      READ (CON,1090) (IWD(JW),  JW=1,NWD)
      READ (CON,1090) (KWD(JW),  JW=1,NWD)
      READ (CON,1090)  NTR
      READ (CON,1100) (TRC(JT),  JT=1,NTR)
      READ (CON,1090) (ITR(JT),  JT=1,NTR)
      READ (CON,1030) (ELTRT(JT),JT=1,NTR)
      READ (CON,1030) (ELTRB(JT),JT=1,NTR)
      READ (CON,1060) (DTRC(JB), JB=1,NBP)

!**** Output control cards (excluding constituents)

      READ (CON,1150)  SCRC, NSCR
      READ (CON,1030) (SCRD(J),J=1,NSCR)
      READ (CON,1030) (SCRF(J),J=1,NSCR)
      READ (CON,1140)  LJPC, (HPRC(J),J=1,4)
      READ (CON,1160)  SNPC, NSNP, NISNP
      READ (CON,1030) (SNPD(J),J=1,NSNP)
      READ (CON,1030) (SNPF(J),J=1,NSNP)
      READ (CON,1090) (IPRI(I),I=1,NISNP)
      READ (CON,1160)  PRFC, NPRF, NIPRF
      READ (CON,1030) (PRFD(J),J=1,NPRF)
      READ (CON,1030) (PRFF(J),J=1,NPRF)
      READ (CON,1090) (IPRF(J),J=1,NIPRF)
      READ (CON,1160)  SPRC, NSPR, NISPR
      READ (CON,1030) (SPRD(J),J=1,NSPR)
      READ (CON,1030) (SPRF(J),J=1,NSPR)
      READ (CON,1090) (ISPR(J),J=1,NISPR)
      READ (CON,1150)  TSRC, NTSR
      READ (CON,1030) (TSRD(J),J=1,NTSR)
      READ (CON,1030) (TSRF(J),J=1,NTSR)
      READ (CON,1150)  VPLC, NVPL
      READ (CON,1030) (VPLD(J),J=1,NVPL)
      READ (CON,1030) (VPLF(J),J=1,NVPL)
      READ (CON,1150)  CPLC, NCPL
      READ (CON,1030) (CPLD(J),J=1,NCPL)
      READ (CON,1030) (CPLF(J),J=1,NCPL)
      READ (CON,1150)  RSOC, NRSO, RSIC
      READ (CON,1030) (RSOD(J),J=1,NRSO)
      READ (CON,1030) (RSOF(J),J=1,NRSO)

!**** Constituent control cards

      READ (CON,1170)  CCC, LIMC, SDC, FREQUK
      READ (CON,1060) (ACC(JC),  JC=1,NCP)
      READ (CON,1030) (C2I(JC),  JC=1,NCP)
      READ (CON,1060) (CPRC(JC), JC=1,NCP)
      READ (CON,1060) (INACC(JC),JC=1,NCP)
      READ (CON,1060) (TRACC(JC),JC=1,NCP)
      READ (CON,1060) (DTACC(JC),JC=1,NCP)
      READ (CON,1060) (PRACC(JC),JC=1,NCP)

!**** Kinetics coefficients

      READ (CON,1030)  EXH2O,  EXSS,   EXOM,   BETA
      READ (CON,1030)  COLQ10, COLDK
      READ (CON,1030)  SSS
      READ (CON,1030)  AG,     AM,     AE,     AR,     AS,    ASAT,  &
                       APOM
      READ (CON,1030)  AT1,    AT2,    AT3,    AT4,    AK1,   AK2,  &
                       AK3,    AK4
      READ (CON,1030)  LDOMDK, LRDDK,  RDOMDK
      READ (CON,1030)  LPOMDK, POMS
      READ (CON,1030)  OMT1,   OMT2,   OMK1,   OMK2
      READ (CON,1030)  SDK,    FSOD
      READ (CON,1030) (SOD(I),I=1,IMP)
      READ (CON,1030)  KBOD,   TBOD,   RBOD
      READ (CON,1030)  PO4R,   PARTP,  AHSP
      READ (CON,1030)  NH4R,   NH4DK,  AHSN
      READ (CON,1030)  NH4T1,  NH4T2,  NH4K1,  NH4K2
      READ (CON,1030)  NO3DK
      READ (CON,1030)  NO3T1,  NO3T2,  NO3K1,  NO3K2
      READ (CON,1030)  CO2R
      READ (CON,1030)  FER,    FES
      READ (CON,1030)  O2NH4,  O2OM,   O2AR,   O2AG,  BIOP,  BION,  &
                       BIOC
      READ (CON,1030)  O2LIM

!**** Input filenames

      READ (CON,1000)  BTHFN
write(*,*) 'BTHFN',BTHFN
      READ (CON,1000)  VPRFN
write(*,*) 'VPRFN ',VPRFN
      READ (CON,1000)  LPRFN
write(*,*) 'QINFN ',NBP,QINFN(1)
      READ (CON,1000)  RSIFN
write(*,*) 'RSIFN ',NBP,RSIFN
      READ (CON,1000)  METFN
write(*,*) 'METFN ',NBP,METFN
      READ (CON,1000)  QWDFN
!write(*,*) 'QWDFN ',NBP,QWDFN
!      READ (CON,1000)  ELOFN
write(*,*) 'ELOFN ',NBP,ELOFN
      READ (CON,1000) (QINFN(JB),JB=1,NBP)
write(*,*) 'QINFN ',NBP,QINFN(1)
      READ (CON,1000) (TINFN(JB),JB=1,NBP)
write(*,*) 'TINFN ',NBP,TINFN(1)
      READ (CON,1000) (CINFN(JB),JB=1,NBP)
      READ (CON,1000) (QOTFN(JB),JB=1,NBP)
write(*,*) 'QINFN ',NBP,QOTFN(1)
      READ (CON,1000) (QTRFN(JT),JT=1,NTP)
      READ (CON,1000) (TTRFN(JT),JT=1,NTP)
      READ (CON,1000) (CTRFN(JT),JT=1,NTP)
      READ (CON,1000) (QDTFN(JB),JB=1,NBP)
      READ (CON,1000) (TDTFN(JB),JB=1,NBP)
      READ (CON,1000) (CDTFN(JB),JB=1,NBP)
      READ (CON,1000) (PREFN(JB),JB=1,NBP)
      READ (CON,1000) (TPRFN(JB),JB=1,NBP)
      READ (CON,1000) (CPRFN(JB),JB=1,NBP)
      READ (CON,1000) (EUHFN(JB),JB=1,NBP)
      READ (CON,1000) (TUHFN(JB),JB=1,NBP)
      READ (CON,1000) (CUHFN(JB),JB=1,NBP)
      READ (CON,1000) (EDHFN(JB),JB=1,NBP)
      READ (CON,1000) (TDHFN(JB),JB=1,NBP)
      READ (CON,1000) (CDHFN(JB),JB=1,NBP)
 
!**** Output filenames

      READ (CON,1000)  SNPFN
write(*,*) 'SNPFN ',SNPFN
      READ (CON,1000)  TSRFN
      READ (CON,1000)  PRFFN
      READ (CON,1000)  VPLFN
      READ (CON,1000)  CPLFN
      READ (CON,1000)  SPRFN

!**** Bathymetry file

      OPEN (BTH,FILE=BTHFN,STATUS='OLD')
      READ (BTH,*)
      READ (BTH,1180) (DLX(I), I=1,IMP)
      READ (BTH,1180) (ELWS(I),I=1,IMP)
      READ (BTH,1180) (PHI0(I),I=1,IMP)
      READ (BTH,1180) (H(K),   K=1,KMP)
      DO I=1,IMP
        READ (BTH,1180) (B(K,I),K=1,KMP)
        ELWS2(I) = ELWS(I)
      END DO

!**** Initialize logical control variables

      CONSTITUENTS = CCC.EQ.' ON'
      RESTART_IN   = RSIC.EQ.' ON'
      ISO_TEMP     = INT(T2I).GE.0
      VERT_TEMP    = INT(T2I).EQ.-1
      LONG_TEMP    = INT(T2I).LT.-1
      OPEN_VPR     = VERT_TEMP.AND..NOT.RESTART_IN
      OPEN_LPR     = LONG_TEMP.AND..NOT.RESTART_IN
      DO JC=1,NCP
        ISO_CONC(JC)  = INT(C2I(JC)).GE.0
        VERT_CONC(JC) = INT(C2I(JC)).EQ.-1
        LONG_CONC(JC) = INT(C2I(JC)).LT.-1
        IF (VERT_CONC(JC).AND..NOT.RESTART_IN) OPEN_VPR = .TRUE.
        IF (LONG_CONC(JC).AND..NOT.RESTART_IN) OPEN_LPR = .TRUE.
      END DO

!**** Close files

      CLOSE (CON)
      CLOSE (BTH)
!      IF (RESTART_IN) CLOSE (RSI)

!**** Input FORMATs

 1000 FORMAT(//(8X,A72))
 1010 FORMAT(//(8X,2F8.0,I8))
 1020 FORMAT(//8X,I8,8F8.0)
 1030 FORMAT(//(8X,9F8.0))
 1040 FORMAT(//(8X,5I8))
 1050 FORMAT(//8X,2F8.0,3X,A5,6(5X,A3))
 1060 FORMAT(//(8X,9(5X,A3)))
 1070 FORMAT(//13X,A3,2A8,6F8.0)
 1080 FORMAT(//8X,A8,F8.0,2A8)
 1090 FORMAT(//(8X,9I8))
 1100 FORMAT(//(8X,9A8))
 1110 FORMAT(:8X,9A8)
 1120 FORMAT(:8X,9F8.0)
 1130 FORMAT(:8X,9I8)
 1140 FORMAT(//11X,A5,8(5X,A3))
 1150 FORMAT(//13X,A3,I8,5X,A3)
 1160 FORMAT(//13X,A3,8I8)
 1170 FORMAT(//8X,3(5X,A3),I8)
 1180 FORMAT(//(10F8.0))
 1190 FORMAT(//(8X,9F8.0))

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                               Task 1.4: Variables                        **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                            Task 1.4.1: Zero Variables                    **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      DO I=1,IMP
        DO K=1,KMP
          CONV(K,I)  = BLANK
          CONV1(K,I) = BLANK1
        END DO
      END DO
      DO I=1,IMC
        DO K=1,KMC
          LFPR(K,I) = BLANK
        END DO
      END DO
      DO I=1,IMP
        ICESW(I) = 1.0
        DO JB=1,NBP
          FETCHU(I,JB) = 0.0
          FETCHD(I,JB) = 0.0
        END DO
!        IF (.NOT.RESTART_IN) THEN
          DO K=1,KMP
            U(K,I) = 0.0
            W(K,I) = 0.0
          END DO
!        END IF
        DO K=1,KMP
          TSS(K,I) = 0.0
          QSS(K,I) = 0.0
          IF (CONSTITUENTS) THEN
            DO JC=1,NCP
!             IF (.NOT.RESTART_IN) THEN
                CSSB(K,I,JC) = 0.0
                CSSK(K,I,JC) = 0.0
!             END IF
            END DO
          END IF
        END DO
      END DO
      DO JB=1,NBP
          VOLEV(JB)  = 0.0
          VOLPR(JB)  = 0.0
          VOLTR(JB)  = 0.0
          VOLDT(JB)  = 0.0
          VOLWD(JB)  = 0.0
          VOLUH(JB)  = 0.0
          VOLDH(JB)  = 0.0
          VOLIN(JB)  = 0.0
          VOLOUT(JB) = 0.0
          VOLSBR(JB) = 0.0
          TSSS(JB)   = 0.0
          TSSB(JB)   = 0.0
          TSSEV(JB)  = 0.0
          TSSPR(JB)  = 0.0
          TSSTR(JB)  = 0.0
          TSSDT(JB)  = 0.0
          TSSWD(JB)  = 0.0
          TSSUH(JB)  = 0.0
          TSSDH(JB)  = 0.0
          TSSIN(JB)  = 0.0
          TSSOUT(JB) = 0.0
          TSSICE(JB) = 0.0
          ESBR(JB)   = 0.0
          ETBR(JB)   = 0.0
          EIBR(JB)   = 0.0
        DO K=1,KMP
          AKBR(K,JB) = 0.0
        END DO
        DO JC=1,NCP
          CMBRT(JC,JB) = 0.0
        END DO
      END DO
      DLMR     = 0.0
      DLVR     = 0.0
      KBPR     = 0
      NAC      = 0
      NACIN    = 0
      NACTR    = 0
      NACDT    = 0
      NACPR    = 0
      NSPRF    = 0
      NTAC     = 0
      KBMAX    = 0
      HMAX     = 0
      DLXMAX   = 0
      HMIN     = 1.0E10
      DLXMIN   = 1.0E10
      TITLE(7) = ' '
        KLOC = 1
        ILOC = 1

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                        Task 1.4.2: Miscellaneous Variables               **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

!**** Logical controls

      OPEN_FILES      = .TRUE.
      VOLUME_WARNING  = .TRUE.
      UPDATE_KINETICS = .TRUE.
      END_RUN         = .FALSE.
      WARNING_OPEN    = .FALSE.
      TRIBUTARIES     = NTR.GT.0
      WITHDRAWALS     = NWD.GT.0
      PLACE_QIN       = PQC.EQ.' ON'
      EVAPORATION     = EVC.EQ.' ON'
      VOLUME_BALANCE  = VBC.EQ.' ON'
      ENERGY_BALANCE  = EBC.EQ.' ON'
      WATER_BALANCE   = WBC.EQ.' ON'
      PRECIPITATION   = PRC.EQ.' ON'
      SCREEN_OUTPUT   = SCRC.EQ.' ON'
      SNAPSHOT        = SNPC.EQ.' ON'
      CONTOUR         = CPLC.EQ.' ON'
      VECTOR          = VPLC.EQ.' ON'
      PROFILE         = PRFC.EQ.' ON'
      SPREADSHEET     = SPRC.EQ.' ON'
      ICE_CALC        = ICEC.EQ.' ON'
      TIME_SERIES     = TSRC.EQ.' ON'
      RESTART_OUT     = RSOC.EQ.' ON'
      INTERP_TRIBS    = TRIC.EQ.' ON'
      INTERP_HEAD     = HDIC.EQ.' ON'
      INTERP_WITHDRWL = WDIC.EQ.' ON'
      INTERP_MET      = METIC.EQ.' ON'
      INTERP_INFLOW   = INFIC.EQ.' ON'
      INTERP_DTRIBS   = DTRIC.EQ.' ON'
      INTERP_OUTFLOW  = OUTIC.EQ.' ON'
      LIMITING_DLT    = HPRC(4).EQ.' ON'
      NO_INFLOW       = QINC.EQ.'OFF'
      NO_HEAT         = HEATC.EQ.'OFF'
      NO_WIND         = WINDC.EQ.'OFF'
      NO_OUTFLOW      = QOUTC.EQ.'OFF'
      LASERJET_II     = LJPC.EQ.'   II'
      LASERJET_III    = LJPC.EQ.'  III'
      LASERJET_IV     = LJPC.EQ.'   IV'
      UPWIND          = SLTR.EQ.'  UPWIND'
      TERM_BY_TERM    = SLHEAT.NE.'      ET'
      DETAILED_ICE    = SLICE.EQ.'  DETAIL'
      DETAILED_ICE    = ICE_CALC.AND.DETAILED_ICE
      MASS_BALANCE    = CONSTITUENTS.AND.MBC.EQ.' ON'
      SUSP_SOLIDS     = CONSTITUENTS.AND.ACC(2).EQ.' ON'
      OXYGEN_DEMAND   = CONSTITUENTS.AND.ACC(12).EQ.' ON'
      LIMITING_FACTOR = CONSTITUENTS.AND.ACC(7).EQ.' ON'  &
                        .AND.LIMC.EQ.' ON'
      FRESH_WATER     = CONSTITUENTS.AND.ACC(4).EQ.' ON'  &
                        .AND.WTYPE.EQ.'FRESH'
      SALT_WATER      = CONSTITUENTS.AND.ACC(4).EQ.' ON'  &
                        .AND.WTYPE.EQ.' SALT'
      SHIFT_DEMAND    = CONSTITUENTS.AND.ACC(12).EQ.' ON'  &
                        .AND.SDC.EQ.' ON'
      INTERPOLATE     = INTERP_INFLOW.OR.INTERP_TRIBS.OR.INTERP_DTRIBS  &
                        .OR.INTERP_OUTFLOW.OR.INTERP_WITHDRWL  &
                        .OR.INTERP_HEAD.OR.INTERP_MET
      LEAP_YEAR       = MOD(YEAR,4).EQ.0
        WINTER = .FALSE.
        IF (TMSTRT.GT.300.0.OR.TMSTRT.LT.40.0) WINTER = .TRUE.
      DO JT=1,NTP
        PLACE_QTR(JT)   = TRC(JT).EQ.' DENSITY'
        SPECIFY_QTR(JT) = TRC(JT).EQ.' SPECIFY'
      END DO
      DO JC=1,NCP
        TRANSPORT(JC) = .TRUE.
        IF (JC.EQ.13.OR.(JC.GT.15.AND.JC.LT.20)) TRANSPORT(JC) = .FALSE.
      END DO
      DO JC=5,13
        IF (ACC(JC).EQ.' ON') UPDATE_RATES = .TRUE.
      END DO
      DO JB=1,NBP
        UP_FLOW(JB)        = UHS(JB).EQ.0
        DN_FLOW(JB)        = DHS(JB).EQ.0
        UP_HEAD(JB)        = .FALSE.
        DN_HEAD(JB)        = .FALSE.
        UH_INTERNAL(JB)    = .FALSE.
        DH_INTERNAL(JB)    = .FALSE.
        UH_EXTERNAL(JB)    = .FALSE.
        DH_EXTERNAL(JB)    = .FALSE.
        DIST_TRIBS(JB)     = .FALSE.
        SEL_WITHDRAWAL(JB) = SWC(JB).EQ.' ON'
        DO JS=1,NSTR(JB)
          POINT_SINK(JS,JB) = SINK(JS,JB).EQ.'   POINT'
        END DO
      END DO

!**** Convert rates from per-day to per-second

      AE     = AE/86400.0
      AM     = AM/86400.0
      AR     = AR/86400.0
      AG     = AG/86400.0
      AS     = AS/86400.0
      FES    = FES/86400.0
      SSS    = SSS/86400.0
      POMS   = POMS/86400.0
      SDK    = SDK/86400.0
      COLDK  = COLDK/86400.0
      LRDDK  = LRDDK/86400.0
      NH4DK  = NH4DK/86400.0
      NO3DK  = NO3DK/86400.0
      LDOMDK = LDOMDK/86400.0
      RDOMDK = RDOMDK/86400.0
      LPOMDK = LPOMDK/86400.0
      KBOD   = KBOD/86400.0
      IF (CONSTITUENTS) THEN
        DO I=1,IMP
          SOD(I)  = SOD(I)/86400.0*FSOD
          SODS(I) = SOD(I)
        END DO
      END IF


!**** Time and printout control variables

        JDAY   = TMSTRT
        JDAYG  = JDAY
        ELTM   = TMSTRT*86400.0
        DLT    = DLTMAX(1)
        DLTS   = DLT
        MINDLT = DLT
        NIT    = 0
        NV     = 0
        DLTDP  = 1
        WSCDP  = 1
        SNPDP  = 1
        TSRDP  = 1
        VPLDP  = 1
        PRFDP  = 1
        SPRDP  = 1
        CPLDP  = 1
        RSODP  = 1
        SCRDP  = 1
        NXTMSN = SNPD(1)
        NXTMTS = TSRD(1)
        NXTMPR = PRFD(1)
        NXTMSP = SPRD(1)
        NXTMCP = CPLD(1)
        NXTMVP = VPLD(1)
        NXTMRS = RSOD(1)
        NXTMSC = SCRD(1)
      DO J=NWSC+1,NDP
        WSCD(J) = TMEND+1.0
      END DO
      DO J=NDT+1,NDP
        DLTD(J) = TMEND+1.0
      END DO
      DO J=NSNP+1,NDP
        SNPD(J) = TMEND+1.0
      END DO
      DO J=NTSR+1,NDP
        TSRD(J) = TMEND+1.0
      END DO
      DO J=NPRF+1,NDP
        PRFD(J) = TMEND+1.0
      END DO
      DO J=NSPR+1,NDP
        SPRD(J) = TMEND+1.0
      END DO
      DO J=NVPL+1,NDP
        VPLD(J) = TMEND+1.0
      END DO
      DO J=NCPL+1,NDP
        CPLD(J) = TMEND+1.0
      END DO
      DO J=NRSO+1,NDP
        RSOD(J) = TMEND+1.0
      END DO
      DO J=NSCR+1,NDP
        SCRD(J) = TMEND+1.0
      END DO
      JDAYNX = JDAYG+1
      NXTVD  = JDAY
      NXELV2 = TMSTRT
      CURMAX = DLTMAX(DLTDP)/DLTF(DLTDP)

!**** Active constituents

      IF (CONSTITUENTS) THEN
        DO JC=1,NCP
          IF (ACC(JC).EQ.' ON') THEN
            NAC     = NAC+1
            CN(NAC) = JC
          END IF
          IF (INACC(JC).EQ.' ON') THEN
            NACIN       = NACIN+1
            INCN(NACIN) = JC
          END IF
          IF (TRACC(JC).EQ.' ON') THEN
            NACTR       = NACTR+1
            TRCN(NACTR) = JC
          END IF
          IF (DTACC(JC).EQ.' ON') THEN
            NACDT       = NACDT+1
            DTCN(NACDT) = JC
          END IF
          IF (PRACC(JC).EQ.' ON') THEN
            NACPR       = NACPR+1
            PRCN(NACPR) = JC
          END IF
        END DO
      END IF
      IF (CONSTITUENTS) THEN
        DO JC=1,NCP
          IF (ACC(JC).EQ.' ON') THEN
            NAC0      = NAC0+1
            CN0(NAC0) = JC
          END IF
        END DO
      END IF
      DEG = CHAR(248)//'C'
      ESC = CHAR(027)
      CALL DATE_TIME (CDATE,CTIME)
      TITLE(7) = 'Model run at '//CTIME//' on '//CDATE                  !FORTRAN
!
!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                              Task 1.4.3: Geometry                        **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

!**** Layer elevations and mimimum/maximum layer heights

      EL(KMP) = ELBOT
      DO K=KMP-1,1,-1
        EL(K) = EL(K+1)+H(K)
        HMIN  = MIN(H(K),HMIN)
        HMAX  = MAX(H(K),HMAX)
      END DO

!**** Water surface and bottom layers

      DO JB=1,NBP
        DO I=US(JB)-1,DS(JB)+1
            KT     = 2
            KTI(I) = 2
            DO WHILE (EL(KTI(I)).GT.ELWS(I))
              KTI(I) = KTI(I)+1
            END DO
            Z(I)    = EL(KTI(I))-ELWS(I)
            KTMAX   = MAX(2,KTI(I))
            KT      = MAX(KTMAX,KT)
            KTI(I)  = MAX(2,KTI(I)-1)
            SKTI(I) = KTI(I)
            SZ(I)   = Z(I)
          K = 2
          DO WHILE (B(K,I).GT.0.0)
            KB(I) = K
            K     = K+1
          END DO
          KBMAX = MAX(KBMAX,KB(I))
        END DO
        KB(US(JB)-1) = KB(US(JB))
        KB(DS(JB)+1) = KB(DS(JB))
      END DO
      DO JB=1,NBP
        IU = US(JB)
        ID = DS(JB)

!****** Upstream active segment and single layer

        IUT = IU
        DO I=IU,ID
          IF (KB(I)-KT.LT.NL(JB)-1) IUT = I+1
          ONE_LAYER(I) = KT.EQ.KB(I)
        END DO
        IF (IUT.GT.DS(JB)-1) THEN
          OPEN  (ERR,FILE='err.opt',STATUS='UNKNOWN')
          WRITE (ERR,7000) JB,JDAY,KT
          STOP
        END IF
        CUS(JB) = IUT

!****** Boundary bottom layers

        IF (UH_INTERNAL(JB)) KB(IUT-1) = MIN(KB(UHS(JB)),KB(IUT))
        IF (DH_INTERNAL(JB)) KB(ID+1)  = MIN(KB(DHS(JB)),KB(ID))

!****** Boundary segment lengths

        DLX(IU-1) = DLX(IU)
        DLX(ID+1) = DLX(ID)

!****** Minimum bottom layers and average segment lengths 

        DO I=IU-1,ID
          KBMIN(I) = MIN(KB(I),KB(I+1))
          DLXR(I)  = (DLX(I)+DLX(I+1))*0.5
        END DO
        DLXR(ID+1) = DLX(ID)

!****** Minimum/maximum segment lengths

        DO I=IU,ID
          DLXMIN = MIN(DLXMIN,DLX(I))
          DLXMAX = MAX(DLXMAX,DLX(I))
        END DO

!****** Boundary widths

        DO K=1,KB(IU)
          B(K,IU-1) = B(K,IU)
          IF (UH_INTERNAL(JB)) B(K,IU-1) = B(K,UHS(JB))
        END DO
        DO K=1,KB(ID)
          B(K,ID+1) = B(K,ID)
          IF (DH_INTERNAL(JB)) B(K,ID+1) = B(K,DHS(JB))
        END DO
        DO I=IU-1,ID+1
          B(1,I) = B(2,I)
          DO K=KB(I)+1,KMP
            B(K,I) = B(KB(I),I)
          END DO
        END DO

!****** Areas and bottom widths

        DO I=IU-1,ID+1
          DO K=1,KMP-1
            BH(K,I) = B(K,I)*H(K)
            BB(K,I) = (B(K,I)+B(K+1,I))*0.5
          END DO
          BH(KB(I)+1,I) = BH(KB(I),I)
        END DO

!****** Derived geometry

        DO I=IU-1,ID+1
          HKT2(I)  = H(KT)-Z(I)
          AVHKT(I) = (HKT2(I)+H(KT+1))*0.5
          BHKT2(I) = B(KTI(I),I)*(EL(KT)-EL(KTI(I)+1)-Z(I))
          DO K=KTI(I)+1,KT
            BHKT2(I) = BHKT2(I)+BH(K,I)
          END DO
          BKT(I) = BHKT2(I)/HKT2(I)
        END DO
        IDT = ID+1
        IF (JB.EQ.NBP) IDT = ID
        DO I=IU-1,IDT
          BHRKT2(I) = (BHKT2(I)+BHKT2(I+1))*0.5
          DO K=1,KMP-1
            BR(K,I)  = (B(K,I)+B(K,I+1))*0.5
            BHR(K,I) = (BH(K,I)+BH(K,I+1))*0.5
          END DO
        END DO
        DO K=1,KMP-1
          AVH(K) = (H(K)+H(K+1))*0.5
        END DO

!****** Branch numbers corresponding to tributaries, withdrawals, and head

        IF (TRIBUTARIES) THEN
          DO JT=1,NTR
            IF (ITR(JT).GE.US(JB).AND.ITR(JT).LE.DS(JB)) THEN
              JBTR(JT) = JB
            END IF
          END DO
        END IF
        IF (WITHDRAWALS) THEN
          DO JW=1,NWD
            IF (IWD(JW).GE.US(JB).AND.IWD(JW).LE.DS(JB)) THEN
              JBWD(JW) = JB
            END IF
          END DO
        END IF
        IF (UH_INTERNAL(JB)) THEN
          JBUH(JB)     = 0
          BRANCH_FOUND = .FALSE.
          DO WHILE (.NOT.BRANCH_FOUND)
            JBUH(JB) = JBUH(JB)+1
            DO I=US(JBUH(JB)),DS(JBUH(JB))
              IF (I.EQ.UHS(JB)) BRANCH_FOUND = .TRUE.
            END DO
          END DO
        END IF
        IF (DH_INTERNAL(JB)) THEN
          JBDH(JB)     = 0
          BRANCH_FOUND = .FALSE.
          DO WHILE (.NOT.BRANCH_FOUND)
            JBDH(JB) = JBDH(JB)+1
            DO I=US(JBDH(JB)),DS(JBDH(JB))
              IF (I.EQ.DHS(JB)) BRANCH_FOUND = .TRUE.
            END DO
          END DO
        END IF

!****** Branch layer area

        DO K=KMP-1,2,-1
          DO I=IU,ID
            IF (K.LE.KB(I)) AKBR(K,JB) = AKBR(K,JB)+B(K,I)*DLX(I)
          END DO
        END DO

!****** Layer bottom and middle depths

        DEPTHB(KT)   = HKT2(I)
        DEPTHM(KT)   = HKT2(I)*0.5
        DEPTHB(KT+1) = DEPTHB(KT)+H(KT+1)
        DEPTHM(KT+1) = DEPTHM(KT)+H(KT+1)*0.5
        DO K=KT+2,KMP
          DEPTHB(K) = DEPTHB(K-1)+H(K)
          DEPTHM(K) = DEPTHM(K-1)+(H(K-1)+H(K))*0.5
        END DO

!****** Total active cells

        DO I=CUS(JB),ID
          DO K=KT,KB(I)
            NTAC = NTAC+1
          END DO
        END DO
        NTACMX = NTAC
        NTACMN = NTAC

!****** Wind fetch lengths

        DO I=IU,ID
          FETCHD(I,JB) = FETCHD(I-1,JB)+DLX(I)
        END DO
        DO I=ID,IU,-1
          FETCHU(I,JB) = FETCHU(I+1,JB)+DLX(I)
        END DO
      END DO

!**** Segment heights

      DO I=1,IMP
        DO K=KB(I),2,-1
          HSEG(K,I) = HSEG(K+1,I)+H(K)
        END DO
      END DO

!**** Ending segment and bottom layer for snapshots

      IEPR = MIN(IMP-2,NISNP)
      DO I=1,IEPR
        IPR(I) = IPRI(I)
      END DO
      DO I=1,IEPR
        KBPR = MAX(KB(IPR(I)),KBPR)
      END DO

!**** Transport interpolation multipliers

      DO I=2,IMP-1
        DO K=2,KMP-1

!******** Positive flows

          DLXT = DLX(I-1) 
          IF (K.GT.KB(I-1)) DLXT = DLX(I)
          DLXM         = MIN(DLX(I+1),DLX(I))
          SF1L(K,I)    = (DLX(I+1)+DLX(I))*0.5
          SF2L(K,I,1)  = DLX(I)/(DLX(I)+DLX(I+1))
          SF3L(K,I,1)  = DLX(I)**2
          SF4L(K,I,1)  = DLX(I+1)/(DLX(I)+DLX(I+1))
          SF5L(K,I,1)  = 0.25*(DLXT+2.0*DLX(I)+DLX(I+1))*(DLXT+DLX(I))
          SF6L(K,I,1)  = -0.25*(DLX(I)+DLX(I+1))*(DLXT+DLX(I))
          SF7L(K,I,1)  = 0.25*(DLX(I)+DLX(I+1))*(DLXT+2.0*DLX(I)  &
                         +DLX(I+1))
          SF8L(K,I,1)  = 0.5*(DLX(I)-DLX(I+1))*DLXM
          SF9L(K,I,1)  = 0.5*(DLXT+2.0*DLX(I)-DLX(I+1))*DLXM
          SF10L(K,I,1) = 0.5*(DLXT+3.0*DLX(I))*DLXM
          SF11L(K,I,1) = SF8L(K,I,1)/SF5L(K,I,1)/SF1L(K,I)
          SF12L(K,I,1) = SF9L(K,I,1)/SF6L(K,I,1)/SF1L(K,I)
          SF13L(K,I,1) = SF10L(K,I,1)/SF7L(K,I,1)/SF1L(K,I)
          HTOP         = H(K-1)
          HMID         = H(K)
          HBOT         = H(K+1)
          HMIN         = MIN(HBOT,HMID)
          SF1V(K)      = (HBOT+HMID)*0.5
          SF2V(K,1)    = HMID**2
          SF3V(K,1)    = HMID/(HMID+HBOT)
          SF4V(K,1)    = HBOT/(HMID+HBOT)
          SF5V(K,1)    = 0.25*(HTOP+2.0*HMID+HBOT)*(HTOP+HMID)
          SF6V(K,1)    = -0.25*(HMID+HBOT)*(HTOP+HMID)
          SF7V(K,1)    = 0.25*(HMID+HBOT)*(HTOP+2.0*HMID+HBOT)
          SF8V(K,1)    = 0.5*(HMID-HBOT)*HMIN
          SF9V(K,1)    = 0.5*(HTOP+2.0*HMID-HBOT)*HMIN
          SF10V(K,1)   = 0.5*(HTOP+3.0*HMID)*HMIN

!******** Negative flows

          IF (I.LT.IMP-1) THEN
            DLXT = DLX(I+2)
            IF (K.GT.KB(I+2)) DLXT = DLX(I+1)
            DLXM         = MIN(DLX(I),DLX(I+1))
            SF1L(K,I)    = (DLX(I+1)+DLX(I))*0.5
            SF2L(K,I,2)  = DLX(I+1)/(DLX(I)+DLX(I+1))
            SF3L(K,I,2)  = DLX(I+1)**2
            SF4L(K,I,2)  = DLX(I)/(DLX(I)+DLX(I+1))
            SF5L(K,I,2)  = 0.25*(DLX(I)+2.0*DLX(I+1)+DLXT)*(DLX(I)  &
                           +DLX(I+1))
            SF6L(K,I,2)  = -0.25*(DLX(I+1)+DLXT)*(DLX(I)+DLX(I+1))
            SF7L(K,I,2)  = 0.25*(DLX(I)+2.0*DLX(I+1)+DLXT)*(DLX(I+1)  &
                           +DLXT)
            SF8L(K,I,2)  = -0.5*(3.0*DLX(I+1)+DLXT)*DLXM
            SF9L(K,I,2)  = 0.5*(DLX(I)-2.0*DLX(I+1)-DLXT)*DLXM
            SF10L(K,I,2) = 0.5*(DLX(I)-DLX(I+1))*DLXM
            SF11L(K,I,2) = SF8L(K,I,2)/SF5L(K,I,2)/SF1L(K,I)
            SF12L(K,I,2) = SF9L(K,I,2)/SF6L(K,I,2)/SF1L(K,I)
            SF13L(K,I,2) = SF10L(K,I,2)/SF7L(K,I,2)/SF1L(K,I)
          END IF
          HTOP = H(K)
          HMID = H(K+1)
          IF (K.LT.KB(I)) THEN
            HBOT = H(K+2)
            IF (K.EQ.KB(I)-1) HBOT = H(K+1)
            HMIN       = MIN(HTOP,HMID)
            SF1V(K)    = (HMID+HTOP)*0.5
            SF2V(K,2)  = HMID**2
            SF3V(K,2)  = HMID/(HTOP+HMID)
            SF4V(K,2)  = HTOP/(HTOP+HMID)
            SF5V(K,2)  = 0.25*(HTOP+2.0*HMID+HBOT)*(HTOP+HMID)
            SF6V(K,2)  = -0.25*(HMID+HBOT)*(HTOP+HMID)
            SF7V(K,2)  = 0.25*(HTOP+2.0*HMID+HBOT)*(HMID+HBOT)
            SF8V(K,2)  = -0.5*(3.0*HMID+HBOT)*HMIN
            SF9V(K,2)  = 0.5*(HTOP-2.0*HMID-HBOT)*HMIN
            SF10V(K,2) = 0.5*(HTOP-HMID)*HMIN
          END IF
        END DO
      END DO

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                          Task 1.4.4: Initial Conditions                  **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

!**** Water surface elevations

      IF (WATER_BALANCE) THEN
        OPEN (QWS,FILE='qws.opt',STATUS='UNKNOWN')
        WRITE (QWS,9000)
      END IF

!**** Vertical profiles

      IF (OPEN_VPR) THEN

!****** Temperature

        OPEN (VPR,FILE=VPRFN,STATUS='OLD')
        READ (VPR,*)
        IF (VERT_TEMP) READ (VPR,1030) (TVP(K),K=KT,KBMAX)
        IF (CONSTITUENTS) THEN
          DO JC=1,NCP
            IF (VERT_CONC(JC)) READ (VPR,1030) (CVP(K,JC),K=KT,KBMAX)
          END DO
        END IF
      END IF

!**** Longitudinal/vertical initial profiles

      IF (OPEN_LPR) THEN
        OPEN (LPR,FILE=LPRFN,STATUS='OLD')
        READ (LPR,*)
      END IF
      DO JB=1,NBP
        IU = CUS(JB)
        ID = DS(JB)
!        IF (.NOT.RESTART_IN) THEN

!******** Temperature

          DO I=IU,ID
            IF (LONG_TEMP) READ (LPR,1190) (T1(K,I),K=KT,KB(I))
            DO K=KT,KB(I)
              IF (ISO_TEMP)  T1(K,I) = T2I
              IF (VERT_TEMP) T1(K,I) = TVP(K)
              T2(K,I) = T1(K,I)
            END DO
          END DO

!******** Constituents

          DO JC=1,NAC0
            DO I=IU,ID
              IF (LONG_CONC(CN0(JC))) THEN
                READ (LPR,1190) (C2(K,I,CN0(JC)),K=KT,KB(I))
              END IF
              DO K=KT,KB(I)
                IF (ISO_CONC(CN0(JC)))  C2(K,I,CN0(JC)) = C2I(CN0(JC))
                IF (VERT_CONC(CN0(JC))) C2(K,I,CN0(JC)) = CVP(K,CN0(JC))
                C1(K,I,CN0(JC))  = C2(K,I,CN0(JC))
                C1S(K,I,CN0(JC)) = C1(K,I,CN0(JC))
              END DO
            END DO
          END DO

!******** Energy

          IF (ENERGY_BALANCE) THEN
            DO I=IU,ID
              EIBR(JB) = EIBR(JB)+T2(KT,I)*DLX(I)*BHKT2(I)
              DO K=KT+1,KB(I)
                EIBR(JB) = EIBR(JB)+T2(K,I)*DLX(I)*BH(K,I)
              END DO
            END DO
          END IF

!******** Constituent mass

          DO JC=1,NAC
            JAC = CN(JC)
            DO I=CUS(JB),DS(JB)
              CMBRT(JAC,JB) = CMBRT(JAC,JB)+C2(KT,I,JAC)*DLX(I)  &
                              *BHKT2(I)
              DO K=KT+1,KB(I)
                CMBRT(JAC,JB) = CMBRT(JAC,JB)+C2(K,I,JAC)*DLX(I)  &
                                *BH(K,I)
              END DO
            END DO
          END DO

!******** Ice cover

          IF (ICE_CALC) THEN
            DO I=IU,ID
              ICETH(I) = ICETHI
              ICE(I)   = ICETH(I).GT.0.0
            END DO
          END IF

!******** Vertical eddy viscosity

          IUT = IU
          IDT = ID-1
          DO I=IUT,IDT
            DO K=KT,KBMIN(I)-1
              AZ(K,I)  = AZMIN
              SAZ(K,I) = AZMIN
            END DO
          END DO
!        END IF

!****** Density

        DO I=IU,ID
          DO K=KT,KB(I)
            IF (CONSTITUENTS) THEN
              RHO(K,I) = DENSITY (T2(K,I),SS(K,I),TDS(K,I))
            ELSE
              RHO(K,I) = DENSITY (T2(K,I),0.0,0.0)
            END IF
          END DO
        END DO

!****** Horizontal diffusivities

        DO I=IU,ID-1
          DO K=KT,KBMIN(I)
            DX(K,I) = DXI
          END DO
        END DO

!****** Saved velocities and eddy viscosities

        DO I=IU-1,ID+1
          DO K=KT,MAX(KB(IU),KB(I))
            SU(K,I)  = U(K,I)
            SW(K,I)  = W(K,I)
            SAZ(K,I) = AZ(K,I)
          END DO
        END DO
      END DO

!**** Density related constants

      RHOWCP = RHOW*CP
      RHOICP = RHOI*CP
      RHORL1 = RHOI*RL1
      DO I=2,IMP-1
        DLXRHO(I) = 0.5/(DLXR(I)*RHOW)
      END DO

      IF (OPEN_VPR) CLOSE (VPR)
      IF (OPEN_LPR) CLOSE (LPR)
      CALL GREGORIAN_DATE (YEAR)

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                                Task 1.5: Outputs                         **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

!**** Open output files (contains UNIX FORTRAN specific I/O)

        IF (SNAPSHOT)    OPEN (SNP,FILE=SNPFN,STATUS='UNKNOWN' )
        IF (TIME_SERIES) OPEN (TSR,FILE=TSRFN,STATUS='UNKNOWN')
        IF (VECTOR)      OPEN (VPL,FILE=VPLFN,STATUS='UNKNOWN')
        IF (PROFILE)     OPEN (PRF,FILE=PRFFN,STATUS='UNKNOWN')
        IF (SPREADSHEET) OPEN (SPR,FILE=SPRFN,STATUS='UNKNOWN')
        IF (CONTOUR)     OPEN (CPL,FILE=CPLFN,STATUS='UNKNOWN')
!
      OPEN (WRN,FILE='w2.wrn',STATUS='UNKNOWN')
      OPEN (ERR,FILE='w2.err',STATUS='UNKNOWN')

!**** Output files

      IF (SNAPSHOT) THEN
          WRITE (SNP,'(A80)') ESC//'E'//ESC//'&l6.0c7E'//ESC//  &
                              '(s0p16.67h8.5v0s0b0T'//ESC//'(10U'  &
                              //ESC//'&a8L'
      END IF

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                              Task 2: Calculations                        **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      DO WHILE (.NOT.END_RUN)
        IF (JDAY.GE.NXTVD)  CALL TIME_VARYING_DATA  (JDAY,NXTVD)
        IF (INTERPOLATE)    CALL INTERPOLATE_INPUTS (JDAY)
        DLTTVDS = DLTTVD
        DLTTVD  = INT((NXTVD-JDAY)*86400.0)+1.0

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                       Task 2.1: Hydrodynamic Sources/Sinks               **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

        DO JB=1,NBP
          IU = CUS(JB)
          ID = DS(JB)
          IF (SEL_WITHDRAWAL(JB)) CALL SELECTIVE_WITHDRAWAL
        END DO

!****** Timestep violation entry point

10010   CONTINUE
        DO JB=1,NBP
          IU = CUS(JB)
          ID = DS(JB)
          IF (EVAPORATION) THEN
            EVBR(JB) = 0.0
            FW       = 9.2+0.46*WIND*WIND
            DO I=IU,ID
              TM    = (T2(KT,I)+TDEW)*0.5
              VPTG  = 0.35+0.015*TM+0.0012*TM*TM
              EV(I) = VPTG*(T2(KT,I)-TDEW)*FW*B(KT,I)*DLX(I)/2.45E9
              IF (EV(I).LT.0.0.OR.ICE(I)) EV(I) = 0.0
              QSS(KT,I) = QSS(KT,I)-EV(I)
              EVBR(JB)  = EVBR(JB)+EV(I)
            END DO
          END IF
          IF (PRECIPITATION) THEN
            QPRBR(JB) = 0.0
            DO I=IU,ID
              QPR(I)    = PR(JB)*B(KTI(I),I)*DLX(I)
              QPRBR(JB) = QPRBR(JB)+QPR(I)
              QSS(KT,I) = QSS(KT,I)+QPR(I)
            END DO
          END IF
          IF (TRIBUTARIES) THEN
            DO JT=1,NTR

!****!******* Inflow fractions

              IF (JB.EQ.JBTR(JT)) THEN
                I = MAX(ITR(JT),IU)
                DO K=KT,KB(I)
                  QTRF(K,JT) = 0.0
                END DO
                IF (PLACE_QTR(JT)) THEN

!****!****!****** Inflow layer

                  K     = KT
                  RHOTR = DENSITY (TTR(JT),CTR(2,JT),CTR(4,JT))
                  DO WHILE (RHOTR.GT.RHO(K,I).AND.K.LT.KB(I))
                    K = K+1
                  END DO
                  KTTR(JT) = K
                  KBTR(JT) = K

!****!****!****** Layer inflows

                  VQTR  = QTR(JT)*DLT
                  VQTRI = VQTR
                  QTRFR = 1.0
                  INCR  = -1
                  DO WHILE (QTRFR.GT.0.0)
                    IF (K.LE.KB(I)) THEN
                      VOL = BH(K,I)*DLX(I)
                      IF (K.EQ.KT) VOL = BHKT2(I)*DLX(I)
                      IF (VQTR.GT.0.5*VOL) THEN
                        QTRF(K,JT) = 0.5*VOL/VQTRI
                        QTRFR      = QTRFR-QTRF(K,JT)
                        VQTR       = VQTR-QTRF(K,JT)*VQTRI
                        IF (K.EQ.KT) THEN
                          K    = KBTR(JT)
                          INCR = 1
                        END IF
                      ELSE
                        QTRF(K,JT) = QTRFR
                        QTRFR      = 0.0
                      END IF
                      IF (INCR.LT.0) KTTR(JT) = K
                      IF (INCR.GT.0) KBTR(JT) = MIN(KB(I),K)
                      K = K+INCR
                    ELSE
                      QTRF(KT,JT) = QTRFR
                      QTRFR       = 0.0
                    END IF
                  END DO
                ELSE
                  IF (SPECIFY_QTR(JT)) THEN
                    KTTR(JT) = 2
                    DO WHILE (EL(KTTR(JT)).GT.ELTRT(JT))
                      KTTR(JT) = KTTR(JT)+1
                    END DO
                    KTTR(JT) = KTTR(JT)-1
                    KBTR(JT) = KB(I)
                    DO WHILE (EL(KBTR(JT)).LT.ELTRB(JT))
                      KBTR(JT) = KBTR(JT)-1
                    END DO
                  ELSE
                    KTTR(JT) = KT
                    KBTR(JT) = KB(I)
                  END IF
                  KTTR(JT) = MAX(KT,KTTR(JT))
                  KBTR(JT) = MIN(KB(I),KBTR(JT))
                  BHSUM    = 0.0
                  DO K=KTTR(JT),KBTR(JT)
                    BHT = BH(K,I)
                    IF (K.EQ.KT) BHT = BHKT2(I)
                    BHSUM = BHSUM+BHT
                  END DO
                  DO K=KTTR(JT),KBTR(JT)
                    BHT = BH(K,I)
                    IF (K.EQ.KT) BHT = BHKT2(I)
                    QTRF(K,JT) = BHT/BHSUM
                  END DO
                END IF
                DO K=KTTR(JT),KBTR(JT)
                  QSS(K,I) = QSS(K,I)+QTR(JT)*QTRF(K,JT)
                END DO
              END IF
            END DO
          END IF
          IF (WITHDRAWALS) THEN
            DO JW=1,NWD
              IF (JB.EQ.JBWD(JW)) THEN
                I        = MAX(CUS(JBWD(JW)),IWD(JW))
                K        = MAX(KT,KWD(JW))
                QSS(K,I) = QSS(K,I)-QWD(JW)
              END IF
            END DO
          END IF
          IF (UH_INTERNAL(JB)) THEN
            DO K=KT,KB(IU-1)
              QSS(K,UHS(JB)) = QSS(K,UHS(JB))-QUH2(K,JB)/DLT
            END DO
          END IF
          IF (DH_INTERNAL(JB)) THEN
            DO K=KT,KB(ID+1)
              QSS(K,DHS(JB)) = QSS(K,DHS(JB))+QDH2(K,JB)/DLT
            END DO
          END IF
        END DO

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                     Task 2.2: Hydrodynamic Calculations                  **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

        DO JB=1,NBP
          IU = CUS(JB)
          ID = DS(JB)

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*        Task 2.2.1: Boundary Concentrations, Temperatures, and Densities  **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

          IUT = IU
          IDT = ID
          IF (UP_FLOW(JB)) THEN
            DO K=KT,KB(IU)
              DO JC=1,NAC
                C1S(K,IU-1,CN(JC)) = CIN(CN(JC),JB)
              END DO
              IF (QIN(JB).GT.0.0) THEN
                T1(K,IU-1) = TIN(JB)
                T2(K,IU-1) = TIN(JB)
              ELSE
                T1(K,IU-1) = T1(K,IU)
                T2(K,IU-1) = T2(K,IU)
              END IF
            END DO
          END IF
          IF (DN_FLOW(JB)) THEN
            DO K=KT,KB(IU)
              DO JC=1,NAC
                C1S(K,ID+1,CN(JC)) = C1S(K,ID,CN(JC))
              END DO
              T1(K,ID+1) = T2(K,ID)
              T2(K,ID+1) = T2(K,ID)
            END DO
          END IF
            IUT = IU-1
            IF (UH_INTERNAL(JB)) THEN
              DO K=KT,KB(IUT)
                DO JC=1,NAC
                  C1S(K,IUT,CN(JC)) = C1S(K,UHS(JB),CN(JC))
                  C1(K,IUT,CN(JC))  = C1S(K,UHS(JB),CN(JC))
                  C2(K,IUT,CN(JC))  = C1S(K,UHS(JB),CN(JC))
                END DO
                T1(K,IUT)  = T2(K,UHS(JB))
                T2(K,IUT)  = T2(K,UHS(JB))
                RHO(K,IUT) = RHO(K,UHS(JB))
              END DO
            END IF

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                          Task 2.2.2: Momentum Terms                      **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

!******** Densities

          DO I=IU,ID
            DO K=KT,KB(I)
              IF (CONSTITUENTS) THEN
                RHO(K,I) = DENSITY (T2(K,I),SS(K,I),TDS(K,I))
              ELSE
                RHO(K,I) = DENSITY (T2(K,I),0.0,0.0)
              END IF
            END DO
          END DO

!******** Density pressures

          DO I=IUT,IDT
            P(KT,I) = RHO(KT,I)*G*H(KT)
            DO K=KT+1,KB(I)
              P(K,I) = P(K-1,I)+RHO(K,I)*G*H(K)
            END DO
          END DO

!******** Horizontal density gradients

          DO I=IUT,IDT-1
            HDG(KT,I) = DLXRHO(I)*0.5*(BKT(I)+BKT(I+1))*((HKT2(I+1)  &
                        *P(KT,I+1))-(HKT2(I)*P(KT,I)))
            DO K=KT+1,KBMIN(I)
              HDG(K,I) = DLXRHO(I)*BHR(K,I)*((P(K-1,I+1)-P(K-1,I))  &
                         +(P(K,I+1)-P(K,I)))
            END DO
          END DO

!******** Wind drag coefficient and surface shear stress

          CZ = 0.0
          IF (WIND.GE.1.0)  CZ = 0.0005*SQRT(WIND)
          IF (WIND.GE.15.0) CZ = 0.0026
          SSC = RHOA*CZ*WIND**2/RHOW

!******** Vertical eddy viscosities

          DO I=IUT,IDT-1
            FETCH  = FETCHD(I,JB)
            SSCCOS = SSC*COS(PHI-PHI0(I))*ICESW(I)
            SSCSIN = SSC*ABS(SIN(PHI-PHI0(I)))*ICESW(I)
            IF (COS(PHI-PHI0(I)).LT.0.0) FETCH = FETCHU(I,JB)
            WWT = 0.0
            IF (WIND.NE.0.0) THEN
              WWT = 6.95E-2*(FETCH**0.233)*ABS(WIND)**0.534
            END IF
            ST(KT,I) = SSCCOS*BR(KTI(I),I)
            IF (.NOT.ONE_LAYER(I)) THEN
              DFC        = (-8.0*3.14159*3.14159)/(G*WWT*WWT+NONZERO)
              EXPDF      = EXP(MAX(DFC*DEPTHB(KT),-30.0))
              ST(KT+1,I) = SSCCOS*EXPDF*BR(KT+1,I)
              SHEARS     = ((U(KT+1,I)-U(KT,I))/((AVHKT(I)+AVHKT(I+1))  &
                           *0.5))**2
              AZ0        = 0.4*HMAX*HMAX*SQRT(SHEARS+(SSCSIN  &
                           *EXPDF/AZ(KT,I))**2)+AZMIN
              RIAZ0      = LOG(AZ0/AZMAX)/1.5  
              BUOY       = (RHO(KT+1,I)-RHO(KT,I)+RHO(KT+1,I+1)  &
                           -RHO(KT,I+1))/(2.0*AVHKT(I))
              RI         = G*BUOY/(RHOW*SHEARS+NONZERO)
              RIAZ1      = MAX(RI,RIAZ0)
              RIAZ1      = MIN(RIAZ1,10.0)
              EXPRAZ     = EXP(-1.5*RIAZ1)
              AZ(KT,I)   = MAX(AZMIN,AZ0*EXPRAZ+AZMIN*(1.0-EXPRAZ))
              DZ(KT,I)   = MAX(DZMIN,FRAZDZ*(AZ0*EXPRAZ+DZMIN  &
                           *(1.0-EXPRAZ)))
              KBT        = KBMIN(I)
              DO K=KT+2,KBT
                EXPDF     = EXP(MAX(DFC*DEPTHB(K-1),-30.0))
                ST(K,I)   = SSCCOS*EXPDF*BR(K,I)
                SHEARS    = ((U(K,I)-U(K-1,I))/AVH(K-1))**2
                AZ0       = 0.4*HMAX*HMAX*SQRT(SHEARS+(SSCSIN  &
                            *EXPDF/AZ(K-1,I))**2)+AZMIN
                RIAZ0     = LOG(AZ0/AZMAX)/1.5
                BUOY      = (RHO(K,I)-RHO(K-1,I)+RHO(K,I+1)  &
                            -RHO(K-1,I+1))/(2.0*AVH(K-1))
                RI        = G*BUOY/(RHOW*SHEARS+NONZERO)
                RIAZ1     = MAX(RI,RIAZ0)
                RIAZ1     = MIN(RIAZ1,10.0)
                EXPRAZ    = EXP(-1.5*RIAZ1)
                AZ(K-1,I) = MAX(AZMIN,AZ0*EXPRAZ+AZMIN*(1.0-EXPRAZ))
                DZ(K-1,I) = MAX(DZMIN,FRAZDZ*(AZ0*EXPRAZ+DZMIN  &
                            *(1.0-EXPRAZ)))
              END DO
              SB(KBT,I) = SSCCOS*EXPDF*(BR(KBT-1,I)+BR(KBT,I))*0.5
            END IF
          END DO

!******** Average eddy diffusivities

          DO K=KT,KB(IDT)-1
            DZ(K,IDT) = DZ(K,IDT-1)
          END DO
          DO I=IUT,IDT-1
            DO K=KT,KB(I)-1
              IF (K.GE.KBMIN(I)) THEN
                IF (KB(I-1).GE.KB(I).AND.I.NE.IUT) THEN
                  DZ(K,I) = DZ(K,I-1)
                ELSE
                  DZ(K,I) = DZMIN
                END IF
              ELSE
                DZ(K,I) = (DZ(K,I)+DZ(K+1,I))*0.5
              END IF
            END DO
          END DO

!******** Density inversions

          DO I=IUT,IDT
            DO K=KT,KB(I)-1
              DZQ(K,I) = DZ(K,I)
              IF (RHO(K,I).GT.RHO(K+1,I)) DZ(K,I) = DZMAX
            END DO
          END DO

!******** Shear stresses

          GC2 = G/(CHEZY*CHEZY)
          DO I=IUT,IDT-1
            KBT = KBMIN(I)
            IF (.NOT.ONE_LAYER(I)) THEN
              ST(KT+1,I) = ST(KT+1,I)+AZ(KT,I)*(BR(KT,I)+BR(KT+1,I))  &
                           *0.5*(U(KT,I)-U(KT+1,I))/((AVHKT(I)  &
                           +AVHKT(I+1))*0.5)
            END IF
            DO K=KT+2,KBT
              ST(K,I) = ST(K,I)+AZ(K-1,I)*(BR(K-1,I)+BR(K,I))*0.5  &
                        *(U(K-1,I)-U(K,I))/AVH(K-1)
            END DO
            DO K=KT,KBT-1
              SB(K,I) = ST(K+1,I)+GC2*(BR(K,I)-BR(K+1,I))*U(K,I)  &
                        *ABS(U(K,I))
            END DO
            SB(KBT,I) = SB(KBT,I)+GC2*BR(KBT,I)*U(KBT,I)*ABS(U(KBT,I))
          END DO

!******** Horizontal momentum

          DO I=IU,ID-1
            UDR        = (1.0+SIGN(1.0,(U(KT,I)+U(KT,I+1))*0.5))*0.5
            UDL        = (1.0+SIGN(1.0,(U(KT,I)+U(KT,I-1))*0.5))*0.5
            ADMX(KT,I) = (BHKT2(I+1)*(U(KT,I+1)+U(KT,I))*0.5  &
                         *(UDR*U(KT,I)+(1.0-UDR)*U(KT,I+1))-BHKT2(I)  &
                         *(U(KT,I)+U(KT,I-1))*0.5*(UDL*U(KT,I-1)  &
                         +(1.0-UDL)*U(KT,I)))/DLXR(I)
            DM(KT,I)   = AX*(BHKT2(I+1)*(U(KT,I+1)-U(KT,I))/DLX(I+1)  &
                         -BHKT2(I)*(U(KT,I)-U(KT,I-1))/DLX(I))/DLXR(I)
            DO K=KT+1,KBMIN(I)  
              UDR       = (1.0+SIGN(1.0,(U(K,I)+U(K,I+1))*0.5))*0.5
              UDL       = (1.0+SIGN(1.0,(U(K,I)+U(K,I-1))*0.5))*0.5
              ADMX(K,I) = (BH(K,I+1)*(U(K,I+1)+U(K,I))*0.5*(UDR*U(K,I)  &
                          +(1.0-UDR)*U(K,I+1))-BH(K,I)*(U(K,I)  &
                          +U(K,I-1))*0.5*(UDL*U(K,I-1)+(1.0-UDL)  &
                          *U(K,I)))/DLXR(I)
              DM(K,I)   = AX*(BH(K,I+1)*(U(K,I+1)-U(K,I))/DLX(I+1)  &
                          -BH(K,I)*(U(K,I)-U(K,I-1))/DLX(I))/DLXR(I)
            END DO
          END DO

!******** Vertical momentum

          DO I=IU,ID-1
            DO K=KT,KB(I)-1
              UD        = (1.0+SIGN(1.0,(W(K,I+1)+W(K,I))*0.5))*0.5
              ADMZ(K,I) = (BR(K,I)+BR(K+1,I))*0.5*(W(K,I+1)+W(K,I))  &
                          *0.5*(UD*U(K,I)+(1.0-UD)*U(K+1,I))
            END DO
          END DO

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                      Task 2.2.3: Water Surface Elevation                 **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

!******** Tridiagonal coefficients

          DO I=IU,ID-1
            BHRHO(I) = BHKT2(I+1)/RHO(KT,I+1)+BHKT2(I)/RHO(KT,I)
            DO K=KT+1,KBMIN(I)
              BHRHO(I) = BHRHO(I)+(BH(K,I+1)/RHO(K,I+1)+BH(K,I)  &
                         /RHO(K,I))
            END DO
            D(I) = U(KT,I)*BHRKT2(I)-U(KT,I-1)*BHRKT2(I-1)-QSS(KT,I)
            F(I) = -SB(KT,I)+ST(KT,I)-ADMX(KT,I)+DM(KT,I)-HDG(KT,I)
            DO K=KT+1,KB(I)
              D(I) = D(I)+(U(K,I)*BHR(K,I)-U(K,I-1)*BHR(K,I-1)  &
                     -QSS(K,I))
              F(I) = F(I)+(-SB(K,I)+ST(K,I)-ADMX(K,I)+DM(K,I)-HDG(K,I))
            END DO
          END DO
          D(IU) = U(KT,IU)*BHRKT2(IU)-QSS(KT,IU)
          DO K=KT+1,KB(IU)
            D(IU) = D(IU)+(U(K,IU)*BHR(K,IU)-QSS(K,IU))
          END DO

!******** Boundary tridiagonal coefficients

          IF (UP_FLOW(JB)) D(IU) = D(IU)-QIN(JB)
          IF (DN_FLOW(JB)) THEN
            D(ID) = -U(KT,ID-1)*BHRKT2(ID-1)-QSS(KT,ID)
            DO K=KT+1,KB(ID)
              D(ID) = D(ID)-(U(K,ID-1)*BHR(K,ID-1)+QSS(K,ID))
            END DO
            DO K=KT,KB(ID)
              D(ID) = D(ID)+QOUT(K,JB)
            END DO
          END IF
        END DO
        DO JB=1,NBP
          IU = CUS(JB)
          ID = DS(JB)

!******** Boundary surface elevations

          IF (UH_INTERNAL(JB)) Z(IU-1) = Z(UHS(JB))
          IF (DH_INTERNAL(JB)) Z(ID+1) = Z(DHS(JB))

!******** Tridiagonal coefficients

          DO I=IU,ID
            A(I) = -RHO(KT,I-1)*G*DLT**2*BHRHO(I-1)*0.5/DLXR(I-1)
            C(I) = -RHO(KT,I+1)*G*DLT**2*BHRHO(I)*0.5/DLXR(I)
            V(I) = RHO(KT,I)*G*DLT**2*(BHRHO(I)*0.5/DLXR(I)  &
                   +BHRHO(I-1)*0.5/DLXR(I-1))+DLX(I)*B(KTI(I),I)
            D(I) = DLT*(D(I)+DLT*(F(I)-F(I-1)))+DLX(I)*B(KTI(I),I)*Z(I)
          END DO

!******** Implicit water surface elevation

          BTA(IU) = V(IU)
          GMA(IU) = D(IU)/BTA(IU)
          DO I=IU+1,ID
            BTA(I) = V(I)-A(I)*C(I-1)/BTA(I-1)
            GMA(I) = (D(I)-A(I)*GMA(I-1))/BTA(I)
          END DO
          Z(ID) = GMA(ID)
          DO K=1,ID-IU
            I    = ID-K
            Z(I) = GMA(I)-C(I)*Z(I+1)/BTA(I)
          END DO
          IF (UP_FLOW(JB)) Z(IU-1) = Z(IU)
          IF (DN_FLOW(JB)) Z(ID+1) = Z(ID)

!******** Updated surface layer and related variables

          DO I=IU-1,ID+1
            IF (EL(KT)-Z(I).GT.EL(KTI(I))) THEN
              DO WHILE (EL(KT)-Z(I).GT.EL(KTI(I)).AND.KTI(I).NE.2)
                Z(I)   = EL(KT)-EL(KTI(I))-(EL(KT)-EL(KTI(I))-Z(I))  &
                         *(B(KTI(I),I)/B(KTI(I)-1,I))
                KTI(I) = KTI(I)-1
              END DO
            ELSE IF (EL(KT)-Z(I).LT.EL(KTI(I)+1)) THEN
              DO WHILE (EL(KT)-Z(I).LT.EL(KTI(I)+1))
                Z(I)   = EL(KT)-EL(KTI(I)+1)-(EL(KT)-EL(KTI(I)+1)  &
                         -Z(I))*(B(KTI(I),I)/B(KTI(I)+1,I))
                KTI(I) = KTI(I)+1
              END DO
            END IF
            HKT1(I)  = H(KT)-Z(I)
            AVHKT(I) = (HKT1(I)+H(KT+1))*0.5
            BHKT1(I) = B(KTI(I),I)*(EL(KT)-Z(I)-EL(KTI(I)+1))
            DO K=KTI(I)+1,KT
              BHKT1(I) = BHKT1(I)+BH(K,I)
            END DO
            BKT(I) = BHKT1(I)/HKT1(I)
          END DO
          DO I=IU-1,ID
            BHRKT1(I) = (BHKT1(I)+BHKT1(I+1))*0.5
          END DO
          BHRKT1(ID+1) = BHKT1(ID+1)
          DLVOL(JB)    = 0.0
          DO I=IU,ID
            DLVOL(JB) = DLVOL(JB)+(BHKT1(I)-BHKT2(I))*DLX(I)
          END DO

!******** Layer bottom and middle depths

          DEPTHB(KT)   = HKT1(I)
          DEPTHM(KT)   = HKT1(I)*0.5
          DEPTHB(KT+1) = DEPTHB(KT)+H(KT+1)
          DEPTHM(KT+1) = DEPTHM(KT)+H(KT+1)*0.5
          DO K=KT+2,KMP-1
            DEPTHB(K) = DEPTHB(K-1)+H(K)
            DEPTHM(K) = DEPTHM(K-1)+(H(K-1)+H(K))*0.5
          END DO

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                      Task 2.2.4: Longitudinal Velocities                 **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

          IUT = IU
          IDT = ID

!******** Pressures

          DO I=IUT,IDT
            P(KT,I) = RHO(KT,I)*G*HKT1(I)
            DO K=KT+1,KB(I)
              P(K,I) = P(K-1,I)+RHO(K,I)*G*H(K)
            END DO
          END DO

!******** Horizontal pressure gradients

          DO I=IUT,IDT-1
            HPG(KT,I) = (DLXRHO(I)*0.5)*(BKT(I)+BKT(I+1))*((HKT1(I+1)  &
                        *P(KT,I+1))-(HKT1(I)*P(KT,I)))
            DO K=KT+1,KBMIN(I)
              HPG(K,I) = DLXRHO(I)*BHR(K,I)*((P(K-1,I+1)-P(K-1,I))  &
                         +(P(K,I+1)-P(K,I)))
            END DO
          END DO

!******** Boundary horizontal velocities

          IF (UP_FLOW(JB)) THEN
            DO K=KT,KB(IU)
              QINF(K,JB) = 0.0
            END DO
            IF (PLACE_QIN) THEN

!****!******* Inflow layer

              K     = KT
              RHOIN = DENSITY (TIN(JB),CIN(2,JB),CIN(4,JB))
              DO WHILE (RHOIN.GT.RHO(K,IU).AND.K.LT.KB(IU))
                K = K+1
              END DO
              KTQIN(JB) = K
              KBQIN(JB) = K

!****!******* Layer inflows

              VQIN  = QIN(JB)*DLT
              VQINI = VQIN
              QINFR = 1.0
              INCR  = -1
              DO WHILE (QINFR.GT.0.0)
                VOL  = BH(K,IU)*DLX(IU)
                IF (K.EQ.KT) VOL = BHKT1(IU)*DLX(IU)
                IF (K.LE.KB(IU)) THEN
                  IF (VQIN.GT.0.5*VOL) THEN
                    QINF(K,JB) = 0.5*VOL/VQINI
                    QINFR      = QINFR-QINF(K,JB)
                    VQIN       = VQIN-QINF(K,JB)*VQINI
                    IF (K.EQ.KT) THEN
                      K    = KBQIN(JB)
                      INCR = 1
                    END IF
                  ELSE
                    QINF(K,JB) = QINFR
                    QINFR      = 0.0
                  END IF
                  IF (INCR.LT.0) KTQIN(JB) = K
                  IF (INCR.GT.0) KBQIN(JB) = MIN(KB(IU),K)
                  K = K+INCR
                ELSE
                  QINF(KT,JB) = QINFR
                  QINFR       = 0.0
                END IF
              END DO
            ELSE
              BHSUM = BHKT1(IU)
              DO K=KT+1,KB(IU)
                BHSUM = BHSUM+BH(K,IU)
              END DO
              QINF(KT,JB) = BHKT1(IU)/BHSUM
              DO K=KT+1,KB(IU)
                QINF(K,JB) = BH(K,IU)/BHSUM
              END DO
              KTQIN(JB) = KT
              KBQIN(JB) = KB(IU)
            END IF
            U(KT,IU-1) = QINF(KT,JB)*QIN(JB)/BHRKT1(IU-1)
            DO K=KT+1,KB(IU)
              U(K,IU-1) = QINF(K,JB)*QIN(JB)/BHR(K,IU-1)
            END DO
          END IF
          IF (DN_FLOW(JB)) THEN
            DO K=KT,KB(ID)
              BHRT = BHR(K,ID)
              IF (K.EQ.KT) BHRT = BHRKT1(ID)
              U(K,ID) = QOUT(K,JB)/BHRT
            END DO
          END IF

!******** Horizontal velocities

          DO I=IU,ID-1
            U(KT,I) = (BHRKT2(I)*U(KT,I)+DLT*(-SB(KT,I)+ST(KT,I)  &
                      -ADMZ(KT,I)+DM(KT,I)-ADMX(KT,I)-HPG(KT,I)))  &
                      /BHRKT1(I)
            DO K=KT+1,KBMIN(I)
              U(K,I) = U(K,I)+DLT/BHR(K,I)*(-SB(K,I)+ST(K,I)-ADMZ(K,I)  &
                       +ADMZ(K-1,I)-ADMX(K,I)+DM(K,I)-HPG(K,I))
            END DO
          END DO

!******** Corrected horizontal velocities

            IS   = IU-1
            IE   = ID
            INCR = 1
            IF (DN_FLOW(JB)) IE = ID-1
            Q(IS) = U(KT,IS)*BHRKT1(IS)
            DO K=KT+1,KB(IU)
              Q(IS) = Q(IS)+U(K,IS)*BHR(K,IS)
            END DO
          QC(IS) = Q(IS)
          DO I=IS+INCR,IE,INCR
            QSSUM(I) = QSS(KT,I)
            DO K=KT+1,KB(I)
              QSSUM(I) = QSSUM(I)+QSS(K,I)
            END DO
            BHRSUM = BHRKT1(I)
            Q(I)   = U(KT,I)*BHRKT1(I)
            DO K=KT+1,KBMIN(I)
              BHRSUM = BHRSUM+BHR(K,I)
              Q(I)   = Q(I)+U(K,I)*BHR(K,I)
            END DO
!            IF (UP_HEAD(JB)) THEN
!              QC(I) = QC(I+1)+(BHKT1(I+1)-BHKT2(I+1))*DLX(I+1)/DLT  &
!                      -QSSUM(I+1)
!            ELSE
              QC(I) = QC(I-1)-(BHKT1(I)-BHKT2(I))*DLX(I)/DLT+QSSUM(I)
!            END IF
            DO K=KT,KBMIN(I)
              U(K,I) = U(K,I)+(QC(I)-Q(I))/BHRSUM
            END DO
          END DO


!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                        Task 2.2.5: Vertical Velocities                   **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

          DO I=IU,ID
            DO K=KB(I)-1,KT,-1
              WT1    = W(K+1,I)*BB(K+1,I)
              WT2    = (BHR(K+1,I)*U(K+1,I)-BHR(K+1,I-1)*U(K+1,I-1)  &
                       -QSS(K+1,I))/DLX(I)
              W(K,I) = (WT1+WT2)/BB(K,I)
            END DO
          END DO
        END DO

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                           Task 2.2.6: Autostepping                       **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

        DO JB=1,NBP
          DO I=CUS(JB),DS(JB)
            IF (HKT1(I).LT.0.0) THEN
              WRITE (ERR,6020) JDAY,I,Z(I),HSEG(KT,I)-Z(I)
              STOP '***SEVERE*** computational error - see "w2.err"'
            END IF
            TAU1   = 2.0*MAX(AX,DXI)/(DLX(I)*DLX(I))
            TAU2   = 2.0*AZ(KT,I)/(HKT1(I)*HKT1(I))
            CELRTY = SQRT((ABS(RHO(KB(I),I)-RHO(KT,I)))/1000.0*G  &
                     *DEPTHB(KB(I))*0.5)
            QTOT   = (ABS(U(KT,I))*BHRKT1(I)+ABS(U(KT,I-1))*BHRKT1(I-1)  &
                     +ABS(W(KT,I))*BB(KT,I)*DLX(I)+DLX(I)*ABS(BHKT2(I)  &
                     -BHKT1(I))/DLT+ABS(QSS(KT,I)))*0.5
            DLTCAL = 1.0/((QTOT/BHKT1(I)+CELRTY)/DLX(I)+TAU1+TAU2)
            IF (DLTCAL.LT.CURMAX) THEN
              KLOC   = KT
              ILOC   = I
              CURMAX = INT(DLTCAL)
              IF (DLTF(DLTDP)*CURMAX.LT.MINDLT) THEN
                KMIN = KT
                IMIN = I
              END IF
            END IF
            DO K=KT+1,KB(I)
              TAU2   = 2.0*AZ(K,I)/(H(K)*H(K))
              QTOT   = (ABS(U(K,I))*BHR(K,I)+ABS(U(K,I-1))*BHR(K,I-1)  &
                       +(ABS(W(K,I))*BB(K,I)+ABS(W(K-1,I))*BB(K-1,I))  &
                       *DLX(I)+ABS(QSS(K,I)))*0.5
              DLTCAL = 1.0/((QTOT/BH(K,I)+CELRTY)/DLX(I)+TAU1+TAU2)
              IF (DLTCAL.LT.CURMAX) THEN
                KLOC   = K
                ILOC   = I
                CURMAX = INT(DLTCAL)
                IF (DLTF(DLTDP)*CURMAX.LT.MINDLT) THEN
                  KMIN = K
                  IMIN = I
                END IF
              END IF
            END DO
          END DO
        END DO
              
!****** Limiting location

        IF (LIMITING_DLT) THEN
          NDLT(KLOC,ILOC) = NDLT(KLOC,ILOC)+1
          DO I=IU,ID
            DO K=KT,KB(I)
              IF (NDLT(KLOC,ILOC).GT.LIMDLT) THEN
                KLIM   = KLOC
                ILIM   = ILOC
                LIMDLT = NDLT(KLOC,ILOC)
              END IF
            END DO
          END DO
        END IF
        
!****** Restore timestep dependent variables and restart calculations

        IF (CURMAX.LT.DLT.AND.DLT.GT.DLTMIN) THEN
          DLT = INT(DLTF(DLTDP)*CURMAX)
          IF (DLT.LE.DLTMIN) THEN
            WARNING_OPEN = .TRUE.
            WRITE (WRN,6000) JDAY,DLT
            DLT = DLTMIN
          END IF
          CURMAX = DLTMAX(DLTDP)/DLTF(DLTDP)
          IF (DLT.LT.MINDLT) THEN
            MINDLT = DLT
            JDMIN  = JDAY+DLT/86400.0
          END IF
          DO JB=1,NBP
            DO I=CUS(JB)-1,DS(JB)+1
              Z(I)   = SZ(I)
              KTI(I) = SKTI(I)
              DO K=KT,MAX(KB(IU),KB(I))
                U(K,I)   = SU(K,I)
                W(K,I)   = SW(K,I)
                AZ(K,I)  = SAZ(K,I)
                QSS(K,I) = 0.0
              END DO
            END DO
          END DO
          NV = NV+1
          GO TO 10010
        END IF

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*               Task 2.3: Temporal Balance Terms and Temperatures          **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

        IF (.NOT.NO_HEAT) THEN
          CALL HEAT_EXCHANGE (JDAY)
          RS = SRO*RHOWCP
          IF (TERM_BY_TERM) CALL RADIATION (RS,RSN,RAN)
        END IF
        DO JB=1,NBP
          IU = CUS(JB)
          ID = DS(JB)

!******** Heat exchange

          IF (.NOT.NO_HEAT) THEN
            DO I=IU,ID

!****!******* Surface

              IF (.NOT.ICE(I)) THEN

!****!****!**** Surface exchange

                IF (TERM_BY_TERM) THEN
                  CALL SURFACE_TERMS (T2(KT,I),RB,RE,RC)
                  RN(I) = RSN+RAN-RB-RE-RC
                  TFLUX = RN(I)/RHOWCP*B(KTI(I),I)*DLX(I)
                ELSE
                  TFLUX = CSHE*(ET-T2(KT,I))*B(KTI(I),I)*DLX(I)
                END IF
                TSS(KT,I) = TSS(KT,I)+TFLUX
                TSSS(JB)  = TSSS(JB)+TFLUX*DLT

!****!****!**** Solar radiation

                IF (CONSTITUENTS) THEN
                  GAMMA = EXH2O+EXSS*SS(KT,I)+EXOM*(ALGAE(KT,I)  &
                          +LPOM(KT,I))
                ELSE
                  GAMMA = EXH2O
                END IF
                TFLUX     = (1.0-BETA)*SRON*EXP(-GAMMA*DEPTHB(KT))  &
                            *B(KTI(I),I)*DLX(I)
                TSS(KT,I) = TSS(KT,I)-TFLUX
                TSSS(JB)  = TSSS(JB)-TFLUX*DLT
                DO K=KT+1,KB(I)
                  IF (CONSTITUENTS) THEN
                    GAMMA = EXH2O+EXSS*SS(K,I)+EXOM*(ALGAE(K,I)  &
                            +LPOM(K,I))
                  ELSE
                    GAMMA = EXH2O
                  END IF
                  TFLUX    = (1.0-BETA)*SRON*(EXP(-GAMMA*DEPTHB(K))  &
                             -EXP(-GAMMA*DEPTHB(K+1)))*B(K,I)*DLX(I)
                  TSS(K,I) = TSS(K,I)+TFLUX
                  TSSS(JB) = TSSS(JB)+TFLUX*DLT
                END DO
              END IF

!****!******* Bottom

              TFLUX     = CBHE*(TSED-T2(KT,I))*(B(KTI(I),I)  &
                          -B(KT+1,I))*DLX(I)
              TSS(KT,I) = TSS(KT,I)+TFLUX
              TSSB(JB)  = TSSB(JB)+TFLUX*DLT
              DO K=KT+1,KB(I)-1
                TFLUX    = CBHE*(TSED-T2(K,I))*(B(K,I)  &
                           -B(K+1,I))*DLX(I)
                TSS(K,I) = TSS(K,I)+TFLUX
                TSSB(JB) = TSSB(JB)+TFLUX*DLT
              END DO
              TFLUX        = CBHE*(TSED-T2(KB(I),I))  &
                             *(B(KB(I),I))*DLX(I)
              TSS(KB(I),I) = TSS(KB(I),I)+TFLUX
              TSSB(JB)     = TSSB(JB)+TFLUX*DLT
            END DO

!****!***** Ice cover

            IF (ICE_CALC) THEN
              HIA = 0.2367*CSHE/5.65E-8
              DO I=IU,ID
                ALLOW_ICE(I) = .TRUE.
                DO K=KT,KB(I)
                  IF (T2(K,I).GT.ICET2) ALLOW_ICE(I) = .FALSE.
                END DO
              END DO
              ICE_IN(JB) = .TRUE.
              DO I=IU,ID
                IF (ICETH(I).LT.ICEMIN) ICE_IN(JB) = .FALSE.
              END DO
              DO I=IU,ID
                IF (DETAILED_ICE) THEN
                  IF (T2(KT,I).LT.0.0) THEN
                    IF (.NOT.ICE(I)) THEN
                      ICETH2 = -T2(KT,I)*RHO(KT,I)*CP*HKT2(I)/RHORL1
                      IF (ICETH2.LT.ICE_TOL) THEN
                        ICETH2 = 0.0
                      ELSE
                        HEAT       = T2(KT,I)*RHO(KT,I)*CP*HKT2(I)  &
                                     *B(KTI(I),I)/(4.186E6*DLT)*DLX(I)
                        TSS(KT,I)  = TSS(KT,I)-HEAT
                        TSSICE(JB) = TSSICE(JB)-HEAT*DLT
                      END IF
                    END IF
                  END IF

!****!****!****** Ice balance

                  IF (ICE(I)) THEN
                    TICE = TAIR
                    DEL  = 2.0
                    J    = 1
                    DO WHILE (DEL.GT.1.0.AND.J.LT.500)
                      CALL SURFACE_TERMS (TICE,RB,RE,RC)
                      RN(I) = RS*(1.0-ALBEDO)*BETAI+RAN-RB-RE-RC

!****!****!****!***** Heat balance

                      DEL = RN(I)+RK1*(RIMT-TICE)/ICETH(I)
                      IF (ABS(DEL).GT.1.0) TICE = TICE+DEL/500.0
                      J = J+1
                    END DO

!****!****!******** Solar radiation attenuation

                    HEAT       = DLX(I)*SRON*(1.0-ALBEDO)*(1.0-BETAI)  &
                                 *EXP(-GAMMAI*ICETH(I))*B(KTI(I),I)
                    TSS(KT,I)  = TSS(KT,I)+HEAT
                    TSSICE(JB) = TSSICE(JB)+HEAT*DLT
                    IF (TICE.GT.0.0) THEN
                      HICE   = RHOICP*0.5*TICE*0.5*ICETH(I)  &
                               *B(KTI(I),I)/4.186E6/DLT
                      ICETHU = -DLT*HICE/B(KTI(I),I)*4.186E6/RHORL1
                      TICE   = 0.0
                    END IF

!****!****!******** Ice growth
                    
                    IF (TICE.LT.0.0) ICETH1 = DLT*(RK1*(RIMT-TICE)  &
                                              /ICETH(I))/RHORL1
!****!****!******** Ice melt

                    IF (T2(KT,I).GT.0.0) THEN
                      ICETH2     = -DLT*HWI*(T2(KT,I)-RIMT)/RHORL1
                      HEAT       = 2.392E-7*HWI*(RIMT-T2(KT,I))  &
                                   *B(KTI(I),I)*DLX(I)
                      TSS(KT,I)  = TSS(KT,I)+HEAT
                      TSSICE(JB) = TSSICE(JB)+HEAT*DLT
                    END IF
                  END IF

!****!****!****** Ice thickness

                  ICETH(I) = ICETH(I)+ICETHU+ICETH1+ICETH2
                  IF (ICETH(I).LT.ICE_TOL) ICETH(I) = 0.0
                  IF (WINTER.AND.(.NOT.ICE_IN(JB))) THEN
                    IF (.NOT.ALLOW_ICE(I)) ICETH(I) = 0.0
                  END IF
                  ICE(I)   = ICETH(I).GT.0.0
                  ICESW(I) = 1.0
                  IF (ICE(I)) ICESW(I) = 0.0
                  ICETHU = 0.0
                  ICETH1 = 0.0
                  ICETH2 = 0.0
                  IF (ICETH(I).LT.ICE_TOL.AND.ICETH(I).GT.0.0) THEN
                    ICETH(I) = ICE_TOL
                  END IF
                ELSE
                  HIA      = 0.2367*CSHE/5.65E-8
                  ICETH(I) = ICETH(I)+DLT*((RIMT-ET)/(ICETH(I)/RK1+1.0  &
                             /HIA)-(T2(KT,I)-RIMT))/RHORL1
                  ICETH(I) = MAX(ICETH(I),0.0)
                  ICE(I)   = ICETH(I).GT.0.0
                  ICESW(I) = 1.0
                  IF (ICE(I)) THEN
                    HEAT       = 2.392E-7*(RIMT-T2(KT,I))*B(KTI(I),I)  &
                                 *DLX(I)
                    TSS(KT,I)  = TSS(KT,I)+HEAT
                    TSSICE(JB) = TSSICE(JB)+HEAT*DLT
                  END IF
                END IF
              END DO
            END IF
          END IF

!******** Heat sources/sinks & total inflow/outflow

          IF (EVAPORATION) THEN
            DO I=IU,ID
              TSS(KT,I) = TSS(KT,I)-EV(I)*T2(KT,I)
              TSSEV(JB) = TSSEV(JB)-EV(I)*T2(KT,I)*DLT
              VOLEV(JB) = VOLEV(JB)-EV(I)*DLT
            END DO
          END IF
          IF (PRECIPITATION) THEN
            DO I=IU,ID
              TSS(KT,I) = TSS(KT,I)+TPR(JB)*QPR(I)
              TSSPR(JB) = TSSPR(JB)+TPR(JB)*QPR(I)*DLT
              VOLPR(JB) = VOLPR(JB)+QPR(I)*DLT
            END DO
          END IF
          IF (TRIBUTARIES) THEN
            DO JT=1,NTR
              IF (JB.EQ.JBTR(JT)) THEN
                I = ITR(JT)
                IF (I.LT.CUS(JB)) I = CUS(JB)
                DO K=KTTR(JT),KBTR(JT)
                  IF (QTR(JT).LT.0.0) THEN
                    TSS(K,I)  = TSS(K,I)+T2(K,I)*QTR(JT)*QTRF(K,JT)
                    TSSTR(JB) = TSSTR(JB)+T2(K,I)*QTR(JT)*QTRF(K,JT)*DLT
                  ELSE
                    TSS(K,I)  = TSS(K,I)+TTR(JT)*QTR(JT)*QTRF(K,JT)
                    TSSTR(JB) = TSSTR(JB)+TTR(JT)*QTR(JT)*QTRF(K,JT)*DLT
                  END IF
                END DO
                VOLTR(JB) = VOLTR(JB)+QTR(JT)*DLT
              END IF
            END DO
          END IF
          IF (WITHDRAWALS) THEN
            DO JW=1,NWD
              IF (JB.EQ.JBWD(JW)) THEN
                I         = MAX(CUS(JBWD(JW)),IWD(JW))
                K         = MAX(KT,KWD(JW))
                TSS(K,I)  = TSS(K,I)-T2(K,I)*QWD(JW)
                TSSWD(JB) = TSSWD(JB)-T2(K,I)*QWD(JW)*DLT
                VOLWD(JB) = VOLWD(JB)-QWD(JW)*DLT
              END IF
            END DO
          END IF
          IF (UP_FLOW(JB)) THEN
            VOLIN(JB) = VOLIN(JB)+QIN(JB)*DLT
            DO K=KT,KB(IU)
              TSS(K,IU)  = TSS(K,IU)+QINF(K,JB)*QIN(JB)*TIN(JB)
              TSSIN(JB)  = TSSIN(JB)+QINF(K,JB)*QIN(JB)*TIN(JB)*DLT
            END DO
          END IF
          IF (DN_FLOW(JB)) THEN
            DO K=KT,KB(ID)
              TSS(K,ID)  = TSS(K,ID)-QOUT(K,JB)*T2(K,ID)
              TSSOUT(JB) = TSSOUT(JB)-QOUT(K,JB)*T2(K,ID)*DLT
              VOLOUT(JB) = VOLOUT(JB)-QOUT(K,JB)*DLT
            END DO
          END IF
          IF (UH_INTERNAL(JB)) THEN
            DO K=KT,KB(IU-1)
              TSS(K,UHS(JB))  = TSS(K,UHS(JB))-TSSUH2(K,JB)/DLT
              TSSUH(JBUH(JB)) = TSSUH(JBUH(JB))-TSSUH2(K,JB)
              VOLUH(JBUH(JB)) = VOLUH(JBUH(JB))-QUH2(K,JB)
            END DO
          END IF

!******** Horizontal advection and diffusion multipliers

          DO I=IU,ID-1
            DO K=KT,KB(I)
              COUR = U(K,I)*DLT/DLXR(I)
              IF (U(K,I).GE.0.0) THEN
                T1L = T2(K,I-1)
                T2L = T2(K,I)
                T3L = T2(K,I+1)
                IF (U(K,I-1).LE.0.0) T1L = T2(K,I)
                IF (UPWIND) THEN
                  DX1(K,I)  = 0.0
                  DX2(K,I)  = -DX(K,I)/SF1L(K,I)
                  DX3(K,I)  = DX(K,I)/SF1L(K,I)
                  AD1L(K,I) = 0.0
                  AD2L(K,I) = 1.0
                  AD3L(K,I) = 0.0
                ELSE
                  DX1(K,I)  = DX(K,I)*SF11L(K,I,1)
                  DX2(K,I)  = DX(K,I)*SF12L(K,I,1)
                  DX3(K,I)  = DX(K,I)*SF13L(K,I,1)  
                  ALFA      = 2.0*(DX(K,I)*DLT/(SF1L(K,I)*SF1L(K,I))  &
                              -(1.0-COUR*COUR)/6.0)*SF3L(K,I,1)
                  AD1L(K,I) = (ALFA-COUR*SF8L(K,I,1)*0.5)/SF5L(K,I,1)
                  AD2L(K,I) = SF4L(K,I,1)+(ALFA-COUR*SF9L(K,I,1)*0.5)  &
                              /SF6L(K,I,1)
                  AD3L(K,I) = SF2L(K,I,1)+(ALFA-COUR*SF10L(K,I,1)*0.5)  &
                              /SF7L(K,I,1)
                END IF
              ELSE
                T1L = T2(K,I)
                T2L = T2(K,I+1)
                T3L = T2(K,I+2)
                IF (U(K,I+2).GE.0.0) T3L = T2(K,I+1)
                IF (UPWIND) THEN
                  DX1(K,I)  = -DX(K,I)/SF1L(K,I)
                  DX2(K,I)  = DX(K,I)/SF1L(K,I)
                  DX3(K,I)  = 0.0
                  AD1L(K,I) = 0.0
                  AD2L(K,I) = 1.0
                  AD3L(K,I) = 0.0
                ELSE
                  DX1(K,I)  = DX(K,I)*SF11L(K,I,2)
                  DX2(K,I)  = DX(K,I)*SF12L(K,I,2)
                  DX3(K,I)  = DX(K,I)*SF13L(K,I,2)
                  ALFA      = 2.0*(DX(K,I)*DLT/(SF1L(K,I)*SF1L(K,I))  &
                              -(1.0-COUR*COUR)/6.0)*SF3L(K,I,2)
                  AD1L(K,I) = SF2L(K,I,2)+(ALFA-COUR*SF8L(K,I,2)*0.5)  &
                              /SF5L(K,I,2)
                  AD2L(K,I) = SF4L(K,I,2)+(ALFA-COUR*SF9L(K,I,2)*0.5)  &
                              /SF6L(K,I,2)
                  AD3L(K,I) = (ALFA-COUR*SF10L(K,I,2)*0.5)/SF7L(K,I,2)
                END IF
              END IF
              TADL(K,I) =  (DX1(K,I)-U(K,I)*AD1L(K,I))*T1L  &
                          +(DX2(K,I)-U(K,I)*AD2L(K,I))*T2L  &
                          +(DX3(K,I)-U(K,I)*AD3L(K,I))*T3L
            END DO
          END DO

!******** Vertical advection multipliers

          DO I=IU,ID
            DO K=KT,KB(I)-1
              IF (W(K,I).GE.0.0) THEN
                T1V = T2(K-1,I)
                T2V = T2(K,I)
                T3V = T2(K+1,I)
                IF (W(K-1,I).LE.0.0) T1V = T2(K,I)
                IF (K.LE.KT+1) THEN
                  T1V  = T2(KT,I)
                  HTOP = HKT1(I)
                  HBOT = H(K+1)
                  HMID = H(K)
                  IF (K.EQ.KT) HMID = HKT1(I)
                  HMIN       = MIN(HBOT,HMID)
                  SF1V(K)    = (HBOT+HMID)*0.5
                  SF2V(K,1)  = HMID**2
                  SF3V(K,1)  = HMID/(HMID+HBOT)
                  SF4V(K,1)  = HBOT/(HMID+HBOT)
                  SF5V(K,1)  = 0.25*(HTOP+2.0*HMID+HBOT)*(HTOP+HMID)
                  SF6V(K,1)  = -0.25*(HMID+HBOT)*(HTOP+HMID)
                  SF7V(K,1)  = 0.25*(HMID+HBOT)*(HTOP+2.0*HMID+HBOT)
                  SF8V(K,1)  = 0.5*(HMID-HBOT)*HMIN
                  SF9V(K,1)  = 0.5*(HTOP+2.0*HMID-HBOT)*HMIN
                  SF10V(K,1) = 0.5*(HTOP+3.0*HMID)*HMIN
                END IF
                IF (UPWIND) THEN
                  AD1V(K,I) = 0.0
                  AD2V(K,I) = 1.0
                  AD3V(K,I) = 0.0
                ELSE
                  COUR      = W(K,I)*DLT/SF1V(K)
                  ALFA      = 2.0*(DZQ(K,I)*DLT/(SF1V(K)*SF1V(K))  &
                              -(1.0-COUR*COUR)/6.0)*SF2V(K,1)
                  AD1V(K,I) = (ALFA-COUR*SF8V(K,1)*0.5)/SF5V(K,1)
                  AD2V(K,I) = SF4V(K,1)+(ALFA-COUR*SF9V(K,1)*0.5)  &
                              /SF6V(K,1)
                  AD3V(K,I) = SF3V(K,1)+(ALFA-COUR*SF10V(K,1)*0.5)  &
                              /SF7V(K,1)
                END IF
              ELSE
                T1V = T2(K,I)
                T2V = T2(K+1,I)
                T3V = T2(K+2,I)
                IF (W(K+2,I).GE.0.0) T3V = T2(K+1,I)
                IF (K.EQ.KB(I)-1)    T3V = T2(K+1,I)
                IF (K.EQ.KT) THEN
                  HTOP       = HKT1(I)
                  HMID       = H(KT+1)
                  HBOT       = H(KT+2)
                  HMIN       = MIN(HTOP,HMID)
                  SF1V(K)    = (HMID+HTOP)*0.5
                  SF2V(K,2)  = HMID**2
                  SF3V(K,2)  = HMID/(HTOP+HMID)
                  SF4V(K,2)  = HTOP/(HTOP+HMID)
                  SF5V(K,2)  = 0.25*(HTOP+2.0*HMID+HBOT)*(HTOP+HMID)
                  SF6V(K,2)  = -0.25*(HMID+HBOT)*(HTOP+HMID)
                  SF7V(K,2)  = 0.25*(HTOP+2.0*HMID+HBOT)*(HMID+HBOT)
                  SF8V(K,2)  = -0.5*(3.0*HMID+HBOT)*HMIN
                  SF9V(K,2)  = 0.5*(HTOP-2.0*HMID-HBOT)*HMIN
                  SF10V(K,2) = 0.5*(HTOP-HMID)*HMIN
                END IF
                IF (UPWIND) THEN
                  AD1V(K,I) = 0.0
                  AD2V(K,I) = 1.0
                  AD3V(K,I) = 0.0
                ELSE
                  COUR      = W(K,I)*DLT/SF1V(K)
                  ALFA      = 2.0*(DZQ(K,I)*DLT/(SF1V(K)*SF1V(K))  &
                              -(1.0-COUR*COUR)/6.0)*SF2V(K,2)
                  AD1V(K,I) = SF3V(K,2)+(ALFA-COUR*SF8V(K,2)*0.5)  &
                              /SF5V(K,2)
                  AD2V(K,I) = SF4V(K,2)+(ALFA-COUR*SF9V(K,2)*0.5)  &
                              /SF6V(K,2)
                  AD3V(K,I) = (ALFA-COUR*SF10V(K,2)*0.5)/SF7V(K,2)
                END IF
              END IF
              TADV(K,I) = -W(K,I)*(AD1V(K,I)*T1V+AD2V(K,I)*T2V  &
                          +AD3V(K,I)*T3V)
            END DO
          END DO
        END DO

!****** Heat transport

        DO JB=1,NBP
          IU = CUS(JB)
          ID = DS(JB)
          DO I=IU,ID
            THETA = THETAI
            IF (ONE_LAYER(I)) THETA = 0.0
            T1(KT,I) = (T2(KT,I)*BHKT2(I)/DLT+(TADL(KT,I)*BHRKT1(I)  &
                       -TADL(KT,I-1)*BHRKT1(I-1))/DLX(I)+(1.0-THETA)  &
                       *TADV(KT,I)*BB(KT,I)+TSS(KT,I)/DLX(I))*DLT  &
                       /BHKT1(I)
            DO K=KT+1,KB(I)
              T1(K,I) = (T2(K,I)*BH(K,I)/DLT+(TADL(K,I)*BHR(K,I)  &
                        -TADL(K,I-1)*BHR(K,I-1))/DLX(I)+(1.0-THETA)  &
                        *(TADV(K,I)*BB(K,I)-TADV(K-1,I)*BB(K-1,I))  &
                        +TSS(K,I)/DLX(I))*DLT/BH(K,I)
            END DO

!****!***** Vertical advection and implicit diffusion

            IF (.NOT.ONE_LAYER(I)) THEN
              K       = KT
              AT(K,I) = 0.0
              CT(K,I) = DLT/BHKT1(I)*(BB(K,I)*(THETA*0.5*W(K,I)-DZ(K,I)  &
                        /AVHKT(I)))
              VT(K)   = 1.0+DLT/BHKT1(I)*(BB(K,I)*(DZ(K,I)/AVHKT(I)  &
                        +THETA*0.5*W(K,I)))
              DT(K)   = T1(K,I)
              K       = KT+1
              AT(K,I) = -DLT/BH(K,I)*(BB(K-1,I)*(DZ(K-1,I)/AVHKT(I)  &
                        +THETA*0.5*W(K-1,I)))
              CT(K,I) = DLT/BH(K,I)*(BB(K,I)*(THETA*0.5*W(K,I)-DZ(K,I)  &
                        /AVH(K)))
              VT(K)   = 1.0+DLT/BH(K,I)*(BB(K,I)*(DZ(K,I)/AVH(K)+THETA  &
                        *0.5*W(K,I))+BB(K-1,I)*(DZ(K-1,I)/AVHKT(I)  &
                        -THETA*0.5*W(K-1,I)))
              DT(K)   = T1(K,I)
              DO K=KT+2,KB(I)-1
                AT(K,I) = -DLT/BH(K,I)*(BB(K-1,I)*(DZ(K-1,I)/AVH(K-1)  &
                          +THETA*0.5*W(K-1,I)))
                CT(K,I) = DLT/BH(K,I)*(BB(K,I)*(THETA*0.5*W(K,I)-DZ(K,I)  &
                          /AVH(K)))
                VT(K)   = 1.0+DLT/BH(K,I)*(BB(K,I)*(DZ(K,I)/AVH(K)+THETA  &
                          *0.5*W(K,I))+BB(K-1,I)*(DZ(K-1,I)/AVH(K-1)  &
                          -THETA*0.5*W(K-1,I)))
                DT(K)   = T1(K,I)
              END DO
              K = KB(I)
              IF (KB(I)-KT.GT.1) THEN
                AT(K,I) = -DLT/BH(K,I)*(BB(K-1,I)*(DZ(K-1,I)/AVH(K-1)  &
                          +THETA*0.5*W(K-1,I)))
                CT(K,I) = 0.0
                VT(K)   = 1.0+DLT/BH(K,I)*(BB(K-1,I)*(DZ(K-1,I)/AVH(K-1)  &
                          -THETA*0.5*W(K-1,I)))
                DT(K)   = T1(K,I)
              ELSE
                AT(K,I) = -DLT/BH(K,I)*(BB(K-1,I)*(DZ(K-1,I)/AVHKT(I)  &
                          +THETA*0.5*W(K-1,I)))
                CT(K,I) = 0.0
                VT(K)   = 1.0+DLT/BH(K,I)*(BB(K-1,I)*(DZ(K-1,I)/AVHKT(I)  &
                          -THETA*0.5*W(K-1,I)))  
                DT(K)   = T1(K,I)
              END IF

!****!******* Tridiagonal solution

              BTAT(KT,I) = VT(KT)
              DO K=KT+1,KB(I)
                BTAT(K,I) = VT(K)-AT(K,I)/BTAT(K-1,I)*CT(K-1,I)
              END DO
              GMAT(KT) = DT(KT)
              DO K=KT+1,KB(I)
                GMAT(K) = DT(K)-AT(K,I)/BTAT(K-1,I)*GMAT(K-1)
              END DO
              T1(KB(I),I) = GMAT(KB(I))/BTAT(KB(I),I)
              DO K=KB(I)-1,KT,-1
                T1(K,I) = (GMAT(K)-CT(K,I)*T1(K+1,I))/BTAT(K,I)
              END DO
            END IF
          END DO
        END DO

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                            Task 2.4:  Constituents                       **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

        IF (CONSTITUENTS) THEN
          PALT = (1.-((EL(KT)-Z(DS(1)))/1000.0)/44.3)**5.25
          DO JB=1,NBP
            IU = CUS(JB)
            ID = DS(JB)

!****!***** Internal sources/sinks

            IF (UPDATE_KINETICS) THEN
              IF (UPDATE_RATES) THEN
                CALL RATE_MULTIPLIERS
                CALL DECAY_CONSTANTS
              END IF
              DO JAC=1,NAC
                JC = CN(JAC)
                IF (JC.EQ.2)  CALL SUSPENDED_SOLIDS
                IF (JC.EQ.3)  CALL COLIFORM
                IF (JC.EQ.5)  CALL LABILE_DOM
                IF (JC.EQ.6)  CALL REFRACTORY_DOM
                IF (JC.EQ.7)  CALL PHYTOPLANKTON
                IF (JC.EQ.8)  CALL LABILE_POM
                IF (JC.EQ.9)  CALL PHOSPHORUS
                IF (JC.EQ.10) CALL AMMONIUM
                IF (JC.EQ.11) CALL NITRATE
                IF (JC.EQ.12) CALL DISSOLVED_OXYGEN
                IF (JC.EQ.13) CALL SEDIMENT
                IF (JC.EQ.14) CALL INORGANIC_CARBON
                IF (JC.EQ.16) CALL PH_CO2
                IF (JC.EQ.20) CALL IRON
                IF (JC.EQ.21) CALL BIOCHEMICAL_O2_DEMAND
              END DO
            END IF
            DO JAC=1,NAC
              JC = CN(JAC)

!****!******* External sources/sinks

              IF (TRIBUTARIES) THEN
                DO JT=1,NTR
                  IF (JB.EQ.JBTR(JT)) THEN
                    I = ITR(JT)
                    IF (I.LT.CUS(JB)) I = CUS(JB)
                    DO K=KTTR(JT),KBTR(JT)
                      IF (QTR(JT).LT.0.0) THEN
                        CSSB(K,I,JC) = CSSB(K,I,JC)+C1(K,I,JC)*QTR(JT)  &
                                       *QTRF(K,JT)
                      ELSE
                        CSSB(K,I,JC) = CSSB(K,I,JC)+CTR(JC,JT)*QTR(JT)  &
                                       *QTRF(K,JT)
                      END IF
                    END DO
                  END IF
                END DO
              END IF
              IF (WITHDRAWALS) THEN
                DO JW=1,NWD
                  IF (JB.EQ.JBWD(JW)) THEN
                    I            = MAX(CUS(JBWD(JW)),IWD(JW))  
                    K            = MAX(KT,KWD(JW))
                    CSSB(K,I,JC) = CSSB(K,I,JC)-C1S(K,I,JC)*QWD(JW)
                  END IF
                END DO
              END IF
              IF (PRECIPITATION) THEN
                DO I=IU,ID
                  CSSB(KT,I,JC) = CSSB(KT,I,JC)+CPR(JC,JB)*QPR(I)
                END DO
              END IF
              IF (UP_FLOW(JB)) THEN
                DO K=KT,KB(IU)
                  CSSB(K,IU,JC) = CSSB(K,IU,JC)+QINF(K,JB)*QIN(JB)  &
                                  *CIN(JC,JB)
                END DO
              END IF
              IF (DN_FLOW(JB)) THEN
                DO K=KT,KB(ID)
                  CSSB(K,ID,JC) = CSSB(K,ID,JC)-QOUT(K,JB)  &
                                  *C1S(K,ID,JC)
                END DO
              END IF
            END DO
          END DO

!******** Transport constituents

          DO JB=1,NBP
            IU = CUS(JB)
            ID = DS(JB)
            DO JAC=1,NAC
              JC = CN(JAC)
              IF (TRANSPORT(JC)) THEN

!****!****!**** Horizontal advection and diffusion terms

                DO I=IU,ID-1
                  DO K=KT,KB(I)
                    IF (U(K,I).GE.0.0) THEN
                      C1L = C1S(K,I-1,JC)
                      C2L = C1S(K,I,JC)
                      C3L = C1S(K,I+1,JC)
                      IF (U(K,I-1).LE.0.0) C1L = C1S(K,I,JC)
                    ELSE
                      C1L = C1S(K,I,JC)
                      C2L = C1S(K,I+1,JC)
                      C3L = C1S(K,I+2,JC)
                      IF (U(K,I+2).GE.0.0) C3L = C1S(K,I+1,JC)
                    END IF
                    CADL(K,I,JC) =  (DX1(K,I)-U(K,I)*AD1L(K,I))*C1L  &
                                   +(DX2(K,I)-U(K,I)*AD2L(K,I))*C2L  &
                                   +(DX3(K,I)-U(K,I)*AD3L(K,I))*C3L
                  END DO
                END DO

!****!****!**** Vertical advection terms

                DO I=IU,ID
                  DO K=KT,KB(I)-1
                    IF (W(K,I).GE.0.0) THEN
                      C1V = C1S(K-1,I,JC)
                      C2V = C1S(K,I,JC)
                      C3V = C1S(K+1,I,JC)
                      IF (W(K-1,I).LE.0.0) C1V = C1S(K,I,JC)
                      IF (K.LE.KT+1) C1V = C1S(KT,I,JC)
                    ELSE
                      C1V = C1S(K,I,JC)
                      C2V = C1S(K+1,I,JC)
                      C3V = C1S(K+2,I,JC)
                      IF (W(K+2,I).GE.0.0) C3V = C1S(K+1,I,JC)
                      IF (K.EQ.KB(I)-1)    C3V = C1S(K+1,I,JC)
                    END IF
                    CADV(K,I,JC) = -W(K,I)*(AD1V(K,I)*C1V+AD2V(K,I)*C2V  &
                                   +AD3V(K,I)*C3V)
                  END DO
                END DO
              END IF
            END DO
          END DO

!******** Constituent transport

          DO JB=1,NBP
            IU = CUS(JB)
            ID = DS(JB)
            DO JAC=1,NAC
              JC = CN(JAC)
              IF (TRANSPORT(JC)) THEN
                DO I=IU,ID
                  C1(KT,I,JC) = ((C1S(KT,I,JC)*BHKT2(I)/DLT  &
                                +(CADL(KT,I,JC)*BHRKT1(I)  &
                                -CADL(KT,I-1,JC)*BHRKT1(I-1))/DLX(I)  &
                                +(1.0-THETA)*CADV(KT,I,JC)*BB(KT,I))  &
                                /BHKT1(I)+CSSB(KT,I,JC)/(DLX(I)  &
                                *BHKT1(I))+CSSK(KT,I,JC))*DLT
                  DO K=KT+1,KB(I)
                    C1(K,I,JC) = ((C1S(K,I,JC)*BH(K,I)/DLT  &
                                 +(CADL(K,I,JC)*BHR(K,I)-CADL(K,I-1,JC)  &
                                 *BHR(K,I-1))/DLX(I)+(1.0-THETA)  &
                                 *(CADV(K,I,JC)*BB(K,I)-CADV(K-1,I,JC)  &
                                 *BB(K-1,I)))/BH(K,I)+CSSB(K,I,JC)  &
                                 /(DLX(I)*BH(K,I))+CSSK(K,I,JC))*DLT
                  END DO

!****!****!****** Time-weighted vertical advection and implicit diffusion

                  IF (.NOT.ONE_LAYER(I)) THEN
                    DT(KB(I)) = C1(KB(I),I,JC)
                    DO K=KT,KB(I)-1
                      DT(K) = C1(K,I,JC)
                    END DO
                    GMAT(KT) = DT(KT)
                    DO K=KT+1,KB(I)
                      GMAT(K) = DT(K)-AT(K,I)/BTAT(K-1,I)*GMAT(K-1)
                    END DO
                    C1(KB(I),I,JC) = GMAT(KB(I))/BTAT(KB(I),I)
                    DO K=KB(I)-1,KT,-1
                      C1(K,I,JC) = (GMAT(K)-CT(K,I)*C1(K+1,I,JC))  &
                                   /BTAT(K,I)
                    END DO
                  END IF
                END DO
              END IF
            END DO
          END DO
        END IF

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*              Task 2.5: Layer - Segment Additions and Subtractions        **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

!****** Water surface minimum thickness

        ZMIN = -1000.0
        DO JB=1,NBP
          DO I=CUS(JB),DS(JB)
            ZMIN = MAX(ZMIN,Z(I))
          END DO
        END DO
        ADD_LAYER = ZMIN.LT.-0.80*H(KT-1).AND.KT.NE.2
        SUB_LAYER = ZMIN.GT.0.60*H(KT)

!****** Add layers

        DO WHILE (ADD_LAYER)
          IF (SNAPSHOT) WRITE (SNP,4000) KT-1,JDAY

!******** Variable initialization

          KT = KT-1
          DO JB=1,NBP
            IU = CUS(JB)
            ID = DS(JB)
            DO I=IU-1,ID+1
              Z(I)     = H(KT)+Z(I)
              U(KT,I)  = U(KT+1,I)
              W(KT,I)  = W(KT+1,I)*(BB(KT+1,I)/BB(KT,I))
              HKT1(I)  = H(KT)-Z(I)
              BHKT1(I) = BHKT1(I)-BH(KT+1,I)
              BKT(I)   = BHKT1(I)/HKT1(I)
              T1(KT,I) = T1(KT+1,I)
              DO JC=1,NAC
                C1(KT,I,CN(JC))   = C1(KT+1,I,CN(JC))
                CSSK(KT,I,CN(JC)) = CSSK(KT+1,I,CN(JC))
              END DO
              IF (CONSTITUENTS) THEN
                RHO(KT,I) = DENSITY (T2(KT,I),SS(KT,I),TDS(KT,I))
              ELSE
                RHO(KT,I) = DENSITY (T2(KT,I),0.0,0.0)
              END IF
            END DO
            DO I=IU-1,ID
              BHRKT1(I) = (BHKT1(I)+BHKT1(I+1))*0.5
            END DO
            DO I=IU,ID-1
              DX(KT,I) = DXI
            END DO
            IUT = IU
            IDT = ID-1
            DO I=IUT,IDT
              AZ(KT,I)  = AZMIN
              SAZ(KT,I) = AZMIN
            END DO

!****!***** Upstream active segment

            IUT = US(JB)
            DO I=US(JB),DS(JB)
              IF (KB(I)-KT.LT.NL(JB)-1) IUT = I+1
            END DO

!****!***** Segment addition

            IF (IUT.NE.IU) THEN
              IF (SNAPSHOT) WRITE (SNP,4010) IUT,IU-1,JDAY
              DO I=IUT-1,IU-1
                Z(I)     = Z(IU)
                KTI(I)   = KTI(IU)
                HKT1(I)  = H(KT)-Z(IU)
                BHKT1(I) = B(KTI(I),I)*(EL(KT)-EL(KTI(I)+1)-Z(I))
                BKT(I)   = BHKT1(I)/HKT1(I)
              END DO
              DO I=IUT-1,IU-1
                BHRKT1(I) = (BHKT1(I)+BHKT1(I+1))*0.5
              END DO
              DO I=IUT,IU-1
                ICE(I)   = ICE(IU)
                ICETH(I) = ICETH(IU)
                DO K=KT,KB(I)
                  DX(K,I) = DXI
                  T1(K,I) = T1(K,IU)
                  T2(K,I) = T1(K,IU)
                  U(K,I)  = U(K,IU-1)
                  SU(K,I) = U(K,IU-1)
                  DO JC=1,NAC
                    BHT = BH(K,I)
                    IF (K.EQ.KT) BHT = BHKT1(I)
                    C1(K,I,CN(JC))   = C1(K,IU,CN(JC))
                    C2(K,I,CN(JC))   = C1(K,IU,CN(JC))
                    CMBRT(CN(JC),JB) = CMBRT(CN(JC),JB)+C1(K,IU,CN(JC))  &
                                       *DLX(I)*BHT
                  END DO
                END DO
                DO K=KT,KB(I)-1
                  AZ(K,I)  = AZ(K,IU)
                  SAZ(K,I) = AZ(K,IU)
                END DO
              END DO
              DO K=KB(IUT),KB(IU)
                U(K,IU-1)  = 0.0
                SU(K,IU-1) = 0.0
                TADL(K,IU) = 0.0
                IF (CONSTITUENTS) THEN
                  DO JAC=1,NAC
                    CADL(K,IU,CN(JAC)) = 0.0
                  END DO
                END IF
              END DO
              IF (SHIFT_DEMAND) THEN
                IDIFF = IUT-US(JB)
                DO I=IUT,ID
                  SOD(I) = SODS(I-IDIFF)
                END DO
              END IF
              IU      = IUT
              CUS(JB) = IU
              IF (UH_INTERNAL(JB)) KB(IU-1) = MIN(KB(UHS(JB)),KB(IU))
            END IF

!****!***** Total active cells and single layers

            DO I=IU,ID
              NTAC         = NTAC+1
              ONE_LAYER(I) = KT.EQ.KB(I)
            END DO
            NTACMX = MAX(NTAC,NTACMX)
          END DO

!******** Additional layers

          ZMIN = -1000.0
          DO JB=1,NBP
            DO I=CUS(JB),DS(JB)
              ZMIN = MAX(ZMIN,Z(I))
            END DO
          END DO
          ADD_LAYER = ZMIN.LT.-0.80*H(KT-1).AND.KT.NE.2
        END DO
        
!****** Subtract layers

        DO WHILE (SUB_LAYER)
          IF (SNAPSHOT) WRITE (SNP,4020) KT,JDAY

!******** Variable initialization

          KT = KT+1
          DO JB=1,NBP
            IU = CUS(JB)
            ID = DS(JB)
            DO I=IU-1,ID+1
              A1(KT-1,I) = 0.0
            END DO
            DO I=IU-1,ID+1
              Z(I)     = Z(I)-H(KT-1)
              HKT1(I)  = H(KT)-Z(I)
              BHKT1(I) = BHKT1(I)+BH(KT,I)
              BKT(I)   = BHKT1(I)/HKT1(I)
              T1(KT,I) = (T1(KT-1,I)*(BHKT1(I)-BH(KT,I))  &
                         +T1(KT,I)*BH(KT,I))/BHKT1(I)
              DO JC=1,NAC
                JAC              = CN(JC)
                C1(KT,I,JAC)     = (C1(KT-1,I,JAC)*(BHKT1(I)-BH(KT,I))  &
                                   +C1(KT,I,JAC)*BH(KT,I))/BHKT1(I)
                CSSB(KT,I,JAC)   = CSSB(KT-1,I,JAC)+CSSB(KT,I,JAC)
                CSSK(KT,I,JAC)   = (CSSK(KT-1,I,JAC)*(BHKT1(I)-BH(KT,I))  &
                                   +CSSK(KT,I,JAC)*BH(KT,I))/BHKT1(I)
                CSSB(KT-1,I,JAC) = 0.0
                CSSK(KT-1,I,JAC) = 0.0
              END DO
            END DO
            DO I=IU-1,ID
              BHRKT1(I) = (BHKT1(I)+BHKT1(I+1))*0.5
            END DO
!****!***** Upstream active segment

            IUT = US(JB)
            DO I=US(JB),DS(JB)
              IF (KB(I)-KT.LT.NL(JB)-1) IUT = I+1
              ONE_LAYER(I) = KT.EQ.KB(I)
            END DO
            IF (IUT.GT.DS(JB)-1) THEN
              OPEN  (ERR,FILE='err.opt',STATUS='UNKNOWN')
              WRITE (ERR,7000) JB,JDAY,KT
              STOP
            END IF

!****!***** Segment subtraction

            IF (IUT.NE.IU) THEN
              IF (SNAPSHOT) WRITE (SNP,4030) IU,IUT-1,JDAY
              DO I=IU,IUT-1
                DO JC=1,NAC
                  VOL              = BHKT1(I)*DLX(I)
                  CMBRT(CN(JC),JB) = CMBRT(CN(JC),JB)-C1(KT,I,CN(JC))  &
                                     *VOL+(CSSB(KT,I,CN(JC))  &
                                     +CSSK(KT,I,CN(JC))*VOL)*DLT
                  DO K=KT+1,KB(I)
                    VOL              = BH(K,I)*DLX(I)
                    CMBRT(CN(JC),JB) = CMBRT(CN(JC),JB)-C1(K,I,CN(JC))  &
                                       *VOL+(CSSB(K,I,CN(JC))  &
                                       +CSSK(K,I,CN(JC))*VOL)*DLT
                  END DO
                END DO
              END DO
              DO I=IU-1,IUT-1
                F(I)     = 0.0
                Z(I)     = 0.0
                ICETH(I) = 0.0
                BHRHO(I) = 0.0
                ICE(I)   = .FALSE.
                DO K=KT,KB(I)
                  AZ(K,I)   = 0.0
                  SAZ(K,I)  = 0.0
                  DX(K,I)   = 0.0
                  U(K,I)    = 0.0
                  T1(K,I)   = 0.0
                  TADL(K,I) = 0.0
                  QSS(K,I)  = 0.0
                  DO JC=1,NAC
                    C1(K,I,CN(JC))   = 0.0
                    C2(K,I,CN(JC))   = 0.0
                    C1S(K,I,CN(JC))  = 0.0
                    CADL(K,I,CN(JC)) = 0.0
                    CSSB(K,I,CN(JC)) = 0.0
                    CSSK(K,I,CN(JC)) = 0.0
                  END DO
                END DO
              END DO
              IF (SHIFT_DEMAND) THEN
                IDIFF = IUT-US(JB)
                DO I=IUT,ID
                  SOD(I) = SODS(I-IDIFF)
                END DO
              END IF
              IU           = IUT
              CUS(JB)      = IU
              Z(IU-1)      = Z(IU)
              KTI(IU-1)    = KTI(IU)
              HKT1(IU-1)   = HKT1(IU)
              BHKT1(IU-1)  = B(KTI(IU),IU-1)*EL(KT)-(EL(KTI(I)+1)-Z(IU))
              BKT(IU-1)    = BHKT1(IU-1)/HKT1(IU-1)
              BHRKT1(IU-1) = (BHKT1(IU-1)+BHKT1(IU))*0.5
              IF (UH_INTERNAL(JB)) KB(IU-1) = MIN(KB(UHS(JB)),KB(IU))
            END IF

!****!***** Total active cells

            DO I=IU,ID
              NTAC = NTAC-1
            END DO
            NTACMN = MIN(NTAC,NTACMN)
          END DO

!******** Additional layer subtractions

          ZMIN = -1000.0
          DO JB=1,NBP
            DO I=CUS(JB),DS(JB)
              ZMIN = MAX(ZMIN,Z(I))
            END DO
          END DO
          SUB_LAYER = ZMIN.GT.0.60*H(KT)
        END DO

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                              Task 2.6: Balances                          **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

        IF (VOLUME_BALANCE) THEN
          DO JB=1,NBP
            VOLSBR(JB) = VOLSBR(JB)+DLVOL(JB)
            VOLTBR(JB) = VOLEV(JB)+VOLPR(JB)+VOLTR(JB)+VOLDT(JB)  &
                         +VOLWD(JB)+VOLUH(JB)+VOLDH(JB)+VOLIN(JB)  &
                         +VOLOUT(JB)
            IF (ABS(VOLSBR(JB)-VOLTBR(JB)).GT.VTOL) THEN
              IF (VOLUME_WARNING) THEN
                WARNING_OPEN = .TRUE.
                WRITE (WRN,6010) JDAY,VOLSBR(JB),VOLTBR(JB),VOLSBR(JB)  &
                                 -VOLTBR(JB)
                VOLUME_WARNING = .FALSE.
              END IF
            END IF
          END DO
        END IF
        IF (ENERGY_BALANCE) THEN
          DO JB=1,NBP
            ETBR(JB) = EIBR(JB)+TSSEV(JB)+TSSPR(JB)+TSSTR(JB)  &
                       +TSSDT(JB)+TSSWD(JB)+TSSUH(JB)+TSSDH(JB)  &
                       +TSSIN(JB)+TSSOUT(JB)+TSSS(JB)+TSSB(JB)  &
                       +TSSICE(JB)
            ESBR(JB) = 0.0
            DO I=CUS(JB),DS(JB)
              ESBR(JB) = ESBR(JB)+T1(KT,I)*DLX(I)*BHKT1(I)
              DO K=KT+1,KB(I)
                ESBR(JB) = ESBR(JB)+T1(K,I)*DLX(I)*BH(K,I)
              END DO
            END DO
          END DO
        END IF
        IF (MASS_BALANCE) THEN
          DO JB=1,NBP
            DO JC=1,NAC
              JAC = CN(JC)
              IF (TRANSPORT(JAC)) THEN
                CMBRS(JAC,JB) = 0.0
                DO I=CUS(JB),DS(JB)
                  CMBRS(JAC,JB) = CMBRS(JAC,JB)+C1(KT,I,JAC)*DLX(I)  &
                                  *BHKT1(I)
                  CMBRT(JAC,JB) = CMBRT(JAC,JB)+(CSSB(KT,I,JAC)  &
                                  +CSSK(KT,I,JAC)*BHKT1(I)*DLX(I))*DLT
                  DO K=KT+1,KB(I)
                    CMBRS(JAC,JB) = CMBRS(JAC,JB)+C1(K,I,JAC)*DLX(I)  &
                                    *BH(K,I)
                    CMBRT(JAC,JB) = CMBRT(JAC,JB)+(CSSB(K,I,JAC)  &
                                    +CSSK(K,I,JAC)*BH(K,I)*DLX(I))*DLT
                  END DO
                END DO
              END IF
            END DO
          END DO
        END IF
        IF (WATER_BALANCE) THEN
          IF (NIT.GT.0.AND.JDAY+DLT/86400.0.GE.NXELV1) THEN
            SAREA = 0.0
            DO JB=1,NBP
              DO I=CUS(JB),DS(JB)
                KTO = 2
                DO WHILE (EL(KTO).GT.ELWSO1)
                  KTO = KTO+1
                END DO
                KTO   = MAX(2,KTO-1)
                SAREA = SAREA+B(KTO,I)*DLX(I)
              END DO
            END DO
            DLVOLO = (ELWSO1-ELWSO2)*SAREA
            DLVOLC = 0.0
            DO JB=1,NBP
              DO I=CUS(JB),DS(JB)
                BKTC = B(KT,I)
                IF (Z(I).LT.0.0) BKTC = B(KT-1,I)
                DLVOLC = DLVOLC+(ELWS(I)-ELWS2(I))*BKTC*DLX(I)
              END DO
            END DO
            DO I=1,IMP
              ELWS2(I) = ELWS(I)
            END DO
            QADJ = (DLVOLO-DLVOLC)/(JDAY-NXELV2)/86400.0
            WRITE (QWS,9010) NXELV2,QADJ
            NXELV2 = NXELV1
          END IF
        END IF

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                 Task 2.7: Variable Updates for Next Timestep             **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

        DO JB=1,NBP
          IU = CUS(JB)
          ID = DS(JB)
          DO I=IU-1,ID+1
            SZ(I)     = Z(I)
            SKTI(I)   = KTI(I)
            HKT2(I)   = HKT1(I)
            AVHKT(I)  = (HKT2(I)+H(KT+1))*0.5
            BHKT2(I)  = BHKT1(I)
            BHRKT2(I) = BHRKT1(I)
            IF (HKT2(I).LT.0.0) THEN
              OPEN  (ERR,FILE='err.opt',STATUS='UNKNOWN')
              WRITE (ERR,7010) JDAY,I,HKT2(I)
              STOP
            END IF
            DO K=KT,MAX(KB(IU),KB(I))
              T2(K,I)  = T1(K,I)
              SU(K,I)  = U(K,I)
              SW(K,I)  = W(K,I)
              SAZ(K,I) = AZ(K,I)
              QSS(K,I) = 0.0
              TSS(K,I) = 0.0
            END DO
          END DO
          IF (UH_INTERNAL(JB)) THEN
            DO K=KT,KB(IU-1)
              QUH2(K,JB)   = QUH1(K,JB)*DLT
              TSSUH2(K,JB) = TSSUH1(K,JB)*DLT
            END DO
          END IF
          IF (DH_INTERNAL(JB)) THEN
            DO K=KT,KB(ID+1)
              QDH2(K,JB)   = QDH1(K,JB)*DLT
              TSSDH2(K,JB) = TSSDH1(K,JB)*DLT
            END DO
          END IF
          DO JC=1,NAC
            DO I=IU-1,ID+1
              DO K=KT,KB(I)
                C1S(K,I,CN(JC))  = C1(K,I,CN(JC))
                C2(K,I,CN(JC))   = MAX(C1(K,I,CN(JC)),0.0)
                CSSB(K,I,CN(JC)) = 0.0
              END DO
            END DO
            IF (UH_INTERNAL(JB)) THEN
              DO K=KT,KB(IU-1)
                CSSUH2(K,CN(JC),JB) = CSSUH1(K,CN(JC),JB)*DLT
              END DO
            END IF
            IF (DH_INTERNAL(JB)) THEN
              DO K=KT,KB(ID+1)
                CSSDH2(K,CN(JC),JB) = CSSDH1(K,CN(JC),JB)*DLT
              END DO
            END IF
          END DO
        END DO
        ELKT = EL(KT)-Z(DS(1))
        DO I=1,IMP
          ELWS(I) = EL(KT)-Z(I)
        END DO

!****** Time variable updates

        NIT     = NIT+1
        ELTM    = ELTM+DLT
        JDAY    = ELTM/86400.0
        ELTMJD  = JDAY-TMSTRT
        END_RUN = JDAY.GE.TMEND 
        DLTS    = DLT
        DLT     = INT(MAX(DLTMIN,DLTF(DLTDP)*CURMAX))
        DLT     = MIN(DLT,DLTS+1.0*DLTS)
        DLT     = MIN(DLT,DLTTVD)
        DLTAV   = (ELTM-TMSTRT*86400.0)/NIT
        IF (JDAY.GE.WSCD(WSCDP+1)) WSCDP = WSCDP+1
        IF (JDAY.GE.DLTD(DLTDP+1)) DLTDP = DLTDP+1
        IF (DLT.GT.DLTMAX(DLTDP))  DLT   = DLTMAX(DLTDP)
        IF (DLTS.LT.MINDLT.AND.DLTS.NE.DLTTVDS) THEN
          MINDLT = DLTS
          JDMIN  = JDAY
        END IF
        CURMAX = DLTMAX(DLTDP)/DLTF(DLTDP)
        IF (INT(JDAY).EQ.JDAYNX) THEN
          JDAYG  = JDAYG+1
          JDAYNX = JDAYNX+1
        END IF
        IF (JDAYG.GT.300) WINTER = .TRUE.
        IF (JDAYG.LT.40)  WINTER = .FALSE.
        WRITE (GDCH,'(I3)') GDAY
        CALL GREGORIAN_DATE (YEAR)
        IF (CONSTITUENTS) THEN
          UPDATE_KINETICS = .FALSE.
          IF (MOD(NIT,FREQUK).EQ.0) UPDATE_KINETICS = .TRUE.
        END IF

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                           Task 2.8: Output Results                       **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

!****** Snapshots

        IF (SNAPSHOT) THEN
          IF (JDAY.GE.NXTMSN.OR.JDAY.GE.SNPD(SNPDP+1)) THEN
            IF (JDAY.GE.SNPD(SNPDP+1)) THEN
              SNPDP  = SNPDP+1
              NXTMSN = SNPD(SNPDP)
            END IF
            NXTMSN = NXTMSN+SNPF(SNPDP)
            WRITE (SNP,3000) (TITLE(I),I=1,7)
            WRITE (SNP,3010) 'Time Parameters',MONTH,GDAY,YEAR,  &
                              INT(JDAY),(JDAY-INT(JDAY))*24.0,  &
                              INT(ELTMJD),(ELTMJD-INT(ELTMJD))*24.0,  &
                              INT(DLTS),KLOC,ILOC,INT(MINDLT),  &
                              INT(JDMIN),(JDMIN-INT(JDMIN))*24.0,KMIN,  &
                              IMIN
            IF (LIMITING_DLT) WRITE (SNP,3020) KLIM,ILIM
            WRITE (SNP,3030)  INT(DLTAV),NIT,NV
            WRITE (SNP,3040) 'Meteorlogical Parameters',TAIR,DEG,TDEW,  &
                              DEG,WIND/(WSC(WSCDP)+NONZERO),PHI,  &
                              WSC(WSCDP),CLOUD,ET,DEG,CSHE,SRO,DEG
            WRITE (SNP,3050) 'Inflows','Upstream inflows'
            DO JB=1,NBP
              IF (UP_FLOW(JB)) THEN
                WRITE (SNP,3060) JB,KTQIN(JB),KBQIN(JB),QIN(JB),TIN(JB),  &
                                 DEG
              END IF
            END DO
            IF (TRIBUTARIES) THEN
              WRITE (SNP,3090) (ITR(JT),JT=1,NTR)
              WRITE (SNP,3100) (KTTR(JT),KBTR(JT),JT=1,NTR)
              WRITE (SNP,3110) (QTR(JT),JT=1,NTR)
              WRITE (SNP,3120) (TTR(JT),JT=1,NTR)
            END IF
            WRITE (SNP,3130) 'Outflows'
            DO JB=1,NBP
              IF (DN_FLOW(JB)) THEN
                IF (SEL_WITHDRAWAL(JB)) THEN
                  WRITE (SNP,3140) JB,(QSTR(JS,JB),JS=1,NSTR(JB))
                END IF
                WRITE (SNP,3150) QSUM(JB),(K,K=KT,KB(DS(JB)))
                WRITE (SNP,3160) (QOUT(K,JB),K=KT,KB(DS(JB)))
              END IF
            END DO
            IF (WITHDRAWALS) THEN
              WRITE (SNP,3170) (MAX(CUS(JBWD(JW)),IWD(JW)),JW=1,NWD)
              WRITE (SNP,3180) (KWD(JW),JW=1,NWD)
              WRITE (SNP,3190) (QWD(JW),JW=1,NWD)
            END IF
            IF (CONSTITUENTS) THEN
              WRITE (SNP,3200) 'Constituent Inflow Concentrations'
              DO JB=1,NBP
                IF (UP_FLOW(JB)) THEN
                  WRITE (SNP,3210) JB,(CNAME1(INCN(JC)),  &
                                   CIN(INCN(JC),JB),CUNIT2(INCN(JC)),  &
                                   JC=1,NACIN)
                END IF
              END DO
              DO JT=1,NTR
                WRITE (SNP,3220) JT,(CNAME1(TRCN(JC)),CTR(TRCN(JC),JT),  &
                                 CUNIT2(JC),JC=1,NACTR)
              END DO
            END IF
            IF (EVAPORATION.OR.PRECIPITATION) WRITE(SNP,3240) 'Surface',  &
                                              ' Calculations'
            IF (EVAPORATION)   WRITE (SNP,3250) (JB,EVBR(JB),JB=1,NBP)
            IF (PRECIPITATION) WRITE (SNP,3260) (JB,QPR(JB),JB=1,NBP)
            IF (VOLUME_BALANCE) THEN
              WRITE (SNP,3300)
              DO JB=1,NBP
                IF (VOLSBR(JB).NE.0.0) DLVR = (VOLTBR(JB)-VOLSBR(JB))  &
                                              /VOLSBR(JB)
                WRITE (SNP,3310) JB,VOLSBR(JB),VOLTBR(JB),VOLTBR(JB)  &
                                 -VOLSBR(JB),DLVR*100.0
              END DO
            END IF
            IF (ENERGY_BALANCE) THEN
              WRITE (SNP,3320)
              DO JB=1,NBP
                WRITE (SNP,3330) JB
                IF (ESBR(JB).NE.0.0) DLE = (ESBR(JB)-ETBR(JB))/ESBR(JB)
                WRITE (SNP,3340) ESBR(JB)*4.184E3,ETBR(JB)*4.1843E3,  &
                                 (ESBR(JB)-ETBR(JB))*4.1843E3,DLE*100.0
              END DO
            END IF
            IF (MASS_BALANCE) THEN
              WRITE (SNP,3350)
              DO JB=1,NBP
                WRITE (SNP,3330) JB
                DO JC=1,NAC
                  JAC = CN(JC)
                  IF (TRANSPORT(JAC)) THEN
                    IF (CMBRS(JAC,JB).NE.0.0) THEN 
                      DLMR = (CMBRT(JAC,JB)-CMBRS(JAC,JB))  &
                             /(CMBRS(JAC,JB)+NONZERO)
                    END IF
                    WRITE (SNP,3360) CNAME1(JAC),CMBRS(JAC,JB),  &
                                     CUNIT3(JAC),CMBRT(JAC,JB),  &
                                     CUNIT3(JAC),(CMBRT(JAC,JB)  &
                                     -CMBRS(JAC,JB)),CUNIT3(JAC),  &
                                     DLMR*100.0
                  END IF
                END DO
              END DO
            END IF
            WRITE (SNP,3370) 'Geometry',KT,ELKT,(JB,CUS(JB),JB=1,NBP)
!
            CALL PRINT_GRID (JDAY,GDAY,MONTH,YEAR)
!
          END IF
        END IF

!****** Vertical profiles

        IF (PROFILE) THEN
          IF (JDAY.GE.NXTMPR.OR.JDAY.GE.PRFD(PRFDP+1)) THEN
            IF (JDAY.GE.PRFD(PRFDP+1)) THEN
              PRFDP  = PRFDP+1
              NXTMPR = PRFD(PRFDP)
            END IF
            NXTMPR = NXTMPR+PRFF(PRFDP)
            NSPRF  = NSPRF+1
            IF (M.EQ.9)                       MON = MONTH(1:3)
            IF (M.EQ.2.OR.M.EQ.11.OR.M.EQ.12) MON = MONTH(2:4)
            IF (M.EQ.1.OR.M.EQ.10)            MON = MONTH(3:5)
            IF (M.EQ.8)                       MON = MONTH(4:6)
            IF (M.EQ.3.OR.M.EQ.4)             MON = MONTH(5:7)
            IF (M.EQ.6.OR.M.EQ.7)             MON = MONTH(6:8)
            IF (M.EQ.5)                       MON = MONTH(7:9)
            WRITE (PRF,2590) JDAY,MON,GDAY,YEAR,KT,SNGL(Z(DS(1))),  &
                             NSPRF
            DO JC=1,NAC
              IF (CPRC(CN(JC)).EQ.' ON') THEN
                DO JPRF=1,NIPRF
                  I   = IPRF(JPRF)
                  NRS = KB(I)-KT+1
                  WRITE (PRF,2560) CN(JC),NRS,(C2(K,I,CN(JC)),  &
                                   K=KT,KB(I))
                END DO
              END IF
            END DO
            DO JPRF=1,NIPRF
              I   = IPRF(JPRF)
              NRS = KB(I)-KT+1
              WRITE (PRF,2560) 22,NRS,(T2(K,I),K=KT,KB(I))
            END DO
          END IF
        END IF

!****** Spreadsheet

        IF (SPREADSHEET) THEN
          IF (JDAY.GE.NXTMSP.OR.JDAY.GE.SPRD(SPRDP+1)) THEN
            IF (JDAY.GE.SPRD(SPRDP+1)) THEN
              SPRDP  = SPRDP+1
              NXTMSP = SPRD(SPRDP)
            END IF
            NXTMSP = NXTMSP+SPRF(SPRDP)
            DO J=1,NISPR
              KBMAX = MAX(KB(ISPR(J)),KBMAX)
              DO K=KT,KB(ISPR(J))
                WRITE (CONV1(K,J),'(F9.2)') T2(K,ISPR(J))
              END DO
            END DO
            DO K=KT,KBMAX
              WRITE (SPR,2570) 'Temperature      ',JDAY,-DEPTHM(K),  &
                               (CONV1(K,J),J=1,NISPR)
            END DO
            DO JC=1,NAC
              IF (CPRC(CN(JC)).EQ.' ON') THEN
                DO J=1,NISPR
                  DO K=KT,KB(ISPR(J))
                    WRITE (CONV1(K,J),'(1PE9.2)') C2(K,ISPR(J),CN(JC))
                  END DO
                END DO
                DO K=KT,KBMAX
                  WRITE (SPR,2570) CNAME1(CN(JC)),JDAY,-DEPTHM(K),  &
                                   (CONV1(K,J),J=1,NISPR)
                END DO
              END IF
            END DO
          END IF
        END IF

!****** Velocity vectors

        IF (VECTOR) THEN
          IF (JDAY.GE.NXTMVP.OR.JDAY.GE.VPLD(VPLDP+1)) THEN
            IF (JDAY.GE.VPLD(VPLDP+1)) THEN
              VPLDP  = VPLDP+1
              NXTMVP = VPLD(VPLDP)
            END IF
            NXTMVP = NXTMVP+VPLF(VPLDP)
            WRITE (VPL,*) JDAY,MONTH//GDCH//', ',YEAR,KT,US
            WRITE (VPL,*) Z
            WRITE (VPL,*) U
            WRITE (VPL,*) W
          END IF
        END IF

!****** Contours

        IF (CONTOUR) THEN
          IF (JDAY.GE.NXTMCP.OR.JDAY.GE.CPLD(CPLDP+1)) THEN
            IF (JDAY.GE.CPLD(CPLDP+1)) THEN
              CPLDP  = CPLDP+1
              NXTMCP = CPLD(CPLDP)
            END IF
            NXTMCP = NXTMCP+CPLF(CPLDP)
            WRITE (CPL,8040) JDAY,MONTH,GDAY,YEAR
            DO JB=1,NBP
              WRITE(CPL,8010) CUS(JB)
              WRITE(CPL,8010) KT
              WRITE(CPL,8050) QIN(JB),QSUM(JB)
              WRITE(CPL,8050) (Z(I),I=CUS(JB),DS(JB))
              DO I=CUS(JB),DS(JB)
                WRITE (CPL,8050) (T2(K,I),K=KT,KB(I))
              END DO
              DO JC=1,NAC
                DO I=CUS(JB),DS(JB)
                  IF (CPRC(CN(JC)).EQ.' ON') THEN
                    WRITE (CPL,8050) (C2(K,I,CN(JC)),K=KT,KB(I))
                  END IF
                END DO
              END DO
            END DO
          END IF
        END IF

!****** Time series

        IF (TIME_SERIES) THEN
          IF (JDAY.GE.NXTMTS.OR.JDAY.GE.TSRD(TSRDP+1)) THEN
            IF (JDAY.GE.TSRD(TSRDP+1)) THEN
              TSRDP  = TSRDP+1
              NXTMTS = TSRD(TSRDP)
            END IF
            NXTMTS = NXTMTS+TSRF(TSRDP)
            DO JB=1,NBP
              IF (DN_FLOW(JB)) THEN
                AVTOUT = 0.0
                DO JC=1,NAC
                  AVCOUT(CN(JC)) = 0.0
                END DO
                DO K=KT,KB(DS(JB))
                  AVTOUT = AVTOUT+T2(K,DS(JB))*QOUT(K,JB)
                  DO JC=1,NAC
                    AVCOUT(CN(JC)) = AVCOUT(CN(JC))+QOUT(K,JB)  &
                                     *C2(K,DS(JB),CN(JC))
                  END DO
                END DO
                IF (QSUM(JB).GT.0.0) THEN
                  AVTOUT = AVTOUT/QSUM(JB)
                  DO JC=1,NAC
                    AVCOUT(CN(JC)) = AVCOUT(CN(JC))/QSUM(JB)
                  END DO
                END IF
              END IF
            END DO
            IF (WITHDRAWALS) THEN
              DO JW=1,NWD
                TWD(JW) = T2(MAX(KT,KWD(JW)),MAX(CUS(JBWD(JW)),IWD(JW)))
                DO JC=1,NAC
                  CWD(CN(JC),JW) = C2(MAX(KT,KWD(JW)),MAX(CUS(JBWD(JW)),  &
                                   IWD(JW)),CN(JC))
                END DO
              END DO
            END IF
            IF (WITHDRAWALS) THEN
              WRITE (TSR,5030) JDAY,DLT,ELKT,AVTOUT,(AVCOUT(CN(JC)),  &
                               JC=1,NAC),(TWD(JW),JW=1,NWD),  &
                               ((CWD(CN(JC),JW),JC=1,NAC),JW=1,NWD)
            ELSE
              WRITE (TSR,5030) JDAY,DLT,ELKT,AVTOUT,(AVCOUT(CN(JC)),  &
                               JC=1,NAC)
            END IF
          END IF
        END IF


!****** Screen output

        IF (SCREEN_OUTPUT) THEN
          IF (JDAY.GE.NXTMSC.OR.JDAY.GE.SCRD(SCRDP+1)) THEN
            IF (JDAY.GE.SCRD(SCRDP+1)) THEN
              SCRDP  = SCRDP+1
              NXTMSC = SCRD(SCRDP)
            END IF
            NXTMSC = NXTMSC+SCRF(SCRDP)
            WRITE(*,9100) JDAY,INT(DLT),INT(DLTAV),NIT,NV
          END IF
        END IF
      END DO

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                            Task 3: End Simulation                        **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      IF (SNAPSHOT) THEN
!        CALL DATE_TIME (CDATE,CTIME)                                    !FORTRAN
        WRITE (SNP,*)                                                   !FORTRAN
        WRITE (SNP,4040) 'Normal termination at '//CTIME//' on '//CDATE !FORTRAN
        WRITE (SNP,*)
        WRITE (SNP,4050) 'Runtime statistics'
        WRITE (SNP,4060) '  Grid',IMP,KMP,NTACMX,NTACMN
        WRITE (SNP,4070) '    Segment lengths      =',DLXMIN,DLXMAX
        WRITE (SNP,4080) '    Layer heights        =',HMIN,HMAX
        WRITE (SNP,4090) '  Timestep'
        WRITE (SNP,4100) '    Total iterations     =',NIT
        WRITE (SNP,4100) '    # of violations      =',NV
        WRITE (SNP,4110) '    Average timestep     =',INT(DLTAV)
        WRITE (SNP,4120) '    Simulation time      =',INT(ELTMJD),  &
                         (ELTMJD-INT(ELTMJD))*24.0
        WRITE (SNP,'(A3)')  ESC//'E'
        CLOSE (SNP)
      END IF
      IF (WATER_BALANCE) THEN
        WRITE (QWS,9010) TMEND,0.0
        CLOSE (QWS)
      END IF
      IF (TIME_SERIES)  CLOSE (TSR)
      IF (VECTOR)       CLOSE (VPL)
      IF (PROFILE)      CLOSE (PRF)
      IF (SPREADSHEET)  CLOSE (SPR)
      IF (CONTOUR)      CLOSE (CPL)
      IF (.NOT.WARNING_OPEN) THEN
        CLOSE (WRN,STATUS='DELETE')
      ELSE
        CLOSE (WRN)
      END IF
      CLOSE (ERR,STATUS='DELETE')

!**** Profile FORMATs

 2500 FORMAT(A72)
 2510 FORMAT(20(1X,A3))
 2520 FORMAT(3(1X,A17,1X,A6))
 2530 FORMAT(L2,(:/20I4))
 2540 FORMAT(20I4)
 2550 FORMAT(10F8.2)
 2560 FORMAT(2I4/(8(1PE10.2E2)))
 2570 FORMAT(A17,2F10.3,17A9)
 2580 FORMAT(17X,F10.3)
 2590 FORMAT(F8.3,1X,A3,I3,', ',2I4,F8.4,I4)
 2600 FORMAT('Constituent      Julian_day     Depth',17(2X,A7))
 2610 FORMAT(37X,17(2X,A7))

!**** Snapshot FORMATs

 3000 FORMAT('1','CE-QUAL-W2 V2.10 - June, 1995'/7(1X,A72/))
 3010 FORMAT(1X,A15/                                           &
             '+',4('_'),1X,10('_')//  &
             3X,'Gregorian date      [GDAY] =',A16,I3,',',I5/  &
             3X,'Julian date         [JDAY] =',I8,' days',F6.2,  &
               ' hours'/  &
             3X,'Elapsed time      [ELTMJD] =',I8,' days',F6.2,  &
               ' hours'/  &
             3X,'Timestep             [DLT] =',I8,' sec'/  &
             5X,'at location  [KLOC,ILOC] = (',I2,',',I2,')'/  &
             3X,'Minimum timestep  [MINDLT] =',I8,' sec '/  &
             5X,'at Julian day    [JDMIN] =',I8,' days',F6.2,' hours'/  &
             5X,'at location  [KMIN,IMIN] = (',I2,',',I2,')')
 3020 FORMAT(3X,'Limiting timestep'/  &
             5X,'at location  [KLIM,ILIM] = (',I2,',',I2,')')  
 3030 FORMAT(3X,'Average timestep   [DLTAV] =',I8,' sec'/  &
             3X,'Number of iterations [NIT] =',I8/  &
             3X,'Number of violations  [NV] =',I8/)
 3040 FORMAT(1X,A24/  &
             '+',13('_'),1X,10('_')//  &
             3X,'Input'/  &
             5X,'Air temperature       [TAIR] =',F9.2,1X,A2/  &
             5X,'Dewpoint temperature  [TDEW] =',F9.2,1X,A2/  &
             5X,'Wind speed            [WIND] =',F9.2,' m/sec'/  &
             5X,'Wind direction         [PHI] =',F9.2,' rad'/  &
             5X,'Wind sheltering        [WSC] =',F9.2/  &
             5X,'Cloud cover          [CLOUD] =',F9.2/  &
             3X,'Calculated'/  &
             5X,'Equilibrium temperature [ET] =',F9.2,1X,A2/  &
             5X,'Surface heat exchange [CSHE] =',1PE9.2,' m/sec'/  &
             5X,'Short wave radiation   [SRO] =',1PE9.2,1X,A2,' m/sec'/)
 3050 FORMAT(1X,A7/  &
             '+',7('_')//  &
             3X,A16)
 3060 FORMAT(5X,'Branch',I3/  &
             7X,'Layer       [KQIN] =',I5,'-',I2/  &
             7X,'Inflow       [QIN] =',F8.2,' m^3/sec'/  &
             7X,'Temperature  [TIN] =',F8.2,1X,A2)
 3070 FORMAT(/3X,'Distributed Tributaries'/)
 3080 FORMAT(5X,'Branch',I3/  &
             7X,'Inflow      [QDTR] =',F8.2,' m^3/sec'/  &
             7X,'Temperature [TDTR] =',F8.2,1X,A2)
 3090 FORMAT('+'://3X,'Tributaries'/  &
             5X,'Segment     [ITR] =',8I8/(T25,8I8))
 3100 FORMAT('+':/5X,'Layer       [KTR] =',8(I5,'-',I2)/  &
            (T25,8(I5,'-',I2)))
 3110 FORMAT('+':/5X,'Inflow      [QTR] =',8F8.1/  &
            (T25,8F8.1))
 3120 FORMAT('+':/5X,'Temperature [TTR] =',8F8.1/  &
            (T25,8F8.1))
 3130 FORMAT(/1X,A8/  &
            '+',8('_'))
 3140 FORMAT(/3X,'Structure outflows [QSTR]'/  &
             5X,'Branch ',I2,' = ',8F8.2/  &
            (T16,8F8.2))
 3150 FORMAT('+'://  &
             3X,'Total outflow [QOUT] =',F8.2,' m^3/s'//  &
             3X,'Outlets'/  &
             5X,'Layer             [KOUT] =',12I7/(31X,12I7))
 3160 FORMAT('+':/  &
             5X,'Outflow (m^3/sec) [QOUT] =',12F7.1/(31X,12F7.1))
 3170 FORMAT('+'://  &
             3X,'Withdrawals'/  &
             '+',2X,11('_')/  &
             5X,'Segment           [IWD] =',12I7/(30X,12I7))
 3180 FORMAT(:5X,'Layer             [KWD] =',12I7/(30X,12I7))
 3190 FORMAT(:5X,'Outflow (m^3/sec) [QWD] =',12F7.1/(30X,12F7.1))
 3200 FORMAT('1',A33/  &
             '+',11('_'),1X,6('_'),1X,14('_')/)
 3210 FORMAT(3X,'Branch',I3,' [CIN]'/  &
            (5X,A17,T24,'=',F8.3,1X,A6))
 3220 FORMAT(3X,'Tributary',I3,' [CTR]'/  &
            (5X,A17,T24,'=',F8.3,1X,A6))
 3230 FORMAT(3X,'Distributed tributary',I3,' [CDT]'/  &
            (5X,A17,T24,'=',F8.3,1X,A6))
 3240 FORMAT(1X,A7,A13/  &
             '+',7('_'),1X,12('_')/)
 3250 FORMAT(3X,'Evaporation [EV]'/  &
             '+',2X,11('_')/  &
            (:/5X,'Branch',I3,' =',F8.2))
 3260 FORMAT(3X,'Precipitation [PR]'/  &
             '+',2X,13('_')//  &
            (5X,'Branch',I3,' =',F8.2/))
 3270 FORMAT(/1X,'External head boundary elevations'/  &
             '+',8('_'),1X,4('_'),1X,8('_'),1X,10('_')/)
 3280 FORMAT(3X,'Branch',I3/  &
             5X,'Upstream elevation   [ELUH] =',F8.3,' m')
 3290 FORMAT(3X,'Branch',I3/  &
             5X,'Downstream elevation [ELDH] =',F8.3,' m')
 3300 FORMAT(/1X,'Water Balance'/  &
             '+',5('_'),1X,7('_')/)
 3310 FORMAT(3X,'Branch',I3/  &
             5X,'Spatial change  [VOLSBR] = ',1PE15.8E2,' m^3'/  &
             5X,'Temporal change [VOLTBR] = ',1PE15.8E2,' m^3'/  &
             5X,'Volume error             = ',1PE15.8E2,' m^3'/  &
             5X,'Percent error            = ',1PE15.8E2,' %')
 3320 FORMAT(/1X,'Energy Balance'/  &
             '+',6('_'),1X,7('_')/)
 3330 FORMAT(3X,'Branch',I3)
 3340 FORMAT(5X,'Spatially integrated energy  [ESBR] = ',1PE15.8E2,  &
               ' kJ'/  &
             5X,'Temporally integrated energy [ETBR] = ',1PE15.8E2,  &
               ' kJ'/  &
             5X,'Energy error                        = ',1PE15.8E2,  &
               ' kJ'/  &
             5X,'Percent error                       = ',1PE15.8E2,  &
               ' %')  
 3350 FORMAT(/1X,'Mass Balance'/  &
             '+',4('_'),1X,7('_')/)
 3360 FORMAT(5X,A17/  &
             7X,'Spatially integrated mass  [CMBRS] = ',1PE15.8E2,1X,A2/  &
             7X,'Temporally integrated mass [CMBRT] = ',1PE15.8E2,1X,A2/  &
             7X,'Mass error                         = ',1PE15.8E2,1X,A2/  &
             7X,'Percent error                      = ',1PE15.8E2,' %')
 3370 FORMAT(/1X,A8/  &
             '+',8('_')//  &
             3X,'Surface layer [KT] =',I8/  &
             3X,'Elevation   [ELKT] =',F8.3,' m'//  &
             3X,'Current upstream segment [CUS]'/  &
             '+',2X,7('_'),1X,8('_'),1X,7('_')//  &
             (5X,'Branch',I3,' =',I3))

!**** Run time information FORMATs

 4000 FORMAT(/1X,20('*'),' Add layer',I3,' at Julian day =',F9.3,  &
               1X,20('*'))
 4010 FORMAT(/1X,11('*'),' Add segments ',I3,' through ',I3,  &
               ' at Julian day =',F9.3,1X,10('*'))
 4020 FORMAT(/1X,18('*'),' Subtract layer',I3,' at Julian day =',F9.3,  &
               1X,17('*'))
 4030 FORMAT(/1X,8('*'),' Subtract segments ',I3,' through ',I3,  &
               ' at Julian day =',F9.3,1X,8('*'))
 4040 FORMAT(//17('*'),2X,A42,2X,17('*')/)
 4050 FORMAT(/1X,A18)
 4060 FORMAT(1X,A6,' =',I3,' x',I3/  &
             '     Maximum active cells = ',I4/  &
             '     Minimum active cells = ',I4)
 4070 FORMAT(1X,A24,F8.1,' -',F8.1,' m')
 4080 FORMAT(1X,A24,F8.2,' -',F8.2,' m')
 4090 FORMAT(1X,A10)
 4100 FORMAT(1X,A24,I8)
 4110 FORMAT(1X,A24,I8,' sec')
 4120 FORMAT(1X,A24,I8,' days',F6.2,' hours')

!**** Time series FORMAT's

 5000 FORMAT(A72)
 5010 FORMAT(30I3)
 5020 FORMAT(4(A17,1X,A6))
 5030 FORMAT(F10.3,2F10.2,100(1PE10.3E2))

!**** Run time warning FORMATs

 6000 FORMAT('Computational warning at Julian day = ',F10.3/  &
             '  timestep        = ',F4.0,' sec')
 6010 FORMAT('Computational warning at Julian day = ',F10.3/  &
             '  spatial change  =',1PE15.8E2,' m^3'/  &
             '  temporal change =',1PE15.8E2,' m^3'/  &
             '  volume error    =',1PE15.8E2,' m^3')
 6020 FORMAT('**SEVERE** computational warning at Julian day = ',F10.3/  &
             '  at segment ',I3/  &
             '    water surface elevation [Z]  =',F10.3,' m'/  &
             '    layer thickness              =',F10.3,' m')

!**** Run time error FORMATs

 7000 FORMAT('Fatal error'/  &
             '  Insufficient segments in branch',I3/  &
             '  Julian day = ',F9.2,'  water surface layer =',I3)
 7010 FORMAT('Fatal error'/  &
             '  Negative surface layer height'/  &
             '  Julian day = ',F9.2,'  segment =',I3,'  height = ',  &
             F8.2,' m')

!**** Contour plot FORMATs

 8000 FORMAT(8(I8,2X))
 8010 FORMAT(9(I8,2X))
 8020 FORMAT(8(E13.6,2X))
 8030 FORMAT(A17)
 8040 FORMAT(F12.4,5X,A9,5X,I2,5X,I4)
 8050 FORMAT(9(E13.6,2X))

!**** Water balance FORMATs

 9000 FORMAT('Computed flows to match observed elevations'//  &
             '    JDAY     QDT')
 9010 FORMAT(F8.3,F8.1)

!**** Screen FORMAT's

 9100 FORMAT('JDAY =',F9.3,' DLT =',I5,' DLTAV =',I5,' NIT =',I6,  &
             ' NV =',I5)
      STOP
      END

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                        F U N C T I O N   D E N S I T Y                   **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      FUNCTION   DENSITY (T,SS,DS)

!****** Type declaration

        LOGICAL FRESH_WATER, SALT_WATER, SUSP_SOLIDS

!****** Common declaration

        COMMON  /DNSPHC/ FRESH_WATER, SALT_WATER, SUSP_SOLIDS

        DS0     = MAX(DS,0.0)
        DENSITY = ((((6.536332E-9*T-1.120083E-6)*T+1.001685E-4)  &
                  *T-9.09529E-3)*T+6.793952E-2)*T+0.842594
        IF (SUSP_SOLIDS) DENSITY = DENSITY+6.2E-4*SS
        IF (FRESH_WATER) DENSITY = DENSITY+DS0*((4.99E-8*T-3.87E-6)  &
                                   *T+8.221E-4)
        IF (SALT_WATER)  DENSITY = DENSITY+DS0*((((5.3875E-9  &
                                   *T-8.2467E-7)*T+7.6438E-5)  &
                                   *T-4.0899E-3)*T+0.824493)  &
                                   +((-1.6546E-6*T+1.0227E-4)  &
                                   *T-5.72466E-3)*DS0**1.5+4.8314E-4  &
                                   *DS0*DS0
        DENSITY = DENSITY+999.0
      END

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*               S U B R O U T I N E    G R E G O R I A N   D A T E         **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      SUBROUTINE GREGORIAN_DATE (YEAR)

!****** Type declaration

        INTEGER   YEAR, GDAY
        LOGICAL   LEAP_YEAR
        CHARACTER MONTH*9
      
!****** Common declaration

        COMMON /GDAYC1/ GDAY, DAYM, JDAYG, LEAP_YEAR, M
        COMMON /GDAYC2/ MONTH

!****** Determine if new year (regular or leap) and increment year

        IF (.NOT.LEAP_YEAR.AND.JDAYG.GE.366) THEN
          JDAYG     = JDAYG-365
          YEAR      = YEAR+1
          LEAP_YEAR = MOD(YEAR,4).EQ.0
        ELSE IF (JDAYG.GE.367) THEN
          JDAYG     = JDAYG-366
          YEAR      = YEAR+1
          LEAP_YEAR = .FALSE.
        END IF

!****** Determine month and day of year

        IF (LEAP_YEAR) THEN
          IF (JDAYG.GE.1.AND.JDAYG.LT.32) THEN
            GDAY  = JDAYG
            DAYM  = 31.0
            MONTH = '  January'
          ELSE IF (JDAYG.GE.32.AND.JDAYG.LT.61) THEN
            GDAY  = JDAYG-31
            DAYM  = 29.0
            MONTH = ' February'
          ELSE IF (JDAYG.GE.61.AND.JDAYG.LT.92) THEN
            GDAY  = JDAYG-60
            DAYM  = 31.0
            MONTH = '    March'
          ELSE IF (JDAYG.GE.92.AND.JDAYG.LT.122) THEN
            GDAY  = JDAYG-91
            DAYM  = 30.0
            MONTH = '    April'
          ELSE IF (JDAYG.GE.122.AND.JDAYG.LT.153) THEN
            GDAY  = JDAYG-121
            DAYM  = 31.0
            MONTH = '      May'
          ELSE IF (JDAYG.GE.153.AND.JDAYG.LT.183) THEN
            GDAY  = JDAYG-152
            DAYM  = 30.0
            MONTH = '     June'
          ELSE IF (JDAYG.GE.181.AND.JDAYG.LT.214) THEN
            GDAY  = JDAYG-182
            DAYM  = 31.0
            MONTH = '     July'
          ELSE IF (JDAYG.GE.214.AND.JDAYG.LT.245) THEN
            GDAY  = JDAYG-213
            DAYM  = 31.0
            MONTH = '   August'
          ELSE IF (JDAYG.GE.245.AND.JDAYG.LT.275) THEN
            GDAY  = JDAYG-244
            DAYM  = 30.0
            MONTH = 'September'
          ELSE IF (JDAYG.GE.275.AND.JDAYG.LT.306) THEN
            GDAY  = JDAYG-274
            DAYM  = 31.0
            MONTH = '  October'
          ELSE IF (JDAYG.GE.306.AND.JDAYG.LT.336) THEN
            GDAY  = JDAYG-305
            DAYM  = 30.0
            MONTH = ' November'
          ELSE IF (JDAYG.GE.336.AND.JDAYG.LT.367) THEN
            GDAY  = JDAYG-335
            DAYM  = 31.0
            MONTH = ' December'
          END IF
        ELSE
          IF (JDAYG.GE.1.AND.JDAYG.LT.32) THEN
            GDAY  = JDAYG
            DAYM  = 31.0
            MONTH = '  January'
          ELSE IF (JDAYG.GE.32.AND.JDAYG.LT.60) THEN
            GDAY  = JDAYG-31
            DAYM  = 29.0
            MONTH = ' February'
          ELSE IF (JDAYG.GE.60.AND.JDAYG.LT.91) THEN
            GDAY  = JDAYG-59
            DAYM  = 31.0
            MONTH = '    March'
          ELSE IF (JDAYG.GE.91.AND.JDAYG.LT.121) THEN
            GDAY  = JDAYG-90
            DAYM  = 30.0
            MONTH = '    April'
          ELSE IF (JDAYG.GE.121.AND.JDAYG.LT.152) THEN
            GDAY  = JDAYG-120
            DAYM  = 31.0
            MONTH = '      May'
          ELSE IF (JDAYG.GE.152.AND.JDAYG.LT.182) THEN
            GDAY  = JDAYG-151
            DAYM  = 30.0
            MONTH = '     June'
          ELSE IF (JDAYG.GE.182.AND.JDAYG.LT.213) THEN
            GDAY  = JDAYG-181
            DAYM  = 31.0
            MONTH = '     July'
          ELSE IF (JDAYG.GE.213.AND.JDAYG.LT.244) THEN
            GDAY  = JDAYG-212
            DAYM  = 31.0
            MONTH = '   August'
          ELSE IF (JDAYG.GE.244.AND.JDAYG.LT.274) THEN
            GDAY  = JDAYG-243
            DAYM  = 30.0
            MONTH = 'September'
          ELSE IF (JDAYG.GE.274.AND.JDAYG.LT.305) THEN
            GDAY  = JDAYG-273
            DAYM  = 31.0
            MONTH = '  October'
          ELSE IF (JDAYG.GE.305.AND.JDAYG.LT.335) THEN
            GDAY  = JDAYG-304
            DAYM  = 30.0
            MONTH = ' November'
          ELSE IF (JDAYG.GE.335.AND.JDAYG.LT.366) THEN
            GDAY  = JDAYG-334
            DAYM  = 31.0
            MONTH = ' December'
          END IF
        END IF
      END

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*             S U B R O U T I N E   T I M E   V A R Y I N G   D A T A       **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      SUBROUTINE TIME_VARYING_DATA (JDAY,NXTVD)
        INCLUDE 'w2.inc'
        SAVE

!****** Type declaration

        INTEGER  OTQ,    TRQ,    WDQ,    TRT,    TRC,    DTQ,    DTT,  &
                 DTC,    PRE,    PRT,    PRC,    UHE,    UHT,    UHC,  &
                 DHE,    DHT,    DHC,    ELO
        INTEGER  CN,     INCN,   TRCN,   DTCN,   PRCN,   UHCN,   DHCN
        INTEGER  UNIT,   WSCDP
        INTEGER  CUS,    DS
        REAL     JDAY,   NXTVD,  LAT,    LONG
        REAL     NXMET1, NXQOT1, NXQIN1, NXTIN1, NXCIN1, NXQWD1, NXQTR1,  &
                 NXTTR1, NXCTR1, NXQDT1, NXTDT1, NXCDT1, NXPR1,  NXTPR1,  &
                 NXCPR1, NXEUH1, NXTUH1, NXCUH1, NXEDH1, NXTDH1, NXCDH1   
        REAL     NXMET2, NXQOT2, NXQIN2, NXTIN2, NXCIN2, NXQWD2, NXQTR2,  &
                 NXTTR2, NXCTR2, NXQDT2, NXTDT2, NXCDT2, NXEUH2, NXTUH2,  &
                 NXCUH2, NXEDH2, NXTDH2, NXCDH2, NXELV1, NXELV2
        LOGICAL  OPEN_FILES,     CONSTITUENTS, WITHDRAWALS, TRIBUTARIES,  &
                 PRECIPITATION,  DIST_TRIBS,   UH_EXTERNAL, DH_EXTERNAL,  &
                 NO_WIND,        NO_INFLOW,    NO_OUTFLOW,  NO_HEAT,  &
                 UP_FLOW,        DN_FLOW,      POINT_SINK   
        LOGICAL  INTERP_INFLOW,  INTERP_OUTFLOW, INTERP_MET,  &
                 INTERP_DTRIBS,  INTERP_TRIBS,   INTERP_HEAD,  &
                 INTERP_WITHDRWL
        LOGICAL  SEL_WITHDRAWAL, TERM_BY_TERM, WATER_BALANCE
        CHARACTER*72 METFN, QINFN, TINFN, CINFN, QOTFN, QWDFN, QTRFN,  &
                     TTRFN, CTRFN, QDTFN, TDTFN, CDTFN, PREFN, TPRFN,  &
                     CPRFN, EUHFN, TUHFN, CUHFN, EDHFN, TDHFN, CDHFN,  &
                     ELOFN
        DOUBLE PRECISION Z

!****** Dimension declaration

        DIMENSION TRQ(NTP), TRT(NTP), TRC(NTP)
        DIMENSION INQ(NBP), DTQ(NBP), PRE(NBP), UHE(NBP), DHE(NBP)
        DIMENSION INT(NBP), DTT(NBP), PRT(NBP), UHT(NBP), DHT(NBP)
        DIMENSION INC(NBP), DTC(NBP), PRC(NBP), UHC(NBP), DHC(NBP)
        DIMENSION OTQ(NBP)
        DIMENSION NXQTR1(NTP),        NXTTR1(NTP),   NXCTR1(NTP),  &
                  NXQIN1(NBP),        NXTIN1(NBP),   NXCIN1(NBP),  &
                  NXQDT1(NBP),        NXTDT1(NBP),   NXCDT1(NBP),  &
                  NXPR1(NBP),         NXTPR1(NBP),   NXCPR1(NBP),  &
                  NXEUH1(NBP),        NXTUH1(NBP),   NXCUH1(NBP),  &
                  NXEDH1(NBP),        NXTDH1(NBP),   NXCDH1(NBP),  &
                  NXQOT1(NBP)
        DIMENSION NXQTR2(NTP),        NXTTR2(NTP),   NXCTR2(NTP),  &
                  NXQIN2(NBP),        NXTIN2(NBP),   NXCIN2(NBP),  &
                  NXQDT2(NBP),        NXTDT2(NBP),   NXCDT2(NBP),  &
                  NXEUH2(NBP),        NXTUH2(NBP),   NXCUH2(NBP),  &
                  NXEDH2(NBP),        NXTDH2(NBP),   NXCDH2(NBP),  &
                  NXQOT2(NBP)  
        DIMENSION CTRNX(NCP,NTP),     CINNX(NCP,NBP),  &
                  QOUTNX(KMP,NBP),    CDTRNX(NCP,NBP),  &
                  CPRNX(NCP,NBP),     TUHNX(KMP,NBP),  &
                  TDHNX(KMP,NBP),     QSTRNX(NSP,NBP),  &
                  CUHNX(KMP,NCP,NBP), CDHNX(KMP,NCP,NBP),  &
                  QDTRNX(NBP),        TDTRNX(NBP),  &
                  PRNX(NBP),          TPRNX(NBP),  &
                  ELUHNX(NBP),        ELDHNX(NBP),  &
                  QWDNX(NWP),         QTRNX(NTP),  &
                  TTRNX(NTP),         QINNX(NBP),  &
                  TINNX(NBP)
        DIMENSION CTRO(NCP,NTP),      CINO(NCP,NBP),  &
                  QOUTO(KMP,NBP),     CDTRO(NCP,NBP),  &
                  TUHO(KMP,NBP),      TDHO(KMP,NBP),  &
                  QSTRO(NSP,NBP),     CUHO(KMP,NCP,NBP),  &
                  CDHO(KMP,NCP,NBP),  QDTRO(NBP),  &
                  TDTRO(NBP),         ELUHO(NBP),  &
                  ELDHO(NBP),         QWDO(NWP),  &
                  QTRO(NTP),          TTRO(NTP),  &
                  QINO(NBP),          TINO(NBP)

!****** Common declaration

        COMMON /SCRNC2/ NWD,    NTR
        COMMON /TVDMTC/ TAIR,   TDEW,    CLOUD,   PHI,   ET,    CSHE,  &
                        SRO,    SRON,    LAT,     LONG,  WINDH, RHOWCP
        COMMON /GLOBLC/ JB,     JC,      IU,      ID,    KT,    ELKT,  &
                        DLT,    KB(IMP), KTI(IMP)
        COMMON /GLBLCC/ PALT,   APOM,    O2LIM,   WIND,  WSCDP, WSC(NDP)
        COMMON /HYDRC2/ Z(IMP)
        COMMON /TVDDSC/ NO_WIND,          NO_INFLOW,     NO_OUTFLOW,  &
                        NO_HEAT
        COMMON /TVDLC1/ PRECIPITATION,    WITHDRAWALS,   TRIBUTARIES,  &
                        DIST_TRIBS(NBP)
        COMMON /TVDLC3/ OPEN_FILES,       TERM_BY_TERM,   WATER_BALANCE
        COMMON /INTERC/ INTERP_INFLOW,    INTERP_OUTFLOW, INTERP_MET,  &
                        INTERP_DTRIBS,    INTERP_TRIBS,   INTERP_HEAD,  &
                        INTERP_WITHDRWL
        COMMON /GRTVDC/ CONSTITUENTS,     CN(NCP),       NAC
        COMMON /SELWC/  NSTR(NBP),        QSTR(NSP,NBP), ESTR(NSP,NBP),  &
                        WSTR(NSP,NBP),    KBSW(NSP,NBP), NOUT(NBP),  &
                        KOUT(KMP,NBP)
        COMMON /TVDFNC/ METFN,            QWDFN,         QOTFN(NBP),  &
                        QINFN(NBP),       TINFN(NBP),    CINFN(NBP),  &
                        QTRFN(NTP),       TTRFN(NTP),    CTRFN(NTP),  &
                        QDTFN(NBP),       TDTFN(NBP),    CDTFN(NBP),  &
                        PREFN(NBP),       TPRFN(NBP),    CPRFN(NBP),  &
                        EUHFN(NBP),       TUHFN(NBP),    CUHFN(NBP),  &
                        EDHFN(NBP),       TDHFN(NBP),    CDHFN(NBP),  &
                        ELOFN
        COMMON /TVDQC/  QIN(NBP),         QTR(NTP),      QDTR(NBP),  &
                        PR(NBP),          ELUH(NBP),     ELDH(NBP),  &
                        QOUT(KMP,NBP),    QWD(NWP),      QSUM(NBP)  
        COMMON /TVDTC/  TIN(NBP),         TTR(NTP),      TDTR(NBP),  &
                        TPR(NBP),         TUH(KMP,NBP),  TDH(KMP,NBP)
        COMMON /TVDCC1/ CIN(NCP,NBP),     CTR(NCP,NTP),  &
                        CDTR(NCP,NBP),    CPR(NCP,NBP),  &
                        CUH(KMC,NCP,NBP), CDH(KMC,NCP,NBP)
        COMMON /TVDCC2/ INCN(NCP),        TRCN(NCP),     DTCN(NCP),  &
                        PRCN(NCP),        UHCN(NCP),     DHCN(NCP)
        COMMON /TVDCC3/ NACIN,            NACTR,         NACPR,  &
                        NACDT
        COMMON /TVDLC2/ UP_FLOW(NBP),        DN_FLOW(NBP),  &
                        UH_INTERNAL(NBP),    UH_EXTERNAL(NBP),  &
                        DH_INTERNAL(NBP),    DH_EXTERNAL(NBP)
        COMMON /TVDSWC/ SEL_WITHDRAWAL(NBP), POINT_SINK(NSP,NBP)
        COMMON /ELBALC/ NXELV1,        NXELV2,         ELWSO1,  &
                        ELWSO2
        COMMON /PRNTC3/ CUS(NBP),      DS(NBP)
        COMMON /GEOMHC/ EL(KMP),       ELWS(IMP),      ELWS2(IMP),  &
                        H(KMP),        HKT1(IMP),      HKT2(IMP),  &
                        DLX(IMP)
        COMMON /GEOMBC/ B(KMP,IMP),    BKT(IMP),       BH(KMP,IMP),  &
                        BHKT1(IMP),    BHKT2(IMP),     BHRKT1(IMP)

!****** Data declaration

        DATA UNIT /30/

!****** Open input files

        IF (OPEN_FILES) THEN
          MET  = UNIT
          UNIT = UNIT+1
          OPEN (MET,FILE=METFN,STATUS='OLD')
          READ (MET,1000)
          IF (WITHDRAWALS) THEN
            WDQ  = UNIT
            UNIT = UNIT+1
            OPEN (WDQ,FILE=QWDFN,STATUS='OLD')
            READ (WDQ,1000)
          END IF
          IF (TRIBUTARIES) THEN
            DO JT=1,NTR
              TRQ(JT) = UNIT
              UNIT    = UNIT+1
              TRT(JT) = UNIT
              UNIT    = UNIT+1
              OPEN (TRQ(JT),FILE=QTRFN(JT),STATUS='OLD')
              OPEN (TRT(JT),FILE=TTRFN(JT),STATUS='OLD')
              READ (TRQ(JT),1000)
              READ (TRT(JT),1000)
              IF (CONSTITUENTS) THEN
                TRC(JT) = UNIT
                UNIT    = UNIT+1
                OPEN (TRC(JT),FILE=CTRFN(JT),STATUS='OLD')
                READ (TRC(JT),1000)
              END IF
            END DO
          END IF
          DO JB=1,NBP
            IF (UP_FLOW(JB)) THEN
              INQ(JB) = UNIT
              UNIT    = UNIT+1
              INT(JB) = UNIT
              UNIT    = UNIT+1
              OPEN (INQ(JB),FILE=QINFN(JB),STATUS='OLD')
              OPEN (INT(JB),FILE=TINFN(JB),STATUS='OLD')
              READ (INQ(JB),1000)
              READ (INT(JB),1000)
              IF (CONSTITUENTS) THEN
                INC(JB) = UNIT
                UNIT    = UNIT+1
                OPEN (INC(JB),FILE=CINFN(JB),STATUS='OLD')
                READ (INC(JB),1000)
              END IF
            END IF
            IF (DN_FLOW(JB)) THEN
              OTQ(JB) = UNIT
              UNIT    = UNIT+1
              OPEN (OTQ(JB),FILE=QOTFN(JB),STATUS='OLD')
              READ (OTQ(JB),1000)
            END IF
            IF (PRECIPITATION) THEN
              PRE(JB) = UNIT
              UNIT    = UNIT+1
              PRT(JB) = UNIT
              UNIT    = UNIT+1
              OPEN (PRE(JB),FILE=PREFN(JB),STATUS='OLD')
              OPEN (PRT(JB),FILE=TPRFN(JB),STATUS='OLD')
              READ (PRE(JB),1000)
              READ (PRT(JB),1000)
              IF (CONSTITUENTS) THEN
                PRC(JB) = UNIT
                UNIT    = UNIT+1
                OPEN (PRC(JB),FILE=CPRFN(JB),STATUS='OLD')
                READ (PRC(JB),1000)
              END IF
            END IF
          END DO
          IF (WATER_BALANCE) THEN
            ELO  = UNIT
            UNIT = UNIT+1
            OPEN (ELO,FILE=ELOFN,STATUS='OLD')
            READ (ELO,1000)
          END IF
          OPEN_FILES = .FALSE.
        END IF
        NXTVD = 1.0E10

!****** Meteorlogical data

        IF (JDAY.GE.NXMET1) THEN
          DO WHILE (JDAY.GE.NXMET1)
            NXMET2 = NXMET1
            TDEW   = TDEWNX
            TDEWO  = TDEWNX
            WIND   = WINDNX
            WINDO  = WINDNX
            PHI    = PHINX
            PHIO   = PHINX
            TAIR   = TAIRNX
            TAIRO  = TAIRNX
            CLOUD  = CLOUDNX
            CLOUDO = CLOUDNX
            READ (MET,1010) NXMET1,TAIRNX,TDEWNX,WINDNX,PHINX,CLOUDNX
            WINDNX = WINDNX*WSC(WSCDP)
          END DO
        END IF
        NXTVD = MIN(NXTVD,NXMET1)

!****** Water surface elevations

        IF (WATER_BALANCE) THEN
          IF (JDAY.GE.NXELV1) THEN
            DO WHILE (JDAY.GE.NXELV1)
              ELWSO2 = ELWSO1
              READ (ELO,1020) NXELV1,ELWSO1
            END DO
          END IF
          NXTVD = MIN(NXTVD,NXELV1)
        END IF

!****** Withdrawals

        IF (WITHDRAWALS) THEN
          IF (JDAY.GE.NXQWD1) THEN
            DO WHILE (JDAY.GE.NXQWD1)
              NXQWD2 = NXQWD1
              DO JW=1,NWD
                QWD(JW)  = QWDNX(JW)
                QWDO(JW) = QWDNX(JW)
              END DO
              READ (WDQ,1020) NXQWD1,(QWDNX(JW),JW=1,NWD)
            END DO
          END IF
          NXTVD = MIN(NXTVD,NXQWD1)
        END IF

!****** Tributaries

        IF (TRIBUTARIES) THEN
          DO JT=1,NTR

!****!***** Inflow

            IF (JDAY.GE.NXQTR1(JT)) THEN
              DO WHILE (JDAY.GE.NXQTR1(JT))
                NXQTR2(JT) = NXQTR1(JT)
                QTR(JT)    = QTRNX(JT)
                QTRO(JT)   = QTRNX(JT)
                READ (TRQ(JT),1020) NXQTR1(JT),QTRNX(JT)
              END DO
            END IF
            NXTVD = MIN(NXTVD,NXQTR1(JT))

!****!***** Inflow temperatures

            IF (JDAY.GE.NXTTR1(JT)) THEN
              DO WHILE (JDAY.GE.NXTTR1(JT))
                NXTTR2(JT) = NXTTR1(JT)
                TTR(JT)    = TTRNX(JT)
                TTRO(JT)   = TTRNX(JT)
                READ (TRT(JT),1020) NXTTR1(JT),TTRNX(JT)
              END DO
            END IF
            NXTVD = MIN(NXTVD,NXTTR1(JT))

!****!***** Inflow constituent concentrations

            IF (CONSTITUENTS) THEN
              IF (JDAY.GE.NXCTR1(JT)) THEN
                DO WHILE (JDAY.GE.NXCTR1(JT))
                  NXCTR2(JT) = NXCTR1(JT)
                  DO JC=1,NACTR
                    CTR(TRCN(JC),JT)  = CTRNX(TRCN(JC),JT)
                    CTRO(TRCN(JC),JT) = CTRNX(TRCN(JC),JT)
                  END DO
                  READ (TRC(JT),1030) NXCTR1(JT),(CTRNX(TRCN(JC),JT),  &
                                      JC=1,NACTR)
                END DO
              END IF
              NXTVD = MIN(NXTVD,NXCTR1(JT))
            END IF
          END DO
        END IF

!****** Branch related inputs

        DO JB=1,NBP

!******** Inflow

          IF (UP_FLOW(JB)) THEN
            IF (JDAY.GE.NXQIN1(JB)) THEN
              DO WHILE (JDAY.GE.NXQIN1(JB))
                NXQIN2(JB) = NXQIN1(JB)
                QIN(JB)    = QINNX(JB)
                QINO(JB)   = QINNX(JB)
                READ (INQ(JB),1020) NXQIN1(JB),QINNX(JB)
              END DO
            END IF
            NXTVD = MIN(NXTVD,NXQIN1(JB))

!****!***** Inflow temperature

            IF (JDAY.GE.NXTIN1(JB)) THEN
              DO WHILE (JDAY.GE.NXTIN1(JB))
                NXTIN2(JB) = NXTIN1(JB)
                TIN(JB)    = TINNX(JB)
                TINO(JB)   = TINNX(JB)
                READ (INT(JB),1020) NXTIN1(JB),TINNX(JB)
              END DO
            END IF
            NXTVD = MIN(NXTVD,NXTIN1(JB))

!****!***** Inflow constituent concentrations

            IF (CONSTITUENTS) THEN
              IF (JDAY.GE.NXCIN1(JB)) THEN
                DO WHILE (JDAY.GE.NXCIN1(JB))
                  NXCIN2(JB) = NXCIN1(JB)
                  DO JC=1,NACIN
                    CIN(INCN(JC),JB)  = CINNX(INCN(JC),JB)
                    CINO(INCN(JC),JB) = CINNX(INCN(JC),JB)
                  END DO
                  READ (INC(JB),1030) NXCIN1(JB),(CINNX(INCN(JC),JB),  &
                                      JC=1,NACIN)
                END DO
              END IF
              NXTVD = MIN(NXTVD,NXCIN1(JB))
            END IF
          END IF

!******** Outflow

          IF (DN_FLOW(JB)) THEN
            IF (JDAY.GE.NXQOT1(JB)) THEN
              DO WHILE (JDAY.GE.NXQOT1(JB))
                NXQOT2(JB) = NXQOT1(JB)
                IF (SEL_WITHDRAWAL(JB)) THEN
                  DO JS=1,NSTR(JB)
                    QSTR(JS,JB)  = QSTRNX(JS,JB)
                    QSTRO(JS,JB) = QSTRNX(JS,JB)
                  END DO
                  READ (OTQ(JB),1020) NXQOT1(JB),(QSTRNX(JS,JB),  &
                                      JS=1,NSTR(JB))
                ELSE
                  DO JO=1,NOUT(JB)
                    QOUT(KOUT(JO,JB),JB)  = QOUTNX(KOUT(JO,JB),JB)
                    QOUTO(KOUT(JO,JB),JB) = QOUTNX(KOUT(JO,JB),JB)
                  END DO
                  READ (OTQ(JB),1020) NXQOT1(JB),(QOUTNX(KOUT(JO,JB),  &
                                      JB),JO=1,NOUT(JB))
                END IF
              END DO
            END IF
            NXTVD = MIN(NXTVD,NXQOT1(JB))
          END IF

!******** Distributed tributaries


!******** Precipitation

          IF (PRECIPITATION) THEN
            IF (JDAY.GE.NXPR1(JB)) THEN
             DO WHILE (JDAY.GE.NXPR1(JB))
                PR(JB) = PRNX(JB)
                READ (PRE(JB),1020) NXPR1(JB),PRNX(JB)
              END DO
            END IF
            NXTVD = MIN(NXTVD,NXPR1(JB))

!****!***** Temperature

            IF (JDAY.GE.NXTPR1(JB)) THEN
              DO WHILE (JDAY.GE.NXTPR1(JB))
                TPR(JB) = TPRNX(JB)
                READ (PRT(JB),1020) NXTPR1(JB),TPRNX(JB)
              END DO
            END IF
            NXTVD = MIN(NXTVD,NXTPR1(JB))

!****!***** Constituent concentrations

            IF (CONSTITUENTS) THEN
              IF (JDAY.GE.NXCPR1(JB)) THEN
                DO WHILE (JDAY.GE.NXCPR1(JB))
                  DO JC=1,NACPR
                    CPR(PRCN(JC),JB) = CPRNX(PRCN(JC),JB)
                  END DO
                  READ (PRC(JB),1030) NXCPR1(JB),(CPRNX(PRCN(JC),JB),  &
                                      JC=1,NACPR)
                END DO
              END IF
              NXTVD = MIN(NXTVD,NXCPR1(JB))
            END IF
          END IF
!
        END DO

!****** Dead sea case

        IF (NO_WIND) THEN
          WIND   = 0.0
          WINDO  = 0.0
          WINDNX = 0.0
        END IF
        IF (NO_INFLOW) THEN
          DO JB=1,NBP
            QIN(JB)    = 0.0
            QINO(JB)   = 0.0
            QINNX(JB)  = 0.0
            QDTR(JB)   = 0.0
            QDTRO(JB)  = 0.0
            QDTRNX(JB) = 0.0
            PR(JB)     = 0.0
            PRNX(JB)   = 0.0
          END DO
          DO JT=1,NTR
            QTR(JT)   = 0.0
            QTRO(JT)  = 0.0
            QTRNX(JT) = 0.0
          END DO
        END IF
        IF (NO_OUTFLOW) THEN
          DO JB=1,NBP
            IF (SEL_WITHDRAWAL(JB)) THEN
              DO JS=1,NSTR(JB)
                QSTR(JS,JB)   = 0.0
                QSTRO(JS,JB)  = 0.0
                QSTRNX(JS,JB) = 0.0
              END DO
            ELSE
              DO K=KT,KB(ID)
                QOUT(K,JB)   = 0.0
                QOUTO(K,JB)  = 0.0
                QOUTNX(K,JB) = 0.0
              END DO
            END IF
          END DO
          DO JW=1,NWD
            QWD(JW)   = 0.0
            QWDO(JW)  = 0.0
            QWDNX(JW) = 0.0
          END DO
        END IF

!****** Input FORMATs

 1000   FORMAT(//)
 1010   FORMAT(10F8.0)
 1020   FORMAT(10F8.0:/(8X,9F8.0))
 1030   FORMAT(21F8.0)
 1040   FORMAT('Computed flows to match observed elevations'//  &
               '    JDAY     QDT')
 1050 FORMAT(F8.3,F8.1)
      RETURN
      
!****!****!****!****!****!****!****!****!****!****!****!****!****!*****
!*                         I N T E R P O L A T E  I N P U T S              **
!****!****!****!****!****!****!****!****!****!****!****!****!****!*****

      ENTRY INTERPOLATE_INPUTS (JDAY)

!****** Meteorlogical data

        IF (INTERP_MET) THEN
          RATIO = (NXMET1-JDAY)/(NXMET1-NXMET2)
          TDEW  = (1.0-RATIO)*TDEWNX+RATIO*TDEWO
          WIND  = (1.0-RATIO)*WINDNX+RATIO*WINDO
          PHI   = (1.0-RATIO)*PHINX+RATIO*PHIO
          TAIR  = (1.0-RATIO)*TAIRNX+RATIO*TAIRO
          CLOUD = (1.0-RATIO)*CLOUDNX+RATIO*CLOUDO
        END IF

!****** Withdrawals

        IF (WITHDRAWALS) THEN
          IF (INTERP_WITHDRWL) THEN
            QRATIO = (NXQWD1-JDAY)/(NXQWD1-NXQWD2)
            DO JW=1,NWD
              QWD(JW) = (1.0-QRATIO)*QWDNX(JW)+QRATIO*QWDO(JW)
            END DO
          END IF
        END IF

!****** Tributaries

        IF (TRIBUTARIES) THEN
          IF (INTERP_TRIBS) THEN
            DO JT=1,NTR
              QRATIO = (NXQTR1(JT)-JDAY)/(NXQTR1(JT)-NXQTR2(JT))
              TRATIO = (NXTTR1(JT)-JDAY)/(NXTTR1(JT)-NXTTR2(JT))
              IF (CONSTITUENTS) THEN
                CRATIO = (NXCTR1(JT)-JDAY)/(NXCTR1(JT)-NXCTR2(JT))
              END IF
              QTR(JT) = (1.0-QRATIO)*QTRNX(JT)+QRATIO*QTRO(JT)
              TTR(JT) = (1.0-TRATIO)*TTRNX(JT)+TRATIO*TTRO(JT)
              DO JC=1,NACTR
                CTR(TRCN(JC),JT) = (1.0-CRATIO)*CTRNX(TRCN(JC),JT)  &
                                   +CRATIO*CTRO(TRCN(JC),JT)
              END DO
            END DO
          END IF
        END IF

!****** Branch related inputs

        DO JB=1,NBP

!******** Inflow

          IF (UP_FLOW(JB)) THEN
            IF (INTERP_INFLOW) THEN
              QRATIO = (NXQIN1(JB)-JDAY)/(NXQIN1(JB)-NXQIN2(JB))
              TRATIO = (NXTIN1(JB)-JDAY)/(NXTIN1(JB)-NXTIN2(JB))
              IF (CONSTITUENTS) THEN
                CRATIO = (NXCIN1(JB)-JDAY)/(NXCIN1(JB)-NXCIN2(JB))
              END IF
              QIN(JB) = (1.0-QRATIO)*QINNX(JB)+QRATIO*QINO(JB)
              TIN(JB) = (1.0-TRATIO)*TINNX(JB)+TRATIO*TINO(JB)
              DO JC=1,NACIN
                CIN(INCN(JC),JB) = (1.0-CRATIO)*CINNX(INCN(JC),JB)  &
                                   +CRATIO*CINO(INCN(JC),JB)
              END DO
            END IF
          END IF

!******** Outflow

          IF (DN_FLOW(JB)) THEN
            IF (INTERP_OUTFLOW) THEN
              QRATIO = (NXQOT1(JB)-JDAY)/(NXQOT1(JB)-NXQOT2(JB))
              IF (SEL_WITHDRAWAL(JB)) THEN
                DO JS=1,NSTR(JB)
                  QSTR(JS,JB) = (1.0-QRATIO)*QSTRNX(JS,JB)  &
                                +QRATIO*QSTRO(JS,JB)
                END DO
              ELSE
                QSUM(JB) = 0.0
                DO JO=1,NOUT(JB)
                  K          = KOUT(JO,JB)
                  QOUT(K,JB) = (1.0-QRATIO)*QOUTNX(K,JB)+QRATIO  &
                               *QOUTO(K,JB)
                  QSUM(JB)   = QSUM(JB)+QOUT(K,JB)
                END DO
              END IF
            END IF
          END IF

        END DO
      RETURN
      END

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                S U B R O U T I N E   H E A T   E X C H A N G E           **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      SUBROUTINE HEAT_EXCHANGE (JDAY)
        INCLUDE 'w2.inc'

!****** Type declaration

        REAL      JDAY, LAT, LONG, LONG0
        INTEGER   GDAY
        LOGICAL   LEAP_YEAR
        CHARACTER MONTH*9

!****** Dimension declaration

        DIMENSION EQT(12)

!****** Common declaration

        COMMON /GLBLCC/ PALT, APOM, O2LIM, WIND, WSCDP, WSC(NDP)
        COMMON /TVDMTC/ TAIR, TDEW, CLOUD, PHI,  ET,    CSHE,    &  
                        SRO,  SRON, LAT,   LONG, WINDH, RHOWCP
        COMMON /GDAYC1/ GDAY, DAYM, JDAYG, LEAP_YEAR,   M
        COMMON /GDAYC2/ MONTH

!****** Data declaration

        DATA EQT /-0.13, -0.23, -0.16, -0.02, 0.06, 0.00, -0.09,  &
                  -0.08,  0.06,  0.22,  0.25, 0.10/
        DATA SBC /2.0411E-7/

!****** Month

        IF (MONTH.EQ.'  January') M = 1
        IF (MONTH.EQ.' February') M = 2
        IF (MONTH.EQ.'    March') M = 3
        IF (MONTH.EQ.'    April') M = 4
        IF (MONTH.EQ.'      May') M = 5
        IF (MONTH.EQ.'     June') M = 6
        IF (MONTH.EQ.'     July') M = 7
        IF (MONTH.EQ.'   August') M = 8
        IF (MONTH.EQ.'September') M = 9
        IF (MONTH.EQ.'  October') M = 10
        IF (MONTH.EQ.' November') M = 11
        IF (MONTH.EQ.' December') M = 12

!****** English units

        WINDT = WIND*2.23714
        TDEWT = TDEW*9.0/5.0+32.0
        TAIRT = TAIR*9.0/5.0+32.0
        WINDT = WINDT*LOG(2.0/0.003)/LOG(WINDH/0.003)

!****** Solar radiation

        LONG0 = 15.0*INT(LONG/15.0)
        D     = 0.409280*COS(0.017214*(172.0-JDAYG))
        X     = (JDAY-JDAYG)*24.0
        H     = 0.261799*(X-(LONG-LONG0)*0.066667+EQT(M)-12.0)
        SINAL = SIN(LAT*.017453)*SIN(D)+COS(LAT*.017453)*COS(D)*COS(H)
        A0    = 57.2985*ASIN(SINAL)
        SRO   = 2.044*A0+0.1296*A0**2-0.001941*A0**3+7.591E-6*A0**4
        SRO   = (1.0-0.0065*CLOUD*CLOUD)*SRO*24.0
        IF (A0.LT.0.0) SRO = 0.0

!****** Equilibrium temperature and heat exchange coefficient

        ET    = TDEWT
        FW    = 70.0+0.7*WINDT*WINDT
        HA    = 3.1872E-08*(TAIRT+459.67)**4
        TSTAR = (ET+TDEWT)*0.5
        BETA  = 0.255-(0.0085*TSTAR)+(0.000204*TSTAR**2)
        CSHE  = 15.7+(0.26+BETA)*FW
        ETP   = (SRO+HA-1801.0)/CSHE+(CSHE-15.7)*(0.26*TAIRT+BETA*TDEWT)  &
                /(CSHE*(0.26+BETA))
        J     = 0
        DO WHILE (ABS(ETP-ET).GT.0.05.AND.J.LT.50)
          ET    = ETP
          TSTAR = (ET+TDEWT)*0.5
          BETA  = 0.255-(0.0085*TSTAR)+(2.04E-4*TSTAR**2)
          CSHE  = 15.7+(0.26+BETA)*FW
          ETP   = (SRO+HA-1801.0)/CSHE+(CSHE-15.7)*(0.26*TAIRT+BETA  &
                  *TDEWT)/(CSHE*(0.26+BETA))
          J     = J+1
        END DO

!****** SI units

        ET   = (ET-32.0)*5.0/9.0
        SRO  = SRO*3.14E-8
        SRON = SRO*0.94
        CSHE = CSHE*5.65E-8
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                                  R A D I A T I O N                       **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY RADIATION (RS,RSN,RAN)

!****** Net solar radiation

        RSN = 0.94*RS

!****** Net atmospheric radiation

        T2K = 273.2+TAIR
        FAC = 1.0+0.0017*CLOUD**2
        RAN = 1000.0/3600.0*9.37E-6*SBC*T2K**6*FAC*0.97
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                             S U R F A C E   T E R M S                    **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY SURFACE_TERMS (TS,RB,RE,RC)

!****** Partial water vapor pressure of the air

        IF (TDEW.GT.0.0) THEN
          EA = EXP(2.3026*(7.5*TDEW/(TDEW+237.3)+0.6609))
        ELSE
          EA = EXP(2.3026*(9.5*TDEW/(TDEW+265.5)+0.6609))
        END IF

!****** Partial water vapor pressure at the water surface temperature

        IF (TS.GT.0.0) THEN
          ES = EXP(2.3026*(7.5*TS/(TS+237.3)+0.6609))
        ELSE
          ES = EXP(2.3026*(9.5*TS/(TS+265.5)+0.6609))
        END IF

!****** Longwave back radiation

        RB = 5.443E-8*(TS+273.2)**4

!****** Evaporation

        FW = 9.2+0.46*WIND**2
        RE = FW*(ES-EA)

!****** Conduction

        RC = 0.47*FW*(TS-TAIR)
      END

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                  S U B R O U T I N E   P R I N T   G R I D               **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      SUBROUTINE PRINT_GRID (JDAY,GDAY,MONTH,YEAR)
        INCLUDE 'w2.inc'

!****** Type declaration

        REAL      JDAY,     ICETH,    MULT
        INTEGER   SNP,      CN,       YEAR,    GDAY,    BL,  &
                  CUS,      DS
        CHARACTER CPRC*3,   HPRC*3,   LFAC*4,  CONV*7,  CUNIT1*6,  &
                  LFPR*7,   BLANK*7,  MONTH*9, CNAME2*24,  &
                  HNAME*24, TITLE*72
        LOGICAL   CONSTITUENTS, ICE, ICE_CALC, OXYGEN_DEMAND,  &
                  LIMITING_FACTOR
        LOGICAL   NEW_PAGE
        DOUBLE PRECISION Z

!****** Common declaration

        COMMON /SNPUC/  SNP
        COMMON /GLOBLC/ JB, JC, IU, ID, KT, ELKT, DLT, KB(IMP), KTI(IMP)
        COMMON /DEPTHC/ DEPTHM(KMP)
        COMMON /TEMPC/  T1(KMP,IMP),    T2(KMP,IMP)
        COMMON /HYDRC1/ U(KMP,IMP),     W(KMP,IMP),    AZ(KMP,IMP),  &
                        RHO(KMP,IMP),   NDLT(KMP,IMP)
        COMMON /HYDRC2/ Z(IMP)
        COMMON /PRNTC1/ IPR(IMP),       IEPR,          KBPR
        COMMON /PRNTC2/ TITLE(7),       CPRC(NCP),     HPRC(4),  &
                        CONV(KMP,IMP)
        COMMON /PRNTC3/ CUS(NBP),       DS(NBP)
        COMMON /LFACC/  LFAC(KMC,IMC),  LFPR(KMC,IMC)
        COMMON /DKSEDC/ SEDD(KMC,IMC),  SODD(KMC,IMC), SOD(IMP)
        COMMON /GRDLGC/ LIMITING_FACTOR,OXYGEN_DEMAND
        COMMON /GRTVDC/ CONSTITUENTS,   CN(NCP),       NAC
        COMMON /CBODC/  KBOD,           TBOD,          RBOD
        COMMON /NAMESC/ HNAME(4),       CNAME2(NCP),   CUNIT1(NCP)
        COMMON /ICEC/   ICE(IMP),       ICETH(IMP),    ICE_CALC

!****** Dimension declarations

        DIMENSION BL(IMP)

!****** Data declaration
!*kevin's mod
        DATA BLANK /'       '/

!****** Variable initialization

        NEW_PAGE = .TRUE.

!****** Inactive segments

        NBL  = 1
        JB   = 1
        IBPR = 1
        DO I=1,IEPR-1
          IF (CUS(JB).GT.IPR(I)) THEN
            BL(NBL) = I
            NBL     = NBL+1
            IF (JB.EQ.1) IBPR = I+1
          END IF
          IF (IPR(I+1).GT.DS(JB)) JB = JB+1
        END DO
        NBL = NBL-1

!****** Water surface and ice cover

        DO I=IBPR,IEPR
          WRITE (CONV(1,I),'(F7.4)') Z(IPR(I))
        END DO
        DO JBL=1,NBL
          CONV(1,BL(JBL)) = BLANK
        END DO
        WRITE (SNP,3000) 'Water Surface Deviation from Top Layer [KT=',  &
                          KT,'], m',(IPR(I),I=IBPR,IEPR)
        WRITE (SNP,3010) (CONV(1,I),I=IBPR,IEPR)
        IF (ICE_CALC) THEN
          DO I=IBPR,IEPR
            WRITE (CONV(1,I),'(F7.4)') ICETH(IPR(I))
          END DO
          DO JBL=1,NBL
            CONV(1,BL(JBL)) = BLANK
          END DO
          WRITE (SNP,3020) 'Ice Thickness, m',(CONV(1,I),I=IBPR,IEPR)
        END IF
        IF (CONSTITUENTS) THEN
          DO I=IBPR,IEPR
            WRITE (CONV(1,I),'(F7.2)') SOD(IPR(I))*86400.0
          END DO
          DO JBL=1,NBL
            CONV(1,BL(JBL)) = BLANK
          END DO
          IF (OXYGEN_DEMAND) THEN
            WRITE (SNP,3030) 'Sediment Oxygen Demand, g/m^2/day',  &
                              (CONV(1,I),I=IBPR,IEPR)
          END IF
        END IF

!****** Velocities, temperatures, and vertical timesteps

        DO J=1,4
          IF (HPRC(J).EQ.' ON') THEN
!*kevin's mod: Another horrifying kluge to initialize the 'CONV' array to include the BLANK
!*values on all cells without water. A definite candidate to clean up later when more
!*time is available!!!
          IF (J.EQ.3) THEN                                                      !9/1/98
             DO I=1,IMP                                                         !9/1/98
                DO K=1,KMP                                                      !9/1/98
                   CONV(K,I) = BLANK                                            !9/1/98
                END DO                                                          !9/1/98
             END DO                                                             !9/1/98
          END IF                                                                !9/1/98
!*end mod

            DO I=IBPR,IEPR
              DO K=KT,KB(IPR(I))
                IF (J.EQ.1) WRITE (CONV(K,I),'(F7.4)') U(K,IPR(I))
                IF (J.EQ.2) WRITE (CONV(K,I),'(F7.4)') W(K,IPR(I))*1000.
                IF (J.EQ.3) WRITE (CONV(K,I),'(F7.2)') T2(K,IPR(I))
                IF (J.EQ.4) WRITE (CONV(K,I),'(I7)')   NDLT(K,IPR(I))
              END DO
            END DO
            DO JBL=1,NBL
              DO K=KT,KB(IPR(BL(JBL)))
                CONV(K,BL(JBL)) = BLANK
              END DO
            END DO
            IF (NEW_PAGE) THEN
              WRITE (SNP,3040) (TITLE(I),I=1,7)
              NLINES = KMP-KT+14
            END IF
            NLINES   = NLINES+KMP-KT+6
            NEW_PAGE = NLINES.GT.72
            WRITE (SNP,3050)  MONTH,GDAY,YEAR,INT(JDAY),  &
                             (JDAY-INT(JDAY))*24.0,HNAME(J)
!*kevin's mod:
!           WRITE (8,3085) JDAY,(IPR(I),I=1,(IMP-5))                    !8/28/98 
            WRITE (SNP,3080) (IPR(I),I=IBPR,IEPR)
            DO K=2,KMP                                                  !8/31/98
!*kevin's mod: write temps to matlab file ('matsnap.dat')
!              WRITE (8,3090) K,DEPTHM(K),(CONV(K,I),I=1,(IMP-5))        !8/28/98            WRITE (SNP,3080) (IPR(I),I=IBPR,IEPR)
              WRITE (SNP,3090) K,DEPTHM(K),(CONV(K,I),I=IBPR,IEPR)
            END DO
          END IF
        END DO

!****** Constituent concentrations

        DO J=1,NAC
          IF (CPRC(CN(J)).EQ.' ON') THEN
            MULT = 1.0
            IF (CN(J).GE.9.AND.CN(J).LE.11) MULT = 1000.0
            DO I=IBPR,IEPR
              DO K=KT,KB(IPR(I))
                WRITE (CONV(K,I),'(F7.2)') C2(K,IPR(I),CN(J))*MULT
              END DO
            END DO
            DO JBL=1,NBL
              DO K=KT,KB(IPR(BL(JBL)))
                CONV(K,BL(JBL)) = BLANK
              END DO
            END DO
            IF (NEW_PAGE) THEN
              WRITE (SNP,3040) (TITLE(I),I=1,7)
              NLINES = KMP-KT+14
            END IF
            NLINES   = NLINES+KMP-KT+6
            NEW_PAGE = NLINES.GT.72
            WRITE (SNP,3050) MONTH,GDAY,YEAR,INT(JDAY),(JDAY-INT(JDAY))  &
                             *24.0,CNAME2(CN(J))
            WRITE (SNP,3080) (IPR(I),I=IBPR,IEPR)
            DO K=KT,KBPR
              WRITE (SNP,3090) K,DEPTHM(K),(CONV(K,I),I=IBPR,IEPR)
            END DO
          END IF
        END DO
        IF (LIMITING_FACTOR) THEN
          IF (NEW_PAGE) THEN
            WRITE (SNP,3040) (TITLE(I),I=1,7)
            NLINES = KMP-KT+14
          END IF
          NLINES   = NLINES+KMP-KT+6
          NEW_PAGE = NLINES.GT.72
          WRITE (SNP,3050)  MONTH,GDAY,YEAR,INT(JDAY),(JDAY-INT(JDAY))  &
                            *24.0,'Limiting Factor    '
          WRITE (SNP,3080) (IPR(I),I=IBPR,IEPR)
          DO K=KT,KBPR
            WRITE (SNP,3090) K,DEPTHM(K),(LFPR(K,IPR(I)),I=IBPR,IEPR)
          END DO
        END IF
      
!****** Snapshot FORMATs
      
 3000   FORMAT(/52X,A43,I2,A4/  &
               '+',51X,5('_'),1X,7('_'),1X,9('_'),1X,4('_'),1X,3('_'),  &
                 1X,5('_')//  &
               1000I7)
 3010   FORMAT(2X,1000A7/)
 3020   FORMAT(/70X,A16/  &
               '+',69X,3('_'),1X,9('_')//  &
               3X,1000A7)
 3030   FORMAT(/63X,A33/  &
               '+',62X,8('_'),1X,6('_'),1X,6('_')//  &
               3X,1000A7/)
 3040   FORMAT('1',7(A72/1X))
 3050   FORMAT(/1X,A9,I3,', ',I4,'       Julian Date =',I6,' days',F6.2,  &
                 ' hours',10X,A24/)
 3060   FORMAT(70X,A24,1X,A6/)
 3080   FORMAT(1X,'Layer  Depth',1000I7)
!*kevin's mod: new format for matlab output 
 3085   FORMAT(1X,'-999',F8.1,1000I7)                          !8/28/98
 3090   FORMAT(1X,I4,F8.2,1000A7)
      END

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*         S U B R O U T I N E   S E L E C T I V E   W I T H D R A W A L    **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      SUBROUTINE SELECTIVE_WITHDRAWAL
        INCLUDE 'w2.inc'

!****** Type declaration

        REAL    NONZERO
        LOGICAL SEL_WITHDRAWAL, POINT_SINK
        DOUBLE PRECISION Z

!****** Dimension declaration

        DIMENSION VNORM(KMP)

!****** Common declaration

        COMMON /GLOBLC/ JB, JC, IU, ID, KT, ELKT, DLT, KB(IMP), KTI(IMP)
        COMMON /NONZC/  NONZERO
        COMMON /GEOMHC/ EL(KMP),       ELWS(IMP),      ELWS2(IMP),  &
                        H(KMP),        HKT1(IMP),      HKT2(IMP),  &
                        DLX(IMP)
        COMMON /TVDQC/  QIN(NBP),      QTR(NTP),       QDTR(NBP),  &
                        PR(NBP),       ELUH(NBP),      ELDH(NBP),  &
                        QOUT(KMP,NBP), QWD(NWP),       QSUM(NBP)
        COMMON /SELWC/  NSTR(NBP),     QSTR(NSP,NBP),  ESTR(NSP,NBP),  &
                        WSTR(NSP,NBP), KBSW(NSP,NBP),  NOUT(NBP),  &
                        KOUT(KMP,NBP)
        COMMON /HYDRC1/ U(KMP,IMP),    W(KMP,IMP),     AZ(KMP,IMP),  &
                        RHO(KMP,IMP),  NDLT(KMP,IMP)
        COMMON /HYDRC2/ Z(IMP)
        COMMON /TVDSWC/ SEL_WITHDRAWAL(NBP), POINT_SINK(NSP,NBP)

!****** Data declaration

        DATA G /9.81/

!****** Variable initialization

        DO K=1,KMP
          QOUT(K,JB) = 0.0
          VNORM(K)   = 0.0
        END DO

!****** Outflows

        DO JS=1,NSTR(JB)

!******** Initial withdrawal limits

          KTOP = KT
          KBOT = MIN(KBSW(JS,JB),KB(ID))
          IF (KBOT.LT.KT+1) KBOT = KT+1

!******** Structure layer

          K = KT
          DO WHILE (EL(K).GE.ESTR(JS,JB))
            K = K+1
          END DO
          KSTR  = MAX(K-1,KT)
          KSTR  = MIN(KSTR,KB(ID))
          ELSTR = ESTR(JS,JB)
          ELKT  = EL(KT)-Z(ID)
          IF (ESTR(JS,JB).LE.EL(KB(ID)+1)) THEN
            KSTR  = KT
            ELSTR = EL(KB(ID))
          END IF
          IF (ESTR(JS,JB).GT.ELKT) ELSTR = ELKT
          IF (KBSW(JS,JB).LT.KSTR) THEN
            KSTR  = KT
            ELSTR = ELKT-(H(KT)-Z(ID))/2.0
          END IF

!******** Boundary interference

          RATIO = (ELSTR-EL(KBOT))/(ELKT-EL(KBOT))
          COEF  = 1.0
          IF (RATIO.LT.0.10.OR.RATIO.GT.0.90) COEF = 2.0

!******** Withdrawal zone above structure

          DO K=KSTR-1,KT,-1

!****!***** Density frequency

            HT    = EL(K)-ELSTR
            RHOFT = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/(HT  &
                    *RHO(KSTR,ID)+NONZERO)*G),NONZERO)

!****!***** Thickness

            IF (POINT_SINK(JS,JB)) THEN
              HSWT = (COEF*QSTR(JS,JB)/RHOFT)**0.333333
            ELSE
              HSWT = SQRT(2.0*COEF*QSTR(JS,JB)/(WSTR(JS,JB)*RHOFT))
            END IF
            IF (HT.GE.HSWT) THEN
              KTOP = K
              GO TO 10000
            END IF
          END DO
10000     CONTINUE

!******** Reference density

          IF ((ELSTR+HSWT).LT.ELKT) THEN
            DLRHOT = ABS(RHO(KSTR,ID)-RHO(KTOP,ID))
          ELSE IF (ELKT.EQ.ELSTR) THEN
            DLRHOT = NONZERO
          ELSE
            DLRHOT = ABS(RHO(KSTR,ID)-RHO(KT,ID))*HSWT/(ELKT-ELSTR)
          END IF
          DLRHOT = MAX(DLRHOT,NONZERO)

!******** Withdrawal zone below structure

          DO K=KSTR+1,KBOT

!****!***** Density frequency

            HB    = ELSTR-EL(K)
            RHOFB = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/(HB  &
                    *RHO(KSTR,ID)+NONZERO)*G),NONZERO)

!****!***** Thickness

            IF (POINT_SINK(JS,JB)) THEN
              HSWB = (COEF*QSTR(JS,JB)/RHOFB)**0.333333
            ELSE
              HSWB = SQRT(2.0*COEF*QSTR(JS,JB)/(WSTR(JS,JB)*RHOFB))
            END IF
            IF (HB.GE.HSWB) THEN
              KBOT = K
              GO TO 10010
            END IF
          END DO
10010     CONTINUE

!******** Reference density

          IF ((ELSTR-HSWB).GT.EL(KBOT+1)) THEN
            DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))
          ELSE
            DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))*HSWB  &
                     /(ELSTR-EL(KBOT+1))
          END IF
          DLRHOB = MAX(DLRHOB,NONZERO)

!******** Velocity profile

          VSUM     = 0.0
          DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)
          DO K=KTOP,KBOT
            VNORM(K) = ABS(1.0-((RHO(K,ID)-RHO(KSTR,ID))/DLRHOMAX)**2)
            VSUM     = VSUM+VNORM(K)
          END DO

!******** Outflows

          DO K=KTOP,KBOT
            QOUT(K,JB) = QOUT(K,JB)+(VNORM(K)/VSUM)*QSTR(JS,JB)
          END DO
        END DO

!****** Boundary velocities and total outflow

        QSUM(JB) = 0.0
        DO K=KT,KB(ID)
          IF (QOUT(K,JB).EQ.0.0) U(K,ID) = 0.0
          QSUM(JB) = QSUM(JB)+QOUT(K,JB)
        END DO
      END

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                    S U B R O U T I N E   D A T E   T I M E               **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

!**** This subroutine contains UNIX specific FORTRAN subroutine calls

      SUBROUTINE DATE_TIME (CDATE,CTIME)

!****** Type declaration

        CHARACTER*8  CTIME                                              !FORTRAN
!        CTIME = '12:00:00'
        CHARACTER*10  CDATE                                              !FORTRAN
!        CDATE ='06/10/2020'
!        CALL DATE (CDATE)                                               !FORTRAN
!        CALL TIME (CTIME)                                               !FORTRAN
      END


!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                      S U B R O U T I N E   K I N E T I C S               **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      SUBROUTINE KINETICS
        INCLUDE "w2.inc"

!****** Type declarations

        REAL    NH4RM, NO3RM
        REAL    NH4K1, NH4K2,  NO3K1,  NO3K2
        REAL    NH4T1, NH4T2,  NO3T1,  NO3T2
        REAL    NH4D,  NO3D,   NH4DK,  NO3DK, NH4R
        REAL    LAM1,  LAM2,   LTCOEF
        REAL    LLIM,  NLIM,   LIMIT
        REAL    LAT,   LONG
        REAL    LPOMD, LRDDK,  LDOMDK, LPOMDK, KBOD
        REAL    K1,    K2,     K3,     K4
        REAL    ICETH, KW,     NETS,   NONZERO
        LOGICAL FRESH_WATER, SALT_WATER, SUSP_SOLIDS
        LOGICAL ICE,         ICE_CALC
        CHARACTER LF*3,      LFAC*4,     LFPR*7

!****** Common declarations

        COMMON /IRONC/  FER
        COMMON /GLOBLC/ JB,     JC,     IU,     ID,     KT,    ELKT,  &
                        DLT,    KB(IMP), KTI(IMP)
        COMMON /SETLC2/ SSS,    POMS,   AS,     FES
        COMMON /ORGDKC/ SDK,    LPOMDK, LDOMDK, RDOMDK, LRDDK
        COMMON /CLFRMC/ COLQ10, COLDK
        COMMON /RTMLTC/ OMT1,   OMT2,   NH4T1,  NH4T2,  NO3T1, NO3T2,  &
                        AT1,    AT2,    AT3,    AT4,    OMK1,  OMK2,  &
                        NH4K1,  NH4K2,  NO3K1,  NO3K2,  AK1,   AK2,  &
                        AK3,    AK4
        COMMON /PHYTC1/ AE,     AM,     AG,     AR,     ASAT,  AHSN,  &
                        AHSP
        COMMON /PHYTC2/ BETA,   EXH2O,  EXSS,   EXOM
        COMMON /TVDMTC/ TAIR,   TDEW,   CLOUD,  PHI,    ET,    CSHE,  &
                        SRO,    SRON,   LAT,    LONG,   WINDH, RHOWCP
        COMMON /NITROC/ BION,   NH4DK,  NH4R,   NO3DK
        COMMON /PHOSPC/ PO4R,   BIOP,   PARTP
        COMMON /OXYGNC/ O2OM,   O2AG,   O2NH4,  O2AR
        COMMON /CBODC/  KBOD,   TBOD,   RBOD
        COMMON /CARBNC/ CO2R,   BIOC
        COMMON /GLBLCC/ PALT,   APOM,   O2LIM,  WIND,   WSCDP, WSC(NDP)
        COMMON /DKBODC/ CBODD(KMC,IMC)
        COMMON /GEOMHC/ EL(KMP),        ELWS(IMP),      ELWS2(IMP),  &
                        H(KMP),         HKT1(IMP),      HKT2(IMP),  &
                        DLX(IMP)
        COMMON /DKSEDC/ SEDD(KMC,IMC),  SODD(KMC,IMC),  SOD(IMP)
        COMMON /DKMLTC/ A1(KMC,IMC),    A2(KMC,IMC),    A3(KMC,IMC)
        COMMON /TEMPC/  T1(KMP,IMP),    T2(KMP,IMP)
        COMMON /ICEC/   ICE(IMP),       ICETH(IMP),     ICE_CALC
        COMMON /DKORGC/ LPOMD(KMC,IMC), DOMD(KMC,IMC)
        COMMON /PHYTGC/ AGR(KMC,IMC),   ARR(KMC,IMC),   AMR(KMC,IMC),  &
                        AER(KMC,IMC)
        COMMON /GEOMBC/ B(KMP,IMP),     BKT(IMP),       BH(KMP,IMP),  &
                        BHKT1(IMP),     BHKT2(IMP),     BHRKT1(IMP)
        COMMON /DKNITC/ NH4D(KMC,IMC),  NO3D(KMC,IMC)
        COMMON /SETLC1/ SETIN(KMC,IMC), SETOUT(KMC,IMC)
        COMMON /DEPTHC/ DEPTHM(KMP)
        COMMON /LFACC/  LFAC(KMC,IMC),  LFPR(KMC,IMC)
        COMMON /GBLRTC/ OMRM(KMC,IMC),  NH4RM(KMC,IMC), NO3RM(KMC,IMC),  &
                        ARMR(KMC,IMC),  ARMF(KMC,IMC)
        COMMON /PO4C1/  FPSS(KMC,IMC),  FPFE(KMC,IMC)
        COMMON /DNSPHC/ FRESH_WATER,    SALT_WATER,     SUSP_SOLIDS
        COMMON /NONZC/  NONZERO

!****** Rising and falling temperature rate functions

        FR(TT,TT1,TT2,K1,K2) = K1*EXP(LOG(K2*(1.0-K1)/(K1*(1.0-K2)))  &
                               /(TT2-TT1)*(TT-TT1))
        FF(TT,TT3,TT4,K3,K4) = K4*EXP(LOG(K3*(1.0-K4)/(K4*(1.0-K3)))  &
                               /(TT4-TT3)*(TT4-TT))

!****** Dissolved oxygen saturation

        SATO(X) = EXP(7.7117-1.31403*(LOG(X+45.93)))*PALT
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                         R A T E   M U L T I P L I E R S                  **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY RATE_MULTIPLIERS
        DO I=IU,ID
          DO K=KT,KB(I)
            LAM1       = FR(T1(K,I),NH4T1,NH4T2,NH4K1,NH4K2)
            NH4RM(K,I) = LAM1/(1.0+LAM1-NH4K1)
            LAM1       = FR(T1(K,I),NO3T1,NO3T2,NO3K1,NO3K2)
            NO3RM(K,I) = LAM1/(1.0+LAM1-NO3K1)
            LAM1       = FR(T1(K,I),OMT1,OMT2,OMK1,OMK2)
            OMRM(K,I)  = LAM1/(1.0+LAM1-OMK1)
            LAM1       = FR(T1(K,I),AT1,AT2,AK1,AK2)
            LAM2       = FF(T1(K,I),AT3,AT4,AK3,AK4)
            ARMR(K,I)  = LAM1/(1.0+LAM1-AK1)
            ARMF(K,I)  = LAM2/(1.0+LAM2-AK4)
          END DO
        END DO
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                         D E C A Y   C O N S T A N T S                    **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY DECAY_CONSTANTS
        DO I=IU,ID
          DO K=KT,KB(I)
            A1(K,I)    = (1.0+SIGN(1.0,DO(K,I)-O2LIM))*0.5
            A2(K,I)    = (1.0+SIGN(1.0,O2LIM-DO(K,I)))*0.5
            A3(K,I)    = (1.0+SIGN(1.0,DO(K,I)-1.E-10))*0.5
            DOMD(K,I)  = OMRM(K,I)*(LDOMDK*LDOM(K,I)+RDOMDK  &
                         *RDOM(K,I))*A3(K,I)
            NH4D(K,I)  = NH4DK*NH4RM(K,I)*NH4(K,I)*A1(K,I)
            NO3D(K,I)  = NO3DK*NO3RM(K,I)*NO3(K,I)*A2(K,I)
            CBODD(K,I) = KBOD*TBOD**(T1(K,I)-20.0)*A3(K,I)
            LPOMD(K,I) = LPOMDK*OMRM(K,I)*LPOM(K,I)*A3(K,I)
            SEDD(K,I)  = SDK*OMRM(K,I)*SED(K,I)*A3(K,I)
          END DO
        END DO
        DO I=IU,ID
          FPSS(KT,I)   = PARTP*SS(KT,I)/(PARTP*(SS(KT,I)+FE(KT,I))+1.0)
          FPFE(KT,I)   = PARTP*FE(KT,I)/(PARTP*(SS(KT,I)+FE(KT,I))+1.0)
          SETOUT(KT,I) = (SSS*FPSS(KT,I)+FES*FPFE(KT,I))/HKT2(I)  &
                         *A1(KT,I)
          DO K=KT+1,KB(I)
            FPSS(K,I)   = PARTP*SS(K,I)/(PARTP*(SS(K,I)+FE(K,I))+1.0)
            FPFE(K,I)   = PARTP*FE(K,I)/(PARTP*(SS(K,I)+FE(K,I))+1.0)
            SETIN(K,I)  = SETOUT(K-1,I)
            SETOUT(K,I) = (SSS*FPSS(K,I)+FES*FPFE(K,I))/H(K)*A1(K,I)
          END DO
        END DO
        DO I=IU,ID
          SODD(KT,I) = SOD(I)/BHKT2(I)*OMRM(KT,I)*(B(KTI(I),I)  &
                       -B(KT+1,I))
          DO K=KT+1,KB(I)-1
            SODD(K,I) = SOD(I)/BH(K,I)*OMRM(K,I)*(B(K,I)-B(K+1,I))
          END DO
          SODD(KB(I),I) = SOD(I)/BH(KB(I),I)*OMRM(KB(I),I)*B(KB(I),I)
        END DO
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                        S U S P E N D E D   S O L I D S                   **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY SUSPENDED_SOLIDS
        DO I=IU,ID
          SSSS(KT,I) = -SSS*SS(KT,I)/HKT2(I)
          DO K=KT+1,KB(I)
            SSSS(K,I) = SSS*(SS(K-1,I)-SS(K,I))/H(K)
          END DO
        END DO
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                           C O L I F O R M                          **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY COLIFORM
        DO I=IU,ID
          DO K=KT,KB(I)
            COLSS(K,I) = -COLDK*COLQ10**(T1(K,I)-20.0)*COL(K,I)
          END DO
        END DO
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                              L A B I L E   D O M                         **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY LABILE_DOM
        DO I=IU,ID
          DO K=KT,KB(I)
            DECAY       = OMRM(K,I)*A3(K,I)*(LDOMDK+LRDDK)*LDOM(K,I)  
            APR         = (AER(K,I)+(1.0-APOM)*AMR(K,I))*ALGAE(K,I)
            LDOMSS(K,I) = APR-DECAY
          END DO
        END DO
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                           R E F R A C T O R Y   D O M                    **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY REFRACTORY_DOM
        DO I=IU,ID
          DO K=KT,KB(I)
            RDOMSS(K,I) = OMRM(K,I)*(LRDDK*LDOM(K,I)-RDOMDK  &
                          *RDOM(K,I))*A3(K,I)
          END DO
        END DO
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                      P H Y T O P L A N K T O N                     **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY PHYTOPLANKTON
        LTCOEF = (1.0-BETA)*SRO*4.186E6/ASAT
        DO I=IU,ID

!******** Limiting factor

          GAMMA = EXH2O+EXSS*SS(KT,I)+EXOM*(ALGAE(KT,I)+LPOM(KT,I))
          LAM1  = LTCOEF
          LAM2  = LTCOEF*EXP(-GAMMA*DEPTHM(KT))
          LLIM  = 2.718282*(EXP(-LAM2)-EXP(-LAM1))/(GAMMA*HKT2(I))
          FDPO4 = 1.0-FPSS(KT,I)-FPFE(KT,I)
          PLIM  = FDPO4*PO4(KT,I)/(FDPO4*PO4(KT,I)+AHSP)
          NLIM  = (NH4(KT,I)+NO3(KT,I))/(NH4(KT,I)+NO3(KT,I)+AHSN)
          LIMIT = MIN(PLIM,NLIM,LLIM)
          IF (LIMIT.EQ.PLIM) THEN
            WRITE (LFAC(KT,I),'(F4.3)') PLIM
            LF         = ' P '
            LFPR(KT,I) = LF//LFAC(KT,I)
          ELSE IF (LIMIT.EQ.NLIM) THEN
            WRITE (LFAC(KT,I),'(F4.3)') NLIM
            LF         = ' N '
            LFPR(KT,I) = LF//LFAC(KT,I)
          ELSE IF (LIMIT.EQ.LLIM) THEN
            WRITE (LFAC(KT,I),'(F4.3)') LLIM
            LF         = ' L '
            LFPR(KT,I) = LF//LFAC(KT,I)
          END IF

!******** Sources/sinks

          ARR(KT,I)   = ARMR(KT,I)*ARMF(KT,I)*AR*A3(KT,I)
          AMR(KT,I)   = (ARMR(KT,I)+1.0-ARMF(KT,I))*AM
          AGR(KT,I)   = ARMR(KT,I)*ARMF(KT,I)*AG*LIMIT
          AGR(KT,I)   = MIN(AGR(KT,I),PO4(KT,I)/(BIOP*DLT*ALGAE(KT,I)  &
                        +NONZERO),(NH4(KT,I)+NO3(KT,I))/(BION*DLT  &
                        *ALGAE(KT,I)+NONZERO))
          AER(KT,I)   = MIN((1.0-LLIM)*AE,AGR(KT,I))
          GROWTH      = (AGR(KT,I)-ARR(KT,I)-AER(KT,I)-AMR(KT,I))  &
                        *ALGAE(KT,I)  
          NETS        = -AS*ALGAE(KT,I)/HKT2(I)
          ALGSS(KT,I) = GROWTH+NETS
          DO K=KT+1,KB(I)

!****!***** Limiting factor

            GAMMA = EXH2O+EXSS*SS(K,I)+EXOM*(ALGAE(K,I)+LPOM(K,I))
            LAM1  = LTCOEF*EXP(-GAMMA*DEPTHM(K))
            LAM2  = LTCOEF*EXP(-GAMMA*DEPTHM(K+1))
            LLIM  = 2.718282*(EXP(-LAM2)-EXP(-LAM1))/(GAMMA*H(K))
            FDPO4 = 1.0-FPSS(K,I)-FPFE(K,I)
            PLIM  = FDPO4*PO4(K,I)/(FDPO4*PO4(K,I)+AHSP)
            NLIM  = (NH4(K,I)+NO3(K,I))/(NH4(K,I)+NO3(K,I)+AHSN)
            LIMIT = MIN(PLIM,NLIM,LLIM)
            IF (LIMIT.EQ.PLIM) THEN
              WRITE (LFAC(K,I),'(F4.3)') PLIM
              LF        = ' P '
              LFPR(K,I) = LF//LFAC(K,I)
            ELSE IF (LIMIT.EQ.NLIM) THEN
              WRITE (LFAC(K,I),'(F4.3)') NLIM
              LF        = ' N '
              LFPR(K,I) = LF//LFAC(K,I)
            ELSE IF (LIMIT.EQ.LLIM) THEN
              WRITE (LFAC(K,I),'(F4.3)') LLIM
              LF        = ' L '
              LFPR(K,I) = LF//LFAC(K,I)
            END IF

!****!***** Sources/sinks

            ARR(K,I)   = ARMR(K,I)*ARMF(K,I)*AR*A3(K,I)
            AMR(K,I)   = (ARMR(K,I)+1.0-ARMF(K,I))*AM
            AGR(K,I)   = ARMR(K,I)*ARMF(K,I)*AG*LIMIT
            AGR(K,I)   = MIN(AGR(K,I),PO4(K,I)/(BIOP*DLT*ALGAE(K,I)  &
                         +NONZERO),(NH4(K,I)+NO3(K,I))/(BION*DLT  &
                         *ALGAE(K,I)+NONZERO))
            AER(K,I)   = MIN((1.0-LLIM)*AE,AGR(K,I))  
            GROWTH     = (AGR(K,I)-ARR(K,I)-AER(K,I)-AMR(K,I))  &
                         *ALGAE(K,I)
            NETS       = AS*(ALGAE(K-1,I)-ALGAE(K,I))/H(K)
            ALGSS(K,I) = GROWTH+NETS
          END DO
        END DO
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                               L A B I L E   P O M                        **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY LABILE_POM
        DO I=IU,ID
          APR          = APOM*AMR(KT,I)*ALGAE(KT,I)
          NETS         = -POMS*LPOM(KT,I)/HKT2(I)
          LPOMSS(KT,I) = APR-LPOMD(KT,I)+NETS
          DO K=KT+1,KB(I)
            APR         = APOM*AMR(K,I)*ALGAE(K,I)
            NETS        = POMS*(LPOM(K-1,I)-LPOM(K,I))/H(K)
            LPOMSS(K,I) = APR-LPOMD(K,I)+NETS
          END DO
        END DO
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                               P H O S P H O R U S                        **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY PHOSPHORUS
        DO I=IU,ID
          DO K=KT,KB(I)
            APR        = (ARR(K,I)-AGR(K,I))*ALGAE(K,I)
            PO4SS(K,I) = BIOP*(APR+LPOMD(K,I)+DOMD(K,I)+SEDD(K,I))  &
                         +PO4R*SODD(K,I)*A2(K,I)+SETIN(K,I)  &
                         *PO4(K-1,I)-SETOUT(K,I)*PO4(K,I)
          END DO
        END DO
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                                A M M O N I U M                           **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY AMMONIUM
        DO I=IU,ID
          DO K=KT,KB(I)
            APR        = (ARR(K,I)-AGR(K,I)*NH4(K,I)/(NH4(K,I)  &
                         +NO3(K,I)+NONZERO))*ALGAE(K,I)
            NH4SS(K,I) = BION*(APR+LPOMD(K,I)+DOMD(K,I)+SEDD(K,I))  &
                         +NH4R*SODD(K,I)*A2(K,I)-NH4D(K,I)
          END DO
        END DO
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                                  N I T R A T E                           **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY NITRATE
        DO I=IU,ID
          DO K=KT,KB(I)
            ALGC       = BION*(1.0-NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO))  &
                         *AGR(K,I)*ALGAE(K,I)
            NO3SS(K,I) = NH4D(K,I)-NO3D(K,I)-ALGC 
          END DO
        END DO
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                   D I S S O L V E D   O X Y G E N                  **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY DISSOLVED_OXYGEN
        O2EX = (0.5+0.05*WIND*WIND)/86400.0
        DO I=IU,ID
          DOSS(KT,I) = 0.0
          DO K=KT,KB(I)
            APR       = (O2AG*AGR(K,I)-O2AR*ARR(K,I))*ALGAE(K,I)
            DOSS(K,I) = APR-O2NH4*NH4D(K,I)-O2OM*(LPOMD(K,I)  &
                        +SEDD(K,I))-SODD(K,I)*A3(K,I)-O2OM*DOMD(K,I)  &
                        -CBODD(K,I)*CBOD(K,I)*RBOD
          END DO
          SATDO = SATO(T1(KT,I))
          IF (.NOT.ICE(I)) DOSS(KT,I) = DOSS(KT,I)+(SATDO-DO(KT,I))*O2EX  &
                                        /HKT2(I)
        END DO
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                                 S E D I M E N T                          **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY SEDIMENT
        DO I=IU,ID
          SETTLE    = (AS*ALGAE(KT,I)+POMS*LPOM(KT,I))*DLT  &
                      /HKT2(I)*(1.0-B(KT+1,I)/BKT(I))
          SED(KT,I) = MAX(SED(KT,I)+SETTLE-SEDD(KT,I)*DLT,0.0)
          DO K=KT+1,KB(I)-1
            SETTLE   = (AS*ALGAE(K,I)+POMS*LPOM(K,I))*DLT/H(K)  &
                       *(1.0-B(K+1,I)/B(K,I))
            SED(K,I) = MAX(SED(K,I)+SETTLE-SEDD(K,I)*DLT,0.0)
          END DO
          SETTLE       = (AS*ALGAE(KB(I),I)+POMS*LPOM(KB(I),I))  &
                         *DLT/H(KB(I))
          SED(KB(I),I) = MAX(SED(KB(I),I)+SETTLE-SEDD(KB(I),I)  &
                         *DLT,0.0)
        END DO
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                         I N O R G A N I C   C A R B O N                  **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY INORGANIC_CARBON
        CO2EX = (0.5+0.05*WIND*WIND)/86400.0
        DO I=IU,ID
          TICSS(KT,I) = 0.0
          DO K=KT,KB(I)
            APR        = (ARR(K,I)-AGR(K,I))*ALGAE(K,I)
            TICSS(K,I) = BIOC*(APR+DOMD(K,I)+LPOMD(K,I)+SEDD(K,I))  &
                         +CO2R*SODD(K,I)*A3(K,I)
          END DO
          IF (.NOT.ICE(I)) TICSS(KT,I) = TICSS(KT,I)+CO2EX*(0.286  &
                                         *EXP(-0.0314*(T2(KT,I))*PALT)  &
                                         -CO2(KT,I))/HKT2(I)
        END DO
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                                   P H   C O 2                            **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY PH_CO2
        DO I=IU,ID
          DO K=KT,KB(I)
            CARB = TIC(K,I)/1.2E4
            ALK  = ALKAL(K,I)/5.0E4
            T1K  = T1(K,I)+273.15

!****!***** Activity equilibrium constants

            KW = 10.0**(35.3944-5242.39/T1K-0.00835*T1K-11.8261  &
                 *LOG10(T1K))
            K1 = 10.0**(14.8435-3404.71/T1K-0.032786*T1K)
            K2 = 10.0**(6.498-2902.39/T1K-0.02379*T1K)

!****!***** Ionic strength

            IF (FRESH_WATER) THEN
              S2 = 2.5E-05*TDS(K,I)
            ELSE
              S2 = 0.00147+0.019885*TDS(K,I)+0.000038*TDS(K,I)**2
            END IF

!****!***** Debye-Huckel terms and activity coefficients

            SQRS2  = SQRT(S2)
            HCO3T  = 10.0**(-0.5085*SQRS2/(1.0+1.3124*SQRS2)+4.745694E-03  &
                     +4.160762E-02*S2-9.284843E-03*S2**2)
            CO3T   = 10.0**(-2.0340*SQRS2/(1.0+1.4765*SQRS2)+1.205665E-02  &
                     +9.715745E-02*S2-2.067746E-02*S2**2)
            H2CO3T = 0.0755*S2
            KW     = KW/HCO3T
            K1     = K1*H2CO3T/HCO3T
            K2     = K2*HCO3T/CO3T

!****!***** pH evaluation

            PH_LEFT  = -14.0
            PH_RIGHT = 0.0
            HION     = 10.0**PH_LEFT
            BICARB   = CARB*K1*HION/(K1*HION+K1*K2+HION**2)
            F_LEFT   = BICARB*(HION+2.0*K2)/HION+KW/HION-ALK-HION
            IF (F_LEFT.LT.0) THEN
             PH_START = PH_LEFT
             PH_DIFF  = PH_RIGHT-PH_LEFT
            ELSE
             PH_START = PH_RIGHT
             PH_DIFF  = PH_LEFT-PH_RIGHT
            ENDIF
            J = 0
            DO WHILE (J.LT.50)
              J       = J+1
              PH_DIFF = PH_DIFF*0.5
              PH_MID  = PH_START+PH_DIFF
              HION    = 10.0**PH_MID
              BICARB  = CARB*K1*HION/(K1*HION+K1*K2+HION**2)
              FMID    = BICARB*(HION+2.0*K2)/HION+KW/HION-ALK-HION
              IF (FMID.LT.0) PH_START = PH_MID
              IF (ABS(PH_DIFF).LT.0.01.OR.FMID.EQ.0.) J = 51
            ENDDO

!****!***** pH, carbon dioxide, bicarbonate, and carbonate concentrations

            PH(K,I)   = -PH_MID
            CO2(K,I)  = TIC(K,I)/(1.0+K1/HION+K1*K2/HION**2)
            HCO3(K,I) = TIC(K,I)/(1.0+HION/K1+K2/HION)
            CO3(K,I)  = TIC(K,I)/((HION*HION)/(K1*K2)+HION/K2+1.0)
          END DO
        END DO
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                                     I R O N                              **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY IRON
        DO I=IU,ID
          NETS       = -FES*FE(KT,I)*A1(KT,I)/HKT2(I)
          SEDR       = FER*SODD(KT,I)*A2(KT,I)
          FESS(KT,I) = SEDR+NETS
          DO K=KT+1,KB(I)
            NETS      = FES*(FE(K-1,I)*A1(K-1,I)-FE(K,I)  &
                        *A1(K,I))/H(K)
            SEDR      = FER*SODD(K,I)*A2(K,I)
            FESS(K,I) = SEDR+NETS
          END DO
        END DO
      RETURN

!****!****!****!****!****!****!****!****!****!****!****!****!****!******
!*                     B I O C H E M I C A L   O 2   D E M A N D            **
!****!****!****!****!****!****!****!****!****!****!****!****!****!******

      ENTRY BIOCHEMICAL_O2_DEMAND
        DO I=IU,ID
          DO K=KT,KB(I)
            CBODSS(K,I) = -CBODD(K,I)*CBOD(K,I)
          END DO
        END DO
      END
