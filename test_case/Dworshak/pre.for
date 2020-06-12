************************************************************************
**                                                                    **
**                         CE-QUAL-W2-PRE                             **
**                                                                    **
**                      A Preprocessor Code for                       **
**                          for CE-QUAL-W2                            **
**                                                                    **
**                          Version 2.10                              **
**                            May, 1997                               **
**                                                                    **
**                    Developed by: Thomas M. Cole                    **
**                    Water Quality Modeling Group                    **
**                    U.S. Army Corps of Engineers                    **
**                    Waterways Experiment Station                    **
**                    Vicksburg, Mississippi 39180                    **
**                    Phone number: (601) 634-3283                    **
**                                                                    **
************************************************************************


************************************************************************
**                    Task 1:  Block data setup                       **
************************************************************************

      BLOCK DATA COMMON_BLOCK
        INCLUDE 'w2.inc'
        CHARACTER*8  CUNIT
        CHARACTER*16 CNAME1
        CHARACTER*23 HNAME1, CNAME2
        CHARACTER*29 HNAME2
        CHARACTER*72 TITLE
        COMMON /NAMESC/ CNAME1(NCP), CNAME2(NCP), HNAME1(3), HNAME2(3),
     .                  CUNIT(NCP),  TITLE(6)
        COMMON /UNITNC/ INI
        DATA HNAME1 /'Horizontal velocity [U]',
     .               'Vertical velocity   [W]',
     .               'Temperature        [T2]'/
        DATA CNAME1 /'Tracer',          'Suspended solids',
     .               'Coliform',        'Dissolved solids',
     .               'Labile DOM',      'Refractory DOM',
     .               'Algae',           'Detritus',
     .               'Phosphorous',     'Ammonium',
     .               'Nitrate-Nitrite', 'Dissolved oxygen',
     .               'Sediment',        'Inorganic carbon',
     .               'Alkalinity',      'pH',
     .               'Carbon dioxide',  'Bicarbonate',
     .               'Carbonate',       'Iron',
     .               'CBOD'/
        DATA HNAME2 /' Horizontal velocity [U], m/s',
     .               ' Vertical velocity [W], mm/s ',
     .               '    Temperature [T2],øC,  '/
        DATA CNAME2 /'    Tracer, g/m^3', 
     .               'Suspended solids, g/m^3',
     .               '    Coliform, g/m^3    ',
     .               'Dissolved solids, g/m^3',
     .               '   Labile DOM, g/m^3   ',
     .               ' Refractory DOM, g/m^3 ',
     .               '      Algae, g/m^3     ',
     .               '    Detritus, g/m^3    ',
     .               '  Phosphorous, mg/m^3  ',
     .               '    Ammonium, mg/m^3   ',
     .               'Nitrate-Nitrite, mg/m^3',
     .               'Dissolved oxygen, g/m^3',
     .               '    Sediment, g/m^3    ',
     .               'Inorganic carbon, g/m^3',
     .               '   Alkalinity, g/m^3   ',
     .               '           pH          ',
     .               ' Carbon dioxide, g/m^3 ',
     .               '   Bicarbonate, g/m^3  ',
     .               '    Carbonate, g/m^3   ',
     .               '      Iron, g/m^3      ',
     .               '      CBOD, g/m^3      '/
        DATA CUNIT /8*'g/m^3', 3*'mg/m^3', 4*'g/m^3', ' ',3*'g/m^3',
     .              'mg/m^3',  'gm/m^3'/
        DATA INI   /15/
      END

************************************************************************
**                       Task 2:  Program setup                       **
************************************************************************

      PROGRAM  W2_PRE
      INCLUDE 'w2.inc'

***** Variable declarations

      INTEGER      CON,    BTH,    VPR,    WRN,    TST,    ERR,
     .             UHS,    DHS,    US,     DS,     CUS
      INTEGER      IPR,    CN,     FREQUK, YEAR
      REAL         NH4R,   NH4DK,  NH4K1,  NH4K2,  NH4T1,  NH4T2,
     .             NO3DK,  NO3T1,  NO3T2,  NO3K1,  NO3K2
      REAL         LDOMDK, LRDDK,  KBOD
      REAL         ICEMIN, ICET2,  ICETHI
      REAL         JDAY,   JDAYO
      REAL         LAT,    LONG
      CHARACTER*1  ESC
      CHARACTER*2  DEG
      CHARACTER*3  EVC,    ICEC,   PRC
      CHARACTER*3  SNPC,   VPLC,   TSRC,   PRFC,   CPLC,   RSIC,
     .             RSOC,   SPRC,   SCRC
      CHARACTER*3  CCC,    ACC,    INACC,  TRACC,  DTACC,  PRACC
      CHARACTER*3  WINDC,  QINC,   QOUTC,  HEATC
      CHARACTER*3  DTRC,   VBC,    EBC,    MBC,    WBC,    PQC,
     .             SWC
      CHARACTER*3  CPRN,   HPRN
      CHARACTER*3  SDC,    LIMC,   LJPC
      CHARACTER*3  INFIC,  TRIC,   DTRIC,  HDIC,   OUTIC,  WDIC,
     .             METIC
      CHARACTER*5  WTYPEC
      CHARACTER*7  ZERO,   BK,     BLANK,  CONV
      CHARACTER*8  AID,    SLTRC,  SLICEC, SLHTC,  PSC,    TRPC,
     .             CUNIT
      CHARACTER*16 CNAME1
      CHARACTER*23 HNAME1, CNAME2
      CHARACTER*29 HNAME2
      CHARACTER*72 TITLE,  METFN,  QOTFN,  QWDFN,  RSIFN,  QINFN,
     .             TINFN,  CINFN,  QTRFN,  TTRFN,  CTRFN,  QDTFN,
     .             TDTFN,  CDTFN,  PREFN,  TPRFN,  CPRFN,  EUHFN,
     .             TUHFN,  CUHFN,  EDHFN,  TDHFN,  CDHFN,  VPRFN,
     .             TSRFN,  PRFFN,  VPLFN,  CPLFN,  SNPFN,  LPRFN,
     .             BTHFN,  INIFN,  WRNFN,  ERRFN,  SPRFN,  ELOFN
      CHARACTER*80 ERROR
      LOGICAL      CONSTITUENTS, VERT_TEMP, LONG_TEMP,   VERT_CONC,
     .             ISO_CONC,     LONG_CONC, ISO_TEMP,    ICE_CALC,
     .             OPEN_VPR,     OPEN_LPR,  UP_FLOW,     DN_FLOW,
     .             UH_EXTERNAL,  DH_EXTERNAL
      LOGICAL      TRIBUTARIES,  WITHDRAWALS,  PRECIPITATION,
     .             DIST_TRIB,    RESTART_IN,   ISOK,
     .             DELETE_ERR,   DELETE_WRN,   VALID_FILENAME
      LOGICAL      DISSOLVED_SOLIDS, COLIFORM,         SUSPENDED_SOLIDS,
     .             LABILE_DOM,       REFRACTORY_DOM,   PHYTOPLANKTON,
     .             PARTICULATE_OM,   PHOSPHORUS,       AMMONIUM,
     .             NITRATE,          DISSOLVED_OXYGEN, SEDIMENT,
     .             TOTIC,            IRON,             CBO_DEMAND
      LOGICAL      LASERJET_II,      LASERJET_III,     LASERJET_IV
      DOUBLE PRECISION VOLB, VOLG

***** DIMENSION declarations

      DIMENSION T2(KMP,IMP),    U(KMP,IMP),    W(KMP,IMP),
     .          B(KMP,IMP)
      DIMENSION NCCLB(KMP,NBP), CVLB(KMP,NBP), ALB(KMP,NBP)
      DIMENSION ESTR(NSP,NBP),  WSTR(NSP,NBP), PSC(NSP,NBP)
      DIMENSION KBSW(NSP,NBP)
      DIMENSION Z(IMP),       DLX(IMP),    KB(IMP),     BK(IMP),
     .          XLOC(IMP)
      DIMENSION ICETH(IMP),   ICE(IMP),    KTI(IMP),    ELWS(IMP)
      DIMENSION IPR(IMP),     IPRI(IMP),   IPRF(IMP),   ISPR(IMP)
      DIMENSION SOD(IMP),     PHI0(IMP)
      DIMENSION TVP(KMP),     H(KMP),      GAL(KMP),    NCCLG(KMP),
     .          CVLG(KMP),    EL(KMP)
      DIMENSION US(NBP),      DS(NBP),     CUS(NBP),    UHS(NBP),
     .          DHS(NBP)
      DIMENSION VOLB(NBP),    SWC(NBP),    NSTR(NBP)
      DIMENSION QINFN(NBP),   TINFN(NBP),  CINFN(NBP),  CDHFN(NBP),
     .          QDTFN(NBP),   TDTFN(NBP),  CDTFN(NBP),  PREFN(NBP),
     .          TPRFN(NBP),   CPRFN(NBP),  EUHFN(NBP),  TUHFN(NBP),
     .          CUHFN(NBP),   EDHFN(NBP),  TDHFN(NBP),  QOTFN(NBP)
      DIMENSION UH_EXTERNAL(NBP), DH_EXTERNAL(NBP), UP_FLOW(NBP),
     .          DN_FLOW(NBP),     DIST_TRIB(NBP)
      DIMENSION CONV(KMP,IMP)
      DIMENSION NOUT(NBP),     KOUT(IMP,NBP)
      DIMENSION IWD(NWP),      KWD(NWP)
      DIMENSION ITR(NTP),      DTRC(NBP),   QTRFN(NTP),  TTRFN(NTP),
     .          CTRFN(NTP),    ELTRT(NTP),  ELTRB(NTP),  TRPC(NTP)      !4/03/95
      DIMENSION ACC(NCP),      C2I(NCP),    CPRN(NCP),   VERT_CONC(NCP),
     .          CVP(KMP,NCP),  CNAME1(NCP), CNAME2(NCP), CUNIT(NCP),
     .          ISO_CONC(NCP), INACC(NCP),  TRACC(NCP),  DTACC(NCP),
     .          PRACC(NCP),    CN(NCP),     LONG_CONC(NCP)
      DIMENSION SNPD(NDP), SNPF(NDP), VPLD(NDP),  VPLF(NDP), PRFD(NDP),
     .          PRFF(NDP), CPLD(NDP), CPLF(NDP),  RSOD(NDP), RSOF(NDP),
     .          TSRD(NDP), TSRF(NDP), WSCD(NDP),  WSC(NDP),  DLTD(NDP),
     .          DLTF(NDP), SPRD(NDP), SPRF(NDP),  SCRD(NDP), SCRF(NDP)
      DIMENSION DLTMAX(NDP)
      DIMENSION TITLE(6),  HPRN(4),   HNAME1(3),  HNAME2(3), ERROR(3)

***** COMMON declarations

      COMMON /GRIDC1/ T2,     U,      W,      Z,      ICETH, IPR,
     .                CN,     ICE
      COMMON /GRIDC2/ NISNP,  KEPR,   NAC
      COMMON /GRIDC3/ CONV,   CPRN,   HPRN
      COMMON /GRIDC4/ US,     DS,     CUS,    KB
      COMMON /NAMESC/ CNAME1, CNAME2, HNAME1, HNAME2, CUNIT, TITLE
      COMMON /UNITNC/ INI
      COMMON /GRDLGC/ CONSTITUENTS,   ICE_CALC

***** DATA declarations

      DATA DMO2  / 2.04E-9/,  DMCO2 /1.63E-9/
      DATA ZERO  /'     0.'/, BLANK /'       '/
      DATA CON   /10/, BTH /11/, VPR /12/, LPR /13/
      DATA ERR   /20/, WRN /21/, TST /22/
      DATA INIFN /'pre.opt'/
      DATA WRNFN /'pre.wrn'/
      DATA ERRFN /'pre.err'/

***** Variable initializations

      NB   = NBP
      KT   = 2
      NWRN = 0
      NERR = 0
      ESC  = CHAR(027)
      DEG  = CHAR(248)//'C'
      DELETE_ERR = .TRUE.
      DELETE_WRN = .TRUE.

************************************************************************
**                   Task 3:  Input and error check                   **
************************************************************************

***** Initialize input, warning, and error files

      OPEN (CON,FILE=CONFN,STATUS='OLD')
      OPEN (ERR,FILE=ERRFN,STATUS='UNKNOWN')
      OPEN (WRN,FILE=WRNFN,STATUS='UNKNOWN')

***** Title cards

      WRITE (*,*) 'Reading control file'
      WRITE (*,*) '  title cards'
      READ (CON,1005,ERR=10040) AID,(TITLE(J),J=1,6)
      IF (AID.NE.'TITLE C ')    GO TO 10050

***** Time control cards

      WRITE (*,*) '  time control cards'
      READ (CON,1006,ERR=10040) AID,TMSTRT,TMEND,YEAR
      IF (AID.NE.'TIME CON')    GO TO 10050
      READ (CON,1100,ERR=10040) AID,NDLT,DLTMIN
      IF (AID.NE.'DLT CON ')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(DLTD(J),J=1,NDLT)
      IF (AID.NE.'DLT DATE')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(DLTMAX(J),J=1,NDLT)
      IF (AID.NE.'DLT MAX ')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(DLTF(J),J=1,NDLT)
      IF (AID.NE.'DLT FRN ')    GO TO 10050

***** Grid definition cards

      WRITE (*,*) '  grid definition cards'
      READ (CON,1040,ERR=10040) AID,(US(JB),DS(JB),UHS(JB),DHS(JB),
     .                          JB=1,NB)
      READ (CON,1010,ERR=10040) AID,LAT,LONG,DATUM
      IF (AID.NE.'LOCATION')    GO TO 10050

C     ***** Initial condition cards

      WRITE (*,*) '  initial condition cards'
      READ (CON,1020,ERR=10040) AID,T2I,ICETHI,WTYPEC
      IF (AID.NE.'INIT CND')    GO TO 10050
      READ (CON,1030,ERR=10040) AID,VBC,EBC,MBC,WBC,PQC,EVC,PRC
      IF (AID.NE.'CALCULAT')    GO TO 10050
      READ (CON,1030,ERR=10040) AID,INFIC,TRIC,DTRIC,HDIC,OUTIC,
     .                          WDIC,METIC
      IF (AID.NE.'INTERPOL')    GO TO 10050
      READ (CON,1030,ERR=10040) AID,WINDC,QINC,QOUTC,HEATC
      IF (AID.NE.'DEAD SEA')    GO TO 10050
      READ (CON,1035,ERR=10040) AID,ICEC,SLICEC,SLHTC,ALBEDO,HWI,BETAI,
     .                          GAMMAI,ICEMIN,ICET2
      IF (AID.NE.'ICE COVE')    GO TO 10050
      READ (CON,1110,ERR=10040) AID,SLTRC,THETA
      IF (AID.NE.'TRANSPOR')    GO TO 10050
      READ (CON,1060,ERR=10040) AID,NWSD
      IF (AID.NE.'WSC NUMB')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(WSCD(J),J=1,NWSD)
      IF (AID.NE.'WSC DATE')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(WSC(J),J=1,NWSD)
      IF (AID.NE.'WSC COEF')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,AX,DXI,CHEZY,CBHE,TSED
      IF (AID.NE.'HYD COEF')    GO TO 10050

***** Inflow-outflow cards

      WRITE (*,*) '  inflow/outflow cards'
      READ (CON,1030,ERR=10040) AID,(SWC(JB),JB=1,NB)
      IF (AID.NE.'SEL WITH')    GO TO 10050
      READ (CON,1060,ERR=10040) AID,(NSTR(JB),JB=1,NB)
      IF (AID.NE.'N STRUC ')    GO TO 10050
      READ (CON,1060,ERR=10040)   AID,(KBSW(JS,1),JS=1,NSTR(1))
      IF (AID.NE.'K BOTTOM')      GO TO 10050
      DO JB=2,NB
        READ (CON,1065,ERR=10040) (KBSW(JS,JB),JS=1,NSTR(JB))
      END DO
      READ (CON,1115,ERR=10040)   AID,(PSC(JS,1),JS=1,NSTR(1))
      IF (AID.NE.'SINK TYP')      GO TO 10050
      DO JB=2,NB
        READ (CON,1116,ERR=10040) (PSC(JS,JB),JS=1,NSTR(JB))
      END DO
      READ (CON,1010,ERR=10040)   AID,(ESTR(JS,1),JS=1,NSTR(1))
      IF (AID.NE.'E STRUC ')      GO TO 10050
      DO JB=2,NB
        READ (CON,1011,ERR=10040) (ESTR(JS,JB),JS=1,NSTR(JB))
      END DO
      READ (CON,1010,ERR=10040)   AID,(WSTR(JS,1),JS=1,NSTR(1))
      IF (AID.NE.'W STRUC ')      GO TO 10050
      DO JB=2,NB
        READ (CON,1011,ERR=10040) (WSTR(JS,JB),JS=1,NSTR(JB))
      END DO
      READ (CON,1060,ERR=10040) AID,(NOUT(JB),JB=1,NB)
      IF (AID.NE.'N OUTLET')    GO TO 10050
      READ (CON,1060,ERR=10040) AID,(KOUT(JO,1),JO=1,NOUT(1))
      IF (AID.NE.'O LAYER ')    GO TO 10050
      DO JB=2,NB
        READ (CON,1065) (KOUT(JO,JB),JO=1,NOUT(JB))
      END DO
      READ (CON,1060,ERR=10040) AID,NWD
      IF (AID.NE.'N WDRWAL')    GO TO 10050
      READ (CON,1060,ERR=10040) AID,(IWD(JWD),JWD=1,NWD)
      IF (AID.NE.'W SEGMNT')    GO TO 10050
      READ (CON,1060,ERR=10040) AID,(KWD(JWD),JWD=1,NWD)
      IF (AID.NE.'W LAYER ')    GO TO 10050
      READ (CON,1060,ERR=10040) AID,NTR
      IF (AID.NE.'N TRIBS ')    GO TO 10050
      READ (CON,1063,ERR=10040) AID,(TRPC(JT),JT=1,NTR)                 !4/03/95
      IF (AID.NE.'TRIB PLA')    GO TO 10050                             !4/03/95
      READ (CON,1060,ERR=10040) AID,(ITR(JT),JT=1,NTR)
      IF (AID.NE.'TRIB SEG')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(ELTRT(JT),JT=1,NTR)                !4/03/95
      IF (AID.NE.'TRIB TOP')    GO TO 10050                             !4/03/95
      READ (CON,1010,ERR=10040) AID,(ELTRB(JT),JT=1,NTR)                !4/03/95
      IF (AID.NE.'TRIB BOT')    GO TO 10050                             !4/03/95
      READ (CON,1030,ERR=10040) AID,(DTRC(JB),JB=1,NB)
      IF (AID.NE.'DST TRIB')    GO TO 10050

***** Output control cards (excluding constituents)

      WRITE (*,*) '  output control cards'
      READ (CON,1080,ERR=10040) AID,SCRC,NSCR
      IF (AID.NE.'SCR PRIN')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(SCRD(J),J=1,NSCR)
      IF (AID.NE.'SCR DATE')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(SCRF(J),J=1,NSCR)
      IF (AID.NE.'SCR FREQ')    GO TO 10050
      READ (CON,1070,ERR=10040) AID,LJPC,(HPRN(J),J=1,4)
      IF (AID.NE.'SNAPSHOT')    GO TO 10050
      READ (CON,1090,ERR=10040) AID,SNPC,NSNP,NISNP
      IF (AID.NE.'SNP PRIN')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(SNPD(J),J=1,NSNP)
      IF (AID.NE.'SNP DATE')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(SNPF(J),J=1,NSNP)
      IF (AID.NE.'SNP FREQ')    GO TO 10050
      READ (CON,1060,ERR=10040) AID,(IPRI(J),J=1,NISNP)
      IF (AID.NE.'SNP SEG ')    GO TO 10050
      READ (CON,1090,ERR=10040) AID,PRFC,NPRF,NIPRF
      IF (AID.NE.'PRF PLOT')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(PRFD(J),J=1,NPRF)
      IF (AID.NE.'PRF DATE')    GO TO 10050
      READ (CON,1010) AID,(PRFF(J),J=1,NPRF)
      IF (AID.NE.'PRF FREQ')    GO TO 10050
      READ (CON,1060,ERR=10040) AID,(IPRF(J),J=1,NIPRF)
      IF (AID.NE.'PRF SEG ')    GO TO 10050
      READ (CON,1090,ERR=10040) AID,SPRC,NSPR,NISPR
      IF (AID.NE.'SPR PLOT')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(SPRD(J),J=1,NSPR)
      IF (AID.NE.'SPR DATE')    GO TO 10050
      READ (CON,1010) AID,(SPRF(J),J=1,NSPR)
      IF (AID.NE.'SPR FREQ')    GO TO 10050
      READ (CON,1060,ERR=10040) AID,(ISPR(J),J=1,NISPR)
      IF (AID.NE.'SPR SEG ')    GO TO 10050
      READ (CON,1080,ERR=10040) AID,TSRC,NTSR
      IF (AID.NE.'TSR PLOT')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(TSRD(J),J=1,NTSR)
      IF (AID.NE.'TSR DATE')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(TSRF(J),J=1,NTSR)
      IF (AID.NE.'TSR FREQ')    GO TO 10050
      READ (CON,1080,ERR=10040) AID,VPLC,NVPL
      IF (AID.NE.'VPL PLOT')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(VPLD(J),J=1,NVPL)
      IF (AID.NE.'VPL DATE')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(VPLF(J),J=1,NVPL)
      IF (AID.NE.'VPL FREQ')    GO TO 10050
      READ (CON,1080,ERR=10040) AID,CPLC,NCPL
      IF (AID.NE.'CPL PLOT')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(CPLD(J),J=1,NCPL)
      IF (AID.NE.'CPL DATE')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(CPLF(J),J=1,NCPL)
      IF (AID.NE.'CPL FREQ')    GO TO 10050
      READ (CON,1080,ERR=10040) AID,RSOC,NRSO,RSIC
      IF (AID.NE.'RESTART ')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(RSOD(J),J=1,NRSO)
      IF (AID.NE.'RSO DATE')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(RSOF(J),J=1,NRSO)
      IF (AID.NE.'RSO FREQ')    GO TO 10050

***** Constituent control cards

      WRITE (*,*) '  constituent control cards'
      READ (CON,1085,ERR=10040) AID,CCC,LIMC,SDC,FREQUK
      IF (AID.NE.'CST COMP')    GO TO 10050
      READ (CON,1030,ERR=10040) AID,(ACC(JC),JC=1,NCP)
      IF (AID.NE.'CST ACT ')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(C2I(JC),JC=1,NCP)
      IF (AID.NE.'CST ICON')    GO TO 10050
      READ (CON,1030,ERR=10040) AID,(CPRN(JC),JC=1,NCP)
      IF (AID.NE.'CST PRIN')    GO TO 10050
      READ (CON,1030,ERR=10040) AID,(INACC(JC),JC=1,NCP)
      IF (AID.NE.'CIN CON ')    GO TO 10050
      READ (CON,1030,ERR=10040) AID,(TRACC(JC),JC=1,NCP)
      IF (AID.NE.'CTR CON ')    GO TO 10050
      READ (CON,1030,ERR=10040) AID,(DTACC(JC),JC=1,NCP)
      IF (AID.NE.'CDT CON ')    GO TO 10050
      READ (CON,1030,ERR=10040) AID,(PRACC(JC),JC=1,NCP)
      IF (AID.NE.'CPR CON ')    GO TO 10050

***** Kinetics coefficients

      WRITE (*,*) '  kinetic cards'
      READ (CON,1010,ERR=10040) AID,EXH2O,EXSS,EXOM,BETA
      IF (AID.NE.'EX COEF ')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,COLQ10,COLDK
      IF (AID.NE.'COLIFORM')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,SSS
      IF (AID.NE.'S SOLIDS')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,AG,AM,AE,AR,AS,ASAT,APOM
      IF (AID.NE.'ALGAE   ')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,AT1,AT2,AT3,AT4,AK1,AK2,AK3,AK4
      IF (AID.NE.'ALG RATE')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,LDOMDK,LRDDK,RDOMDK
      IF (AID.NE.'DOM     ')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,POMDK,POMS
      IF (AID.NE.'POM     ')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,OMT1,OMT2,OMK1,OMK2
      IF (AID.NE.'OM RATE ')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,SDK,FSOD
      IF (AID.NE.'SEDIMENT')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,(SOD(I),I=1,IMP)
      IF (AID.NE.'S DEMAND')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,KBOD,TBOD,RBOD
      IF (AID.NE.'CBOD    ')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,PO4R,PARTP,AHSP
      IF (AID.NE.'PHOSPHOR')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,NH4R,NH4DK,PARTN,AHSN
      IF (AID.NE.'AMMONIUM')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,NH4T1,NH4T2,NH4K1,NH4K2
      IF (AID.NE.'NH4 RATE')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,NO3DK
      IF (AID.NE.'NITRATE ')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,NO3T1,NO3T2,NO3K1,NO3K2
      IF (AID.NE.'NO3 RATE')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,CO2R
      IF (AID.NE.'SED CO2 ')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,FER,FES
      IF (AID.NE.'IRON    ')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,O2NH4,O2OM,O2AR,O2AG,BIOP,BION,
     .                          BIOC
      IF (AID.NE.'STOICHMT')    GO TO 10050
      READ (CON,1010,ERR=10040) AID,O2LIM
      IF (AID.NE.'O2 LIMIT')    GO TO 10050

***** Input filenames

      WRITE (*,*) '  input filenames'
      READ (CON,1000,ERR=10040) AID,BTHFN
      IF (AID.NE.'BTH FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,VPRFN
      IF (AID.NE.'VPR FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,LPRFN
      IF (AID.NE.'LPR FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,RSIFN
      IF (AID.NE.'RSI FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,METFN
      IF (AID.NE.'MET FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,QWDFN
      IF (AID.NE.'QWD FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,ELOFN
      IF (AID.NE.'ELO FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(QINFN(JB),JB=1,NB)
      IF (AID.NE.'QIN FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(TINFN(JB),JB=1,NB)
      IF (AID.NE.'TIN FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(CINFN(JB),JB=1,NB)
      IF (AID.NE.'CIN FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(QOTFN(JB),JB=1,NB)
      IF (AID.NE.'QOT FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(QTRFN(JT),JT=1,NTP)
      IF (AID.NE.'QTR FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(TTRFN(JT),JT=1,NTP)
      IF (AID.NE.'TTR FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(CTRFN(JT),JT=1,NTP)
      IF (AID.NE.'CTR FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(QDTFN(JB),JB=1,NB)
      IF (AID.NE.'QDT FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(TDTFN(JB),JB=1,NB)
      IF (AID.NE.'TDT FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(CDTFN(JB),JB=1,NB)
      IF (AID.NE.'CDT FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(PREFN(JB),JB=1,NB)
      IF (AID.NE.'PRE FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(TPRFN(JB),JB=1,NB)
      IF (AID.NE.'TPR FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(CPRFN(JB),JB=1,NB)
      IF (AID.NE.'CPR FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(EUHFN(JB),JB=1,NB)
      IF (AID.NE.'EUH FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(TUHFN(JB),JB=1,NB)
      IF (AID.NE.'TUH FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(CUHFN(JB),JB=1,NB)
      IF (AID.NE.'CUH FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(EDHFN(JB),JB=1,NB)
      IF (AID.NE.'EDH FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(TDHFN(JB),JB=1,NB)
      IF (AID.NE.'TDH FILE')    GO TO 10050
      READ (CON,1000,ERR=10040) AID,(CDHFN(JB),JB=1,NB)
      IF (AID.NE.'CDH FILE')    GO TO 10050

***** Output filenames

      WRITE (*,*) '  output filenames'
      READ (CON,1000,ERR=10040)  AID,SNPFN
      IF (AID.NE.'SNP FILE')     GO TO 10050
      READ (CON,1000,ERR=10040)  AID,TSRFN
      IF (AID.NE.'TSR FILE')     GO TO 10050
      READ (CON,1000,ERR=10040)  AID,PRFFN
      IF (AID.NE.'PRF FILE')     GO TO 10050
      READ (CON,1000,ERR=10040)  AID,VPLFN
      IF (AID.NE.'VPL FILE')     GO TO 10050
      READ (CON,1000,ERR=10040)  AID,CPLFN
      IF (AID.NE.'CPL FILE')     GO TO 10050
      READ (CON,1000,ERR=10040)  AID,SPRFN
      IF (AID.NE.'SPR FILE')     GO TO 10050
      CLOSE (CON)

***** Bathymetry definition

      WRITE (*,*) 'Reading bathymetry file'
      IF (.NOT.VALID_FILENAME(BTHFN)) THEN
        DELETE_WRN = .FALSE.
        WRITE(WRN,9165) BTHFN, CONFN
        NWRN = NWRN+1
      END IF
      OPEN (BTH,FILE=BTHFN,STATUS='OLD',IOSTAT=IERR)
      IF (IERR.EQ.0) THEN
        READ (BTH,*)
        READ (BTH,1050) (DLX(I),I=1,IMP)
        READ (BTH,1050) (ELWS(I),I=1,IMP)
        READ (BTH,1050) (PHI0(I),I=1,IMP)
        READ (BTH,1050) (H(K),K=1,KMP)
        DO I=1,IMP
          READ (BTH,1050) (B(K,I),K=1,KMP)
        END DO
        CLOSE (BTH)
      ELSE
        DELETE_ERR = .FALSE.
        WRITE(ERR,9160) BTHFN, 'BTHFN', CONFN
        NERR = NERR+1
        CLOSE (BTH)
      END IF

***** Initialize logical control variables

      WITHDRAWALS   = NWD.GT.0
      TRIBUTARIES   = NTR.GT.0
      PRECIPITATION = PRC.EQ.' ON'
      CONSTITUENTS  = CCC.EQ.' ON'
      ISO_TEMP      = T2I.GE.0.
      VERT_TEMP     = T2I.EQ.-1.0
      LONG_TEMP     = T2I.LT.-1.0
      ICE_CALC      = ICEC.EQ.' ON'
      OPEN_VPR      = VERT_TEMP
      OPEN_LPR      = LONG_TEMP
      RESTART_IN    = RSIC.EQ.' ON'
      LASERJET_II   = LJPC.EQ.' II'
      LASERJET_III  = LJPC.EQ.'III'
      LASERJET_IV   = LJPC.EQ.' IV'
      DO JB=1,NBP
        UP_FLOW(JB)     = UHS(JB).EQ.0
        DN_FLOW(JB)     = DHS(JB).EQ.0
        UH_EXTERNAL(JB) = UHS(JB).LT.0
        DH_EXTERNAL(JB) = DHS(JB).LT.0
        DIST_TRIB(JB)   = DTRC(JB).EQ.' ON'
      END DO
      DO JC=1,NCP
        ISO_CONC(JC)  = C2I(JC).GE.0.
        VERT_CONC(JC) = C2I(JC).EQ.-1.0
        LONG_CONC(JC) = C2I(JC).LT.-1.0
        IF (VERT_CONC(JC)) OPEN_VPR = .TRUE.
        IF (LONG_CONC(JC)) OPEN_LPR = .TRUE.
      END DO
      IF (ACC(2).EQ.' ON') THEN
        SUSPENDED_SOLIDS = .TRUE.
      END IF
      IF (ACC(3).EQ.' ON') THEN
        COLIFORM = .TRUE.
      END IF
      IF (ACC(4).EQ.' ON') THEN
        DISSOLVED_SOLIDS = .TRUE.
      END IF
      IF (ACC(5).EQ.' ON') THEN
        LABILE_DOM = .TRUE.
      END IF
      IF (ACC(6).EQ.' ON') THEN
        REFRACTORY_DOM = .TRUE.
      END IF
      IF (ACC(7).EQ.' ON') THEN
        PHYTOPLANKTON = .TRUE.
      END IF
      IF (ACC(8).EQ.' ON') THEN
        PARTICULATE_OM = .TRUE.
      END IF
      IF (ACC(9).EQ.' ON') THEN
        PHOSPHORUS = .TRUE.
      END IF
      IF (ACC(10).EQ.' ON') THEN
        AMMONIUM = .TRUE.
      END IF
      IF (ACC(11).EQ.' ON') THEN
        NITRATE = .TRUE.
      END IF
      IF (ACC(12).EQ.' ON') THEN
        DISSOLVED_OXYGEN = .TRUE.
      END IF
      IF (ACC(13).EQ.' ON') THEN
        SEDIMENT = .TRUE.
      END IF
      IF (ACC(14).EQ.' ON') THEN
        TOTIC = .TRUE.
      END IF
      IF (ACC(20).EQ.' ON') THEN
        IRON = .TRUE.
      END IF
      IF (ACC(21).EQ.' ON') THEN
        CBO_DEMAND = .TRUE.
      END IF
      GO TO 10060

***** Error message for input cards

10040 CONTINUE
        DELETE_ERR = .FALSE.
        BACKSPACE (CON)
        BACKSPACE (CON)
        BACKSPACE (CON)
        READ  (CON,9020) ERROR(1)
        READ  (CON,9020) ERROR(2)
        READ  (CON,9020) ERROR(3)
        WRITE (ERR,9000) ERROR
        WRITE (*,*)      '    ERROR in control file - cannot continue ',
     .                   'processing'
        WRITE (*,*)      '    refer to "pre.err" for more information'
        STOP
10050 CONTINUE
        DELETE_ERR = .FALSE.
        BACKSPACE (CON)
        BACKSPACE (CON)
        BACKSPACE (CON)
        READ  (CON,9020) ERROR(1)
        READ  (CON,9020) ERROR(2)
        READ  (CON,9020) ERROR(3)
        WRITE (ERR,9010) ERROR
        WRITE (*,*)      '    ERROR in control file - cannot continue ',
     .                   'processing'
        WRITE (*,*)      '    refer to "pre.err" for more information'
        STOP
10060 CONTINUE

************************************************************************
**                      Input/Output File Checks                      **
************************************************************************

***** Active constituents

      IF (CONSTITUENTS) THEN
        NAC   = 0
        NACIN = 0
        NACTR = 0
        NACDT = 0
        NACPR = 0
        DO JC=1,NCP
          IF (ACC(JC).EQ.' ON') THEN
            NAC     = NAC+1
            CN(NAC) = JC
          END IF
          IF (INACC(JC).EQ.' ON') THEN
            NACIN = NACIN+1
          END IF
          IF (TRACC(JC).EQ.' ON') THEN
            NACTR = NACTR+1
          END IF
          IF (DTACC(JC).EQ.' ON') THEN
            NACDT = NACDT+1
          END IF
          IF (PRACC(JC).EQ.' ON') THEN
            NACPR = NACPR+1
          END IF
        END DO
      END IF

***** Boundary cells

      WRITE (*,*) 'Checking bathymetry file'
      DO JB=1,NBP
        DO K=1,KMP
          IF (B(K,US(JB)-1).NE.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9070) JB
            NERR = NERR+1
          END IF
          IF (B(K,DS(JB)+1).NE.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9060) K,DS(JB)+1,B(K,DS(JB)+1),JB
            NERR = NERR+1
          END IF
        END DO
        DO I=1,IMP
          IF (B(1,I).NE.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9060) 1,I,B(1,I),JB
            NERR = NERR+1
          ELSE IF (B(KMP,I).NE.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9060) KMP,I,B(KMP,I),JB
            NERR = NERR+1
          END IF
        END DO
      END DO
      DO I=1,IMP-1
        DO K=2,KMP-1
          IF (B(K+1,I).NE.0.0) THEN
            IF (B(K+1,I).GT.B(K,I)) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9061) K+1,I,K,I
              NERR = NERR+1
            END IF
          END IF
        END DO
      END DO

***** Restart input file

      WRITE (*,*) 'Checking time-varying input files'
      IF (RESTART_IN) THEN
        WRITE (*,*) '  restart'
        IF (.NOT.VALID_FILENAME(RSIFN)) THEN
          DELETE_WRN = .FALSE.
          WRITE(WRN,9165) RSIFN,CONFN
          NWRN = NWRN+1
        END IF
        OPEN (TST,FILE=RSIFN,STATUS='OLD',IOSTAT=IERR)
        IF (IERR.NE.0) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9160) RSIFN,'RSIFN',CONFN
          NERR = NERR+1
        ELSE 
          CLOSE(TST)
        END IF
      END IF

***** Meteorologic data

      WRITE (*,*) '  meteorology'
      IF (.NOT.VALID_FILENAME(METFN)) THEN
        DELETE_WRN = .FALSE.
        WRITE(WRN,9165) METFN,CONFN
        NWRN = NWRN+1
      END IF
      OPEN (UNIT=TST,FILE=METFN,STATUS='OLD',IOSTAT=IERR)
      IF (IERR.EQ.0) THEN
        READ(TST,9500)
        J     = 1
        JDAYO = 0.0
        DO WHILE (.TRUE.)
          READ(TST,9520,END=9305) JDAY,TAIR,TDEW,WIND,DUMMY,CLOUD
          IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9090) METFN,JDAY,TMSTRT
            NERR = NERR+1
          ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9095) METFN,JDAY,JDAYO
            NERR = NERR+1
          END IF
          IF (TAIR.GT.50.0) THEN 
            DELETE_WRN = .FALSE.
            WRITE(WRN,9100) 'TAIR',METFN,JDAY
            NWRN = NWRN+1
          END IF
          IF (TDEW.GT.50.0) THEN 
            DELETE_WRN = .FALSE.
            WRITE(WRN,9100) 'TDEW',METFN,JDAY
            NWRN = NWRN+1
          END IF
          IF (WIND.GT.20.0) THEN 
            DELETE_WRN = .FALSE.
            WRITE(WRN,9100) 'WIND',METFN,JDAY
            NWRN = NWRN+1
          END IF
          IF (CLOUD.GT.10.0.OR.CLOUD.LT.0.0) THEN 
            DELETE_ERR = .FALSE.
            WRITE(ERR,9100) 'CLOUD',METFN,JDAY
            NERR = NERR+1
          END IF
          IF (TAIR.EQ.0.0.AND.TDEW.EQ.0.0.AND.ABS(TAIR-TAIRO)
     .        .GT.5.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9101) JDAY
            NWRN = NWRN+1
          END IF
          J     = J+1
          JDAYO = JDAY
          TAIRO = TAIR
        END DO
 9305   CONTINUE
        IF (JDAY.LT.TMEND) THEN 
          DELETE_ERR = .FALSE.
          WRITE(ERR,9096) METFN,JDAY,TMEND
          NERR = NERR+1
        END IF
        CLOSE(TST)
      ELSE 
        DELETE_ERR = .FALSE.
        WRITE(ERR,9160) METFN,'METFN',CONFN
        NERR = NERR+1
      END IF

***** Withdrawals

      IF (WITHDRAWALS) THEN
        WRITE (*,*) '  withdrawals'
        IF (.NOT.VALID_FILENAME(QWDFN)) THEN
          DELETE_WRN = .FALSE.
          WRITE(WRN,9165) QWDFN,CONFN
          NWRN = NWRN+1
        END IF
        OPEN (UNIT=TST,FILE=QWDFN,STATUS='OLD',IOSTAT=IERR)
        IF (IERR.EQ.0) THEN
          READ(TST,9500)
          J     = 1
          JDAYO = 0.0
          DO WHILE (.TRUE.)
            READ(TST,9520,END=9308) JDAY,(DUMMY,JW=1,NWD)
            IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9090) QWDFN,JDAY,TMSTRT
              NERR = NERR+1
            ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9095) QWDFN,JDAY
              NERR = NERR+1
            END IF
            J     = J+1
            JDAYO = JDAY
          END DO
 9308     CONTINUE
          IF (JDAY.LT.TMEND) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9096) QWDFN,JDAY,TMEND
            NERR = NERR+1
          END IF
          CLOSE(TST)
        ELSE
          DELETE_ERR = .FALSE.
          WRITE(ERR,9160) QWDFN,'QWDFN',CONFN
          NERR = NERR+1
        END IF
      END IF

***** Time varying data

      DO JB=1,NBP
        IF (UP_FLOW(JB)) THEN

********* Inflows

          WRITE (*,*) '  inflows for branch',JB
          IF (.NOT.VALID_FILENAME(QINFN(JB))) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9165) QINFN(JB),CONFN
            NWRN = NWRN+1
          END IF
          OPEN (UNIT=TST,FILE=QINFN(JB),STATUS='OLD',IOSTAT=IERR)
          IF (IERR.EQ.0) THEN
            READ(TST,9500)
            J     = 1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
              READ(TST,9520,END=9315) JDAY
              IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9090) QINFN(JB),JDAY,TMSTRT
                NERR = NERR+1
              ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9095) QINFN(JB),JDAY,JDAYO
                NERR = NERR+1
              END IF
              J     = J+1
              JDAYO = JDAY
            END DO
 9315       CONTINUE
            IF (JDAY.LT.TMEND) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9096) QINFN(JB),JDAY,TMEND
              NERR = NERR+1
            END IF
            CLOSE(TST)
          ELSE
            DELETE_ERR = .FALSE.
            WRITE(ERR,9160) QINFN(JB),'QINFN',CONFN
            NERR = NERR+1
          END IF

********* Inflow temperatures

          WRITE (*,*) '  inflow temperatures'
          IF (.NOT.VALID_FILENAME(TINFN(JB))) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9165) TINFN(JB),CONFN
            NWRN = NWRN+1
          END IF
          OPEN (UNIT=TST,FILE=TINFN(JB),STATUS='OLD',IOSTAT=IERR)
          IF (IERR.EQ.0) THEN
            READ(TST,9500)
            J     = 1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
              READ(TST,9520,END=9325) JDAY,DUMMY
              IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9090) TINFN(JB),JDAY,TMSTRT
                NERR = NERR+1
              ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9095) TINFN(JB),JDAY,JDAYO
                NERR = NERR+1
              END IF
              J     = J+1
              JDAYO = JDAY
            END DO
 9325       CONTINUE
            IF (JDAY.LT.TMEND) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9096) TINFN(JB),JDAY,TMEND
              NERR = NERR+1
            END IF
            CLOSE(TST)
          ELSE
            DELETE_ERR = .FALSE.
            WRITE(ERR,9160) TINFN(JB),'TINFN',CONFN
            NERR = NERR+1
          END IF

********* Inflow concentrations

          IF (CONSTITUENTS) THEN
            WRITE (*,*) '  inflow concentrations'
            IF (.NOT.VALID_FILENAME(CINFN(JB))) THEN
              DELETE_WRN = .FALSE.
              WRITE(WRN,9165) CINFN(JB),CONFN
              NWRN = NWRN+1
            END IF
            OPEN (UNIT=TST,FILE=CINFN(JB),IOSTAT=IERR,STATUS='OLD')
            IF (IERR.EQ.0) THEN
              READ(TST,9500)
              J     = 1
              JDAYO = 0.0
              DO WHILE (.TRUE.)
                READ(TST,*,END=9335) JDAY,(DUMMY,JC=1,NACIN)
                IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                  DELETE_ERR = .FALSE.
                  WRITE(ERR,9090) CINFN(JB),JDAY,TMSTRT
                  NERR = NERR+1
                ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                  DELETE_ERR = .FALSE.
                  WRITE(ERR,9095) CINFN(JB),JDAY,JDAYO
                  NERR = NERR+1
                END IF
                J     = J+1
                JDAYO = JDAY
              END DO
 9335         CONTINUE
              IF (JDAY.LT.TMEND) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9096) CINFN(JB),JDAY,TMEND
                NERR = NERR+1
              END IF
              CLOSE(TST)
            ELSE
              DELETE_WRN = .FALSE.
              WRITE(WRN,9110) 'UP_FLOW',CINFN(JB)
              NWRN = NWRN+1
            END IF
          END IF
        END IF

******* Distributed tributaries
        
        IF (DIST_TRIB(JB)) THEN

********* Inflows

          WRITE (*,*) '  distributed tributary'
          WRITE (*,*) '    inflows'
          IF (.NOT.VALID_FILENAME(QDTFN(JB))) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9165) QDTFN(JB),CONFN
            NWRN = NWRN+1
          END IF
          OPEN (UNIT=TST,FILE=QDTFN(JB),STATUS='OLD',IOSTAT=IERR)
          IF (IERR.EQ.0) THEN
            READ(TST,9500)
            J     = 1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
              READ(TST,9520,END=9385) JDAY,DUMMY
              IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9090) QDTFN(JB),JDAY,TMSTRT
                NERR = NERR+1
              ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9095) QDTFN(JB),JDAY,JDAYO
                NERR = NERR+1
              END IF
              J     = J+1
              JDAYO = JDAY
            END DO
 9385       CONTINUE
            IF (JDAY.LT.TMEND) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9096) QDTFN(JB),JDAY,TMEND
              NERR = NERR+1
            END IF
            CLOSE(TST)
          ELSE
            DELETE_ERR = .FALSE.
            WRITE(ERR,9160) QDTFN(JB),'QDTFN',CONFN
            NERR = NERR+1
          END IF

********* Inflow temperatures

          WRITE (*,*) '    inflow temperatures'
          IF (.NOT.VALID_FILENAME(TDTFN(JB))) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9165) TDTFN(JB),CONFN
            NWRN = NWRN+1
          END IF
          OPEN (UNIT=TST,FILE=TDTFN(JB),STATUS='OLD',IOSTAT=IERR)
          IF (IERR.EQ.0) THEN
            READ(TST,9500)
            J     = 1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
              READ(TST,9520,END=9395) JDAY,DUMMY
              IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9090) TDTFN(JB),JDAY,TMSTRT
                NERR = NERR+1
              ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9090) TDTFN(JB),JDAY,JDAYO
                NERR = NERR+1
              END IF
              J     = J+1
              JDAYO = JDAY
            END DO
 9395       CONTINUE
            IF (JDAY.LT.TMEND) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9096) TDTFN(JB),JDAY,TMEND
              NERR = NERR+1
            END IF
            CLOSE(TST)
          ELSE
            DELETE_ERR = .FALSE.
            WRITE(ERR,9160) TDTFN(JB),'TDTFN',CONFN
            NERR = NERR+1
          END IF

********* Inflow concentrations

          IF (CONSTITUENTS) THEN
            WRITE (*,*) '    inflow concentrations'
            IF (.NOT.VALID_FILENAME(CDTFN(JB))) THEN
              DELETE_WRN = .FALSE.
              WRITE(WRN,9165) CDTFN(JB),CONFN
              NWRN = NWRN+1
            END IF
            OPEN (UNIT=TST,FILE=CDTFN(JB),IOSTAT=IERR,STATUS='OLD')
            IF (IERR.EQ.0) THEN
              READ(TST,9500)
              J     = 1
              JDAYO = 0.0
              DO WHILE (.TRUE.)
                READ(TST,*,END=9405) JDAY,(DUMMY,JC=1,NACDT)
                IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                  DELETE_ERR = .FALSE.
                  WRITE(ERR,9090) CDTFN(JB),JDAY,TMSTRT
                  NERR = NERR+1
                ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                  DELETE_ERR = .FALSE.
                  WRITE(ERR,9095) CDTFN(JB),JDAY,JDAYO
                  NERR = NERR+1
                END IF
                J     = J+1
                JDAYO = JDAY
              END DO
 9405         CONTINUE
              IF (JDAY.LT.TMEND) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9096) CDTFN(JB),JDAY,TMEND
                NERR = NERR+1
              END IF
              CLOSE(TST)
            ELSE
              DELETE_WRN = .FALSE.
              WRITE(WRN,9110) 'DIST_TR',CDTFN(JB)
              NWRN = NWRN+1
            END IF
          END IF
        END IF

******* Precipitation
        
        IF (PRECIPITATION) THEN    
          WRITE (*,*) '  precipitation'
          WRITE (*,*) '    inflows'
          IF (.NOT.VALID_FILENAME(PREFN(JB))) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9165) PREFN(JB), CONFN
            NWRN = NWRN+1
          END IF

********* Inflows

          OPEN (UNIT=TST,FILE=PREFN(JB),STATUS='OLD')
          IF (IERR.GT.0) THEN 
            DELETE_ERR = .FALSE.
            WRITE(ERR,9050) PREFN(JB)
            NERR = NERR+1
          END IF
          READ(TST,9500)
          J     = 1
          JDAYO = 0.0
          DO WHILE (.TRUE.)
            READ(TST,9520,END=9415) JDAY,DUMMY
            IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9090) PREFN(JB),JDAY,TMSTRT
              NERR = NERR+1
            ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9095) PREFN(JB),JDAY,JDAYO
              NERR = NERR+1
            END IF
            J     = J+1
            JDAYO = JDAY
          END DO
 9415     CONTINUE
          IF (JDAY.LT.TMEND) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9096) PREFN(JB),JDAY,TMEND
            NERR = NERR+1
          END IF
          CLOSE(TST)

********* Inflow temperatures

          WRITE (*,*) '    inflow temperatures'
          IF (.NOT.VALID_FILENAME(TPRFN(JB))) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9165) TPRFN(JB),CONFN
            NWRN = NWRN+1
          END IF
          OPEN (UNIT=TST,FILE=TPRFN(JB),STATUS='OLD',IOSTAT=IERR)
          IF (IERR.EQ.0) THEN
            READ(TST,9500)
            J     = 1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
              READ(TST,9520,END=9425) JDAY,DUMMY
              IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9090) TPRFN(JB),JDAY,TMSTRT
                NERR = NERR+1
              ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9095) TPRFN(JB),JDAY,JDAYO
                NERR = NERR+1
              END IF
              J     = J+1
              JDAYO = JDAY
            END DO
 9425       CONTINUE
            IF (JDAY.LT.TMEND) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9096) TPRFN(JB),JDAY,TMEND
              NERR = NERR+1
            END IF
            CLOSE(TST)
          ELSE
            DELETE_ERR = .FALSE.
            WRITE(ERR,9160) TPRFN(JB),'TPRFN',CONFN
            NERR = NERR+1
          END IF

********* Inflow concentrations
          
          IF (CONSTITUENTS) THEN
            WRITE (*,*) '    inflow concentrations'
            IF (.NOT.VALID_FILENAME(CPRFN(JB))) THEN
              DELETE_WRN = .FALSE.
              WRITE(WRN,9165) CPRFN(JB),CONFN
              NWRN = NWRN+1
            END IF
            OPEN (UNIT=TST,FILE=CPRFN(JB),IOSTAT=IERR,STATUS='OLD')
            IF (IERR.EQ.0) THEN
              READ(TST,9500)
              J     = 1
              JDAYO = 0.0
              DO WHILE (.TRUE.)
                READ(TST,*,END=9435) JDAY,(DUMMY,JC=1,NACPR)
                IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                  DELETE_ERR = .FALSE.
                  WRITE(ERR,9090) CPRFN(JB),JDAY,TMSTRT
                  NERR = NERR+1
                ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                  DELETE_ERR = .FALSE.
                  WRITE(ERR,9095) CPRFN(JB),JDAY,JDAYO
                  NERR = NERR+1
                END IF
                J     = J+1
                JDAYO = JDAY
              END DO
 9435         CONTINUE
              IF (JDAY.LT.TMEND) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9096) CPRFN(JB),JDAY,TMEND
                NERR = NERR+1
              END IF
              CLOSE(TST)
            ELSE
              DELETE_WRN = .FALSE.
              WRITE(WRN,9110) 'PRECIP',CPRFN(JB)
              NWRN = NWRN+1
            END IF
          END IF
        END IF

******* Upstream head
        
        IF (UH_EXTERNAL(JB)) THEN        

********* Elevations

          WRITE (*,*) '  external upstream head'
          WRITE (*,*) '    elevations'
          IF (.NOT.VALID_FILENAME(EUHFN(JB))) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9165) EUHFN(JB),CONFN
            NWRN = NWRN+1
          END IF
          OPEN (UNIT=TST,FILE=EUHFN(JB),STATUS='OLD',IOSTAT=IERR)
          IF (IERR.EQ.0) THEN
            READ(TST,9500)
            J     = 1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
              READ(TST,9520,END=9445) JDAY,DUMMY
              IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9090) EUHFN(JB),JDAY,TMSTRT
                NERR = NERR+1
              ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9095) EUHFN(JB),JDAY,JDAYO
                NERR = NERR+1
              END IF
              J     = J+1
              JDAYO = JDAY
            END DO
 9445       CONTINUE
            IF (JDAY.LT.TMEND) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9096) EUHFN(JB),JDAY,TMEND
              NERR = NERR+1
            END IF
            CLOSE(TST)
          ELSE
            DELETE_ERR = .FALSE.
            WRITE(ERR,9160) EUHFN(JB),'EUHFN',CONFN
            NERR = NERR+1
          END IF

********* Temperatures

          WRITE (*,*) '    temperatures'
          IF (.NOT.VALID_FILENAME(TUHFN(JB))) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9165) TUHFN(JB),CONFN
            NWRN = NWRN+1
          END IF
          OPEN (UNIT=TST,FILE=TUHFN(JB),STATUS='OLD',IOSTAT=IERR)
          IF (IERR.EQ.0) THEN
            READ(TST,9500)
            J     = J+1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
              READ(TST,9520,END=9455) JDAY,(DUMMY,K=2,KMP-1)
              IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9090) TUHFN(JB),JDAY,TMSTRT
                NERR = NERR+1
              ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9095) TUHFN(JB),JDAY,JDAYO
                NERR = NERR+1
              END IF
              J     = J+1
              JDAYO = JDAY
            END DO
 9455       CONTINUE
            IF (JDAY.LT.TMEND) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9096) TUHFN(JB),JDAY,TMEND
              NERR = NERR+1
            END IF
            CLOSE(TST)
          ELSE
            DELETE_ERR = .FALSE.
            WRITE(ERR,9160) TUHFN(JB),'TUHFN',CONFN
            NERR = NERR+1
          END IF

********* Constituent concentrations

          IF (CONSTITUENTS) THEN
            WRITE (*,*) '    constituent concentrations'
            IF (.NOT.VALID_FILENAME(CUHFN(JB))) THEN
              DELETE_WRN = .FALSE.
              WRITE(WRN,9165) CUHFN(JB),CONFN
              NWRN = NWRN+1
            END IF
            OPEN (UNIT=TST,FILE=CUHFN(JB),STATUS='OLD',IOSTAT=IERR)
            IF (IERR.EQ.0) THEN
              READ(TST,9500)
              J     = 1
              JDAYO = 0.0
              DO WHILE (.TRUE.)
                DO JC=1,NAC
                  READ(TST,*,END=9465) JDAY,(DUMMY,K=2,KMP-1)
                  IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                    DELETE_ERR = .FALSE.
                    WRITE(ERR,9090) CUHFN(JB),JDAY,TMSTRT
                    NERR = NERR+1
                  ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                    DELETE_ERR = .FALSE.
                    WRITE(ERR,9095) CUHFN(JB),JDAY,JDAYO
                    NERR = NERR+1
                  END IF
                END DO
                J     = J+1
                JDAYO = JDAY
              END DO
 9465         CONTINUE
              IF (JDAY.LT.TMEND) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9096) CUHFN(JB),JDAY,TMEND
                NERR = NERR+1
              END IF
              CLOSE(TST)
            ELSE
              DELETE_ERR = .FALSE.
              WRITE(ERR,9160) CUHFN(JB),'CUHFN',CONFN
              NERR = NERR+1
            END IF
          END IF
        END IF

******* Downstream head

        IF (DH_EXTERNAL(JB)) THEN

********* Elevations

          WRITE (*,*) '  downstream head'
          WRITE (*,*) '    elevations'
          IF (.NOT.VALID_FILENAME(EDHFN(JB))) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9165) EDHFN(JB)
            NWRN = NWRN+1
          END IF
          OPEN (UNIT=TST,FILE=EDHFN(JB),STATUS='OLD',IOSTAT=IERR)
          IF (IERR.EQ.0) THEN
            READ(TST,9500)
            J     = 1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
              READ(TST,9520,END=9475) JDAY,DUMMY
              IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9090) EDHFN(JB),JDAY,TMSTRT
                NERR = NERR+1
              ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9095) EDHFN(JB),JDAY,JDAYO
                NERR = NERR+1
              END IF
              J     = J+1
              JDAYO = JDAY
            END DO
 9475       CONTINUE
            IF (JDAY.LT.TMEND) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9096) EDHFN(JB),JDAY,TMEND
              NERR = NERR+1
            END IF
            CLOSE(TST)
          ELSE
            DELETE_ERR = .FALSE.
            WRITE(ERR,9160) EDHFN(JB),'EDHFN',CONFN
            NERR = NERR+1
          END IF

********* Temperatures

          WRITE (*,*) '    temperatures'
          IF (.NOT.VALID_FILENAME(TDHFN(JB))) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9165) TDHFN(JB),CONFN
            NWRN = NWRN+1
          END IF
          OPEN (UNIT=TST,FILE=TDHFN(JB),STATUS='OLD',IOSTAT=IERR)
          IF (IERR.EQ.0) THEN
            READ(TST,9500)
            J     = 1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
              READ(TST,1120,END=9485) JDAY,(DUMMY,K=2,KMP-1)
              IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9090) TDHFN(JB),JDAY,TMSTRT
                NERR = NERR+1
              ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9095) TDHFN(JB),JDAY,JDAYO
                NERR = NERR+1
              END IF
              J     = J+1
              JDAYO = JDAY
            END DO
 9485       CONTINUE
            IF (JDAY.LT.TMEND) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9096) TDHFN(JB),JDAY,TMEND
              NERR = NERR+1
            END IF
            CLOSE(TST)
          ELSE
            DELETE_ERR = .FALSE.
            WRITE(ERR,9160) TDHFN(JB),'TDHFN',CONFN
            NERR = NERR+1
          END IF

********* Constituent concentrations

          IF (CONSTITUENTS) THEN
            WRITE (*,*) '    concentrations'
            IF (.NOT.VALID_FILENAME(CDHFN(JB))) THEN
              DELETE_WRN = .FALSE.
              WRITE(WRN,9165) CDHFN(JB),CONFN
              NWRN = NWRN+1
            END IF
            OPEN (UNIT=TST,FILE=CDHFN(JB),STATUS='OLD',IOSTAT=IERR)
            IF (IERR.EQ.0) THEN
              READ(TST,9500)
              J     = 1
              JDAYO = 0.0
              DO WHILE (.TRUE.)
                DO JC=1,NAC
                  READ(TST,1120,END=9495) JDAY,(DUMMY,K=2,KMP-1)
                  IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                    DELETE_ERR = .FALSE.
                    WRITE(ERR,9090) CDHFN(JB),JDAY,TMSTRT
                    NERR = NERR+1
                  ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                    DELETE_ERR = .FALSE.
                    WRITE(ERR,9095) CDHFN(JB),JDAY,JDAYO
                    NERR = NERR+1
                  END IF
                END DO
                J     = J+1
                JDAYO = JDAY
              END DO
 9495         CONTINUE
              IF (JDAY.LT.TMEND) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9096) CDHFN(JB),JDAY,TMEND
                NERR = NERR+1
              END IF
              CLOSE(TST)
            ELSE
              DELETE_ERR = .FALSE.
              WRITE(ERR,9160) CDHFN(JB),'CDHFN',CONFN
              NERR = NERR+1
            END IF
          END IF
        END IF
      END DO    

***** Tributaries

      IF (TRIBUTARIES) THEN
        DO JT=1,NTR

********* Inflows

          WRITE (*,*) '  tributary',JT
          WRITE (*,*) '    inflows'
          IF (.NOT.VALID_FILENAME(QTRFN(JT))) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9165) QTRFN(JT),CONFN
            NWRN = NWRN+1
          END IF
          OPEN (UNIT=TST,FILE=QTRFN(JT),STATUS='OLD',IOSTAT=IERR)
          IF (IERR.EQ.0) THEN
            READ(TST,9500)
            J     = 1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
              READ(TST,9520,END=9355) JDAY,DUMMY
              IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9090) QTRFN(JT),JDAY,TMSTRT
                NERR = NERR+1
              ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9095) QTRFN(JT),JDAY,JDAYO
                NERR = NERR+1
              END IF
              J     = J+1
              JDAYO = JDAY
            END DO
 9355       CONTINUE
            IF (JDAY.LT.TMEND) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9096) QTRFN(JT),JDAY,TMEND
              NERR = NERR+1
            END IF
            CLOSE(TST)
          ELSE
            DELETE_ERR = .FALSE.
            WRITE(ERR,9160) QTRFN(JT),'QTRFN',CONFN
            NERR = NERR+1
          END IF

********* Inflow temperatures

          WRITE (*,*) '    inflow temperatures'
          IF (.NOT.VALID_FILENAME(TTRFN(JT))) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9165) TTRFN(JT),CONFN
            NWRN = NWRN+1
          END IF
          OPEN (UNIT=TST,FILE=TTRFN(JT),STATUS='OLD',IOSTAT=IERR)
          IF (IERR.EQ.0) THEN
            READ(TST,9500)
            J     = 1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
              READ(TST,9520,END=9365) JDAY,DUMMY
              IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9090) TTRFN(JT),JDAY,TMSTRT
                NERR = NERR+1
              ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9095) TTRFN(JT),JDAY,JDAYO
                NERR = NERR+1
              END IF
              J     = J+1
              JDAYO = JDAY
            END DO
 9365       CONTINUE
            IF (JDAY.LT.TMEND) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9096) TTRFN(JT),JDAY,TMEND
              NERR = NERR+1
            END IF
            CLOSE(TST)
          ELSE
            DELETE_ERR = .FALSE.
            WRITE(ERR,9160) TTRFN(JT),'TTRFN',CONFN
            NERR = NERR+1
          END IF

********* Inflow concentrations

          IF (CONSTITUENTS) THEN
            WRITE (*,*) '    inflow concentrations'
            IF (.NOT.VALID_FILENAME(CTRFN(JT))) THEN
              DELETE_WRN = .FALSE.
              WRITE(WRN,9165) CTRFN(JT),CONFN
              NWRN = NWRN+1
            END IF
            OPEN (UNIT=TST,FILE=CTRFN(JT),IOSTAT=IERR,STATUS='OLD')
            IF (IERR.EQ.0) THEN
              READ(TST,9500)
              J     = 1
              JDAYO = 0.0
              DO WHILE (.TRUE.)
                READ(TST,*,END=9375) JDAY,(DUMMY,JC=1,NACTR)
                IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                  DELETE_ERR = .FALSE.
                  WRITE(ERR,9090) CTRFN(JT),JDAY,TMSTRT
                  NERR = NERR+1
                ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                  DELETE_ERR = .FALSE.
                  WRITE(ERR,9095) CTRFN(JT),JDAY,JDAYO
                  NERR = NERR+1
                END IF
                J     = J+1
                JDAYO = JDAY
              END DO
 9375         CONTINUE
              IF (JDAY.LT.TMEND) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9096) CTRFN(JT),JDAY,TMEND
                NERR = NERR+1
              END IF
              CLOSE(TST)
            ELSE
              DELETE_WRN = .FALSE.
              WRITE(WRN,9110) 'TRIB',CTRFN(JT)
              NWRN = NWRN+1
            END IF
          END IF
        END DO
      END IF

***** Outflows

      DO JB=1,NBP
        IF (DN_FLOW(JB)) THEN
          WRITE (*,*) '  outflows for branch',JB
          IF (.NOT.VALID_FILENAME(QOTFN(JB))) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9165) QOTFN(JB),CONFN
            NWRN = NWRN+1
          END IF
          OPEN (UNIT=TST,FILE=QOTFN(JB),STATUS='OLD',IOSTAT=IERR)
          IF (IERR.EQ.0) THEN
            READ(TST,9500)
            J     = 1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
              READ(TST,*,END=9345) JDAY,(DUMMY,JS=1,NSTR(JB))
              IF (J.EQ.1.AND.JDAY.GT.TMSTRT) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9090) QOTFN(JB),JDAY,TMSTRT
                NERR = NERR+1
              ELSE IF (JDAY.LE.JDAYO.AND.J.NE.1) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9095) QOTFN(JB),JDAY,JDAYO
                NERR = NERR+1
              END IF
              J     = J+1
              JDAYO = JDAY
            END DO
 9345       CONTINUE
            IF (JDAY.LT.TMEND) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9096) QOTFN(JB),JDAY,TMEND
              NERR = NERR+1
            END IF
            CLOSE(TST)
          ELSE
            DELETE_ERR = .FALSE.
            WRITE(ERR,9160) QOTFN(JB),'QOTFN',CONFN
            NERR = NERR+1
          END IF
        END IF
      END DO
      DUMMY = DUMMY

************************************************************************
**                         Control File Inputs                        **
************************************************************************

      WRITE (*,*) 'Control file input checks'
      WRITE (*,*) '  timestep control'
      IF (TMEND.LT.TMSTRT) THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9029) TMSTRT,TMEND
        NERR = NERR+1
      END IF
      IF (DLTD(1).GT.TMSTRT) THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9030) DLTD(1),TMSTRT
        NERR = NERR+1
      END IF
      DO J=1,NDLT-1
        IF (DLTD(J+1).LT.DLTD(J)) THEN
          DELETE_ERR = .FALSE.
          WRITE (ERR,9032) J+1,DLTD(J+1),J,DLTD(J)
          NERR = NERR+1
        END IF
      END DO
      DO J=1,NDLT
        IF (DLTMAX(J).EQ.0.0) THEN
          DELETE_ERR = .FALSE.
          WRITE (ERR,9025) J
          NERR = NERR+1
        END IF
        IF (DLTF(J).EQ.0.0) THEN
          DELETE_ERR = .FALSE.
          WRITE (ERR,9026) J
          NERR = NERR+1
        END IF
      END DO

***** Water body type

      WRITE (*,*) '  waterbody type'
      IF (WTYPEC.EQ.'SALT') THEN
        IF (.NOT. CONSTITUENTS) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,*) 'Since water type is SALT, constituents [CCC] ',
     .                 'must be set to "ON".'    
          NERR = NERR+1
        ELSE IF (.NOT.DISSOLVED_SOLIDS) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,*) 'Since water type is SALT, ACC(4) must be set',
     .                 ' to ON.'
          NERR = NERR+1
        END IF
      END IF

***** Calculations

      IF (VBC.NE.' ON'.AND.VBC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'VBC',CONFN
        NERR = NERR+1
      END IF
      IF (EBC.NE.' ON'.AND.EBC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'EBC',CONFN
        NERR = NERR+1
      END IF
      IF (MBC.NE.' ON'.AND.MBC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'MBC',CONFN
        NERR = NERR+1
      END IF
      IF (WBC.NE.' ON'.AND.WBC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'WBC',CONFN
        NERR = NERR+1
      END IF
      IF (PQC.NE.' ON'.AND.PQC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'PQC',CONFN
        NERR = NERR+1
      END IF
      DO JT=1,NTR
        IF (TRPC(JT).NE.' SPECIFY'.AND.TRpC(JT).NE.'   DISTR'           !4/03/95
     .      .AND.TRPC(JT).NE.' DENSITY') THEN                           !4/03/95
          DELETE_ERR = .FALSE.
          WRITE (ERR,9098) 'TRPC',CONFN                                 !4/03/95
          NERR = NERR+1
        END IF
      END DO
      IF (EVC.NE.' ON'.AND.EVC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'EVC',CONFN
        NERR = NERR+1
      END IF
      IF (PRC.NE.' ON'.AND.PRC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'PRC',CONFN
        NERR = NERR+1
      END IF

***** Interpolation controls

      IF (INFIC.NE.' ON'.AND.INFIC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'INFIC',CONFN
        NERR = NERR+1
      END IF
      IF (TRIC.NE.' ON'.AND.TRIC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'TRIC',CONFN
        NERR = NERR+1
      END IF
      IF (DTRIC.NE.' ON'.AND.DTRIC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'DTRIC',CONFN
        NERR = NERR+1
      END IF
      IF (HDIC.NE.' ON'.AND.HDIC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'HDIC',CONFN
        NERR = NERR+1
      END IF
      IF (OUTIC.NE.' ON'.AND.OUTIC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'OUTIC',CONFN
        NERR = NERR+1
      END IF
      IF (WDIC.NE.' ON'.AND.WDIC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'WDIC',CONFN
        NERR = NERR+1
      END IF
      IF (METIC.NE.' ON'.AND.METIC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'METIC',CONFN
        NERR = NERR+1
      END IF

***** Dead Sea

      IF (WINDC.NE.' ON'.AND.WINDC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'WINDC',CONFN
        NERR = NERR+1
      END IF
      IF (QINC.NE.' ON'.AND.QINC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'QINC',CONFN
        NERR = NERR+1
      END IF
      IF (QOUTC.NE.' ON'.AND.QOUTC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'QOUTC',CONFN
        NERR = NERR+1
      END IF
      IF (HEATC.NE.' ON'.AND.HEATC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'HEATC',CONFN
        NERR = NERR+1
      END IF

***** Ice cover

      IF (ICEC.NE.' ON'.AND.ICEC.NE.'OFF') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9099) 'ICEC',CONFN
        NERR = NERR+1
      END IF
      IF (SLICEC.NE.'  DETAIL'.AND.SLICEC.NE.'  SIMPLE') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9011)
        NERR = NERR+1
      END IF
      IF (SLHTC.NE.'    TERM'.AND.SLHTC.NE.'      ET') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9012)
        NERR = NERR+1
      END IF
      IF (ALBEDO.LT.0.0) THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9155) ALBEDO,CONFN
        NERR = NERR+1
      END IF
      IF (HWI.LT.0.0) THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9155) HWI,CONFN
        NERR = NERR+1
      END IF
      IF (BETAI.LT.0.0) THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9155) BETAI,CONFN
        NERR = NERR+1
      END IF
      IF (BETAI.GT.0.9) THEN
        DELETE_WRN = .FALSE.
        WRITE (WRN,9013)
        NWRN = NWRN+1
      END IF
      IF (GAMMAI.LT.0.0) THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9155) GAMMAI,CONFN
        NERR = NERR+1
      END IF
      IF (ICEMIN.LT.0.0) THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9155) ICEMIN,CONFN
        NERR = NERR+1
      END IF
      IF (ICEMIN.GT.0.1) THEN
        DELETE_WRN = .FALSE.
        WRITE (WRN,9014)
        NWRN = NWRN+1
      END IF
      IF (ICET2.LT.0.0) THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9155) ICET2,CONFN
        NERR = NERR+1
      END IF
      IF (ICET2.GT.4.0) THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9015)
        NERR = NERR+1
      END IF

***** Transport solution

      IF (SLTRC.NE.'QUICKEST'.AND.SLTRC.NE.'  UPWIND') THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9016)
        NERR = NERR+1
      END IF
      IF (THETA.LT.0.0) THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9017)
        NERR = NERR+1
      END IF
      IF (THETA.GT.1.0) THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9018)
        NERR = NERR+1
      END IF

***** Wind sheltering date

      WRITE (*,*) '  wind sheltering'
      IF (NWSD.LE.0) THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9034)
        NERR = NERR+1
      END IF
      IF (WSCD(1).GT.TMSTRT) THEN
        DELETE_ERR = .FALSE.
        WRITE (ERR,9035) WSCD(1),TMSTRT
        NERR = NERR+1
      END IF
      DO J=1,NWSD
        IF (WSC(J).GT.1.OR.WSC(J).LT.0) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9107) 'WSC',J,CONFN
          NERR = NERR+1
        END IF
      END DO
      DO J=1,NWSD-1
        IF (NWSD.GT.1.AND.WSCD(J+1).LT.WSCD(J)) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9036) J+1,WSCD(J+1),J,WSCD(J)
          NERR = NERR+1
        END IF
      END DO

***** Hydraulic coefficients

      WRITE (*,*) '  hydraulic coefficients'
      IF (AX.LT.0.1.OR.AX.GT.10.0) THEN
        DELETE_WRN = .FALSE.
        WRITE(WRN,9105) 'AX',CONFN
        NWRN = NWRN+1
      END IF
      IF (DXI.LT.0.1.OR.DXI.GT.10.0) THEN
        DELETE_WRN = .FALSE.
        WRITE(WRN,9105) 'DX',CONFN
        NWRN = NWRN+1
      END IF
      IF (CHEZY.LT.30.0.OR.CHEZY.GT.100.0) THEN
        DELETE_WRN = .FALSE.
        WRITE(WRN,9105) 'CHEZY',CONFN
        NWRN = NWRN+1
      END IF
      IF (CBHE.GT.1.0E-6) THEN
        DELETE_WRN = .FALSE.
        WRITE(WRN,9105) 'CBHE',CONFN
        NWRN = NWRN+1
      END IF
      IF (TSED.GT.30.OR.TSED.LT.0.0) THEN
        DELETE_ERR = .FALSE.
        WRITE(ERR,9105) 'TSED',CONFN
        NERR = NERR+1
      END IF

***** Light extinction

      WRITE (*,*) '  light extinction'
      IF (EXH2O.GT.0.95) THEN
        DELETE_WRN = .FALSE.
        WRITE(WRN,9150) 'EXH2O',0.95,CONFN
        NWRN = NWRN+1
      ELSE IF (EXH2O.LT.0.1) THEN
        DELETE_WRN = .FALSE.
        WRITE(WRN,9155) 'EXH2O',CONFN
        NWRN = NWRN+1
      END IF
      IF (BETA.GT.0.75) THEN
        DELETE_WRN = .FALSE.
        WRITE(WRN,9150) 'BETA',0.75,CONFN
        NWRN = NWRN+1
      ELSE IF (BETA.LE.0.1) THEN
        DELETE_WRN = .FALSE.
        WRITE(WRN,9155) 'BETA',CONFN
        NWRN = NWRN+1
      END IF

***** Tributaries

      IF (TRIBUTARIES) THEN
        WRITE (*,*) '  tributaries'
        DO JT=1,NTP
          ISOK = .FALSE.
          DO JB=1,NBP
            IF (ITR(JT).GE.US(JB).AND.ITR(JT).LE.DS(JB)) THEN
              ISOK = .TRUE.
            END IF
          END DO
          IF (.NOT.ISOK) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9107) 'ITR',JT,CONFN
            NERR = NERR+1
          END IF
        END DO
      END IF

***** Distributed tributaries

      DO JB=1,NBP
        IF (DIST_TRIB(JB)) THEN
          WRITE (*,*) '  distributed tributaries'
          IF ((DTRC(JB).NE.' ON').AND.(DTRC(JB).NE.'OFF')) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9099) 'DTRC',JB,CONFN
            NERR = NERR+1
          END IF
        END IF
      END DO

***** Screen updates

      IF ((SCRC.NE.' ON').AND.(SCRC.NE.'OFF')) THEN
        DELETE_ERR = .FALSE.
        WRITE(ERR,9099) 'SCRC',CONFN
        NERR = NERR+1
      END IF
      DO J=1,NSCR-1
        IF (SCRD(J+1).LT.SCRD(J)) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9125) 'SCRD',J+1,'SCRD',J,CONFN
          NERR = NERR+1
        END IF
      END DO
      DO J=1,NSCR
        IF (SCRF(J).EQ.0.0) THEN
          DELETE_WRN = .FALSE.
          WRITE(WRN,9107) 'SCRF',J,CONFN
          NWRN = NWRN+1
        END IF
      END DO

***** Snapshot

      WRITE (*,*) 'Outputs'
      IF ((SNPC.NE.' ON').AND.(SNPC.NE.'OFF')) THEN
        DELETE_ERR = .FALSE.
        WRITE(ERR,9099) 'SNPC',CONFN
        NERR = NERR+1
      END IF
      WRITE (*,*) '  snapshot'
      IF (.NOT.VALID_FILENAME(SNPFN)) THEN
        DELETE_WRN = .FALSE.
        WRITE(WRN,9165) SNPFN,CONFN
        NWRN = NWRN+1
      END IF
      IF ((HPRN(1).NE.' ON').AND.(HPRN(1).NE.'OFF')) THEN
        DELETE_ERR = .FALSE.
        WRITE(ERR,9099) 'UPRNC',CONFN
        NERR = NERR+1
      END IF
      IF ((HPRN(2).NE.' ON').AND.(HPRN(2).NE.'OFF')) THEN
        DELETE_ERR = .FALSE.
        WRITE(ERR,9099) 'WPRNC',CONFN
        NERR = NERR+1
      END IF
      IF ((HPRN(3).NE.' ON').AND.(HPRN(3).NE.'OFF')) THEN
        DELETE_ERR = .FALSE.
        WRITE(ERR,9099) 'TPRNC',CONFN
        NERR = NERR+1
      END IF
      IF ((HPRN(4).NE.' ON').AND.(HPRN(4).NE.'OFF')) THEN
        DELETE_ERR = .FALSE.
        WRITE(ERR,9099) 'DTPRNC',CONFN
        NERR = NERR+1
      END IF

***** Output segments

      DO JB=1,NBP
        JBT = JB
        IF (JB.EQ.NBP) JBT = JB-1
        DO J=1,NISNP
          IF (IPRI(J).EQ.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9120) 'IPR',CONFN
            NERR = NERR+1
          ELSE IF ((IPRI(J).GT.DS(JBT)).AND.(IPRI(J)
     .             .LT.US(JBT+1))) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9115) IPRI(J)
            NWRN = NWRN+1
          END IF
        END DO
      END DO

***** Snapshot dates

      DO J=1,NSNP-1
        IF (SNPD(J+1).LT.SNPD(J)) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9125) 'SNPD',J+1,'SNPD',J,CONFN
          NERR = NERR+1
        END IF
      END DO
      DO J=1,NSNP
        IF (SNPF(J).EQ.0.0) THEN
          DELETE_WRN = .FALSE.
          WRITE(WRN,9107) 'SNPF',J,CONFN
          NWRN = NWRN+1
        END IF
      END DO

***** Profile plots

      IF ((PRFC.NE.' ON').AND.(PRFC.NE.'OFF')) THEN
        DELETE_ERR = .FALSE.
        WRITE(ERR,9099) 'PRFC',CONFN
        NERR = NERR+1
      END IF
      IF (PRFC.EQ.' ON') THEN
        IF (.NOT.VALID_FILENAME(PRFFN)) THEN
          DELETE_WRN = .FALSE.
          WRITE(WRN,9165) PRFFN,CONFN
          NWRN = NWRN+1
        END IF
        WRITE (*,*) '  profiles'
        IF (NPRF.LE.0) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9105) 'NPRF',CONFN
          NERR = NERR+1
        END IF
        IF (NIPRF.LE.0) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9105) 'NIPRF',CONFN
          NERR = NERR+1
        END IF
        DO J=2,NPRF
          IF (PRFD(J-1).GE.PRFD(J)) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9125) 'PRFD',J,'PRFD',J-1,CONFN
            NERR = NERR+1
          END IF
        END DO
        DO J=1,NPRF
          IF (PRFF(J).EQ.0.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9107) 'PRFF',J,CONFN
            NWRN = NWRN+1
          END IF
        END DO
        DO JB=1,NBP
          DO J=1,NIPRF
            IF (IPRF(J).LT.2) THEN
              DELETE_WRN = .FALSE.
              WRITE(WRN,9117) 'IPRF',J,CONFN
              NWRN = NWRN+1
            END IF
            IF ((IPRF(J).GT.DS(JB)).AND.(IPRF(J).LT.US(JB+1))) THEN
              DELETE_WRN = .FALSE.
              WRITE(WRN,9117) 'IPRF',J,CONFN
              NWRN = NWRN+1
            END IF
          END DO
        END DO
      END IF

***** Spreadsheet plots

      IF ((SPRC.NE.' ON').AND.(SPRC.NE.'OFF')) THEN
        DELETE_ERR = .FALSE.
        WRITE(ERR,9099) 'SPRC',CONFN
        NERR = NERR+1
      END IF
      IF (SPRC.EQ.' ON') THEN
        IF (.NOT.VALID_FILENAME(SPRFN)) THEN
          DELETE_WRN = .FALSE.
          WRITE(WRN,9165) SPRFN,CONFN
          NWRN = NWRN+1
        END IF
        WRITE (*,*) '  profiles'
        IF (NSPR.LE.0) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9105) 'NSPR',CONFN
          NERR = NERR+1
        END IF
        IF (NISPR.LE.0) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9105) 'NISPR',CONFN
          NERR = NERR+1
        END IF
        DO J=2,NSPR
          IF (SPRD(J-1).GE.SPRD(J)) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9125) 'SPRD',J,'SPRD',J-1,CONFN
            NERR = NERR+1
          END IF
        END DO
        DO J=1,NSPR
          IF (SPRF(J).EQ.0.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9107) 'SPRF',J,CONFN
            NWRN = NWRN+1
          END IF
        END DO
        DO JB=1,NBP
          DO J=1,NISPR
            IF (ISPR(J).LT.2) THEN
              DELETE_WRN = .FALSE.
              WRITE(WRN,9117) 'ISPR',J,CONFN
              NWRN = NWRN+1
            END IF
            IF ((ISPR(J).GT.DS(JB)).AND.(ISPR(J).LT.US(JB+1))) THEN
              DELETE_WRN = .FALSE.
              WRITE(WRN,9117) 'ISPR',J,CONFN
              NWRN = NWRN+1
            END IF
          END DO
        END DO
      END IF

***** Time series

      IF ((TSRC.NE.' ON').AND.(TSRC.NE.'OFF')) THEN
        DELETE_ERR = .FALSE.
        WRITE(ERR,9099) 'TSRC',CONFN
        NERR = NERR+1
      END IF
      IF (TSRC.EQ.' ON') THEN
        WRITE (*,*) '  time series'
        IF (.NOT.VALID_FILENAME(TSRFN)) THEN
          DELETE_WRN = .FALSE.
          WRITE(WRN,9165) TSRFN,CONFN
          NWRN = NWRN+1
        END IF
        DO J=2,NTSR
          IF (TSRD(J-1).GE.TSRD(J)) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9125) 'TSRD',J,'TSRD',J-1,CONFN
            NERR = NERR+1
          END IF
        END DO
        DO J=1,NTSR
          IF (TSRF(J).EQ.0.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9107) 'TSRF',J,CONFN
            NWRN = NWRN+1
          END IF
        END DO
      END IF

***** Vector plots

      IF ((VPLC.NE.' ON').AND.(VPLC.NE.'OFF')) THEN
        DELETE_ERR = .FALSE.
        WRITE(ERR,9099) 'VPLC',CONFN
        NERR = NERR+1
      END IF
      IF (VPLC.EQ.' ON') THEN
        WRITE (*,*) '  vector plot'
        IF (.NOT.VALID_FILENAME(VPLFN)) THEN
          DELETE_WRN = .FALSE.
          WRITE(WRN,9165) VPLFN,CONFN
          NWRN = NWRN+1
        END IF
        WRITE (*,*) '  vector plots'
        DO J=2,NVPL
          IF (VPLD(J-1).GE.VPLD(J)) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9125) 'VPLD',J,'VPLD',J-1,CONFN
            NERR = NERR+1
          END IF
        END DO
        DO J=1,NVPL
          IF (VPLF(J).EQ.0.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9107) 'VPLF',J,CONFN
            NWRN = NWRN+1
          END IF
        END DO
      END IF

***** Contour plots

      IF ((CPLC.NE.' ON').AND.(CPLC.NE.'OFF')) THEN
        DELETE_ERR = .FALSE.
        WRITE(ERR,9099) 'CPLC',CONFN
        NERR = NERR+1
      END IF
      IF (CPLC.EQ.' ON') THEN
        WRITE (*,*) '  contour plot'
        IF (.NOT.VALID_FILENAME(CPLFN)) THEN
          DELETE_WRN = .FALSE.
          WRITE(WRN,9165) CPLFN,CONFN
          NWRN = NWRN+1
        END IF
        WRITE (*,*) '  contour plots'
        DO J=2,NCPL
          IF (CPLD(J-1).GE.CPLD(J)) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9125) 'CPLD',J,'CPLD',J-1,CONFN
            NERR = NERR+1
          END IF
        END DO
        DO J=1,NCPL
          IF (CPLF(J).EQ.0.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9107) 'CPLF',J,CONFN
            NWRN = NWRN+1
          END IF
        END DO
      END IF

***** Spreadsheet output

      IF (.NOT.VALID_FILENAME(SPRFN)) THEN
        DELETE_WRN = .FALSE.
        WRITE(WRN,9165) SPRFN,CONFN
        NWRN = NWRN+1
      END IF

***** Restart

      IF ((RSOC.NE.' ON').AND.(RSOC.NE.'OFF')) THEN
        DELETE_ERR = .FALSE.
        WRITE(ERR,9099) 'RSOC',CONFN
        NERR = NERR+1
      END IF
      IF (RSOC.EQ.' ON') THEN
        WRITE (*,*) '  restarts'
        IF (NRSO.LE.0) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9105) 'NRSO',CONFN
          NERR = NERR+1
        END IF
        DO J=1,NRSO
          IF (RSOF(J).EQ.0.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9107) 'RSOF',J,CONFN
            NWRN = NWRN+1
          END IF
        END DO
      END IF
      
***** Constituent computations

      IF ((CCC.NE.' ON').AND.(CCC.NE.'OFF')) THEN
        DELETE_ERR = .FALSE.
        WRITE(ERR,9099) 'CCC',CONFN
        NERR = NERR+1
      END IF
      IF (CONSTITUENTS) THEN
        WRITE (*,*) 'Constituents'
        IF ((LIMC.NE.' ON').AND.(LIMC.NE.'OFF')) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9099) 'LIMC',CONFN
          NERR = NERR+1
        END IF
        IF ((SDC.NE.' ON').AND.(SDC.NE.'OFF')) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9099) 'SDC',CONFN
          NERR = NERR+1
        END IF
        IF (FREQUK.GT.12) THEN
          DELETE_WRN = .FALSE.
          WRITE(WRN,9130) 'FREQUK',12.0,CONFN
          NWRN = NWRN+1
        END IF
      END IF

***** Active constituents

      IF (CONSTITUENTS) THEN
        WRITE (*,*) '  active concentrations'
        ISOK = .FALSE.
        DO J=1,NCP
          IF ((ACC(J).NE.' ON').AND.(ACC(J).NE.'OFF')) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9099) 'ACC',CONFN
            NERR = NERR+1
          ELSE IF (ACC(J).EQ.' ON') THEN
            ISOK = .TRUE.
          END IF
        END DO
        IF (.NOT.ISOK) THEN
          DELETE_WRN = .FALSE.
          WRITE(WRN,9135) 'CCC','ACC',CONFN
          NWRN = NWRN+1
        END IF

******* Initial constituent concentrations

        WRITE (*,*) '  initial concentrations'
        DO J=1,NCP
          IF (C2I(J).LT.-2.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9107) 'C2I',J,CONFN
            NERR = NERR+1
          END IF
          IF (ACC(J).EQ.' ON') THEN
            IF (C2I(J).EQ.0.0.AND.J.NE.13) THEN
              DELETE_WRN = .FALSE.
              WRITE(WRN,9142) 'C2I',J,0.0,CONFN
              NWRN = NWRN+1
            END IF
          END IF
        END DO

******* Constituent output

        WRITE (*,*) '  constituent output'
        DO J=1,NCP
          IF ((CPRN(J).NE.' ON').AND.(CPRN(J).NE.'OFF')) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9099) 'CPRN',CONFN
            NERR = NERR+1
          END IF
          IF (CPRN(J).EQ.' ON') THEN
            IF (ACC(J).NE.' ON') THEN
              DELETE_WRN = .FALSE.
              WRITE(WRN,9145) 'CPRN',J,'ACC',J
              NWRN = NWRN+1
            END IF
          END IF
        END DO

******* Inflow constituent concentrations

        WRITE (*,*) '  active concentrations'
        WRITE (*,*) '    inflows'
        DO J=1,NCP
          IF ((INACC(J).NE.' ON').AND.(INACC(J).NE.'OFF')) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9099) 'INACC',CONFN
            NERR = NERR+1
          END IF
          IF (ACC(J).EQ.' ON'.AND.J.NE.13) THEN
            IF (INACC(J).NE.' ON') THEN
              DELETE_WRN = .FALSE.
              WRITE(WRN,9145) 'ACC',J,'INACC',J
              NWRN = NWRN+1
            END IF
          END IF
        END DO

******* Tributary constituent concentrations

        IF (NTR.GE.1) THEN
          WRITE (*,*) '    tributaries'
          DO J=1,NCP
            IF ((TRACC(J).NE.' ON').AND.(TRACC(J).NE.'OFF')) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9099) 'TRACC',CONFN
              NERR = NERR+1
            END IF  
            IF (ACC(J).EQ.' ON'.AND.J.NE.13) THEN
              IF (TRACC(J).NE.' ON') THEN
                DELETE_WRN = .FALSE.
                WRITE(WRN,9145) 'ACC',J,'TRACC',J
                NWRN = NWRN+1
              END IF
            END IF
          END DO
        END IF

******* Distributed tributary concentrations

        DO JB=1,NBP
          IF (DIST_TRIB(JB)) THEN
            WRITE (*,*) '    distributed tributaries'
            DO J=1,NCP
              IF ((DTACC(J).NE.' ON').AND.(DTACC(J).NE.'OFF')) THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9099) 'DTACC',CONFN
                NERR = NERR+1
              END IF  
              IF (ACC(J).EQ.' ON'.AND.J.NE.13) THEN
                IF (DTACC(J).NE.' ON') THEN
                  DELETE_WRN = .FALSE.
                  WRITE(WRN,9145) 'ACC',J,'DTACC',J
                  NWRN = NWRN+1
                END IF
              END IF
            END DO
          END IF    
        END DO
      
******* Precipitation concentrations

        IF (PRC.EQ.' ON') THEN
          WRITE (*,*) '    precipitation'
          DO J=1,NCP
            IF ((PRACC(J).NE.' ON').AND.(PRACC(J).NE.'OFF')) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9099) 'PRACC',CONFN
              NERR = NERR+1
            END IF  
            IF (ACC(J).EQ.' ON'.AND.J.NE.13) THEN
              IF (PRACC(J).NE.' ON') THEN
                DELETE_WRN = .FALSE.
                WRITE(WRN,9145) 'ACC',J,'PRACC',J
                NWRN = NWRN+1
              END IF
            END IF
          END DO
        END IF
      END IF

***** Light extinction

      IF (CONSTITUENTS) THEN
        IF (DISSOLVED_SOLIDS) THEN
          WRITE (*,*) '  light extinction (inorganic)'
          IF (EXSS.GT.0.1) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'EXSS',0.1,CONFN
            NWRN = NWRN+1
          ELSE IF (EXSS.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'EXSS',CONFN
            NERR = NERR+1
          END IF
        END IF
        IF (LABILE_DOM.OR.REFRACTORY_DOM.OR.PHYTOPLANKTON.
     .      OR.PARTICULATE_OM) THEN
          WRITE (*,*) '  light extinction (organic)'
          IF (EXOM.GT.0.85) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'EXOM',0.65,CONFN
            NWRN = NWRN+1
          ELSE IF (EXOM.LE.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'EXOM',CONFN
            NERR = NERR+1
          END IF
        END IF

******* Coliform

        IF (COLIFORM) THEN
          WRITE (*,*) '  coliform'
          IF (COLQ10.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'COLQ10',CONFN
            NERR = NERR+1
          END IF
          IF (COLDK.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'COLDK',CONFN
            NERR = NERR+1
          END IF
        END IF

******* Suspended solids

        IF (SUSPENDED_SOLIDS) THEN
          WRITE (*,*) '  suspended solids'
          IF (SSS.GT.10.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'SSS',10.0,CONFN
            NWRN = NWRN+1
          ELSE IF (SSS.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'SSS',CONFN
            NERR = NERR+1
          END IF
        END IF
      
******* Algae

        IF (PHYTOPLANKTON) THEN
          WRITE (*,*) '  algae'
          IF (AG.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'AG',CONFN
            NERR = NERR+1
          ELSE IF (AG.LT.0.5) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9151) 'AG',0.5,CONFN
            NWRN = NWRN+1
          ELSE IF (AG.GT.4.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'AG',4.0,CONFN
            NWRN = NWRN+1
          END IF
          IF (AM.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'AM',CONFN
            NERR = NERR+1
          ELSE IF (AM.LT.0.001) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9151) 'AM',0.001,CONFN
            NWRN = NWRN+1
          ELSE IF (AM.GT.1.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'AM',1.0,CONFN
            NWRN = NWRN+1
          END IF
          IF (AE.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'AE',CONFN
            NERR = NERR+1
          ELSE IF (AE.LT.0.001) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9151) 'AE',0.001,CONFN
            NWRN = NWRN+1
          ELSE IF (AE.GT.1.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'AE',1.0,CONFN
            NWRN = NWRN+1
          END IF
          IF (AR.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'AR',CONFN
            NERR = NERR+1
          ELSE IF (AR.LT.0.001) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9151) 'AR',0.001,CONFN
            NWRN = NWRN+1
          ELSE IF (AR.GT.1.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'AR',1.0,CONFN
            NWRN = NWRN+1
          END IF
          IF (AS.GT.0.5) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'AS',0.5,CONFN
            NWRN = NWRN+1
          ELSE IF (AS.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'AS',CONFN
            NERR = NERR+1
          END IF
          IF (ASAT.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'ASAT',CONFN
            NERR = NERR+1
          END IF
          IF (APOM.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'APOM',CONFN
            NERR = NERR+1
          ELSE IF (APOM.LT.0.5) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9151) 'APOM',0.5,CONFN
            NWRN = NWRN+1
          END IF

********* Algal temperature rate constants

          IF (AT1.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'AT1',CONFN
            NERR = NERR+1
          END IF
          IF (AT2.GT.40.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'AT2',40.0,CONFN
            NWRN = NWRN+1
          ELSE IF (AT2.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'AT2',CONFN
            NERR = NERR+1
          END IF
          IF (AT3.GT.40.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'AT3',40.0,CONFN
            NWRN = NWRN+1
          ELSE IF (AT3.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'AT3',CONFN
            NERR = NERR+1
          END IF
          IF (AT4.GT.50.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'AT4',40.0,CONFN
            NWRN = NWRN+1
          ELSE IF (AT4.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'AT4',CONFN
            NERR = NERR+1
          END IF
          IF (AK1.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'AK1',CONFN
            NERR = NERR+1
          END IF
          IF (AK2.GT.1.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9150) 'AK2',1.0,CONFN
            NERR = NERR+1
          ELSE IF (AK2.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'AK2',CONFN
            NERR = NERR+1
          END IF
          IF (AK3.GT.1.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9150) 'AK3',1.0,CONFN
            NERR = NERR+1
          ELSE IF (AK3.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'AK3',CONFN
            NERR = NERR+1
          END IF
          IF (AK4.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'AK4',CONFN
            NERR = NERR+1
          END IF
        END IF

******* Dissolved organics

        IF (LABILE_DOM) THEN
          WRITE (*,*) '  dissolved organics'
          IF (LDOMDK.GT.0.3) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'LDOMDK',0.3,CONFN
            NWRN = NWRN+1
          ELSE IF (LDOMDK.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'LDOMDK',CONFN
            NERR = NERR+1
          END IF
          IF (LRDDK.GT.0.01) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'LRDDK',0.01,CONFN
            NWRN = NWRN+1
          ELSE IF (LRDDK.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'LRDDK',CONFN
            NERR = NERR+1
          END IF
        END IF
        IF (REFRACTORY_DOM) THEN
          IF (RDOMDK.GT.0.01) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'RDOMDK',0.01,CONFN
            NWRN = NWRN+1
          ELSE IF (RDOMDK.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'RDOMDK',CONFN
            NERR = NERR+1
          END IF
        END IF

******* Particulate organic matter

        IF (PARTICULATE_OM) THEN
          WRITE (*,*) '  POM'
          IF (POMDK.GT.0.5) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'POMDK',0.5,CONFN
            NWRN = NWRN+1
          ELSE IF (POMDK.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'POMDK',CONFN
            NERR = NERR+1
          END IF
          IF (POMS.GT.1.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'POMS',1.0,CONFN
            NWRN = NWRN+1
          ELSE IF (POMS.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'POMS',CONFN
            NERR = NERR+1
          END IF
        END IF

******* Organic matter temperature rate multipliers

        IF (LABILE_DOM.OR.REFRACTORY_DOM.OR.PARTICULATE_OM) THEN
          IF (OMT1.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'OMT1',CONFN
            NERR = NERR+1
          END IF
          IF (OMT2.GT.40.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'OMT2',40.0,CONFN
            NWRN = NWRN+1
          ELSE IF (OMT2.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'OMT2',CONFN
            NERR = NERR+1
          END IF
          IF (OMK1.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'OMK1',CONFN
            NERR = NERR+1
          END IF
          IF (OMK2.GT.1.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9150) 'OMK2',1.0,CONFN
            NERR = NERR+1
          ELSE IF (OMK2.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'OMK2',CONFN
            NERR = NERR+1
          END IF
        END IF

******* Sediment decay

        IF (SEDIMENT) THEN
          WRITE (*,*) '  sediments'
          IF (SDK.GT.0.3) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'SDK',0.3,CONFN
            NWRN = NWRN+1
          ELSE IF (SDK.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'SDK',CONFN
            NERR = NERR+1
          END IF
        END IF

******* Sediment oxygen demand

        IF (DISSOLVED_OXYGEN) THEN
          WRITE (*,*) '  sediment oxygen demand'
          DO I=1,IMP
            IF (SOD(I).GT.5.0) THEN
              DELETE_WRN = .FALSE.
              WRITE(WRN,9153) 'SOD',I,5.0,CONFN
              NWRN = NWRN+1
            ELSE IF (SOD(I).LT.0.0) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9155) 'SOD',CONFN
              NERR = NERR+1
            END IF
          END DO
        END IF

******* Carbonaceous biochemical oxygen demand

        IF (CBO_DEMAND) THEN
          WRITE (*,*) '  CBOD'
          IF (KBOD.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'KBOD',CONFN
            NERR = NERR+1
          END IF
          IF (TBOD.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'TBOD',CONFN
            NERR = NERR+1
          END IF
          IF (RBOD.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'RBOD',CONFN
            NERR = NERR+1
          END IF
        END IF

******* Phosphorus

        IF (PHOSPHORUS) THEN
          WRITE (*,*) '  phosphorus'
          IF (PO4R.GT.0.05) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'PO4R',0.05,CONFN
            NWRN = NWRN+1
          ELSE IF (PO4R.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'PO4R',CONFN
            NERR = NERR+1
          END IF
          IF (PARTP.GT.2.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'PARTP',2.0,CONFN
            NWRN = NWRN+1
          ELSE IF (PARTP.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'PARTP',CONFN
            NERR = NERR+1
          END IF
          IF (AHSP.GT.0.05) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'AHSP',0.05,CONFN
            NWRN = NWRN+1
          ELSE IF (AHSP.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'AHSP',CONFN
            NERR = NERR+1
          END IF
        END IF

******* Ammonium

        IF (AMMONIUM) THEN
          WRITE (*,*) '  ammonium'
          IF (NH4R.GT.0.5) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'NH4R',0.5,CONFN
            NWRN = NWRN+1
          ELSE IF (NH4R.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'NH4R',CONFN
            NERR = NERR+1
          END IF
          IF (NH4DK.GT.0.5) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'NH4DK',0.5,CONFN
            NWRN = NWRN+1
          ELSE IF (NH4DK.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'NH4DK',CONFN
            NERR = NERR+1
          END IF
          IF (PARTN.GT.2.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'PARTN',2.0,CONFN
            NWRN = NWRN+1
          ELSE IF (PARTN.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'PARTN',CONFN
            NERR = NERR+1
          END IF
          IF (AHSN.GT.0.05) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'AHSN',0.05,CONFN
            NWRN = NWRN+1
          ELSE IF (AHSN.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'AHSN',CONFN
            NERR = NERR+1
          END IF

********* Ammonium temperature rate multipliers

          IF (NH4T1.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'NH4T1',CONFN
            NERR = NERR+1
          END IF
          IF (NH4T2.GT.40.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'NH4T2',40.0,CONFN
            NWRN = NWRN+1
          ELSE IF (NH4T2.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'NH4T2',CONFN
            NERR = NERR+1
          END IF
          IF (NH4K1.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'NH4K1',CONFN
            NERR = NERR+1
          END IF
          IF (NH4K2.GT.1.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9150) 'NH4K2',1.0,CONFN
            NERR = NERR+1
          ELSE IF (NH4K2.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'NH4K2',CONFN
            NERR = NERR+1
          END IF
        END IF

******* Nitrate

        IF (NITRATE) THEN
          WRITE (*,*) '  nitrate-nitrite'
          IF (NO3DK.GT.0.5) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'NO3DK',0.5,CONFN
            NWRN = NWRN+1
          ELSE IF (NO3DK.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'NO3DK',CONFN
            NERR = NERR+1
          END IF

********* Nitrate rate multipliers

          IF (NO3T1.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'NO3T1',CONFN
            NERR = NERR+1
          END IF
          IF (NO3T2.GT.40.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'NO3T2',4.0,CONFN
            NWRN = NWRN+1
          ELSE IF (NO3T2.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'NO3T2',CONFN
            NERR = NERR+1
          END IF
          IF (NO3K1.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'NO3K1',CONFN
            NERR = NERR+1
          END IF
          IF (NO3K2.GT.1.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9150) 'NO3K2',1.0,CONFN
            NERR = NERR+1
          ELSE IF (NO3K2.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'NO3K2',CONFN
            NERR = NERR+1
          END IF
        END IF

******* Sediment carbon dioxide release

        IF (TOTIC) THEN
          WRITE (*,*) '  total inorganic carbon'
          IF (CO2R.GT.0.5) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'CO2R',0.5,CONFN
            NWRN = NWRN+1
          ELSE IF (CO2R.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'CO2R',CONFN
            NERR = NERR+1
          END IF
        END IF

******* Iron

        IF (IRON) THEN
          WRITE (*,*) '  iron'
          IF (FER.GT.2.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'FER',2.0,CONFN
            NWRN = NWRN+1
          ELSE IF (FER.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'FER',CONFN
            NERR = NERR+1
          END IF
          IF (FES.GT.10.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'FES',10.0,CONFN
            NWRN = NWRN+1
          ELSE IF (FES.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'FES',CONFN
            NERR = NERR+1
          END IF
        END IF

******* Stoichiometry

        WRITE (*,*) '  stoichiometry'
        IF (O2NH4.LT.0.0) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9155) 'O2NH4',CONFN
          NERR = NERR+1
        ELSE IF (O2NH4.NE.4.57) THEN
          DELETE_WRN = .FALSE.
          WRITE(WRN,9152) 'O2NH4',4.57,CONFN
          NWRN = NWRN+1
        END IF
        IF (O2OM.GT.2.0) THEN
          DELETE_WRN = .FALSE.
          WRITE(WRN,9150) 'O2OM',2.0,CONFN
          NWRN = NWRN+1
        ELSE IF (O2OM.LT.0.0) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9155) 'O2OM',CONFN
          NERR = NERR+1
        END IF
        IF (O2AR.GT.2.0) THEN
          DELETE_WRN = .FALSE.
          WRITE(WRN,9150) 'O2AR',2.0,CONFN
          NWRN = NWRN+1
        ELSE IF (O2AR.LT.0.0) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9155) 'O2AR',CONFN
          NERR = NERR+1
        END IF
        IF (O2AG.LT.0.0) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9155) 'O2AG',CONFN
          NERR = NERR+1
        END IF
        IF (BIOP.LT.0.0) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9155) 'BIOP',CONFN
          NERR = NERR+1
        ELSE IF (BIOP.GT.0.011) THEN
          DELETE_WRN = .FALSE.
          WRITE(WRN,9152) 'BIOP',0.011,CONFN
          NWRN = NWRN+1
        END IF
        IF (BION.LT.0.0) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9155) 'BION',CONFN
          NERR = NERR+1
        ELSE IF (BION.GT.0.08) THEN
          DELETE_WRN = .FALSE.
          WRITE(WRN,9152) 'BION',0.08,CONFN
          NWRN = NWRN+1
        END IF
        IF (BIOC.LT.0.0) THEN
          DELETE_ERR = .FALSE.
          WRITE(ERR,9155) 'BIOC',CONFN
          NERR = NERR+1
        ELSE IF (BIOC.GT.0.45) THEN
          DELETE_WRN = .FALSE.
          WRITE(WRN,9152) 'BIOC',0.45,CONFN
          NWRN = NWRN+1
        END IF

******* Oxygen limit

        IF (DISSOLVED_OXYGEN) THEN
          WRITE (*,*) '  oxygen limit'
          IF (O2LIM.GT.1.0) THEN
            DELETE_WRN = .FALSE.
            WRITE(WRN,9150) 'O2LIM',1.0,CONFN
            NWRN = NWRN+1
          ELSE IF (O2LIM.LT.0.0) THEN
            DELETE_ERR = .FALSE.
            WRITE(ERR,9155) 'O2LIM',CONFN
            NERR = NERR+1
          END IF
        END IF
      END IF

***** Input FORMAT statements

 1000 FORMAT(/A8/(8X,A72))
 1005 FORMAT(//A8/(8X,A72))
 1006 FORMAT(/A8/8X,2F8.0,I8)
 1010 FORMAT(/A8/(8X,9F8.0))
 1011 FORMAT(:8X,9F8.0)
 1020 FORMAT(/A8/8X,2F8.0,3X,A5,6(5X,A3))
 1030 FORMAT(/A8/(8X,9(5X,A3)))
 1035 FORMAT(/A8/13X,A3,2A8,6F8.0)
 1040 FORMAT(/A8/(8X,4I8))
 1050 FORMAT(//(10F8.0))
 1060 FORMAT(/A8/(8X,9I8))
 1063 FORMAT(/A8/(8X,9A8))
 1065 FORMAT(:8X,9I8)
 1070 FORMAT(/A8/13X,A3,8(5X,A3))
 1080 FORMAT(/A8/13X,A3,I8,5X,A3)
 1085 FORMAT(/A8/8X,3(5X,A3),I8)
 1090 FORMAT(/A8/13X,A3,8I8)
 1100 FORMAT(/A8/8X,I8,8F8.0)
 1110 FORMAT(/A8/8X,A8,F8.0)
 1115 FORMAT(/A8/(8X,9A8))
 1116 FORMAT(:8X,9A8)
 1120 FORMAT(10F8.0:/(9F8.0))

***** Error FORMAT statements

 9000 FORMAT('Either illegal value or incorrect card in the',
     .       ' following set '/3(A80/))
 9010 FORMAT('Incorrect card in the following set '/3(A80/))
 9011 FORMAT('SLICEC is not set to either "DETAIL" or "SIMPLE"')
 9012 FORMAT('SLHTC is not set to either "TERM" or "ET"')
 9013 FORMAT ('BETAI is greater than 0.9')
 9014 FORMAT ('ICEMIN is greater than 0.1')
 9015 FORMAT('ICET2 is greater than 4.0')
 9016 FORMAT('SLTRC must be set to either "QUICKEST" or "UPWIND"')
 9017 FORMAT('THETA must be positive')
 9018 FORMAT('THETA must be less than or equal to 1.0')
 9020 FORMAT(A80)
 9025 FORMAT('DLTMAX(',I2,') is set to zero')
 9026 FORMAT('DLTF(',I2,') is set to zero')
 9029 FORMAT('TMEND (=',F8.2,') is less than TMSTRT (=',F8.2,')')
 9030 FORMAT('[DLTD(1)] (=',F8.2,' must be less than [TMSTRT] (=',
     .         F8.1,')')
 9032 FORMAT('DLTD(',I2,') (=',F8.1,'), must be greater than ',
     .       'DLTD(',I2,') (=',F8.1,')')
 9034 FORMAT('NWSD must be greater than 0')
 9035 FORMAT('WSCD(1) = ',F8.1, ', must be less than',
     .         ' TMSTRT (=',F8.1)
 9036 FORMAT('WSCD(',I2,') (=',F8.1,'), must be greater than ',
     .       'WSCD(',I2,') (=',F8.1,')')
 9040 FORMAT('The file, ',A12,' must exist because T2I or C2I = -1.0')
 9045 FORMAT('The file, ',A12,' must exist because T2I or C2I = -2.0')
 9050 FORMAT('The file, ',A12,' must exist because PRC is ON')
 9060 FORMAT('B(',I2,',',I2,') = ',F8.1,' for branch ',I2,' down',
     .       'stream boundary segment')
 9061 FORMAT('B(',I2,',',I2,') is greater than B(',I2,',',I2,')')
 9070 FORMAT('Branch geometry disagrees for upstream segment on branch',
     .        I3)
 9090 FORMAT('Start date in ', A12,', ',F8.3, ', is ',
     .       'greater than TMSTRT, ', F8.3, '.' )
 9095 FORMAT('Julian date in ',A12,', ',F8.3,' is <= previous ',
     .       'date ',F8.3)
 9096 FORMAT('End date in ',A12,', ',F8.3,', is less than TMEND (',
     .        F8.3, ').'   )
 9098 FORMAT(A6,' must be either '' SPECIFY'' or ''   DISTR'' or',
     .         ''' DENSITY'' in ',A12)
 9099 FORMAT(A6,' must be either '' ON'' or ''OFF'' in ',A12)
 9100 FORMAT(A6,' is outside its expected range in ',A12,
     .       ' for JDAY ',F8.3) 
 9101 FORMAT('Possible missing data values included in '/
     .       '  meteorologic file for Julian day',F8.3)
 9105 FORMAT(A6,' is outside its expected range in ',A12) 
 9106 FORMAT(A6,'(',I2,',',I2,') = ',I2,' is below the bottom active',
     .       'cell') 
 9107 FORMAT(A6,'(',I2,') is outside its expected range in ',A12) 
 9108 FORMAT(A6,'(',I2,',',I2,') = ',F8.2,' is outside its expected ',
     .       'range in ',A12) 
 9109 FORMAT(A6,'(',I2,',',I2,') = ',I2,' is outside its expected ',
     .       'range in ',A12) 
 9110 FORMAT('You specifed that some ',A8,' constituents',/
     .       '    be turned on, but the file ',A12,' does not exist.')
 9115 FORMAT('Segment ',I2,' is a boundary segment in the SEG PRINT ',
     .         'card')
 9117 FORMAT(A6,' segment ',I2,
     . ' is a boundary segment in ',A12)
 9120 FORMAT('Not enough values for ',A6,' in ',A12)
 9125 FORMAT(A6,'(',I2,') must be greater than ',A6,'(',I2,') in ',A12)
 9130 FORMAT('Recommend ',A6,' not be greater than ',F8.3,
     .       ' in ',A12)
 9135 FORMAT('You specified ',A6,', as '' ON'', but all values',
     .       ' in ',A6,' are ''OFF'' ',/,'         in ',A12)
 9142 FORMAT(A6,'(',I2,') should have a value other than ',
     .F4.1,' in ',A12)
 9145 FORMAT(A6,'(',I2,') is '' ON'', so ',A6,'(',I2,
     .       ') should be '' ON''')
 9150 FORMAT(A6,' is greater than upper limit of ',F8.3,
     .       ' in ',A12)
 9151 FORMAT(A6,' is less than lower limit of ',F8.3,
     .       ' in ',A12)
 9152 FORMAT(A6,' should be less than ',F8.3,' in ',A12)
 9153 FORMAT(A6,'(',I2,') is greater than upper limit of ',
     .       F8.3,' in ',A12)
 9155 FORMAT(A6,' must be non-negative in ',A12)
 9160 FORMAT('Could not open file ',A12,' specified as ',A6,' in ',A12)
 9165 FORMAT('Filename, ',A18,' in ',A12,
     .       ' does not conform to',/,'   DOS 8.3 character standard.')
 9500 FORMAT(//)
 9520 FORMAT(10F8.0/8X,9F8.0)

************************************************************************
**           Task 4:  Calculate volume-area-elevation table           **
************************************************************************

***** Blank printout grid

      DO I=1,NISNP
        DO K=1,KMP
          CONV(K,I) = BLANK
        END DO
      END DO

***** Grid geometry

      VOLG = 0.
      DO JB=1,NBP
        IU = US(JB)
        ID = DS(JB)

******* Layer elevations

        EL(KMP) = DATUM                                                 !4/06/93
        DO K=KMP-1,1,-1                                                 !4/06/93
          EL(K) = EL(K+1)+H(K)                                          !4/06/93
        END DO                                                          !4/06/93

******* Segment locations

        XLOC(IU) = DLX(IU)*0.5
        DO I=IU+1,ID
          XLOC(I) = XLOC(I-1)+(DLX(I-1)+DLX(I))*0.5
        END DO

******* Water surface and bottom layers

        DO I=IU-1,ID+1
          KTI(I) = 2                                                    !6/15/94
          DO WHILE (EL(KTI(I)).GT.ELWS(I))                              !6/15/94
            KTI(I) = KTI(I)+1                                           !6/15/94
          END DO                                                        !6/15/94
          KT = MAX(KT,KTI(I))                                           !6/15/94
        END DO
        DO I=IU-1,ID+1
          Z(I) = EL(KT)-ELWS(I)                                         !6/15/94
        END DO

******* Bottom layers

        DO I=IU,ID
          K = 2
          DO WHILE (B(K,I).GT.0.)
            KB(I) = K
            K     = K+1
          END DO
          KBMAX = MAX(KBMAX,KB(I))
        END DO

******* Upstream active cell

        I = IU
        DO WHILE (KB(I)-KT.LT.1)
          I = I+1
        END DO
        IUC     = I
        CUS(JB) = IUC

******* Branch and grid total volume

        DO I=IU,ID
          VOLB(JB) = VOLB(JB)+DLX(I)*B(KT,I)*(H(KT)-Z(I))
          VOLG     = VOLG+DLX(I)*B(KT,I)*(H(KT)-Z(I))                   !9/06/94
          DO K=KT+1,KB(I)
            VOLB(JB) = VOLB(JB)+DLX(I)*B(K,I)*H(K)
            VOLG     = VOLG+DLX(I)*B(K,I)*H(K)                          !9/06/94
          END DO
        END DO

******* Branch and grid area and volume by layer

        DO K=KMP-1,2,-1
          NCCLB(K,JB) = NCCLB(K+1,JB)
          CVLB(K,JB)  = CVLB(K+1,JB)
          DO I=IU,ID
            IF (K.LE.KB(I)) THEN
              ALB(K,JB)   = ALB(K,JB)+B(K,I)*DLX(I)
              NCCLB(K,JB) = NCCLB(K,JB)+1
              CVLB(K,JB)  = CVLB(K,JB)+B(K,I)*H(K)*DLX(I)               !9/06/94
            END IF
          END DO
          GAL(K)   = GAL(K)+ALB(K,JB)
          NCCLG(K) = NCCLG(K)+NCCLB(K,JB)
          CVLG(K)  = CVLG(K)+CVLB(K,JB)
        END DO
      END DO

***** Beginning and ending segment and bottom layer for snapshots

      DO I=1,NISNP
        IPR(I) = IPRI(I)
      END DO
      DO I=1,NISNP
        KEPR = MAX(KB(IPR(I)),KEPR)
      END DO

***** Outflows

      DO JB=1,NBP
        IF (DN_FLOW(JB)) THEN
          IF (SWC(JB).EQ.' ON') THEN
            WRITE (*,*) 'Selective withdrawal structures'
            ID = DS(JB)
            IF (NSTR(JB).LE.0) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9105) 'NSTR',CONFN 
              NERR = NERR+1
            END IF
            DO JS=1,NSTR(JB)
              IF (KB(ID).LT.KBSW(JS,JB)) THEN
                DELETE_WRN = .FALSE.
                WRITE(WRN,9109) 'KBSW',JS,JB,KBSW(JS,JB),CONFN
                NWRN = NWRN+1
              END IF
              IF (PSC(JS,JB).NE.'    LINE'.AND.
     .            PSC(JS,JB).NE.'   POINT') THEN
                DELETE_ERR = .FALSE.
                WRITE(ERR,9108) 'SINK',JS,JB,CONFN
                NERR = NERR+1
              END IF
              IF (ESTR(JS,JB).LT.EL(KB(ID)+1)) THEN
                DELETE_WRN = .FALSE.
                WRITE(WRN,9108) 'ESTR',JS,JB,ESTR(JS,JB),CONFN
                NWRN = NWRN+1
              END IF
              IF (PSC(JS,JB).EQ.'    LINE') THEN
                DO K=2,KB(ID)
                  IF (EL(K).GT.ESTR(JS,JB)) KSTR = K
                END DO
                IF (WSTR(JS,JB).LE.0.0.
     .              OR.WSTR(JS,JB).GT.B(KSTR,ID)) THEN
                  DELETE_ERR = .FALSE.
                  WRITE(ERR,9108) 'WSTR',JS,JB,WSTR(JS,JB),CONFN
                  NERR = NERR+1
                END IF
             END IF
            END DO
          ELSE
            WRITE (*,*) 'Outlet structures for branch',JB
            IF (NOUT(JB).LE.0) THEN
              DELETE_ERR = .FALSE.
              WRITE(ERR,9107) 'NOUT',JB,CONFN
              NERR = NERR+1
              DO JO=1,NOUT(JB)
                IF (KOUT(JO,JB).LE.1.OR.KOUT(JO,JB).GT.KB(DS(JB))) THEN
                  DELETE_ERR = .FALSE.
                  WRITE(ERR,9106) 'KOUT',JO,JB,KOUT(JO,JB)
                  NERR = NERR+1
                END IF
              END DO
            END IF
          END IF
        END IF
      END DO

************************************************************************
**             Task 5:  Initial Conditions for Simulation             **
************************************************************************


***** Temperature and constituents

      IF (OPEN_VPR) THEN
        WRITE (*,*) 'Reading vertical profile file'
        OPEN (VPR,FILE=VPRFN,STATUS='OLD',IOSTAT=IERR)
        IF (IERR.EQ.0) THEN
          READ (VPR,*)
          IF (VERT_TEMP) READ (VPR,1010) AID,(TVP(K),K=KT,KBMAX)
          IF (CONSTITUENTS) THEN
            DO JC=1,NCP
              IF (VERT_CONC(JC)) READ (VPR,1010) AID,(CVP(K,JC),
     .                                           K=KT,KBMAX)
            END DO
          ENDIF
        ELSE
          DELETE_ERR = .FALSE.
          WRITE(ERR,9040) VPRFN
          NERR = NERR+1
        END IF
      ENDIF

***** Longitudinal/vertical initial profiles

      IF (OPEN_LPR) THEN
        OPEN (LPR,FILE=LPRFN,STATUS='OLD',IOSTAT=IERR)
        IF (IERR.EQ.0) THEN
          READ (LPR,*)
        ELSE
          DELETE_ERR = .FALSE.
          WRITE(ERR,9045) LPRFN
          NERR = NERR+1
        END IF
      END IF
      DO JB=1,NBP
        IU = CUS(JB)
        ID = DS(JB)

******* Temperature

        DO I=IU,ID
          IF (LONG_TEMP) READ (LPR,1010) AID,(T2(K,I),K=KT,KB(I))
          DO K=KT,KB(I)
            IF (ISO_TEMP)  T2(K,I) = T2I
            IF (VERT_TEMP) T2(K,I) = TVP(K)
          END DO
        END DO

******* Constituents

        DO JC=1,NAC
          DO I=IU,ID
            IF (LONG_CONC(CN(JC))) READ (LPR,1010) AID,(C2(K,I,CN(JC)),
     .                                             K=KT,KB(I))
            DO K=KT,KB(I)
              IF (ISO_CONC(CN(JC)))  C2(K,I,CN(JC)) = C2I(CN(JC))
              IF (VERT_CONC(CN(JC))) C2(K,I,CN(JC)) = CVP(K,CN(JC))
            END DO
          END DO
        END DO
      END DO

************************************************************************
**                      Task 6:  Output section                       **
************************************************************************

***** Open output files

       OPEN (INI,FILE=INIFN,STATUS='UNKNOWN')

***** Initial input and conditions

      WRITE (*,*) 'Preprocessor output'
      IF (LASERJET_II) THEN
        WRITE (INI,'(''+'',A80)') ESC//'E'//ESC//'(s16.66H'//ESC//'(10U'
     .                            //ESC//'&a8L'//ESC//'&l7E'
      ELSE IF (LASERJET_III) THEN
        WRITE (INI,'(''+'',A80)') ESC//'E'//ESC//'&l6.0C'//ESC//
     .                            '(s0p16.67h8.5v0s0b0T'//ESC//'(10U'
     .                            //ESC//'&a8L'//ESC//'&l7E'
      ELSE IF (LASERJET_IV) THEN
        WRITE (INI,'(A80)') ESC//'E'//ESC//'&l6.0c7E'//ESC//
     .                      '(s0p16.67h8.5v0s0b0T'//ESC//'(10U'//ESC//
     .                      '&a8L'
      END IF
      WRITE (INI,2000)
      WRITE (INI,2010) (TITLE(I),I=1,6)
      WRITE (INI,2020) 'Time Control',TMSTRT,TMEND,YEAR
      WRITE (INI,2030)  NDLT,DLTMIN
      WRITE (INI,2040) (DLTD(J),J=1,NDLT)
      WRITE (INI,2050) (DLTMAX(J),J=1,NDLT)
      WRITE (INI,2060) (DLTF(J),J=1,NDLT)
      WRITE (INI,2070) 'Initial Conditions'
      IF (ISO_TEMP)     WRITE (INI,2080) 'Temperature',T2I,DEG
      IF (VERT_TEMP)    WRITE (INI,2090) 'Temperature'
      WRITE (INI,2100) 'Water type',WTYPEC
      WRITE (INI,2110) 'Ice thickness',ICETHI
      WRITE (INI,2120) 'Calculations',EVC,PRC,VBC,EBC,MBC,WBC,PQC
      WRITE (INI,2130)  WINDC,QINC,QOUTC,HEATC
      WRITE (INI,2135) 'Input Interpolations',INFIC,TRIC,DTRIC,
     .                  HDIC,OUTIC,WDIC,METIC
      WRITE (INI,2140) 'Meteorological Parameters'
      WRITE (INI,2150) 'Location',LAT,LONG
      WRITE (INI,2160) (WSCD(I),I=1,NWSD)
      WRITE (INI,2170) (WSC(I),I=1,NWSD)
      IEGR = 1
      WRITE (INI,2175)
      DO WHILE (IEGR.LT.IMP-1)
        IBGR = IEGR+1
        IEGR = IEGR+19
        IF (IEGR.GT.IMP-1) IEGR = IMP-1
        WRITE (INI,2176) (I,I=IBGR,IEGR)
        WRITE (INI,2180) (PHI0(I),I=IBGR,IEGR)
      END DO
      WRITE (INI,2190) 'Transport Solution',SLTRC,THETA
      WRITE (INI,2200) 'Hydraulic Coefficients',AX,DXI,CHEZY,TSED,DEG,
     .                  CBHE,DEG
      WRITE (INI,2210) 'Ice cover',ICEC,SLICEC,SLHTC,ALBEDO,HWI,BETAI,
     .                  GAMMAI
      WRITE (INI,2220) 'Output Control',LJPC,SNPC,(HNAME1(I),HPRN(I),
     .                  I=1,3)
      WRITE (INI,2230)  NSNP
      WRITE (INI,2240) (SNPD(J),J=1,NSNP)
      WRITE (INI,2250) (SNPF(J),J=1,NSNP)
      WRITE (INI,2251)  SCRC,NSCR
      WRITE (INI,2252) (SCRD(J),J=1,NSCR)
      WRITE (INI,2253) (SCRF(J),J=1,NSCR)
      WRITE (INI,2260)  TSRC
      WRITE (INI,2270)  NTSR
      WRITE (INI,2280) (TSRD(J),J=1,NTSR)
      WRITE (INI,2290) (TSRF(J),J=1,NTSR)
      WRITE (INI,2300)  VPLC
      WRITE (INI,2310)  NVPL
      WRITE (INI,2320) (VPLD(J),J=1,NVPL)
      WRITE (INI,2330) (VPLF(J),J=1,NVPL)
      WRITE (INI,2340)  PRFC
      WRITE (INI,2350)  NPRF
      WRITE (INI,2360)  NIPRF
      WRITE (INI,2370) (IPRF(I),I=1,NIPRF)
      WRITE (INI,2380) (PRFD(J),J=1,NPRF)
      WRITE (INI,2390) (PRFF(J),J=1,NPRF)
      WRITE (INI,2391)  SPRC
      WRITE (INI,2392)  NSPR
      WRITE (INI,2393)  NISPR
      WRITE (INI,2394) (ISPR(I),I=1,NISPR)
      WRITE (INI,2395) (SPRD(J),J=1,NSPR)
      WRITE (INI,2396) (SPRF(J),J=1,NSPR)
      WRITE (INI,2400)  CPLC
      WRITE (INI,2410)  NCPL
      WRITE (INI,2420) (CPLD(J),J=1,NCPL)
      WRITE (INI,2430) (CPLF(J),J=1,NCPL)
      WRITE (INI,2440)  RSOC,RSIC
      WRITE (INI,2450)  NRSO
      WRITE (INI,2460) (RSOD(J),J=1,NRSO)
      WRITE (INI,2470) (RSOF(J),J=1,NRSO)
      WRITE (INI,2480) (JB,SWC(JB),NSTR(JB),JB=1,NB)
      DO JB=1,NBP
        WRITE (INI,2490) (JS,PSC(JS,JB),WSTR(JS,JB),ESTR(JS,JB),
     .                    KBSW(JS,JB),JS=1,NSTR(JB))
      END DO
      WRITE (INI,2500) (NOUT(JB),JB=1,NB)
      DO JB=1,NBP
        IF (NOUT(JB).GT.0) WRITE (INI,2510) JB,(KOUT(K,JB),K=1,NOUT(JB))
      END DO
      WRITE (INI,2520)  NWD,(IWD(JW),JW=1,NWD)
      WRITE (INI,2530) (KWD(JW),     JW=1,NWD)
      WRITE (INI,2540)  NTR,(ITR(JT),JT=1,NTR)
      WRITE (INI,2544) (TRPC(JT),    JT=1,NTR)
      WRITE (INI,2545) (ELTRT(JT),   JT=1,NTR)
      WRITE (INI,2546) (ELTRB(JT),   JT=1,NTR)
      WRITE (INI,2550)
      WRITE (INI,2560) (JB,DTRC(JB), JB=1,NB)
      WRITE (INI,2570)  CONFN,BTHFN,VPRFN,LPRFN,RSIFN,METFN,ELOFN,QWDFN
      DO JB=1,NBP
        IF (MOD(JB-1,3).EQ.0.AND.JB.GT.1) THEN
          WRITE (INI,2580) JB,QINFN(JB),TINFN(JB),CINFN(JB),QOTFN(JB),
     .                     QDTFN(JB),TDTFN(JB),CDTFN(JB),PREFN(JB),
     .                     TPRFN(JB),CPRFN(JB),EUHFN(JB),TUHFN(JB),
     .                     CUHFN(JB),EDHFN(JB),TDHFN(JB),CDHFN(JB)
        ELSE
          WRITE (INI,2585) JB,QINFN(JB),TINFN(JB),CINFN(JB),QOTFN(JB),
     .                     QDTFN(JB),TDTFN(JB),CDTFN(JB),PREFN(JB),
     .                     TPRFN(JB),CPRFN(JB),EUHFN(JB),TUHFN(JB),
     .                     CUHFN(JB),EDHFN(JB),TDHFN(JB),CDHFN(JB)
        END IF
      END DO
      WRITE (INI,2590) (JT,QTRFN(JT),TTRFN(JT),CTRFN(JT),JT=1,NTR)
      WRITE (INI,2600)  ERRFN,WRNFN,SNPFN,TSRFN,PRFFN,VPLFN,CPLFN
      IF (CONSTITUENTS) THEN
        WRITE (INI,2610) 'Constituents',CCC
        WRITE (INI,2620)  FREQUK
        WRITE (INI,2625)
        WRITE (INI,2630) (CNAME1(JC),ACC(JC),C2I(JC),CPRN(JC),INACC(JC),
     .                   TRACC(JC),DTACC(JC),PRACC(JC),JC=1,NCP),
     .                   'Limiting factor ',LIMC
        WRITE (INI,2640) 'Constituent Rates'
        WRITE (INI,2650)  SSS
        WRITE (INI,2660)  COLDK
        WRITE (INI,2670)  LDOMDK,LRDDK
        WRITE (INI,2680)  RDOMDK
        WRITE (INI,2690)  AG,AM,AE,AR,AS
        WRITE (INI,2700)  POMDK,POMS
        WRITE (INI,2710)  PO4R
        WRITE (INI,2720)  NH4DK,NH4R
        WRITE (INI,2730)  NO3DK
        WRITE (INI,2740)  SDK
        WRITE (INI,2750)  FES,FER
        WRITE (INI,2760)  SOD
        WRITE (INI,2764)  FSOD
        WRITE (INI,2765)  SDC
        WRITE (INI,2770) 'Lower Temperature Bounds',DEG,DEG
        WRITE (INI,2780)  AT1,AT2
        WRITE (INI,2790)  NH4T1,NH4T2
        WRITE (INI,2800)  NO3T1,NO3T2
        WRITE (INI,2810)  OMT1,OMT2
        WRITE (INI,2820) 'Upper Temperature Bounds',DEG,DEG
        WRITE (INI,2830)  AT3,AT4
        WRITE (INI,2840) 'Stoichiometric Equivalence',O2NH4,O2OM,
     .                    O2AR,O2AG
        WRITE (INI,2850)  BIOP,BION,BIOC
        WRITE (INI,2860)  KBOD,TBOD,RBOD
        WRITE (INI,2870) 'Half Saturation',AHSP,AHSN
        WRITE (INI,2880) 'Light',BETA,EXH2O,EXSS,EXOM,ASAT
        WRITE (INI,2890) 'Diffusion',DMO2,DMCO2
        WRITE (INI,2900) 'Partitioning Coefficients',PARTP,PARTN
        WRITE (INI,2910) 'Miscellaneous Constants',O2LIM,COLQ10,APOM,
     .                    CO2R
      END IF
      WRITE (INI,2920) 'Geometry',IMP,KMP,NBP,DATUM,KT
      DO WHILE (KEGR.LT.KMP)
        KBGR = KEGR+1
        KEGR = KEGR+19
        IF (KEGR.GT.KMP) KEGR = KMP
        WRITE (INI,2930) (K,K=KBGR,KEGR)
        WRITE (INI,2940) (H(K),K=KBGR,KEGR)
      END DO
      DO JB=1,NBP
        WRITE (INI,2950) JB,US(JB),DS(JB),UHS(JB),DHS(JB)
      END DO

***** Volume-area-elevation table

      IF (NB.GT.1) THEN
        NLINES = KMP+27
        DO JB=1,NBP
          WRITE (*,*) 'Volume-area-elevation table for branch',JB
          IF (NLINES.LT.73) THEN
            WRITE (INI,3010) VOLB(JB)
          ELSE
            WRITE (INI,3015) VOLB(JB)
            NLINES = KMP+3
          END IF
          WRITE (INI,3020) JB
          WRITE (INI,3030)
          DO K=2,KMP-1
            WRITE (INI,3040) K,EL(K),ALB(K,JB)/1.E6,CVLB(K,JB)/1.E6,
     .                       NCCLB(K,JB)
            IF (K.EQ.KT) WRITE (INI,3050)
          END DO
          NLINES = NLINES+KMP+3
        END DO
      END IF
      IF (NLINES.LT.73) WRITE (INI,3060) VOLG
      IF (NLINES.GE.73) WRITE (INI,3065) VOLG
      WRITE (INI,3070)
      WRITE (INI,3030)
      write(50,3035)
      vfctr=(3.2808**3)/(43650.e6)
      DO K=2,KMP-1
        WRITE (INI,3040) K,EL(K),GAL(K)/1.E6,CVLG(K)/1.E6,NCCLG(K)
        write(50,3045) k,el(k),el(k)*3.2808,cvlg(k)/1.e6,cvlg(k)*vfctr
        IF (K.EQ.KT) WRITE (INI,3050)
      END DO
      ziltch=0.0
      write(50,3040) kmp,ziltch,ziltch,ziltch

***** Bathymetry

      IEGR = 1
      NLINES = KMP+4
      IF (LASERJET_II) THEN
        WRITE (INI,'(A50)') ESC//'&l1o4.8C'//ESC//'&a8L'
      ELSE IF (LASERJET_III.OR.LASERJET_IV) THEN
        WRITE (INI,'(A50)') ESC//'&l1o4e4.8C'//ESC//'&a8L'
      END IF
      WRITE (INI,2960)
      DO WHILE (IEGR.LT.IMP-1)
        IBGR = IEGR+1
        IEGR = IEGR+21
        IF (IEGR.GT.IMP-1) IEGR = IMP-1
        IF (NLINES.LT.72) THEN
          WRITE (INI,2970) (I,I=IBGR,IEGR)
        ELSE
          WRITE (INI,2972) (I,I=IBGR,IEGR)
          NLINES = KMP+4
        END IF
        WRITE (INI,2975) (INT(XLOC(I)),I=IBGR,IEGR)
        DO K=1,KMP
          DO I=IBGR,IEGR
            WRITE (BK(I),2980) B(K,I)
            IF (BK(I).EQ.ZERO) BK(I) = BLANK
          END DO
          WRITE (INI,2990) K,(BK(I),I=IBGR,IEGR)
          IF (K.EQ.KT) WRITE (INI,3000)
        END DO
        NLINES = NLINES+KMP+4
      END DO
      WRITE (*,*) 'Initial conditions'
      CALL PRINT_INITIAL_GRID (KT)
      CLOSE (INI)

***** Close error/warning files

      IF (DELETE_ERR) THEN
        CLOSE (ERR,STATUS='DELETE')
        WRITE (*,*)
        WRITE (*,*) NERR,' errors'
      ELSE 
        CLOSE (ERR)
        WRITE (*,*)
        WRITE (*,*) NERR,' errors - see "pre.err"',
     .              ' for a more detailed description'
      END IF
      IF (DELETE_WRN) THEN
        CLOSE (WRN,STATUS='DELETE')
        WRITE (*,*) NWRN,' warnings'
      ELSE
        CLOSE (WRN)
        WRITE (*,*) NWRN,' warnings - see "pre.wrn"',
     .              ' for a more detailed description'
      END IF

***** Output formats

 2000 FORMAT(1X,'CE-QUAL-W2 - V2.10 - June, 1995'/)
 2010 FORMAT(1X,A72)
 2020 FORMAT(/1X,A12/
     .       '+',4('_'),1X,7('_')//
     .       3X,'Starting time (Julian day) [TMSTRT] =',F8.2/
     .       3X,'Ending time (Julian day)    [TMEND] =',F8.2/
     .       3X,'Year                         [YEAR] =',I8)
 2030 FORMAT(3X,'# Timestep intervals         [NDLT] =',I8/
     .       3X,'Minimum timestep (sec)     [DLTMIN] =',F8.1)
 2040 FORMAT(3X,'Timestep day (Julian day)    [DLTD] =',7F8.1)
 2050 FORMAT(3X,'Maximum timestep (sec)     [DLTMAX] =',7F8.1,' sec')
 2060 FORMAT(3X,'Fraction of timestep         [DLTF] =',7F8.2)
 2070 FORMAT(/1X,A18/
     .       '+',7('_'),1X,10('_')/)
 2080 FORMAT(3X,A11,6X,'[T2I] = ',F5.1,1X,A2)
 2090 FORMAT(3X,A11,6X,'[T2I] = Downstream vertical profile')
 2100 FORMAT(3X,A10,4X,'[WTYPEC] = ',A5,' water')
 2110 FORMAT(3X,A13,' [ICETHI] = ',F5.1,' m'/)
 2120 FORMAT(3X,A12/
     .       '+',2X,12('_')//
     .       5X,'Evaporation     [EVC] = ',A3/
     .       5X,'Precipitation   [PRC] = ',A3/
     .       5X,'Volume balance  [VBC] = ',A3/
     .       5X,'Energy balance  [EBC] = ',A3/
     .       5X,'Mass balance    [MBC] = ',A3/
     .       5X,'Water balance   [WBC] = ',A3/
     .       5X,'Place inflows   [PQC] = ',A3)
 2130 FORMAT(5X,'Wind          [WINDC] = ',A3/
     .       5X,'Inflow         [QINC] = ',A3/
     .       5X,'Outflow       [QOUTC] = ',A3/
     .       5X,'Heat exchange [HEATC] = ',A3/)
 2135 FORMAT(3X,A20/
     .       '+',2X,5('_'),1X,14('_')//
     .       5X,'Inflows                 [QINIC] = ',A3/
     .       5X,'Tributaries              [TRIC] = ',A3/
     .       5X,'Distributed tributaries [DTRIC] = ',A3/
     .       5X,'Head boundaries          [HDIC] = ',A3/
     .       5X,'Outflows               [QOUTIC] = ',A3/
     .       5X,'Withdrawals              [WDIC] = ',A3/
     .       5X,'Meteorologic data       [METIC] = ',A3/)
 2140 FORMAT(3X,A24/
     .       '+',2X,13('_'),1X,10('_')/)
 2150 FORMAT(5X,A8/
     .       7X,'Latitude   [LAT] =',F8.2,'ø'/
     .       7X,'Longitude [LONG] =',F8.2,'ø'/)
 2160 FORMAT('+':/
     .       5X,'Wind shading date (Julian day) [WSCD] =',6F7.2)
 2170 FORMAT('+':/
     .       5X,'Wind shading coefficient        [WSC] =',6F7.2)
 2175 FORMAT(/5X,'Axis orientation')
 2176 FORMAT(7X,'Segment #    ',(1X,21I5))
 2180 FORMAT(7X,'[PHI0] (rads)',(1X,21F5.2))
 2190 FORMAT(/3X,A18/
     .       '+',2X,9('_'),1X,8('_')//
     .       5X,'Transport [SLTRC] =',A9/
     .       5X,'Theta     [THETA] =',F9.2/)
 2200 FORMAT('1',2X,A22/
     .       '+',2X,9('_'),1X,12('_')//
     .       5X,'Longitudinal eddy viscosity           [AX] =',F9.2,
     .         ' m^2/sec'/
     .       5X,'Longitudinal eddy diffusivity        [DXI] =',F9.2,
     .         ' m^2/sec'/
     .       5X,'Chezy coefficient                  [CHEZY] =',F9.2,
     .         ' m^0.5/sec'/
     .       5X,'Sediment temperature                [TSED] =',F9.2,1X,
     .         A2/
     .       5X,'Coefficient of bottom heat exchange [CBHE] =',1PE9.2,
     .         1X,A2,' m/s'/)
 2210 FORMAT(1X,A9/
     .       '+',3('_'),1X,5('_')//
     .       3X,'Ice calculations          [ICEC] = ',5X,A3/
     .       3X,'Solution                [SLICEC] = ',A8/
     .       3X,'Heat exchange            [SLHTC] = ',A8/
     .       3X,'Albedo                  [ALBEDO] = ',F8.2/
     .       3X,'Ice-water heat exchange    [HWI] = ',F8.2,' øC m/s'/
     .       3X,'Light absorption         [BETAI] = ',F8.2/
     .       3X,'Light decay             [GAMMAI] = ',F8.2,' /m'/)
 2220 FORMAT(1X,A14/
     .       '+',6('_'),1X,7('_')//
     .       3X,'Laserjet [LJPC] = ',A3/
     .       3X,'Snapshot [SNPC] = ',A3/
     .       3(5X,A23,' = ',A3/))
 2230 FORMAT(5X,'Number of time intervals [NSNP] =',I7)
 2240 FORMAT('+':/
     .       5X,'Date  (Julian day)       [SNPD] =',11F7.2)
 2250 FORMAT('+':/
     .       5X,'Frequency  (days)        [SNPF] =',11F7.2)
 2251 FORMAT(/3X,'Screen updates [SCRC] = ',A3/
     .       5X,'Number of time intervals [NSCR] =',I7)
 2252 FORMAT('+':/
     .       5X,'Date  (Julian day)       [SCRD] =',11F7.2)
 2253 FORMAT('+':/
     .       5X,'Frequency  (days)        [SCRF] =',11F7.2)
 2260 FORMAT(/3X,'Time series [TSRC] = ',A3)
 2270 FORMAT(5X,'Number of time intervals [NTSR] =',I7)
 2280 FORMAT('+':/
     .       5X,'Date  (Julian day)       [TSRD] =',11F7.2)
 2290 FORMAT('+':/
     .       5X,'Frequency  (days)        [TSRF] =',11F7.2)
 2300 FORMAT(/3X,'Vector plot [VPLC] = ',A3)
 2310 FORMAT(5X,'Number of time intervals [NVPL] =',I7)
 2320 FORMAT('+':/
     .       5X,'Date  (Julian day)       [VPLD] =',11F7.2)
 2330 FORMAT('+':/
     .       5X,'Frequency  (days)        [VPLF] =',11F7.2)
 2340 FORMAT(/3X,'Profile plot [PRFC] = ',A3)
 2350 FORMAT(5X,'Number of profiles       [NPRF] =',I7)
 2360 FORMAT(5X,'Number of stations      [NIPRF] =',I7)
 2370 FORMAT('+':/
     .       5X,'Segment location         [IPRF] =',11I7)
 2380 FORMAT('+':/
     .       5X,'Date  (Julian day)       [PRFD] =',11F7.2)
 2390 FORMAT('+':/
     .       5X,'Frequency  (days)        [PRFF] =',11F7.2)
 2391 FORMAT(/3X,'Spreadsheet plot [SPRC] = ',A3)
 2392 FORMAT(5X,'Number of profiles       [NSPR] =',I7)
 2393 FORMAT(5X,'Number of stations      [NISPR] =',I7)
 2394 FORMAT('+':/
     .       5X,'Segment location         [ISPR] =',11I7)
 2395 FORMAT('+':/
     .       5X,'Date  (Julian day)       [SPRD] =',11F7.2)
 2396 FORMAT('+':/
     .       5X,'Frequency  (days)        [SPRF] =',11F7.2)
 2400 FORMAT(/3X,'Contour plot [CPLC] = ',A3)
 2410 FORMAT(5X,'Number of time intervals [NCPL] =',I7)
 2420 FORMAT('+':/
     .       5X,'Date (Julian day)        [CPLD] =',11F7.2)
 2430 FORMAT('+':/
     .       5X,'Frequency  (days)        [CPLF] =',11F7.2)
 2440 FORMAT(/3X,'Restart out [RSOC] = ',A3/
     .       3X,'Restart in  [RSIC] = ',A3)
 2450 FORMAT(5X,'Number of time intervals [NRSO] =',I7)
 2460 FORMAT('+':/
     .       5X,'Date (Julian day)        [RSOD] =',11F7.2)
 2470 FORMAT('+':/
     .       5X,'Frequency  (days)        [RSOF] =',11F7.2)
 2480 FORMAT(/1X,'Inflow/Outflow'/
     .       '+',14('_')//
     .       3X,'Selective Withdrawal'/
     .       5X,'Branch ',3X,'Calculations',3X,'# of structures'/
     .      (7X,I2,10X,A3,12X,I3))
 2490 FORMAT('+':/
     .       5X,'Structure',3X,'Type',3X,'Width (m)',3X,
     .         'Elevation (m)',3X,'Bottom Layer'/
     .      (8X,I2,4X,A8,1X,F8.1,5X,F8.1,9X,I5))
 2500 FORMAT(/3X,'Number of outlets    [NOUT] =',10I3,
     .      (:/T40,10I3))
 2510 FORMAT('+':/
     .       3X,'Branch',I3,' location at layer ',T34,'[KOUT] =',10I3,
     .      (:/T40,10I3))
 2520 FORMAT(3X,'Number of withdrawals [NWD] =',I3:/
     .       5X,'at segment [IWD] =',10I3,
     .      (:/T26,10I3))
 2530 FORMAT('+':/
     .       5X,'and layer  [KWD] =',10I3,
     .      (:/T26,10I3))
 2540 FORMAT(3X,'Number of tributaries [NTR] =',I3/
     .       '+':/
     .       5X,'at segment [ITR] =',10I3/
     .      (:T26,40I3/))
 2544 FORMAT('+':/
     .       5X,'Inflow placement  [TRPC] = ',7A8:/
     .      (:T33,7A8))
 2545 FORMAT('+':/
     .       5X,'Top elevation    [ELTRT] = ',7F8.2:/
     .      (:T33,7F8.2))
 2546 FORMAT('+':/
     .       5X,'Bottom elevation [ELTRB] = ',7F8.2:/
     .      (:T33,7F8.2))
 2550 FORMAT(3X,'Distributed tributaries [DTRC]')
 2560 FORMAT(5X,'Branch ',I2,' = ',A3)
 2570 FORMAT('1Input Filenames'/
     .       '+',5('_'),1X,9('_')//
     .       3X,'Control               = ',A72/
     .       3X,'Bathymetry            = ',A72/
     .       3X,'Vertical profiles     = ',A72/
     .       3X,'Longitudinal profiles = ',A72/
     .       3X,'Restart               = ',A72/
     .       3X,'Meteorology           = ',A72/
     .       3X,'Water surface         = ',A72/
     .       3X,'Withdrawal            = ',A72/)
 2580 FORMAT('1',2X,'Branch',I2/
     .       '+',2X,6('_')/
     .       5X,'Inflow                               = ',A72/
     .       5X,'Inflow temperature                   = ',A72/
     .       5X,'Inflow concentrations                = ',A72/
     .       5X,'Outflow                              = ',A72/
     .       5X,'Distributed tributary inflows        = ',A72/
     .       5X,'Distributed tributary temperatures   = ',A72/
     .       5X,'Distributed tributary concentrations = ',A72/
     .       5X,'Precipitation                        = ',A72/
     .       5X,'Precipitation temperatures           = ',A72/
     .       5X,'Precipitation concentrations         = ',A72/
     .       5X,'Upstream head                        = ',A72/
     .       5X,'Upstream head temperatures           = ',A72/
     .       5X,'Upstream head concentrations         = ',A72/
     .       5X,'Downstream head                      = ',A72/
     .       5X,'Downstream head temperatures         = ',A72/
     .       5X,'Downstream head concentrations       = ',A72/)
 2585 FORMAT(3X,'Branch',I2/
     .       '+',2X,6('_')/
     .       5X,'Inflow                               = ',A72/
     .       5X,'Inflow temperature                   = ',A72/
     .       5X,'Inflow concentrations                = ',A72/
     .       5X,'Outflow                              = ',A72/
     .       5X,'Distributed tributary inflows        = ',A72/
     .       5X,'Distributed tributary temperatures   = ',A72/
     .       5X,'Distributed tributary concentrations = ',A72/
     .       5X,'Precipitation                        = ',A72/
     .       5X,'Precipitation temperatures           = ',A72/
     .       5X,'Precipitation concentrations         = ',A72/
     .       5X,'Upstream head                        = ',A72/
     .       5X,'Upstream head temperatures           = ',A72/
     .       5X,'Upstream head concentrations         = ',A72/
     .       5X,'Downstream head                      = ',A72/
     .       5X,'Downstream head temperatures         = ',A72/
     .       5X,'Downstream head concentrations       = ',A72/)
 2590 FORMAT('+':/(3X,'Tributary',I2/
     .       '+',2X,9('_')/
     .       5X,'Inflow               = ',A72/
     .       5X,'Inflow temperature   = ',A72/
     .       5X,'Inflow concentration = ',A72)/)
 2600 FORMAT(1X,'Output Filenames'/
     .       '+',6('_'),1X,9('_')//
     .       3X,'Error        = ',A72/
     .       3X,'Warning      = ',A72/
     .       3X,'Snapshot     = ',A72/
     .       3X,'Time-series  = ',A72/
     .       3X,'Profile      = ',A72/
     .       3X,'Vector plot  = ',A72/
     .       3X,'Contour plot = ',A72)
 2610 FORMAT('1',A12,' [CCC] = ',A3/
     .       '+',12('_')/)
 2620 FORMAT(3X,'Kinetics update frequency [KFU] = ',I3/)
 2625 FORMAT(3X,'Constituent',5X,'Computation',2X,'Initial ',
     .         'Concentration',2X,'Printout',2X,'Inflows',2X,
     .         'Tributaries',2X,'Distr Tributaries',2X,'Precipitation'/
     .       '+',2X,11('_'),5X,11('_'),2X,7('_'),1X,13('_'),2X,8('_'),
     .         2X,7('_'),2X,11('_'),2X,5('_'),1X,11('_'),2X,13('_')/
     .       5X,'[CNAME]',9X,'[ACC]',9X,'[C2I] (g/m^3)',7X,'[CPRNC]',
     .         3X,'[INACC]',4X,'[TRACC]',9X,'[DTACC]',10X,'[PRACC]'/)
 2630 FORMAT(3X,A16,3X,A3,12X,F8.3,12X,A3,7X,A3,8X,A3,13X,A3,14X,A3)
 2640 FORMAT(/3X,A17/
     .       '+',2X,11('_'),1X,5('_')//
     .       5X,'Constituent',27X,'Rate'/
     .       '+',4X,11('_'),27X,4('_')/)
 2650 FORMAT(5X,'Suspended solids',T25,'Settling',T40,'  [SSS] =',F6.3,
     .         ' m/day')
 2660 FORMAT(5X,'Coliform',T25,'Decay',T40,'[COLDK] =',F6.3,' /day')
 2670 FORMAT(5X,'Labile DOM',T25,'Decay',T39,'[LDOMDK] =',F6.3,' /day'/
     .       T25,'to refractory  [LRDDK] =',F6.3,' /day')
 2680 FORMAT(5X,'Refractory DOM',T25,'Decay',T39,'[RDOMDK] =',F6.3,
     .         ' /day')
 2690 FORMAT(5X,'Algae',T25,'Growth',T40,'   [AG] =',F6.3,' /day'/
     .       T25,'Mortality         [AM] =',F6.3,' /day'/
     .       T25,'Excretion         [AE] =',F6.3,' /day'/
     .       T25,'Respiration       [AR] =',F6.3,' /day'/
     .       T25,'Settling          [AS] =',F6.3,' m/day')
 2700 FORMAT(5X,'POM',T25,'Decay',T39,'[LPOMDK] =',F6.3,' /day'/
     .       T25,'Settling',T40,'[LPOMS] =',F6.3,' m/day')
 2710 FORMAT(5X,'Phosphorous',T25,'Release',T41,'[PO4R] =',F6.3,
     .       ' g/m^2/day')
 2720 FORMAT(5X,'Ammonium',T25,'Decay',T40,'[NH4DK] =',F6.3,' /day'/
     .       T25,'Release',T41,'[NH4R] =',F6.3,' g/m^2/day')
 2730 FORMAT(5X,'Nitrate-Nitrite',T25,'Decay',T40,'[NO3DK] =',F6.3,
     .         ' /day')
 2740 FORMAT(5X,'Sediment',T25,'Decay',T42,'[SDK] =',F6.3,' /day')
 2750 FORMAT(5X,'Iron',T25,'Settling',T42,'[FES] =',F6.3,' m/day'/
     .       T25,'Release',T42,'[FER] =',F6.3,' g/m^2/day')
 2760 FORMAT(5X,'Oxygen',T25,'Sediment demand  [SOD] = ',13F5.1/
     .      (T50,13F5.1))
 2764 FORMAT(T25,'Fraction of SOD',T41,'[FSOD] = ',F5.1)
 2765 FORMAT(T25,'Shift demand',T42,'[SDC] = ',2X,A3)
 2770 FORMAT(/3X,A24/
     .       '+',2X,5('_'),1X,11('_'),1X,6('_')//
     .       5X,'Constituent',T20,'Rate',T31,'Lower',T44,'Max Lower'/
     .       '+',4X,11('_'),T20,4('_'),T31,5('_'),T44,3('_'),1X,5('_')/
     .       T32,'(',A2,')',T46,'(',A2,')'/)
 2780 FORMAT(5X,'Algae',T19,'Growth',T29,'[AT1]  =',F5.1,T44,
     .       '[AT2]  =',F5.1)
 2790 FORMAT(5X,'Ammonium',T19,'Decay',T28,'[NH4T1] =',F5.1,T43,
     .       '[NH4T2] =',F5.1)
 2800 FORMAT(5X,'Nitrate',T19,'Decay',T28,'[NO3T1] =',F5.1,T43,
     .       '[NO3T2] =',F5.1)
 2810 FORMAT(5X,'Organic',T19,'Decay',T28,'[OMT1]  =',F5.1,T43,
     .       '[SEDT2] =',F5.1)
 2820 FORMAT(/3X,A24/
     .       '+',2X,5('_'),1X,11('_'),1X,6('_')//
     .       5X,'Constituent',T20,'Rate',T31,'Upper',T44,'Max Upper'/
     .       '+',4X,11('_'),T20,4('_'),T31,5('_'),T44,3('_'),1X,5('_')/
     .       T32,'(',A2,')',T46,'(',A2,')'/)
 2830 FORMAT(5X,'Algae',T19,'Growth',T29,'[AT3] =',F5.1,T44,'[AT4] =',
     .         F5.1)
 2840 FORMAT('1',2X,A26/
     .       '+',2X,14('_'),1X,11('_')//
     .       5X,'Oxygen'/'+',4X,6('_')//
     .       7X,'Ammonium       [O2NH4] =',F5.2/
     .       7X,'Organic matter  [O2OM] =',F5.2/
     .       7X,'Respiration     [O2AR] =',F5.2/
     .       7X,'Algal growth    [O2AG] =',F5.2/)
 2850 FORMAT(5X,'Algae'/
     .       '+',4X,5('_')//
     .       7X,'Phosphorous [BIOP] =',F6.3/
     .       7X,'Nitrogen    [BION] =',F6.3/
     .       7X,'Carbon      [BIOC] =',F6.3/)
 2860 FORMAT(3X,'CBOD'/
     .       '+',2X,4('_')//
     .       5X,'Decay rate                   [KBOD] =',F6.3,' /day'/
     .       5X,'Temperature adjustment       [TBOD] =',F6.3/
     .       5X,'Ultimate CBOD to CBOD5 ratio [RBOD] =',F6.3/)
 2870 FORMAT(3X,A15/
     .       '+',2X,4('_'),1X,10('_')//
     .       5X,'Phosphorous [AHSP] =',F6.3,' g/m^3'/
     .       5X,'Nitrogen    [AHSN] =',F6.3,' g/m^3'/)
 2880 FORMAT(3X,A5/
     .       '+',2X,5('_')//
     .       5X,'Attenuation'/
     .       '+',4X,11('_')//
     .       7X,'Surface layer      [BETA] =',F5.2/
     .       7X,'Water             [EXH2O] =',F5.2,' /m'/
     .       7X,'Inorganic solids   [EXSS] =',F5.2,' /m'/
     .       7X,'Organic solids     [EXOM] =',F5.2,' /m'//
     .       5X,'Saturation Intensity'/
     .       '+',4X,10('_'),1X,9('_')//
     .       7X,'Algae [ASAT] =',F6.1,' W/m^2'/)
 2890 FORMAT(3X,A9/
     .       '+',2X,9('_')//
     .       5X,'Oxygen          [DMO2] =',1PE10.3,' m^2/g'/
     .       5X,'Carbon dioxide [DMCO2] =',1PE10.3,' m^2/g'/)
 2900 FORMAT(3X,A24/
     .       '+',2X,12('_'),1X,11('_')//
     .       5X,'Phosphorous [PARTP] =',F6.3,' m^3/g'/
     .       5X,'Nitrogen    [PARTN] =',F6.3,' m^3/g'/)
 2910 FORMAT(3X,A23/
     .       '+',2X,13('_'),1X,9('_')//
     .       5X,'Aerobic oxygen limit     [O2LIM] =',F5.2,' g/m^3'/
     .       5X,'Coliform Q10            [COLQ10] =',F5.2/
     .       5X,'Fraction algae to POM     [APOM] =',F5.2/
     .       5X,'Sediment release of CO2   [CO2R] =',F5.2,' g/m^2/day')
 2920 FORMAT('1',A8/'+',8('_')//
     .       3X,'Overall Grid'/
     .       '+',2X,7('_'),1X,4('_')//
     .       5X,'Total segments  [IMP] =',I3,3X,'Total layers ',
     .         '      [KMP] =   ',I4/
     .       5X,'Total branches  [NBP] =',I3,3X,'Bottom elevation ',
     .         '[ELBOT] =',F8.2,' m'/
     .       5X,'Surface layer    [KT] =',I3//5X,'Vertical Spacing ',
     .         '[H]'/
     .       '+',4X,8('_'),1X,7('_'))
 2930 FORMAT(7X,'Layer',6X,19I5)
 2940 FORMAT(7X,'Height (m) ',19F5.1)
 2950 FORMAT(/3X,'Branch',I2/
     .       '+',2X,6('_')/
     .       5X,'Upstream segment',T29,'[US] =',I3,T41,'Downstream ',
     .         'segment       [DS] =',I3/
     .       5X,'Upstream head segment [UHS] =',I3,T41,'Downstream ',
     .         'head segment [DHS] =',I3)
 2960 FORMAT(76X,'Bathymetry [B], m'/
     .       '+',75X,10('_')/)
 2970 FORMAT('0',5X,21I7)
 2972 FORMAT('1',5X,21I7)
 2975 FORMAT(6X,21I7)
 2980 FORMAT(F7.0)
 2990 FORMAT(1X,I2,4X,21A7)
 3000 FORMAT('+','   [KT]')
 3010 FORMAT(/1X,'Initial Branch Volume [VOLB] =',F14.0,' m^3'/)
 3015 FORMAT('1','Initial Branch Volume [VOLB] =',F14.0,' m^3'/)
 3020 FORMAT(T14,'Branch ',I2,' Volume-Area-Elevation Table'/
     .       '+',T14,6('_'),4X,6('_'),1X,4('_'),1X,9('_'),1X,5('_')/)
 3030 FORMAT(1X,'Layer',T12,'Elevation',T27,'Area',T39,'Volume',T51,
     .         'Active Cells'/
     .       '+',5('_'),T12,9('_'),T27,4('_'),T39,6('_'),T51,6('_'),
     .         1X,5('_')/
     .       T15,'(m)',T23,'(1.0E6 m^2)',T36,'(1.0E6 m^3)'/)
 3035 FORMAT(1X,'Layer',11x,'Elevation',20x,'Volume'/
     .     6('_'),6x,19('_'),5x,31('_')/
     .     15x,'(m)',8x,'(ft)',7x,'(1.0E6 m^3)',2x
     .     ,'(1.0e6 acre-feet)'/)
 3040 FORMAT(1X,I3,T12,F7.2,T24,F8.3,T37,F9.3,T54,I4)
 3045 FORMAT(1X,I3,8x,F7.2,5x,F8.1,5x,F9.3,6x,f9.3)
 3050 FORMAT('+',T7,'[KT]')
 3060 FORMAT(/1X,'Initial Grid Volume [VOLG] = ',F15.0,' m^3'/)
 3065 FORMAT('1','Initial Grid Volume [VOLG] = ',F15.0,' m^3'/)
 3070 FORMAT(T14,'Grid Volume-Area-Elevation Table'/
     .       '+',T14,4('_'),1X,6('_'),1X,4('_'),1X,9('_'),1X,5('_')/)
      END

************************************************************************* 
**       S U B R O U T I N E   P R I N T  I N I T I A L  G R I D       ** 
************************************************************************* 
 
      SUBROUTINE PRINT_INITIAL_GRID (KT) 
      INCLUDE   'w2.inc' 

***** Type declarations

      REAL         ICETH 
      INTEGER      IPR,  CN, US, DS,  CUS
      LOGICAL      CONSTITUENTS, ICE, ICE_CALC
      LOGICAL      NEW_PAGE
      CHARACTER*1  ESC
      CHARACTER*3  CPRN,   HPRN
      CHARACTER*8  CUNIT
      CHARACTER*7  CONV
      CHARACTER*16 CNAME1
      CHARACTER*23 HNAME1, CNAME2
      CHARACTER*29 HNAME2
      CHARACTER*72 TITLE

***** Common declarations

      COMMON /GRIDC1/ T2(KMP,IMP),   U(KMP,IMP),  W(KMP,IMP), Z(IMP),
     .                ICETH(IMP),    IPR(IMP),    CN(NCP),    ICE(IMP)
      COMMON /GRIDC2/ IEPR,          KEPR,        NAC
      COMMON /GRIDC3/ CONV(KMP,IMP), CPRN(NCP),   HPRN(4)
      COMMON /GRIDC4/ US(NBP),       DS(NBP),     CUS(NBP),   KB(IMP)
      COMMON /NAMESC/ CNAME1(NCP),   CNAME2(NCP), HNAME1(3),  HNAME2(3),
     .                CUNIT(NCP),    TITLE(6)
      COMMON /GRDLGC/ CONSTITUENTS,  ICE_CALC
      COMMON /UNITNC/ INI

***** Variable initialization

      ESC      = CHAR(027)
      NEW_PAGE = .TRUE.

***** Water surface variables 
 
      DO I=1,IEPR
        WRITE (CONV(1,I),'(F7.4)') Z(IPR(I))
      END DO
      WRITE (INI,3000) 'Water Surface [Z], m' 
      WRITE (INI,3010) (IPR(I),I=1,IEPR) 
      WRITE (INI,3020) (CONV(1,I),I=1,IEPR)
      IF (ICE_CALC) WRITE (INI,3030) 'Ice Thickness (m)',(ICETH(IPR(I)), 
     .                                I=1,IEPR) 
 
***** Velocities and temperatures
 
      IF (HPRN(3).EQ.' ON') THEN 
        DO I=1,IEPR 
          DO K=KT,KB(IPR(I)) 
            WRITE (CONV(K,I),4010) T2(K,IPR(I)) 
          END DO
        END DO
        IF (NEW_PAGE) THEN
          WRITE (INI,3070) HNAME2(3)
          NLINES = (KMP-KT+5)
        ELSE
          WRITE (INI,3075) HNAME2(3)
        END IF
        NLINES   = NLINES+(KMP-KT+5)
        NEW_PAGE = NLINES.GT.73
        WRITE (INI,3010) (IPR(I),I=1,IEPR) 
        DO K=KT,KEPR 
          WRITE (INI,3080) K,(CONV(K,I),I=1,IEPR) 
        END DO
      END IF 
 
***** Constituent concentrations 
 
      DO J=1,NAC 
        IF (CPRN(CN(J)).EQ.' ON') THEN 
          MULT = 1.0 
          IF (CN(J).GE.9.AND.CN(J).LE.11) MULT = 1000. 
          DO I=1,IEPR 
            DO K=KT,KB(IPR(I)) 
              WRITE (CONV(K,I),4010) C2(K,IPR(I),CN(J))*MULT 
            END DO
          END DO
          IF (NEW_PAGE) THEN
            WRITE (INI,3100)  CNAME2(CN(J))
            NLINES = (KMP-KT+4)
          ELSE
            WRITE (INI,3110)  CNAME2(CN(J))
          END IF
          NLINES   = NLINES+(KMP-KT+4)
          NEW_PAGE = NLINES.GT.73
          WRITE (INI,3010) (IPR(I),I=1,IEPR) 
          DO K=KT,KEPR 
            WRITE (INI,3080) K,(CONV(K,I),I=1,IEPR) 
          END DO
        END IF 
      END DO
      WRITE (INI,'(1X,A2)') ESC//'E'

***** FORMAT statements

 3000 FORMAT(///75X,A21/
     .       '+',75X,5('_'),1X,7('_')/) 
 3010 FORMAT(2X,1000I7) 
 3020 FORMAT(3X,1000A7/)
 3030 FORMAT(/75X,A17/
     .       '+',74X,3('_'),1X,9('_')// 
     .       3X,24F7.4/)
 3070 FORMAT('1',75X,A23/) 
 3075 FORMAT(/75X,A23/) 
 3080 FORMAT(1X,I2,200A7)
 3100 FORMAT('1',75X,A23/) 
 3110 FORMAT(/75X,A23/) 
 4010 FORMAT(F7.2) 
      RETURN 
      END

***********************************************************************
**           F U N C T I O N   V A L I D   F I L E N A M E           **
***********************************************************************

      LOGICAL FUNCTION VALID_FILENAME(FNAME)
        CHARACTER*30 FNAME
        LOGICAL VALID_NAME, VALID_EXTENSION

        VALID_NAME = .FALSE.
        DO I=1,9
          IF ((FNAME(I:I).EQ.'.').AND.(.NOT.VALID_NAME)) THEN
            VALID_NAME = .TRUE.
            NPER       = I
          END IF
        END DO
        IF (VALID_NAME) THEN
          IF (FNAME(NPER+4:NPER+4).EQ.' ') THEN
            VALID_EXTENSION = .TRUE.
          ELSE 
            VALID_EXTENSION = .FALSE.
          END IF
        END IF
        IF (VALID_NAME.AND.VALID_EXTENSION) THEN
          VALID_FILENAME = .TRUE.
        ELSE
          VALID_FILENAME = .FALSE.
        END IF
      END
