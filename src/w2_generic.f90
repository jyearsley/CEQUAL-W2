!***********************************************************************************************************************************
!**                                                                                                                               **
!**                                                         CE-QUAL-W2                                                            **
!**                                            A Two-dimensional, Laterally Averaged,                                             **
!**                                             Hydrodynamic and Water Quality Model                                              **
!**                                                            for                                                                **
!**                                           Rivers, Lakes, Reservoirs, and Estuaries                                            **
!**                                                                                                                               **
!**                                                        Version 3.5                                                            **
!**                                                                                                                               **
!**                                                       Thomas M. Cole                                                          **
!**                                                Water Quality Modeling Group                                                   **
!**                                                U.S. Army Corps of Engineers                                                   **
!**                                                Waterways Experiment Station                                                   **
!**                                                Vicksburg, Mississippi 39180                                                   **
!**                                                phone number: (601) 634-3283                                                   **
!**                                                 fax number: (601) 634-3129                                                    **
!**                                                 e-mail: colet@wes.army.mil                                                    **
!**                                                                                                                               **
!**                                                        Scott Wells                                                            **
!**                                       Department of Civil and Environmental Engineering                                       **
!**                                                  Portland State University                                                    **
!**                                                         PO Box 751                                                            **
!**                                                 Portland, Oregon  97207-0751                                                  **
!**                                                 phone number: (503) 725-4276                                                  **
!**                                                 fax   number: (503) 725-5950                                                  **
!**                                                   e-mail: scott@cecs.pdx.edu                                                  **
!**                                                                                                                               **
!***********************************************************************************************************************************
 
!***********************************************************************************************************************************
!**                                                                                                                               **
!**                  The long arm of the lawyers has found its way into the water quality modeling arena, so:                     **
!**                                                                                                                               **
!**  This model was developed and is maintained by the U.S. Army Engineer Waterways Experiment Station, Vicksburg, MS.  The US    **
!**  government and its components are not responsible for any damages, including incidental or consequential damages, arising    **
!**  from use or misuse of this model, or from results achieved or conclusions drawn by others.  Distribution of this model is    **
!**  restricted by the Export Administration Act of 1969,  50 app. USC subsections 2401-2420, as amended, and other applicable    **
!**  laws or regulations.                                                                                                         **
!**                                                                                                                               **
!***********************************************************************************************************************************
 
!***********************************************************************************************************************************
!**                                                      Module Declaration                                                       **
!***********************************************************************************************************************************

MODULE PREC
  INTEGER, PARAMETER :: I2=SELECTED_INT_KIND (3)
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(15)
END MODULE PREC
MODULE RSTART
  USE PREC
  REAL                                               :: DLTS,   CURMAX
  INTEGER                                            :: RSODP,  DLTDP,  TSRDP,  WDODP,  CUF,    RSO=31
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: SNPDP,  VPLDP,  CPLDP,  PRFDP,  SCRDP,  SPRDP,  FLXDP, NSPRF
  REAL                                               :: NXTMRS, NXTMWD, NXTMTS
  REAL,              ALLOCATABLE, DIMENSION(:)       :: NXTMSN, NXTMPR, NXTMSP, NXTMCP, NXTMVP, NXTMSC, NXTMFL
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ELWS,   SBKT,   ELTMF
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: TSSUH2, TSSDH2, SAVH2,  SAVHR,  SU,     SW,     SAZ
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: CSSUH2, CSSDH2
  REAL(R8)                                           :: ELTM
  REAL(R8),          ALLOCATABLE, DIMENSION(:)       :: VOLIN,  VOLOUT, VOLUH,  VOLDH,  VOLPR,  VOLTRB, VOLDT, VOLWD,  VOLEV
  REAL(R8),          ALLOCATABLE, DIMENSION(:)       :: VOLSBR, VOLTBR, VOLSR,  VOLTR
  REAL(R8),          ALLOCATABLE, DIMENSION(:)       :: TSSEV,  TSSPR,  TSSTR,  TSSDT,  TSSWD,  TSSUH,  TSSDH, TSSIN,  TSSOUT
  REAL(R8),          ALLOCATABLE, DIMENSION(:)       :: TSSS,   TSSB,   TSSICE
  REAL(R8),          ALLOCATABLE, DIMENSION(:)       :: ESBR,   ETBR,   EBRI,   SZ
  REAL(R8),          ALLOCATABLE, DIMENSION(:,:)     :: CMBRT
END MODULE RSTART
MODULE GLOBAL
  REAL,   PARAMETER                                  :: DAY=86400.0,  NONZERO=1.0E-20, REFL=0.94, FRAZDZ=0.14, DZMIN=1.4E-7
  REAL,   PARAMETER                                  :: AZMIN=1.4E-6, DZMAX=1.0E3,     RHOW=1000.0
  REAL                                               :: DLT,    DLTMIN, DLTTVD
  REAL                                               :: BETABR, START,  HMAX2,  CURRENT
  REAL,   POINTER,                DIMENSION(:,:)     :: U,      W,      T2,     AZ,     RHO,    ST,     SB
  REAL,   POINTER,                DIMENSION(:,:)     :: DLTLIM, VSH,    ADMX,   DM,     ADMZ,   HDG,    HPG,    GRAV
  REAL,   TARGET,    ALLOCATABLE, DIMENSION(:,:)     :: T1,     TSS
  REAL,   TARGET,    ALLOCATABLE, DIMENSION(:,:,:)   :: C1,     C2,     C1S,    CSSB,   CSSK,   HYD,    KF,     CD
  REAL,   TARGET,    ALLOCATABLE, DIMENSION(:,:,:,:) :: AF,     EF
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ICETH,  ELKT,   HMULT,  CMULT,  CDMULT, WIND2,  AZMAX,  PALT
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: QSS,    VOLUH2, VOLDH2, QUH1,   QDH1,   UXBR,   UYBR,   VOL
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: ALLIM,  APLIM,  ANLIM,  ASLIM,  KFS
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: ELLIM,  EPLIM,  ENLIM,  ESLIM
  INTEGER                                            :: IMX,    KMX,    NBR,    NTR,    NWD,    NWB,    NCT,    NBOD
  INTEGER                                            :: NST,    NSP,    NGT,    NPI,    NPU,    NWDO,   NIKTSR, NUNIT
  INTEGER                                            :: JW,     JB,     JC,     IU,     ID,     KT,     I,      JJB
  INTEGER                                            :: NOD,    NDC,    NAL,    NSS,    NHY,    NFL,    NEP,    NEPT
  INTEGER                                            :: NZP,    nzpt, JZ,     NZOOS,  NZOOE,  nmc,   nmct  ! number of zooplankton groups, CONSTIUENT NUMBER FOR ZOOPLANKTON, START AND END
  INTEGER, POINTER,               DIMENSION(:)       :: SNP,    PRF,    VPL,    CPL,    SPR,    FLX
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: BS,     BE,     US,     CUS,    DS,     JBDN
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: KB,     KTI,    SKTI,   KTWB,   KBMIN,  CDHS
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: UHS,    DHS,    UQB,    DQB
  INTEGER, TARGET,   ALLOCATABLE, DIMENSION(:,:)     :: OPT
  LOGICAL,           ALLOCATABLE, DIMENSION(:)       :: ICE,    ICE_CALC
  CHARACTER(10)                                      :: CTIME
  CHARACTER(12)                                      :: CDATE
  CHARACTER(72)                                      :: RSIFN
  REAL,     SAVE, ALLOCATABLE, DIMENSION(:,:)   :: RATZ,   CURZ1,  CURZ2,  CURZ3    ! SW 5/15/06
  DATA                                                  NDC /23/, NHY /15/, NFL /118/
  DATA                                                  G /9.81/, PI/3.14159265359/      
END MODULE GLOBAL
MODULE GEOMC
  USE PREC
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: JBUH,   JBDH,   JWUH,   JWDH
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ALPHA,  SINA,   COSA,   SLOPE,  BKT,    ELWS2,   DLX,    DLXR
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: H,      H1,     H2,     BH1,    BH2,    BHR1,    BHR2,   AVHR
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: B,      BI,     BB,     BH,     BHR,    BR,      EL,     AVH1,  AVH2, bnew ! SW 1/23/06
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: DEPTHB, DEPTHM, FETCHU, FETCHD
  REAL(R8),          ALLOCATABLE, DIMENSION(:)       :: Z
END MODULE GEOMC
MODULE NAMESC
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: LNAME
  CHARACTER(6),      ALLOCATABLE, DIMENSION(:)       :: CUNIT,  CUNIT2
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:)       :: CNAME2, CDNAME2
  CHARACTER(9),      ALLOCATABLE, DIMENSION(:)       :: FMTH,   FMTC,   FMTCD
  CHARACTER(19),     ALLOCATABLE, DIMENSION(:)       :: CNAME1
  CHARACTER(43),     ALLOCATABLE, DIMENSION(:)       :: CNAME,  CNAME3, CDNAME, CDNAME3, HNAME
  CHARACTER(72),     ALLOCATABLE, DIMENSION(:)       :: TITLE
  CHARACTER(10),     ALLOCATABLE, DIMENSION(:,:)     :: CONV
END MODULE NAMESC
MODULE STRUCTURES
  REAL                                               :: DIA,    FMAN,   CLEN,   CLOSS,  UPIE,   DNIE
  REAL,              ALLOCATABLE, DIMENSION(:)       :: QOLD,   QOLDS,  VMAX,   DTP,    DTPS
  REAL,              ALLOCATABLE, DIMENSION(:)       :: EGT,    A1GT,   B1GT,   G1GT,   A2GT,   B2GT,   G2GT
  REAL,              ALLOCATABLE, DIMENSION(:)       :: QGT,    GTA1,   GTB1,   GTA2,   GTB2,   BGT
  REAL,              ALLOCATABLE, DIMENSION(:)       :: QSP,    A1SP,   B1SP,   A2SP,   B2SP,   ESP
  REAL,              ALLOCATABLE, DIMENSION(:)       :: EUPI,   EDPI,   WPI,    DLXPI,  FPI,    FMINPI, QPI
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: YS,     VS,     YSS,    VSS,    YST,    VST,    YSTS,   VSTS
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: IUPI,   IDPI,   JWUPI,  JWDPI,  JBDPI,  JBUPI
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: IUSP,   IDSP,   JWUSP,  JWDSP,  JBUSP,  JBDSP
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: IUGT,   IDGT,   JWUGT,  JWDGT,  JBUGT,  JBDGT
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: IWR,    KTWR,   KBWR
  LOGICAL,           ALLOCATABLE, DIMENSION(:)       :: LATERAL_SPILLWAY, LATERAL_PIPE, LATERAL_GATE, LATERAL_PUMP, BEGIN, WLFLAG
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:)       :: LATGTC, LATSPC, LATPIC, LATPUC, DYNGTC 
  DATA                                                  THR/0.01/, NN/19/, OMEGA/0.8/, EPS2/0.0001/
  DATA                                                  NNPIPE /19/, NC/7/
END MODULE STRUCTURES
MODULE TRANS
  USE PREC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: THETA
  REAL,    POINTER,               DIMENSION(:,:)     :: COLD,   CNEW,   SSB,    SSK
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: DX,     DZ,     DZQ
  REAL(R8),          ALLOCATABLE, DIMENSION(:,:)     :: ADX,    ADZ,    AT,     VT,     CT,     DT
END MODULE TRANS
MODULE SURFHE
  REAL                                               :: RHOWCP, PHISET
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ET,     CSHE,   LAT,    LONGIT, SHADE,  RB,     RE,     RC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: WIND,   WINDH,  WSC,    AFW,    BFW,    CFW,    PHI0
  LOGICAL,           ALLOCATABLE, DIMENSION(:)       :: RH_EVAP
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: IWIND  !MLM 08/12/05
END MODULE SURFHE
MODULE TVDC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: QIN,    QTR,    QDTR,   PR,     ELUH,   ELDH,   QWD,    QSUM
  REAL,              ALLOCATABLE, DIMENSION(:)       :: TIN,    TTR,    TDTR,   TPR,    TOUT,   TWDO,   TIND,   QIND   
  REAL,              ALLOCATABLE, DIMENSION(:)       :: TAIR,   TDEW,   CLOUD,  PHI,    SRON   
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: CIN,    CTR,    CDTR,   CPR,    CIND,   TUH,    TDH,    QOUT
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: CUH,    CDH
  INTEGER                                            :: NAC,    NOPEN
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: NACPR,  NACIN,  NACDT,  NACTR,  NACD,   CN
  INTEGER,           ALLOCATABLE, DIMENSION(:,:)     :: TRCN,   INCN,   DTCN,   PRCN
  LOGICAL                                            :: CONSTITUENTS
  CHARACTER(72)                                      :: QGTFN,  QWDFN,  WSCFN,  SHDFN
  CHARACTER(72),     ALLOCATABLE, DIMENSION(:)       :: METFN,  QOTFN,  QINFN,  TINFN,  CINFN,  QTRFN,  TTRFN,  CTRFN,  QDTFN
  CHARACTER(72),     ALLOCATABLE, DIMENSION(:)       :: TDTFN,  CDTFN,  PREFN,  TPRFN,  CPRFN,  EUHFN,  TUHFN,  CUHFN,  EDHFN
  CHARACTER(72),     ALLOCATABLE, DIMENSION(:)       :: EXTFN,  CDHFN,  TDHFN
END MODULE TVDC
MODULE KINETIC
  REAL                                               :: kdo                        !v3.5
  REAL,    POINTER,               DIMENSION(:,:)     :: TDS,    COL,    NH4,    NO3,    PO4,    FE,     DSI,    PSI,    LDOM
  REAL,    POINTER,               DIMENSION(:,:)     :: RDOM,   LPOM,   RPOM,   O2,     TIC,    ALK
  REAL,    POINTER,               DIMENSION(:,:)     :: COLSS,  NH4SS,  NO3SS,  PO4SS,  FESS,   DSISS,  PSISS,  LDOMSS
  REAL,    POINTER,               DIMENSION(:,:)     :: RDOMSS, LPOMSS, RPOMSS, DOSS,   TICSS,  CASS
  REAL,    POINTER,               DIMENSION(:,:)     :: PH,     CO2,    HCO3,   CO3
  REAL,    POINTER,               DIMENSION(:,:)     :: TN,     TP,     TKN
  REAL,    POINTER,               DIMENSION(:,:)     :: DON,    DOP,    DOC
  REAL,    POINTER,               DIMENSION(:,:)     :: PON,    POP,    POC
  REAL,    POINTER,               DIMENSION(:,:)     :: TON,    TOP,    TOC
  REAL,    POINTER,               DIMENSION(:,:)     :: APR,    CHLA,   ATOT
  REAL,    POINTER,               DIMENSION(:,:)     :: O2DG
  REAL,    POINTER,               DIMENSION(:,:)     :: SSSI,   SSSO,   TISS,   TOTSS
  REAL,    POINTER,               DIMENSION(:,:)     :: PO4AR,  PO4AG,  PO4AP,  PO4SD,  PO4SR,  PO4NS,  PO4POM, PO4DOM, PO4OM
  REAL,    POINTER,               DIMENSION(:,:)     :: PO4ER,  PO4EG,  PO4EP,  TICEP,  DOEP,   DOER
  REAL,    POINTER,               DIMENSION(:,:)     :: NH4ER,  NH4EG,  NH4EP,  NO3EG,  DSIEG,  LDOMEP, LPOMEP
  REAL,    POINTER,               DIMENSION(:,:)     :: NH4AR,  NH4AG,  NH4AP,  NH4SD,  NH4SR,  NH4D,   NH4POM, NH4DOM, NH4OM
  REAL,    POINTER,               DIMENSION(:,:)     :: NO3AG,  NO3D,   NO3SED
  REAL,    POINTER,               DIMENSION(:,:)     :: DSIAG,  DSID,   DSISD,  DSISR,  DSIS
  REAL,    POINTER,               DIMENSION(:,:)     :: PSIAM,  PSID,   PSINS
  REAL,    POINTER,               DIMENSION(:,:)     :: FENS,   FESR
  REAL,    POINTER,               DIMENSION(:,:)     :: LDOMAP, LDOMD,  LRDOMD, RDOMD
  REAL,    POINTER,               DIMENSION(:,:)     :: LPOMAP, LPOMD,  LRPOMD, RPOMD,  LPOMNS, RPOMNS
  REAL,    POINTER,               DIMENSION(:,:)     :: DOAP,   DOAR,   DODOM,  DOPOM,  DOOM,   DONIT
  REAL,    POINTER,               DIMENSION(:,:)     :: DOSED,  DOSOD,  DOBOD,  DOAE
  REAL,    POINTER,               DIMENSION(:,:)     :: CBODU,  CBODDK, TICAP
  REAL,    POINTER,               DIMENSION(:,:)     :: SEDD,   SODD,   SEDAS,  SEDOMS, SEDNS
  REAL,    POINTER,               DIMENSION(:,:,:)   :: SS,     SSSS,   ALG,    ASS,    CBOD,   CBODSS, CG,   CGSS
  REAL,    POINTER,               DIMENSION(:,:,:)   :: AGR,    ARR,    AER,    AMR,    ASR
  REAL,    POINTER,               DIMENSION(:,:,:)   :: EGR,    ERR,    EER,    EMR,    EBR
! v3.5 start
  REAL,    POINTER,               DIMENSION(:,:)     :: LDOMP,  RDOMP,  LPOMP,  RPOMP,  LDOMN,  RDOMN,  LPOMN,  RPOMN
  REAL,    POINTER,               DIMENSION(:,:)     :: LDOMPSS,  RDOMPSS, LPOMPSS, RPOMPSS, LDOMNSS, RDOMNSS
  REAL,    POINTER,               DIMENSION(:,:)     :: LPOMNSS,  RPOMNSS
  REAL,    POINTER,               DIMENSION(:,:)     :: LDOMPAP,  LDOMPEP, LPOMPAP, LPOMPNS, RPOMPNS
  REAL,    POINTER,               DIMENSION(:,:)     :: LDOMNAP,  LDOMNEP, LPOMNAP, LPOMNNS, RPOMNNS
  REAL,    POINTER,               DIMENSION(:,:)     :: SEDDP,    SEDASP,  SEDOMSP, SEDNSP,  LPOMEPP
  REAL,    POINTER,               DIMENSION(:,:)     :: SEDDN,    SEDASN,  SEDOMSN, SEDNSN,  LPOMEPN, SEDNO3
  REAL,    POINTER,               DIMENSION(:,:)     :: SEDDC,    SEDASC,  SEDOMSC, SEDNSC,  LPOMEPC  
  REAL,    POINTER,               DIMENSION(:,:)     :: CBODNS,   SEDCB,   SEDCBP,  SEDCBN,  SEDCBC
  REAL,    POINTER,               DIMENSION(:,:)     :: sedbr,    sedbrp,  sedbrc,  sedbrn        !cb 11/30/06
! v3.5 end 
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: EPM,    EPD,    EPC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: CGQ10,  CG0DK,  CG1DK,  CGS
  REAL,              ALLOCATABLE, DIMENSION(:)       :: SOD,    SDK,    LPOMDK, RPOMDK, LDOMDK, RDOMDK, LRDDK,  LRPDK
  REAL,              ALLOCATABLE, DIMENSION(:)       :: SSS,    TAUCR,  POMS,   FES, seds, sedb  !cb 11/27/06
  REAL,              ALLOCATABLE, DIMENSION(:)       :: AG,     AR,     AE,     AM,     AS,     AHSN,   AHSP,   AHSSI,  ASAT
  REAL,              ALLOCATABLE, DIMENSION(:)       :: AP,     AN,     AC,     ASI,    ACHLA,  APOM,   ANPR   
  REAL,              ALLOCATABLE, DIMENSION(:)       :: EG,     ER,     EE,     EM,     EB
  REAL,              ALLOCATABLE, DIMENSION(:)       :: EHSN,   EHSP,   EHSSI,  ESAT,   EHS,    ENPR
  REAL,              ALLOCATABLE, DIMENSION(:)       :: EP,     EN,     EC,     ESI,    ECHLA,  EPOM
  REAL,              ALLOCATABLE, DIMENSION(:)       :: BETA,   EXH2O,  EXSS,   EXOM,   EXA
  REAL,              ALLOCATABLE, DIMENSION(:)       :: DSIR,   PSIS,   PSIDK,  PARTSI
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ORGP,   ORGN,   ORGC,   ORGSI
  REAL,              ALLOCATABLE, DIMENSION(:)       :: BODP,   BODN,   BODC   
  REAL,              ALLOCATABLE, DIMENSION(:)       :: PO4R,   PARTP
  REAL,              ALLOCATABLE, DIMENSION(:)       :: NH4DK,  NH4R,   NO3DK,  NO3S
  REAL,              ALLOCATABLE, DIMENSION(:)       :: O2AG,   O2AR,   O2OM,   O2NH4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: O2EG,   O2ER
  REAL,              ALLOCATABLE, DIMENSION(:)       :: CO2R,   FER
  REAL,              ALLOCATABLE, DIMENSION(:)       :: KBOD,   TBOD,   RBOD
  REAL,              ALLOCATABLE, DIMENSION(:)       :: CAQ10,  CADK,   CAS
  REAL,              ALLOCATABLE, DIMENSION(:)       :: OMT1,   OMT2,   SODT1,  SODT2,  NH4T1,  NH4T2,  NO3T1,  NO3T2
  REAL,              ALLOCATABLE, DIMENSION(:)       :: OMK1,   OMK2,   SODK1,  SODK2,  NH4K1,  NH4K2,  NO3K1,  NO3K2
  REAL,              ALLOCATABLE, DIMENSION(:)       :: AT1,    AT2,    AT3,    AT4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: AK1,    AK2,    AK3,    AK4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ET1,    ET2,    ET3,    ET4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: EK1,    EK2,    EK3,    EK4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: REAER,  WIND10, CZ,     QC,     QERR
  REAL,              ALLOCATABLE, DIMENSION(:)       :: RCOEF1, RCOEF2, RCOEF3, RCOEF4
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: DO1,    DO2,    DO3,    GAMMA
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: SED,    FPSS,   FPFE
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: CBODD
! v3.5 start
  REAL,              ALLOCATABLE, DIMENSION(:)       :: CBODS  
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: ORGPLD,  ORGPRD,   ORGPLP,    ORGPRP,  ORGNLD,  ORGNRD, ORGNLP, ORGNRP
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: LDOMPMP, LDOMNMP,  LPOMPMP,   LPOMNMP, RPOMPMP, RPOMNMP
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: LPZOOINP,LPZOOINN, LPZOOOUTP, LPZOOOUTN
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: SEDC,    SEDN, SEDP
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: SEDVPC,  SEDVPP, SEDVPN
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: SDKV,    SEDDKTOT
  INTEGER                                            :: nldomp,nrdomp,nlpomp,nrpomp,nldomn,nrdomn,nlpomn,nrpomn
! v3.5 end
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: NAF,    NEQN,   ANEQN,  ENEQN
  INTEGER,           ALLOCATABLE, DIMENSION(:,:)     :: KFCN
  LOGICAL,           ALLOCATABLE, DIMENSION(:)       :: SEDIMENT_RESUSPENSION
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:)       :: CAC,    REAERC
  CHARACTER(10),     ALLOCATABLE, DIMENSION(:,:)     :: LFPR
  CONTAINS
    FUNCTION SATO (T,SAL,P,SALT_WATER)
      LOGICAL :: SALT_WATER
      SATO = EXP(7.7117-1.31403*(LOG(T+45.93)))*P
      IF (SALT_WATER) SATO = EXP(LOG(SATO)-SAL*(1.7674E-2-1.0754E1/(T+273.15)+2.1407E3/(T+273.15)**2))
    END FUNCTION SATO
    FUNCTION FR (TT,TT1,TT2,SK1,SK2)
      FR = SK1*EXP(LOG(SK2*(1.0-SK1)/(SK1*(1.0-SK2)))/(TT2-TT1)*(TT-TT1))
    END FUNCTION FR     
    FUNCTION FF (TT,TT3,TT4,SK3,SK4)
      FF = SK4*EXP(LOG(SK3*(1.0-SK4)/(SK4*(1.0-SK3)))/(TT4-TT3)*(TT4-TT))
    END FUNCTION FF
END MODULE KINETIC
MODULE SELWC
  REAL,              ALLOCATABLE, DIMENSION(:)   :: EWD,    VNORM,  QNEW
  REAL,              ALLOCATABLE, DIMENSION(:,:) :: QSTR,   QSW,    ESTR,   WSTR
  INTEGER,           ALLOCATABLE, DIMENSION(:)   :: NSTR,   NOUT,   KTWD,   KBWD,   KTW,   KBW
  INTEGER,           ALLOCATABLE, DIMENSION(:,:) :: KTSW,   KBSW,   KOUT
END MODULE SELWC
MODULE GDAYC
  REAL                                           :: DAYM,   EQTNEW 
  INTEGER                                        :: JDAYG,  M,      YEAR,   GDAY
  LOGICAL                                        :: LEAP_YEAR
  CHARACTER(9)                                   :: MONTH
END MODULE GDAYC
MODULE SCREENC
  USE PREC
  REAL                                           :: JDAY,   DLTS1,  JDMIN,  MINDLT, DLTAV,  ELTMJD
  REAL(R8),          ALLOCATABLE, DIMENSION(:)   :: ZMIN,   CMIN,   CMAX,   HYMIN,  HYMAX,  CDMIN,  CDMAX
  INTEGER                                        :: ILOC,   KLOC,   IMIN,   KMIN,   NIT,    NV,     JTT,     JWW
  INTEGER,           ALLOCATABLE, DIMENSION(:)   :: IZMIN
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:)   :: ACPRC,  AHPRC,  ACDPRC
END MODULE SCREENC
MODULE TDGAS
  REAL,              ALLOCATABLE, DIMENSION(:)   :: AGASSP, BGASSP, CGASSP, AGASGT, BGASGT, CGASGT
  INTEGER,           ALLOCATABLE, DIMENSION(:)   :: EQSP,   EQGT
END MODULE TDGAS
MODULE LOGICC
  LOGICAL                                        :: SUSP_SOLIDS,        OXYGEN_DEMAND,    UPDATE_GRAPH,     INITIALIZE_GRAPH
  LOGICAL                                        :: WITHDRAWALS,        TRIBUTARIES,      GATES
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: NO_WIND,            NO_INFLOW,        NO_OUTFLOW,       NO_HEAT
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: UPWIND,             ULTIMATE,         FRESH_WATER,      SALT_WATER
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: LIMITING_DLT,       TERM_BY_TERM,     MANNINGS_N,       PH_CALC
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: ONE_LAYER,          DIST_TRIBS,       PRECIPITATION
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: PRINT_SEDIMENT,     LIMITING_FACTOR,  READ_EXTINCTION,  READ_RADIATION
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: UH_INTERNAL,        DH_INTERNAL,      UH_EXTERNAL,      DH_EXTERNAL
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: UQ_INTERNAL,        DQ_INTERNAL,      UQ_EXTERNAL,      DQ_EXTERNAL
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: UP_FLOW,            DN_FLOW,          INTERNAL_FLOW
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: DAM_INFLOW,         DAM_OUTFLOW                                    !TC 08/03/04
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: INTERP_METEOROLOGY, INTERP_INFLOW,    INTERP_DTRIBS,    INTERP_TRIBS
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: INTERP_WITHDRAWAL,  INTERP_HEAD,      INTERP_EXTINCTION
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: VISCOSITY_LIMIT,    CELERITY_LIMIT,   IMPLICIT_AZ,      TRAPEZOIDAL !SW 07/16/04
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: HYDRO_PLOT,         CONSTITUENT_PLOT, DERIVED_PLOT
  LOGICAL,           ALLOCATABLE, DIMENSION(:,:) :: PRINT_DERIVED,      PRINT_HYDRO,      PRINT_CONST,      PRINT_EPIPHYTON
  LOGICAL,           ALLOCATABLE, DIMENSION(:,:) :: POINT_SINK,         INTERNAL_WEIR,    INTERP_OUTFLOW
END MODULE LOGICC
MODULE SHADEC 
  PARAMETER (IANG=18)      
  REAL,PARAMETER                                 :: GAMA=(3.1415926*2.)/REAL(IANG)                         ! SW 10/17/05
  REAL,                           DIMENSION(IANG):: ANG                                                    ! SW 10/17/05
  REAL,              ALLOCATABLE, DIMENSION(:)   :: A00,    DECL,   HH,     TTLB,   TTRB,   CLLB,   CLRB   ! SW 10/17/05
  REAL,              ALLOCATABLE, DIMENSION(:)   :: SRLB1,  SRRB1,  SRLB2,  SRRB2,  SRFJD1, SRFJD2, SHADEI
  REAL,              ALLOCATABLE, DIMENSION(:,:) :: TOPO   
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: DYNAMIC_SHADE
  DATA ANG  /0.00000, 0.34907, 0.69813, 1.04720, 1.39626, 1.74533, 2.09440, 2.44346, &
            2.79253, 3.14159, 3.49066, 3.83972, 4.18879, 4.53786, 4.88692, 5.23599, 5.58505, 5.93412/      ! SW 10/17/05
END MODULE SHADEC
MODULE EDDY
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:)      :: AZC
  REAL,              ALLOCATABLE, DIMENSION(:)      :: WSHY,   FRIC                     
  REAL,              ALLOCATABLE, DIMENSION(:,:)    :: AZT,    DZT                              ! SW 10/20/07
  REAL,              ALLOCATABLE, DIMENSION(:,:)    :: FRICBR, DECAY
  REAL,              ALLOCATABLE, DIMENSION (:,:,:) :: TKE
END MODULE EDDY
! v3.5 start
MODULE MACROPHYTEC
  REAL,    POINTER,               DIMENSION(:,:)     :: NH4MR,  NH4MG,  LDOMMAC, RPOMMAC, LPOMMAC, DOMP, DOMR, TICMC
  REAL,    POINTER,               DIMENSION(:,:)     :: PO4MR,  PO4MG
  REAL,              ALLOCATABLE, DIMENSION(:)       :: MG,     MR,     MM, MMAX,   MBMP 
  REAL,              ALLOCATABLE, DIMENSION(:)       :: MT1,    MT2,    MT3,    MT4,    MK1,    MK2,    MK3,    MK4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: MP,     MN,     MC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: PSED,   NSED,   MHSP,   MHSN,   MHSC,   msat,   exm
  REAL,              ALLOCATABLE, DIMENSION(:)       :: CDdrag, dwv,    dwsa,  anorm    ! cb 6/29/06
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ARMAC  
  REAL,              ALLOCATABLE, DIMENSION(:)       :: O2MG,   O2MR,   LRPMAC,  MPOM
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: MACMBRS,MACMBRT,SSMACMB     
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: CW,     BIC, macwbci
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: MACTRMR,MACTRMF,MACTRM
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: MMR,    MRR  
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: MAC,    MACT
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: MPLIM,  MNLIM, MCLIM
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: SMAC,   SMACT  
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: GAMMAJ
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: MGR  
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: MACRC,  MACRM    
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: MLLIM
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: MACSS
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: SMACRC, SMACRM  
  LOGICAL,           ALLOCATABLE, DIMENSION(:)       :: KTICOL
  LOGICAL,           ALLOCATABLE, DIMENSION(:,:)     :: PRINT_MACROPHYTE, MACROPHYTE_CALC
  LOGICAL                                            :: MACROPHYTE_ON  
  CHARACTER(3),      ALLOCATABLE, DIMENSION(:,:)     :: mprwbc, macwbc
  CHARACTER(10),      ALLOCATABLE, DIMENSION(:,:)    :: CONV2
  CHARACTER(10),     ALLOCATABLE, DIMENSION(:,:,:,:) :: MLFPR      
!  DATA                                                  SAVOLRAT /9000.0/, DEN /6.0E4/   !cb 6/30/06
END MODULE MACROPHYTEC
MODULE POROSITYC
    REAL,              ALLOCATABLE, DIMENSION(:)     :: SAREA, VOLKTI
    REAL,              ALLOCATABLE, DIMENSION(:,:)   :: POR,   VOLI,   VSTEMKT
    REAL,              ALLOCATABLE, DIMENSION(:,:,:) :: VSTEM
    LOGICAL,       ALLOCATABLE, DIMENSION(:)         :: HEAD_FLOW
    LOGICAL,       ALLOCATABLE, DIMENSION(:)         :: UP_HEAD
END MODULE POROSITYC
MODULE ZOOPLANKTONC
  LOGICAL                                            :: ZOOPLANKTON_CALC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: zg,zm,zeff,PREFP,zr,ZOOMIN,ZS2P,EXZ
  REAL,              ALLOCATABLE, DIMENSION(:)       :: Zt1,Zt2,Zt3,Zt4,Zk1,Zk2,Zk3,Zk4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ZP,ZN,ZC,o2zr
    REAL,              ALLOCATABLE, DIMENSION(:,:)   :: PREFA, PREFZ ! OMNIVOROUS ZOOPLANKTON
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: po4zr,NH4ZR,DOZR,TICZR,LPZOOOUT,LPZOOIN
  REAL,    POINTER,               DIMENSION(:,:,:)   :: ZOO, ZOOSS 
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: ZMU,TGRAZE,ZRT,ZMT
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: ZOORM,ZOORMR,ZOORMF
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: agzt
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: AGZ, ZGZ ! OMNIVOROUS ZOOPLANKTON
END MODULE ZOOPLANKTONC
! V3.5 END

! CE-QUAL-W2 computations

PROGRAM CE_QUAL_W2

  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  use macrophytec; use porosityc; use zooplanktonc  !v3.5
  EXTERNAL RESTART_OUTPUT

! Variable declaration

  REAL          :: JDAYTS, JDAY1
  REAL          :: NXTVD
  REAL          :: ICETHU, ICETH1, ICETH2, ICE_TOL
  REAL          :: DLTCAL
  INTEGER       :: CON,    RSI,    W2ERR,  WRN,  GRF,  NDG=16
  integer       :: vsf,    sif !v3.5
  LOGICAL       :: ADD_LAYER,      SUB_LAYER,          WARNING_OPEN,    ERROR_OPEN,      VOLUME_WARNING, SURFACE_WARNING
  LOGICAL       :: END_RUN,        BRANCH_FOUND,       NEW_PAGE,        UPDATE_KINETICS, UPDATE_RATES
  LOGICAL       :: WEIR_CALC,      DERIVED_CALC,       RESTART_IN,      RESTART_OUT
  LOGICAL       :: SPILLWAY,        PIPES,          PUMPS
  LOGICAL       :: TIME_SERIES,    DOWNSTREAM_OUTFLOW, ICE_COMPUTATION, WINTER
  CHARACTER(1)  :: ESC
  CHARACTER(2)  :: DEG
  CHARACTER(3)  :: GDCH
  CHARACTER(8)  :: RSOC,   RSIC,   CCC,   LIMC,   WDOC,   TSRC,   EXT
  CHARACTER(10) :: BLANK,  BLANK1, sedch,   sedpch,   sednch,   sedcch  !v3.5
  CHARACTER(72) :: WDOFN,  RSOFN,  TSRFN, SEGNUM, LINE

! Allocatable array declarations

  REAL,          ALLOCATABLE, DIMENSION(:)     :: ETUGT,  EBUGT,  ETDGT,  EBDGT
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ETUSP,  EBUSP,  ETDSP,  EBDSP
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ETUPI,  EBUPI,  ETDPI,  EBDPI,  ETPU,   EBPU,   TSEDF
  REAL,          ALLOCATABLE, DIMENSION(:)     :: CSUM,   CDSUM  
  REAL,          ALLOCATABLE, DIMENSION(:)     :: RSOD,   RSOF,   DLTD,   DLTF,   DLTMAX, QWDO   
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EPU,    STRTPU, ENDPU,  EONPU,  EOFFPU, QPU
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ICETHI, ALBEDO, HWI,    BETAI,  GAMMAI, ICEMIN, ICET2,  CBHE,   TSED
  REAL,          ALLOCATABLE, DIMENSION(:)     :: FI,     SEDCI,  FSOD,   FSED,   AX,     RAN,    T2I,    ELBOT,  DXI
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QINT,   QOUTT
  REAL,          ALLOCATABLE, DIMENSION(:)     :: WSHX,   SROSH,  EV
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QDT,    QPR,    ICESW,  RS,     RN,     DLXRHO, Q,      QSSUM
  REAL,          ALLOCATABLE, DIMENSION(:)     :: XBR,    QPRBR,  EVBR,   TPB    
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ELTRT,  ELTRB
  REAL,          ALLOCATABLE, DIMENSION(:)     :: TSRD,   TSRF,   WDOD,   WDOF
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QOAVR,  QIMXR,  QOMXR,  QTAVB,  QTMXB
  REAL,          ALLOCATABLE, DIMENSION(:)     :: FETCH,  ETSR
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QINSUM, TINSUM 
  REAL,          ALLOCATABLE, DIMENSION(:)     :: CDTOT
  REAL,          ALLOCATABLE, DIMENSION(:)     :: SEDCIp, sedcin, sedcic, sedcis   !v3.5
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: ESTRT,  WSTRT,  CINSUM
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: P,      HSEG,   QTOT
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: CPB,    COUT,   CWDO,   CDWDO
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: C2I,    EPICI
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: QTRF
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: SNPD,   SCRD,   PRFD,   SPRD,   CPLD,   VPLD,   FLXD
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: SNPF,   SCRF,   PRFF,   SPRF,   CPLF,   VPLF,   FLXF
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: TVP,    SEDVP,  QINF
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: TSSUH1, TSSDH1
  REAL,          ALLOCATABLE, DIMENSION(:,:,:) :: CSSUH1, CSSDH1
  REAL,          ALLOCATABLE, DIMENSION(:,:,:) :: EPIVP,  CVP
  REAL(R8),      ALLOCATABLE, DIMENSION(:)     :: VOLB
  REAL(R8),      ALLOCATABLE, DIMENSION(:)     :: DLVOL,  VOLG
  REAL(R8),      ALLOCATABLE, DIMENSION(:)     :: A,      C,      D,      F,      V,      BTA,    GMA,    BHRHO
  REAL(R8),      ALLOCATABLE, DIMENSION(:)     :: DLVR,   ESR,    ETR
  REAL(R8),      ALLOCATABLE, DIMENSION(:,:)   :: CMBRS
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: KTUGT,  KBUGT,  KTDGT,  KBDGT
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: KTUSP,  KBUSP,  KTDSP,  KBDSP
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: KTUPI,  KBUPI,  KTDPI,  KBDPI
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: NSNP,   NSCR,   NSPR,   NVPL,   NFLX,   NCPL,   BTH
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: VPR,    LPR,    NIPRF,  NISPR,  NPRF
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: IUPU,   IDPU,   KTPU,   KBPU,   JWUPU,  JWDPU,  JBUPU,  JBDPU
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: NISNP
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: NBL,    KBMAX,  KBI
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: KBR,    IBPR
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: TSR
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: NPOINT, NL,     KTQIN,  KBQIN, ilayer    ! SW 1/23/06
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: ITR,    KTTR,   KBTR,   JBTR
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: IWD,    KWD,    JBWD
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: IWDO,   ITSR
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: ILAT,   JBDAM,  JSS    
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)   :: KTSWT,  KBSWT  
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)   :: IPRF,   ISPR,   ISNP,   BL,     WDO,    CDN
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: ALLOW_ICE,      ICE_IN,         PUMPON,        FETCH_CALC
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: DN_HEAD,        HEAD_BOUNDARY   
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: PLACE_QIN,      PLACE_QTR,      SPECIFY_QTR
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: OPEN_VPR,       OPEN_LPR
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: ISO_TEMP,       VERT_TEMP,      LONG_TEMP,     VERT_PROFILE,  LONG_PROFILE
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: SEDIMENT_CALC,  DETAILED_ICE,   IMPLICIT_VISC, SNAPSHOT,      PROFILE
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: VECTOR,         CONTOUR,        SPREADSHEET,   SCREEN_OUTPUT
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: FLUX,           EVAPORATION,    ZERO_SLOPE
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: ISO_SEDIMENT,   VERT_SEDIMENT,  LONG_SEDIMENT
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: VOLUME_BALANCE, ENERGY_BALANCE, MASS_BALANCE
  LOGICAL,       ALLOCATABLE, DIMENSION(:,:)   :: ISO_EPIPHYTON,  VERT_EPIPHYTON, LONG_EPIPHYTON, EPIPHYTON_CALC
  LOGICAL,       ALLOCATABLE, DIMENSION(:,:)   :: ISO_CONC,       VERT_CONC,      LONG_CONC,      TDG_SPILLWAY,   TDG_GATE
  CHARACTER(4),  ALLOCATABLE, DIMENSION(:)     :: CUNIT1
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: SEG,    SEDRC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: HPLTC,  CPLTC,  CDPLTC 
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: EXC,    EXIC   
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: GASGTC, GASSPC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: CWDOC,  CDWDOC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: ICEC,   SEDCc,  SEDPRC, SNPC,   SCRC,   SPRC,   PRFC ! v3.5
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: RHEVC,  VPLC,   CPLC,   AZSLC,  FETCHC 
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: DTRC,   SROC,   KFAC,   CDAC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: INCAC,  TRCAC,  DTCAC,  PRCAC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: WTYPEC, GRIDC                                                        !SW 07/16/04
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: PUSPC,  PDSPC,  PUGTC,  PDGTC,  PDPIC,  PUPIC,  PPUC,   TRC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: SLICEC, FLXC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: VBC,    MBC,    EBC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: PQC,    EVC,    PRC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: QINC,   QOUTC,  WINDC,  HEATC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: VISC,   CELC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: SLTRC,  SLHTC,  FRICC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: QINIC,  TRIC,   DTRIC,  WDIC,   HDIC,   METIC
  CHARACTER(10), ALLOCATABLE, DIMENSION(:)     :: C2CH,   CDCH,   EPCH,   macch  !v3.5
  CHARACTER(45), ALLOCATABLE, DIMENSION(:)     :: KFNAME 
  CHARACTER(72), ALLOCATABLE, DIMENSION(:)     :: SNPFN,  PRFFN,  VPLFN,  CPLFN,  SPRFN,  FLXFN,  BTHFN,  VPRFN,  LPRFN
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:,:)   :: SINKC,  SINKCT 
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:,:)   :: CPRBRC, CDTBRC, CPRWBC, CINBRC, CTRTRC, HPRWBC, STRIC,  CDWBC,  KFWBC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:,:)   :: EPIC,   EPIPRC
  CHARACTER(10), ALLOCATABLE, DIMENSION(:,:)   :: CONV1
  CHARACTER(72), PARAMETER                     :: CONFN='w2_con.npt'
  CHARACTER(72)                                :: TEXT

! Data declarations

  DATA RK1   /2.12/,         RL1    /333507.0/, RIMT /0.0/, RHOA /1.25/, RHOI /916.0/, VTOL /1.0E3/, CP /4186.0/, ICE_TOL /0.005/
  DATA BLANK /'          '/, BLANK1 /'     m    '/
  DATA CON   /10/,  RSI /11/,  WRN /32/, W2ERR /33/
  data thrkti /0.10/  !v3.5

!***********************************************************************************************************************************
!**                                                       Task 1: Inputs                                                          **
!***********************************************************************************************************************************
 
! Open control file

  OPEN (CON,FILE=CONFN,STATUS='OLD',IOSTAT=I)
  IF (I /= 0) THEN
    TEXT = 'Could not open w2_con.npt'
    GO TO 240
  END IF

! Title and array dimensions

  ALLOCATE (TITLE(11))
  READ (CON,'(///(8X,A72))') (TITLE(J),J=1,10)
  READ (CON,'(//8X,4I8)')     NWB, NBR, IMX, KMX
  READ (CON,'(//8X,8I8)')     NTR, NST, NIW, NWD, NGT, NSP, NPI, NPU
  READ (CON,'(//8X,7I8,a8)')  NGC, NSS, NAL, NEP, NBOD, nmc, nzp  !v3.5
  READ (CON,'(//8X,I8)')      NOD

! Constituent numbers
  
  NTDS  = 1
  NGCS  = 2
  NGCE  = NGCS+NGC-1
  NSSS  = NGCE+1
  NSSE  = NSSS+NSS-1
  NPO4  = NSSE+1
  NNH4  = NPO4+1
  NNO3  = NNH4+1
  NDSI  = NNO3+1
  NPSI  = NDSI+1
  NFE   = NPSI+1
  NLDOM = NFE+1
  NRDOM = NLDOM+1
  NLPOM = NRDOM+1
  NRPOM = NLPOM+1
  NBODS = NRPOM+1
  NBODE = NBODS+NBOD-1
  NAS   = NBODE+1
  NAE   = NAS+NAL-1
  NDO   = NAE+1
  NTIC  = NDO+1
  NALK  = NTIC+1

! v3.5 start
  NZOOS = NALK + 1
  NZOOE = NZOOS + NZP - 1
  NLDOMP=NZOOE+1
  NRDOMP=nldomp+1
  NLPOMP=nrdomp+1
  NRPOMP=nlpomp+1
  NLDOMn=nrpomp+1
  NRDOMn=nldomn+1
  NLPOMn=nrdomn+1
  NRPOMn=nlpomn+1
  nct=nrpomn  
! v3.5 end

! Constituent, tributary, and widthdrawal totals

  NTRT = NTR+NGT+NSP+NPI+NPU
  NWDT = NWD+NGT+NSP+NPI+NPU
  NEPT = MAX(NEP,1)
  Nmct = MAX(nmc,1)
  nzpt=max(nzp,1)
  ALLOCATE (CDAC(NDC))
  ALLOCATE (WSC(IMX),    KBI(IMX))
  ALLOCATE (VBC(NWB),    EBC(NWB),    MBC(NWB),    PQC(NWB),    EVC(NWB),    PRC(NWB))
  ALLOCATE (WINDC(NWB),  QINC(NWB),   QOUTC(NWB),  HEATC(NWB),  SLHTC(NWB))
  ALLOCATE (QINIC(NBR),  DTRIC(NBR),  TRIC(NTR),   WDIC(NWD),   HDIC(NBR),   METIC(NWB))
  ALLOCATE (EXC(NWB),    EXIC(NWB))
  ALLOCATE (SLTRC(NWB),  THETA(NWB),  FRICC(NWB),  NAF(NWB),    ELTMF(NWB))
  ALLOCATE (ZMIN(NWB),   IZMIN(NWB))
  ALLOCATE (C2CH(NCT),   CDCH(NDC),   EPCH(NEPT),  macch(nmct))  !v3.5
  ALLOCATE (CPLTC(NCT),  HPLTC(NHY),  CDPLTC(NDC))
  ALLOCATE (CMIN(NCT),   CMAX(NCT),   HYMIN(NHY),  HYMAX(NHY),  CDMIN(NDC),  CDMAX(NDC))
  ALLOCATE (JBDAM(NBR),  ILAT(NWDT))
  ALLOCATE (QINSUM(NBR), TINSUM(NBR), TIND(NBR),   JSS(NBR),    QIND(NBR))
  ALLOCATE (QOLD(NPI),   DTP(NPI),    DTPS(NPI),   QOLDS(NPI))
  ALLOCATE (LATGTC(NGT), LATSPC(NSP), LATPIC(NPI), LATPUC(NPU), DYNGTC(NGT))
  ALLOCATE (OPT(NWB,7),         CIND(NCT,NBR),         CINSUM(NCT,NBR))
  ALLOCATE (CDWBC(NDC,NWB),     KFWBC(NFL,NWB),        CPRWBC(NCT,NWB),    CINBRC(NCT,NBR),     CTRTRC(NCT,NTR))
  ALLOCATE (CDTBRC(NCT,NBR),    CPRBRC(NCT,NBR),       SINKCT(NST,NBR))
  ALLOCATE (STRIC(NST,NBR),     ESTRT(NST,NBR),        WSTRT(NST,NBR),     KTSWT(NST,NBR),      KBSWT(NST,NBR))
  ALLOCATE (YSS(NNPIPE,NPI),    VSS(NNPIPE,NPI),       YS(NNPIPE,NPI),     VS(NNPIPE,NPI),      VSTS(NNPIPE,NPI))
  ALLOCATE (YSTS(NNPIPE,NPI),   YST(NNPIPE,NPI),       VST(NNPIPE,NPI))
  ALLOCATE (CBODD(KMX,IMX,NBOD))
  ALLOCATE (ALLIM(KMX,IMX,NAL), APLIM(KMX,IMX,NAL),    ANLIM(KMX,IMX,NAL), ASLIM(KMX,IMX,NAL))
  ALLOCATE (ELLIM(KMX,IMX,NEP), EPLIM(KMX,IMX,NEP),    ENLIM(KMX,IMX,NEP), ESLIM(KMX,IMX,NEP))
  ALLOCATE (CSSK(KMX,IMX,NCT),  C1(KMX,IMX,NCT),       C2(KMX,IMX,NCT),    CD(KMX,IMX,NDC),     KF(KMX,IMX,NFL))
  ALLOCATE (KFS(KMX,IMX,NFL),   AF(KMX,IMX,NAL,5),     EF(KMX,IMX,NEP,5),  HYD(KMX,IMX,NHY))
  ALLOCATE (TKE(KMX,IMX,2), AZT(KMX,IMX),DZT(KMX,IMX))                                   ! SW 10/20/07
  ALLOCATE (HYDRO_PLOT(NHY),    CONSTITUENT_PLOT(NCT), DERIVED_PLOT(NDC))
  ALLOCATE (ZERO_SLOPE(NWB),    DYNAMIC_SHADE(IMX))
  ALLOCATE (AZSLC(NWB))
  ALLOCATE (NSPRF(NWB))
  ALLOCATE (KBMAX(NWB),  ELKT(NWB),   WIND2(IMX))
  ALLOCATE (VISC(NWB),   CELC(NWB),   REAERC(NWB))
  ALLOCATE (QOAVR(NWB),  QIMXR(NWB),  QOMXR(NWB))
  ALLOCATE (LAT(NWB),    LONGIT(NWB), ELBOT(NWB))
  ALLOCATE (BTH(NWB),    VPR(NWB),    LPR(NWB))
  ALLOCATE (NISNP(NWB),  NIPRF(NWB),  NISPR(NWB))
  ALLOCATE (A00(NWB),    HH(NWB),     DECL(NWB))
  ALLOCATE (T2I(NWB),    KTWB(NWB),   KBR(NWB),    IBPR(NWB))
  ALLOCATE (DLVR(NWB),   ESR(NWB),    ETR(NWB),    NBL(NWB))
  ALLOCATE (LPRFN(NWB),  EXTFN(NWB),  BTHFN(NWB),  METFN(NWB),  VPRFN(NWB))
  ALLOCATE (SNPFN(NWB),  PRFFN(NWB),  SPRFN(NWB),  CPLFN(NWB),  VPLFN(NWB),  FLXFN(NWB))
  ALLOCATE (AFW(NWB),    BFW(NWB),    CFW(NWB),    WINDH(NWB),  RHEVC(NWB),  FETCHC(NWB))
  ALLOCATE (SDK(NWB),    FSOD(NWB),   FSED(NWB),   SEDCI(NWB),  SEDCc(NWB),   SEDPRC(NWB), seds(nwb), sedb(nwb))  !cb 11/28/06  
  ALLOCATE (ICEC(NWB),   SLICEC(NWB), ICETHI(NWB), ALBEDO(NWB), HWI(NWB),    BETAI(NWB),  GAMMAI(NWB), ICEMIN(NWB), ICET2(NWB))
  ALLOCATE (EXH2O(NWB),  BETA(NWB),   EXOM(NWB),   EXSS(NWB),   DXI(NWB),    CBHE(NWB),   TSED(NWB),   TSEDF(NWB),  FI(NWB))
  ALLOCATE (AX(NWB),     WTYPEC(NWB), JBDN(NWB),   AZC(NWB),    AZMAX(NWB),  QINT(NWB),   QOUTT(NWB),  GRIDC(NWB))     !SW 07/14/04
  ALLOCATE (TAIR(NWB),   TDEW(NWB),   WIND(NWB),   PHI(NWB),    CLOUD(NWB),  CSHE(IMX),   SRON(NWB),   RAN(NWB))
  ALLOCATE (SNPC(NWB),   SCRC(NWB),   PRFC(NWB),   SPRC(NWB),   CPLC(NWB),   VPLC(NWB),   FLXC(NWB))
  ALLOCATE (NXTMSN(NWB), NXTMSC(NWB), NXTMPR(NWB), NXTMSP(NWB), NXTMCP(NWB), NXTMVP(NWB), NXTMFL(NWB))
  ALLOCATE (SNPDP(NWB),  SCRDP(NWB),  PRFDP(NWB),  SPRDP(NWB),  CPLDP(NWB),  VPLDP(NWB),  FLXDP(NWB))
  ALLOCATE (NSNP(NWB),   NSCR(NWB),   NPRF(NWB),   NSPR(NWB),   NCPL(NWB),   NVPL(NWB),   NFLX(NWB))
  ALLOCATE (NEQN(NWB),   PO4R(NWB),   PARTP(NWB))
  ALLOCATE (NH4DK(NWB),  NH4R(NWB))
  ALLOCATE (NO3DK(NWB),  NO3S(NWB))
  ALLOCATE (FER(NWB),    FES(NWB))
  ALLOCATE (CO2R(NWB),   SROC(NWB))
  ALLOCATE (O2ER(NEPT),  O2EG(NEPT))
  ALLOCATE (CAQ10(NWB),  CADK(NWB),   CAS(NWB))
  ALLOCATE (BODP(NBOD),  BODN(NBOD),  BODC(NBOD))
  ALLOCATE (KBOD(NBOD),  TBOD(NBOD),  RBOD(NBOD))
  ALLOCATE (LDOMDK(NWB), RDOMDK(NWB), LRDDK(NWB))
  ALLOCATE (OMT1(NWB),   OMT2(NWB),   OMK1(NWB),   OMK2(NWB))
  ALLOCATE (LPOMDK(NWB), RPOMDK(NWB), LRPDK(NWB),  POMS(NWB))
  ALLOCATE (ORGP(NWB),   ORGN(NWB),   ORGC(NWB),   ORGSI(NWB))
  ALLOCATE (RCOEF1(NWB), RCOEF2(NWB), RCOEF3(NWB), RCOEF4(NWB))
  ALLOCATE (NH4T1(NWB),  NH4T2(NWB),  NH4K1(NWB),  NH4K2(NWB))
  ALLOCATE (NO3T1(NWB),  NO3T2(NWB),  NO3K1(NWB),  NO3K2(NWB))
  ALLOCATE (DSIR(NWB),   PSIS(NWB),   PSIDK(NWB),  PARTSI(NWB))
  ALLOCATE (SODT1(NWB),  SODT2(NWB),  SODK1(NWB),  SODK2(NWB))
  ALLOCATE (O2NH4(NWB),  O2OM(NWB))
  ALLOCATE (O2AR(NAL),   O2AG(NAL))
  ALLOCATE (CGQ10(NGC),  CG0DK(NGC),  CG1DK(NGC),  CGS(NGC))
  ALLOCATE (CUNIT(NCT),  CUNIT1(NCT), CUNIT2(NCT))
  ALLOCATE (CAC(NCT),    INCAC(NCT),  TRCAC(NCT),  DTCAC(NCT),  PRCAC(NCT))
  ALLOCATE (CNAME(NCT),  CNAME1(NCT), CNAME2(NCT), CNAME3(NCT), CMULT(NCT),  CSUM(NCT))
  ALLOCATE (CN(NCT))
  ALLOCATE (SSS(NSS),    TAUCR(NSS),  SEDRC(NSS))
  ALLOCATE (CDSUM(NDC))
  ALLOCATE (DTRC(NBR))
  ALLOCATE (NSTR(NBR),   XBR(NBR))
  ALLOCATE (QTAVB(NBR),  QTMXB(NBR))
  ALLOCATE (BS(NWB),     BE(NWB),     JBUH(NBR),   JBDH(NBR),   JWUH(NBR),   JWDH(NBR))
  ALLOCATE (TSSS(NBR),   TSSB(NBR),   TSSICE(NBR))
  ALLOCATE (ESBR(NBR),   ETBR(NBR),   EBRI(NBR))
  ALLOCATE (QIN(NBR),    PR(NBR),     QPRBR(NBR),  QDTR(NBR),   EVBR(NBR))
  ALLOCATE (TIN(NBR),    TOUT(NBR),   TPR(NBR),    TDTR(NBR),   TPB(NBR))
  ALLOCATE (NACPR(NBR),  NACIN(NBR),  NACDT(NBR),  NACTR(NTR),  NACD(NWB))
  ALLOCATE (QSUM(NBR),   NOUT(NBR),   KTQIN(NBR),  KBQIN(NBR),  ELUH(NBR),   ELDH(NBR))
  ALLOCATE (NL(NBR),     NPOINT(NBR), SLOPE(NBR),  ALPHA(NBR),  COSA(NBR),   SINA(NBR), ilayer(imx))    ! SW 1/23/06
  ALLOCATE (CPRFN(NBR),  EUHFN(NBR),  TUHFN(NBR),  CUHFN(NBR),  EDHFN(NBR),  TDHFN(NBR),  QOTFN(NBR),  PREFN(NBR))
  ALLOCATE (QINFN(NBR),  TINFN(NBR),  CINFN(NBR),  CDHFN(NBR),  QDTFN(NBR),  TDTFN(NBR),  CDTFN(NBR),  TPRFN(NBR))
  ALLOCATE (VOLWD(NBR),  VOLSBR(NBR), VOLTBR(NBR), DLVOL(NBR),  VOLG(NWB),   VOLSR(NWB),  VOLTR(NWB),  VOLEV(NBR))
  ALLOCATE (VOLB(NBR),   VOLPR(NBR),  VOLTRB(NBR), VOLDT(NBR),  VOLUH(NBR),  VOLDH(NBR),  VOLIN(NBR),  VOLOUT(NBR))
  ALLOCATE (US(NBR),     DS(NBR),     CUS(NBR),    UHS(NBR),    DHS(NBR),    UQB(NBR),    DQB(NBR),    CDHS(NBR))
  ALLOCATE (TSSEV(NBR),  TSSPR(NBR),  TSSTR(NBR),  TSSDT(NBR),  TSSWD(NBR),  TSSUH(NBR),  TSSDH(NBR),  TSSIN(NBR),  TSSOUT(NBR))
  ALLOCATE (ET(IMX),     RS(IMX),     RN(IMX),     RB(IMX),     RC(IMX),     RE(IMX),     SHADE(IMX))
  ALLOCATE (DLTMAX(NOD), QWDO(IMX),   TWDO(IMX))                                                                        ! SW 1/24/05
  ALLOCATE (SOD(IMX),    ELWS(IMX),   BKT(IMX),    REAER(IMX))
  ALLOCATE (ICETH(IMX),  ICE(IMX),    ICESW(IMX))
  ALLOCATE (Q(IMX),      QC(IMX),     QERR(IMX),   QSSUM(IMX))
  ALLOCATE (KTI(IMX),    SKTI(IMX),   SROSH(IMX),  SEG(IMX),    DLXRHO(IMX))
  ALLOCATE (DLX(IMX),    DLXR(IMX))
  ALLOCATE (A(IMX),      C(IMX),      D(IMX),      F(IMX),      V(IMX),      BTA(IMX),    GMA(IMX))
  ALLOCATE (KBMIN(IMX),  EV(IMX),     QDT(IMX),    QPR(IMX),    SBKT(IMX),   BHRHO(IMX))
  ALLOCATE (SZ(IMX),     WSHX(IMX),   WSHY(IMX),   WIND10(IMX), CZ(IMX),     FETCH(IMX),  PHI0(IMX),   FRIC(IMX))
  ALLOCATE (Z(IMX),      KB(IMX),     PALT(IMX))
  ALLOCATE (VNORM(KMX))
  ALLOCATE (ANPR(NAL),   ANEQN(NAL),  APOM(NAL))
  ALLOCATE (AC(NAL),     ASI(NAL),    ACHLA(NAL),  AHSP(NAL),   AHSN(NAL),   AHSSI(NAL))
  ALLOCATE (AT1(NAL),    AT2(NAL),    AT3(NAL),    AT4(NAL),    AK1(NAL),    AK2(NAL),    AK3(NAL),    AK4(NAL))
  ALLOCATE (AG(NAL),     AR(NAL),     AE(NAL),     AM(NAL),     AS(NAL),     EXA(NAL),    ASAT(NAL),   AP(NAL),   AN(NAL))
  ALLOCATE (ENPR(NEPT),  ENEQN(NEPT))     
  ALLOCATE (EG(NEPT),    ER(NEPT),    EE(NEPT),    EM(NEPT),    EB(NEPT),    ESAT(NEPT),  EP(NEPT),    EN(NEPT))
  ALLOCATE (EC(NEPT),    ESI(NEPT),   ECHLA(NEPT), EHSP(NEPT),  EHSN(NEPT),  EHSSI(NEPT), EPOM(NEPT),  EHS(NEPT))
  ALLOCATE (ET1(NEPT),   ET2(NEPT),   ET3(NEPT),   ET4(NEPT),   EK1(NEPT),   EK2(NEPT),   EK3(NEPT),   EK4(NEPT))
  ALLOCATE (HNAME(NHY),  FMTH(NHY),   HMULT(NHY),  FMTC(NCT),   FMTCD(NDC))
  ALLOCATE (KFAC(NFL),   KFNAME(NFL), KFCN(NFL,NWB))
  ALLOCATE (C2I(NCT,NWB),    TRCN(NCT,NTR))
  ALLOCATE (CDN(NDC,NWB),    CDNAME(NDC),     CDNAME2(NDC),    CDNAME3(NDC),    CDMULT(NDC))
  ALLOCATE (CMBRS(NCT,NBR),  CMBRT(NCT,NBR),  INCN(NCT,NBR),   DTCN(NCT,NBR),   PRCN(NCT,NBR))
  ALLOCATE (FETCHU(IMX,NBR), FETCHD(IMX,NBR))
  ALLOCATE (IPRF(IMX,NWB),   ISNP(IMX,NWB),   ISPR(IMX,NWB),   BL(IMX,NWB))
  ALLOCATE (H1(KMX,IMX),     H2(KMX,IMX),     BH1(KMX,IMX),    BH2(KMX,IMX),    BHR1(KMX,IMX),   BHR2(KMX,IMX),   QTOT(KMX,IMX))
  ALLOCATE (SAVH2(KMX,IMX),  AVH1(KMX,IMX),   AVH2(KMX,IMX),   AVHR(KMX,IMX),   SAVHR(KMX,IMX))
  ALLOCATE (LFPR(KMX,IMX),   BI(KMX,IMX), bnew(kmx,imx))        ! SW 1/23/06
  ALLOCATE (ADX(KMX,IMX),    ADZ(KMX,IMX),    DO1(KMX,IMX),    DO2(KMX,IMX),    DO3(KMX,IMX),    SED(KMX,IMX))
  ALLOCATE (B(KMX,IMX),      CONV(KMX,IMX),   CONV1(KMX,IMX),  EL(KMX,IMX),     DZ(KMX,IMX),     DZQ(KMX,IMX),    DX(KMX,IMX))
  ALLOCATE (P(KMX,IMX),      SU(KMX,IMX),     SW(KMX,IMX),     SAZ(KMX,IMX),    T1(KMX,IMX),     TSS(KMX,IMX),    QSS(KMX,IMX))
  ALLOCATE (BB(KMX,IMX),     BR(KMX,IMX),     BH(KMX,IMX),     BHR(KMX,IMX),    VOL(KMX,IMX),    HSEG(KMX,IMX),   DECAY(KMX,IMX))
  ALLOCATE (DEPTHB(KMX,IMX), DEPTHM(KMX,IMX), FPSS(KMX,IMX),   FPFE(KMX,IMX),   FRICBR(KMX,IMX), UXBR(KMX,IMX),   UYBR(KMX,IMX))
  ALLOCATE (QUH1(KMX,NBR),   QDH1(KMX,NBR),   VOLUH2(KMX,NBR), VOLDH2(KMX,NBR), TUH(KMX,NBR),    TDH(KMX,NBR))
  ALLOCATE (TSSUH1(KMX,NBR), TSSUH2(KMX,NBR), TSSDH1(KMX,NBR), TSSDH2(KMX,NBR))
  ALLOCATE (TVP(KMX,NWB),    SEDVP(KMX,NWB),  H(KMX,NWB))
  ALLOCATE (QINF(KMX,NBR),   QOUT(KMX,NBR),   KOUT(KMX,NBR))
  ALLOCATE (CT(KMX,IMX),     AT(KMX,IMX),     VT(KMX,IMX),     DT(KMX,IMX),     GAMMA(KMX,IMX))
  ALLOCATE (CWDO(NCT,NOD),   CDWDO(NDC,NOD),  CWDOC(NCT),      CDWDOC(NDC),     CDTOT(NDC))
  ALLOCATE (CIN(NCT,NBR),    CDTR(NCT,NBR),   CPR(NCT,NBR),    CPB(NCT,NBR),    COUT(NCT,NBR))
  ALLOCATE (RSOD(NOD),       RSOF(NOD),       DLTD(NOD),       DLTF(NOD))
  ALLOCATE (TSRD(NOD),       TSRF(NOD),       WDOD(NOD),       WDOF(NOD))
  ALLOCATE (SNPD(NOD,NWB),   SNPF(NOD,NWB),   SPRD(NOD,NWB),   SPRF(NOD,NWB))
  ALLOCATE (SCRD(NOD,NWB),   SCRF(NOD,NWB),   PRFD(NOD,NWB),   PRFF(NOD,NWB))
  ALLOCATE (CPLD(NOD,NWB),   CPLF(NOD,NWB),   VPLD(NOD,NWB),   VPLF(NOD,NWB),   FLXD(NOD,NWB),   FLXF(NOD,NWB))
  ALLOCATE (EPIC(NWB,NEPT),  EPICI(NWB,NEPT), EPIPRC(NWB,NEPT))     
  ALLOCATE (EPIVP(KMX,NWB,NEP))     
  ALLOCATE (CUH(KMX,NCT,NBR),     CDH(KMX,NCT,NBR))
  ALLOCATE (EPM(KMX,IMX,NEPT),    EPD(KMX,IMX,NEPT),    EPC(KMX,IMX,NEPT))
  ALLOCATE (C1S(KMX,IMX,NCT),     CSSB(KMX,IMX,NCT),    CVP(KMX,NCT,NWB))
  ALLOCATE (CSSUH1(KMX,NCT,NBR),  CSSUH2(KMX,NCT,NBR),  CSSDH2(KMX,NCT,NBR), CSSDH1(KMX,NCT,NBR))
  ALLOCATE (OPEN_VPR(NWB),        OPEN_LPR(NWB))
  ALLOCATE (READ_EXTINCTION(NWB), READ_RADIATION(NWB))
  ALLOCATE (DIST_TRIBS(NBR),      LIMITING_FACTOR(NAL))
  ALLOCATE (UPWIND(NWB),          ULTIMATE(NWB))
  ALLOCATE (FRESH_WATER(NWB),     SALT_WATER(NWB),      TRAPEZOIDAL(NWB))                                              !SW 07/16/04
  ALLOCATE (UH_EXTERNAL(NBR),     DH_EXTERNAL(NBR),     UH_INTERNAL(NBR),    DH_INTERNAL(NBR))
  ALLOCATE (UQ_EXTERNAL(NBR),     DQ_EXTERNAL(NBR),     UQ_INTERNAL(NBR),    DQ_INTERNAL(NBR))
  ALLOCATE (UP_FLOW(NBR),         DN_FLOW(NBR),         UP_HEAD(NBR),        DN_HEAD(NBR))
  ALLOCATE (INTERNAL_FLOW(NBR),   DAM_INFLOW(NBR),      DAM_OUTFLOW(NBR),    HEAD_FLOW(NBR),      HEAD_BOUNDARY(NWB))  !TC 08/03/04
  ALLOCATE (ISO_CONC(NCT,NWB),    VERT_CONC(NCT,NWB),   LONG_CONC(NCT,NWB))
  ALLOCATE (ISO_SEDIMENT(NWB),    VERT_SEDIMENT(NWB),   LONG_SEDIMENT(NWB))
  ALLOCATE (VISCOSITY_LIMIT(NWB), CELERITY_LIMIT(NWB),  IMPLICIT_AZ(NWB))
  ALLOCATE (FETCH_CALC(NWB),      ONE_LAYER(IMX),       IMPLICIT_VISC(NWB))
  ALLOCATE (LIMITING_DLT(NWB),    TERM_BY_TERM(NWB),    MANNINGS_N(NWB))
  ALLOCATE (PLACE_QIN(NWB),       PLACE_QTR(NTRT),      SPECIFY_QTR(NTRT))
  ALLOCATE (PRINT_CONST(NCT,NWB), PRINT_HYDRO(NHY,NWB), PRINT_SEDIMENT(NWB))
  ALLOCATE (VOLUME_BALANCE(NWB),  ENERGY_BALANCE(NWB),  MASS_BALANCE(NWB))
  ALLOCATE (DETAILED_ICE(NWB),    ICE_CALC(NWB),        ICE_IN(NBR),          ALLOW_ICE(IMX))
  ALLOCATE (EVAPORATION(NWB),     PRECIPITATION(NWB),   RH_EVAP(NWB),         PH_CALC(NWB))
  ALLOCATE (NO_INFLOW(NWB),       NO_OUTFLOW(NWB),      NO_HEAT(NWB),         NO_WIND(NWB))
  ALLOCATE (ISO_TEMP(NWB),        VERT_TEMP(NWB),       LONG_TEMP(NWB),       VERT_PROFILE(NWB),  LONG_PROFILE(NWB))
  ALLOCATE (SNAPSHOT(NWB),        PROFILE(NWB),         VECTOR(NWB),          CONTOUR(NWB),       SPREADSHEET(NWB))
  ALLOCATE (SCREEN_OUTPUT(NWB),   FLUX(NWB))
  ALLOCATE (PRINT_DERIVED(NDC,NWB),  PRINT_EPIPHYTON(NWB,NEPT))
  ALLOCATE (SEDIMENT_CALC(NWB),      EPIPHYTON_CALC(NWB,NEPT), SEDIMENT_RESUSPENSION(NSS))
  ALLOCATE (TDG_SPILLWAY(NWDT,NSP),  TDG_GATE(NWDT,NGT),       INTERNAL_WEIR(KMX,IMX))     
  ALLOCATE (ISO_EPIPHYTON(NWB,NEPT), VERT_EPIPHYTON(NWB,NEPT), LONG_EPIPHYTON(NWB,NEPT))     
  ALLOCATE (LATERAL_SPILLWAY(NSP),   LATERAL_GATE(NGT),        LATERAL_PUMP(NPU),        LATERAL_PIPE(NPI))
  ALLOCATE (INTERP_HEAD(NBR),        INTERP_WITHDRAWAL(NWD),   INTERP_EXTINCTION(NWB),   INTERP_DTRIBS(NBR))
  ALLOCATE (INTERP_OUTFLOW(NST,NBR), INTERP_INFLOW(NBR),       INTERP_METEOROLOGY(NWB),  INTERP_TRIBS(NTR))
  ALLOCATE (LNAME(NCT+NHY+NDC))
  ALLOCATE (IWR(NIW),    KTWR(NIW),   KBWR(NIW))
  ALLOCATE (JWUSP(NSP),  JWDSP(NSP),  QSP(NSP))
  ALLOCATE (KTWD(NWDT),  KBWD(NWDT),  JBWD(NWDT))
  ALLOCATE (GTA1(NGT),   GTB1(NGT),   GTA2(NGT),   GTB2(NGT))
  ALLOCATE (BGT(NGT),    IUGT(NGT),   IDGT(NGT),   EGT(NGT))
  ALLOCATE (QTR(NTRT),   TTR(NTRT),   KTTR(NTRT),  KBTR(NTRT))
  ALLOCATE (AGASGT(NGT), BGASGT(NGT), CGASGT(NGT), GASGTC(NGT))
  ALLOCATE (PUGTC(NGT),  ETUGT(NGT),  EBUGT(NGT),  KTUGT(NGT),  KBUGT(NGT))
  ALLOCATE (PDGTC(NGT),  ETDGT(NGT),  EBDGT(NGT),  KTDGT(NGT),  KBDGT(NGT))
  ALLOCATE (A1GT(NGT),   B1GT(NGT),   G1GT(NGT),   A2GT(NGT),   B2GT(NGT),   G2GT(NGT))
  ALLOCATE (EQGT(NGT),   JBUGT(NGT),  JBDGT(NGT),  JWUGT(NGT),  JWDGT(NGT),  QGT(NGT))
  ALLOCATE (JBUPI(NPI),  JBDPI(NPI),  JWUPI(NPI),  JWDPI(NPI),  QPI(NPI))
  ALLOCATE (IUPI(NPI),   IDPI(NPI),   EUPI(NPI),   EDPI(NPI),   WPI(NPI),    DLXPI(NPI),  FPI(NPI),    FMINPI(NPI), PUPIC(NPI))
  ALLOCATE (ETUPI(NPI),  EBUPI(NPI),  KTUPI(NPI),  KBUPI(NPI),  PDPIC(NPI),  ETDPI(NPI),  EBDPI(NPI),  KTDPI(NPI),  KBDPI(NPI))
  ALLOCATE (PUSPC(NSP),  ETUSP(NSP),  EBUSP(NSP),  KTUSP(NSP),  KBUSP(NSP),  PDSPC(NSP),  ETDSP(NSP),  EBDSP(NSP))
  ALLOCATE (KTDSP(NSP),  KBDSP(NSP),  IUSP(NSP),   IDSP(NSP),   ESP(NSP),    A1SP(NSP),   B1SP(NSP),   A2SP(NSP))
  ALLOCATE (B2SP(NSP),   AGASSP(NSP), BGASSP(NSP), CGASSP(NSP), EQSP(NSP),   GASSPC(NSP), JBUSP(NSP),  JBDSP(NSP))
  ALLOCATE (IUPU(NPU),   IDPU(NPU),   EPU(NPU),    STRTPU(NPU), ENDPU(NPU),  EONPU(NPU),  EOFFPU(NPU), QPU(NPU),   PPUC(NPU))
  ALLOCATE (ETPU(NPU),   EBPU(NPU),   KTPU(NPU),   KBPU(NPU),   JWUPU(NPU),  JWDPU(NPU),  JBUPU(NPU),  JBDPU(NPU), PUMPON(NPU))
  ALLOCATE (IWD(NWDT),   KWD(NWDT),   QWD(NWDT),   EWD(NWDT),   KTW(NWDT),   KBW(NWDT))
  ALLOCATE (ITR(NTRT),   QTRFN(NTR),  TTRFN(NTR),  CTRFN(NTR),  ELTRT(NTRT), ELTRB(NTRT), TRC(NTRT),   JBTR(NTRT), QTRF(KMX,NTRT))
  ALLOCATE (TTLB(IMX),   TTRB(IMX),   CLLB(IMX),   CLRB(IMX))     
  ALLOCATE (SRLB1(IMX),  SRRB1(IMX),  SRLB2(IMX),  SRRB2(IMX),  SRFJD1(IMX), SHADEI(IMX), SRFJD2(IMX))
  ALLOCATE (TOPO(IMX,IANG))                                                                                        ! SW 10/17/05
  ALLOCATE (QSW(KMX,NWDT),  CTR(NCT,NTRT), HPRWBC(NHY,NWB))
  ALLOCATE (RATZ(KMX,NWB),   CURZ1(KMX,NWB),  CURZ2(KMX,NWB),   CURZ3(KMX,NWB))   ! SW 5/15/06
! v3.5 start
  ALLOCATE (zg(NZPt),zm(NZPt),zeff(NZPt),prefp(NZPt),zr(NZPt),zoomin(NZPt),zs2p(NZPt),exz(NZPt),PREFZ(NZPt,nzpt))
  ALLOCATE (zt1(NZPt),zt2(NZPt),zt3(NZPt),zt4(NZPt),zk1(NZPt),zk2(NZPt),zk3(NZPt),zk4(NZPt),o2zr(nzpt))
  ALLOCATE (zp(NZPt),zn(NZPt),zc(NZPt))
  allocate (prefa(nal,nzpt))
  allocate (po4zr(kmx,imx),nh4zr(kmx,imx))
  allocate (zmu(kmx,imx,nzp),tgraze(kmx,imx,nzp),zrt(kmx,imx,nzp),zmt(kmx,imx,nzp)) ! MLM POINTERS:,zoo(kmx,imx,NZP),zooss(kmx,imx,NZP))
  allocate (zoorm(kmx,imx,nzp),zoormr(kmx,imx,nzp),zoormf(kmx,imx,nzp))
  allocate (lpzooout(kmx,imx),lpzooin(kmx,imx),dozr(kmx,imx),ticzr(kmx,imx))
  allocate (agz(kmx,imx,nal,nzp), zgz(kmx,imx,nzp,nzp),agzt(kmx,imx,nal)) !omnivorous zooplankton
  allocate (ORGPLD(kmx,imx), ORGPRD(kmx,imx), ORGPLP(kmx,imx), ORGPRP(kmx,imx), ORGNLD(kmx,imx), ORGNRD(kmx,imx), ORGNLP(kmx,imx))
  allocate (ORGNRP(kmx,imx))
  allocate (ldompmp(kmx,imx),ldomnmp(kmx,imx),lpompmp(kmx,imx),lpomnmp(kmx,imx),rpompmp(kmx,imx),rpomnmp(kmx,imx))
  allocate (lpzooinp(kmx,imx),lpzooinn(kmx,imx),lpzoooutp(kmx,imx),lpzoooutn(kmx,imx))
  allocate (SEDVPp(KMX,NWB),SEDVPc(KMX,NWB),SEDVPn(KMX,NWB))
  allocate (sedp(kmx,imx),sedn(kmx,imx),sedc(kmx,imx))
  allocate (sdkv(kmx,imx),seddktot(kmx,imx))
  allocate (sedcip(nwb),sedcin(nwb),sedcic(nwb),sedcis(nwb))
  ALLOCATE (cbods(NBOD), cbodns(kmx,imx), sedcb(kmx,imx), sedcbp(kmx,imx), sedcbn(kmx,imx), sedcbc(kmx,imx))
  allocate  (print_macrophyte(nwb,nmct), macrophyte_calc(nwb,nmct),macwbc(nwb,nmct),conv2(kmx,kmx),mprwbc(nwb,nmct))    
  allocate  (mac(kmx,imx,nmct), macrc(kmx,kmx,imx,nmct),mact(kmx,kmx,imx), macrm(kmx,kmx,imx,nmct), macss(kmx,kmx,imx,nmct))
  allocate  (mgr(kmx,kmx,imx,nmct),mmr(kmx,imx,nmct), mrr(kmx,imx,nmct))
  allocate  (smacrc(kmx,kmx,imx,nmct), smacrm(kmx,kmx,imx,nmct))
  allocate  (smact(kmx,kmx,imx), smac(kmx,imx,nmct))
  allocate  (mt1(nmct),mt2(nmct),mt3(nmct),mt4(nmct),mk1(nmct),mk2(nmct),mk3(nmct),mk4(nmct),mg(nmct),mr(nmct),mm(nmct))
  allocate  (mbmp(nmct), mmax(nmct), cddrag(nmct), dwv(nmct), dwsa(nmct), anorm(nmct))
  allocate  (mp(nmct), mn(nmct), mc(nmct),psed(nmct),nsed(nmct),mhsp(nmct),mhsn(nmct),mhsc(nmct),msat(nmct),exm(nmct))
  allocate  (O2MG(nmct), O2MR(nmct),  LRPMAC(nmct),  MPOM(nmct))
  allocate  (kticol(imx),armac(imx),macwbci(nwb,nmct))
  allocate  (macmbrs(nbr,nmct), macmbrt(nbr,nmct),ssmacmb(nbr,nmct))
  allocate  (cw(kmx,imx), bic(kmx,imx))  
  allocate  (mactrmr(kmx,imx,nmct), mactrmf(kmx,imx,nmct),mactrm(kmx,imx,nmct))     
  allocate  (mlfpr(kmx,kmx,imx,nmct))
  allocate  (mllim(kmx,kmx,imx,nmct), mplim(kmx,imx,nmct),mclim(kmx,imx,nmct),mnlim(kmx,imx,nmct))
  ALLOCATE  (GAMMAj(kmx,KMX,IMX))	  
  allocate (por(kmx,imx),VOLKTi(imx),VOLi(Kmx,Imx),vstem(kmx,imx,nmct),vstemkt(imx,nmct),sarea(nmct))  
  ALLOCATE (IWIND(NWB))
! v3.5 end

! Allocate subroutine variables

  CALL TRANSPORT
  CALL KINETICS
  CALL WATERBODY
  CALL OPEN_CHANNEL_INITIALIZE
  CALL PIPE_FLOW_INITIALIZE

! State variables

  TDS  => C2(:,:,1);         PO4  => C2(:,:,NPO4);      NH4  => C2(:,:,NNH4);        NO3  => C2(:,:,NNO3);   DSI  => C2(:,:,NDSI)
  PSI  => C2(:,:,NPSI);      FE   => C2(:,:,NFE);       LDOM => C2(:,:,NLDOM);       RDOM => C2(:,:,NRDOM);  LPOM => C2(:,:,NLPOM)
  RPOM => C2(:,:,NRPOM);     O2   => C2(:,:,NDO);       TIC  => C2(:,:,NTIC);        ALK  => C2(:,:,NALK)
  CG   => C2(:,:,NGCS:NGCE); SS   => C2(:,:,NSSS:NSSE); CBOD => C2(:,:,NBODS:NBODE); ALG  => C2(:,:,NAS:NAE)
! v3.5 start
  ZOO  => C2(:,:,NZOOS:NZOOE)
  LDOMP  => C2(:,:,nldomp); RDOMP  => C2(:,:,nrdomp); LPOMP  => C2(:,:,nlpomp); RPOMP  => C2(:,:,nrpomp)
  LDOMN  => C2(:,:,nldomN); RDOMN  => C2(:,:,nrdomn); LPOMN  => C2(:,:,nlpomn); RPOMN  => C2(:,:,nrpomn)  
! v3.5 end

! State variable source/sinks

  CGSS   => CSSK(:,:,NGCS:NGCE);   SSSS   => CSSK(:,:,NSSS:NSSE); PO4SS  => CSSK(:,:,NPO4);  NH4SS  => CSSK(:,:,NNH4)
  NO3SS  => CSSK(:,:,NNO3);        DSISS  => CSSK(:,:,NDSI);      PSISS  => CSSK(:,:,NPSI);  FESS   => CSSK(:,:,NFE)
  LDOMSS => CSSK(:,:,NLDOM);       RDOMSS => CSSK(:,:,NRDOM);     LPOMSS => CSSK(:,:,NLPOM); RPOMSS => CSSK(:,:,NRPOM)
  CBODSS => CSSK(:,:,NBODS:NBODE); ASS    => CSSK(:,:,NAS:NAE);   DOSS   => CSSK(:,:,NDO);   TICSS  => CSSK(:,:,NTIC)
! v3.5 start
  zooss  => cssk(:,:,NZOOS:NZOOE)
  LDOMPSS  => cssk(:,:,nldomp); RDOMPSS  => cssk(:,:,nrdomp); LPOMPSS  => cssk(:,:,nlpomp); RPOMPSS  => cssk(:,:,nrpomp)
  LDOMNSS  => cssk(:,:,nldomN); RDOMNSS  => cssk(:,:,nrdomn); LPOMNSS  => cssk(:,:,nlpomn); RPOMNSS  => cssk(:,:,nrpomn)
! v3.5 end

! Derived variables

  DOC   => CD(:,:,1);  POC  => CD(:,:,2);  TOC  => CD(:,:,3);  DON  => CD(:,:,4);  PON   => CD(:,:,5);  TON  => CD(:,:,6)
  TKN   => CD(:,:,7);  TN   => CD(:,:,8);  DOP  => CD(:,:,9);  POP  => CD(:,:,10); TOP   => CD(:,:,11); TP   => CD(:,:,12)
  APR   => CD(:,:,13); CHLA => CD(:,:,14); ATOT => CD(:,:,15); O2DG => CD(:,:,16); TOTSS => CD(:,:,17); TISS => CD(:,:,18)
  CBODU => CD(:,:,19); PH   => CD(:,:,20); CO2  => CD(:,:,21); HCO3 => CD(:,:,22); CO3   => CD(:,:,23)
 
! Kinetic fluxes

  SSSI   => KF(:,:,1);  SSSO   => KF(:,:,2);  PO4AR  => KF(:,:,3);  PO4AG  => KF(:,:,4);  PO4AP  => KF(:,:,5)
  PO4ER  => KF(:,:,6);  PO4EG  => KF(:,:,7);  PO4EP  => KF(:,:,8);  PO4POM => KF(:,:,9);  PO4DOM => KF(:,:,10)
  PO4OM  => KF(:,:,11); PO4SD  => KF(:,:,12); PO4SR  => KF(:,:,13); PO4NS  => KF(:,:,14); NH4D   => KF(:,:,15)
  NH4AR  => KF(:,:,16); NH4AG  => KF(:,:,17); NH4AP  => KF(:,:,18); NH4ER  => KF(:,:,19); NH4EG  => KF(:,:,20)
  NH4EP  => KF(:,:,21); NH4POM => KF(:,:,22); NH4DOM => KF(:,:,23); NH4OM  => KF(:,:,24); NH4SD  => KF(:,:,25)
  NH4SR  => KF(:,:,26); NO3D   => KF(:,:,27); NO3AG  => KF(:,:,28); NO3EG  => KF(:,:,29); NO3SED => KF(:,:,30)
  DSIAG  => KF(:,:,31); DSIEG  => KF(:,:,32); DSID   => KF(:,:,33); DSISD  => KF(:,:,34); DSISR  => KF(:,:,35)
  DSIS   => KF(:,:,36); PSIAM  => KF(:,:,37); PSINS  => KF(:,:,38); PSID   => KF(:,:,39); FENS   => KF(:,:,40)
  FESR   => KF(:,:,41); LDOMD  => KF(:,:,42); LRDOMD => KF(:,:,43); RDOMD  => KF(:,:,44); LDOMAP => KF(:,:,45)
  LDOMEP => KF(:,:,46); LPOMD  => KF(:,:,47); LRPOMD => KF(:,:,48); RPOMD  => KF(:,:,49); LPOMAP => KF(:,:,50)
  LPOMEP => KF(:,:,51); LPOMNS => KF(:,:,52); RPOMNS => KF(:,:,53); CBODDK => KF(:,:,54); DOAP   => KF(:,:,55)
  DOEP   => KF(:,:,56); DOAR   => KF(:,:,57); DOER   => KF(:,:,58); DOPOM  => KF(:,:,59); DODOM  => KF(:,:,60)
  DOOM   => KF(:,:,61); DONIT  => KF(:,:,62); DOBOD  => KF(:,:,63); DOAE   => KF(:,:,64); DOSED  => KF(:,:,65)
  DOSOD  => KF(:,:,66); TICAP  => KF(:,:,67); TICEP  => KF(:,:,68); SEDD   => KF(:,:,69); SEDAS  => KF(:,:,70)
  SEDOMS => KF(:,:,71); SEDNS  => KF(:,:,72); SODD   => KF(:,:,73)
! v3.5 start
  LDOMPAP => KF(:,:,74); LDOMPeP => KF(:,:,75); LPOMpAP => KF(:,:,76); LPOMPNS => KF(:,:,77); RPOMPNS => KF(:,:,78)
  LDOMnAP => KF(:,:,79); LDOMneP => KF(:,:,80); LPOMnAP => KF(:,:,81); LPOMnNS => KF(:,:,82); RPOMnNS => KF(:,:,83) 
  SEDDp   => KF(:,:,84); SEDASp  => KF(:,:,85); SEDOMSp => KF(:,:,86); SEDNSp  => KF(:,:,87); lpomepp => KF(:,:,88) 
  SEDDn   => KF(:,:,89); SEDASn  => KF(:,:,90); SEDOMSn => KF(:,:,91); SEDNSn  => KF(:,:,92); lpomepn => KF(:,:,93) 
  SEDDc   => KF(:,:,94); SEDASc  => KF(:,:,95); SEDOMSc => KF(:,:,96); SEDNSc  => KF(:,:,97); lpomepc => KF(:,:,98)
  sedno3  => KF(:,:,99) 
  po4mr   => KF(:,:,100);po4mg   => KF(:,:,101); nh4mr   => KF(:,:,102); nh4mg => KF(:,:,103); ldommac => KF(:,:,104)
  rpommac => KF(:,:,105);lpommac => KF(:,:,106); domp    => KF(:,:,107); domr  => KF(:,:,108); ticmc   => KF(:,:,109)
  cbodns  => KF(:,:,110);sedcb   => KF(:,:,111); sedcbp  => KF(:,:,112); sedcbn => KF(:,:,113); sedcbc  => KF(:,:,114)
  sedbr   => KF(:,:,115);sedbrp  => KF(:,:,116); sedbrn  => KF(:,:,117); sedbrc  => KF(:,:,118)
! v3.5 end

! Algal rate variables

  AGR => AF(:,:,:,1); ARR => AF(:,:,:,2); AER => AF(:,:,:,3); AMR => AF(:,:,:,4); ASR => AF(:,:,:,5)
  EGR => EF(:,:,:,1); ERR => EF(:,:,:,2); EER => EF(:,:,:,3); EMR => EF(:,:,:,4); EBR => EF(:,:,:,5)

! Hydrodynamic variables

  DLTLIM => HYD(:,:,1);  U   => HYD(:,:,2);  W    => HYD(:,:,3); T2   => HYD(:,:,4);  RHO => HYD(:,:,5);  AZ  => HYD(:,:,6)
  VSH    => HYD(:,:,7);  ST  => HYD(:,:,8);  SB   => HYD(:,:,9); ADMX => HYD(:,:,10); DM  => HYD(:,:,11); HDG => HYD(:,:,12)
  ADMZ   => HYD(:,:,13); HPG => HYD(:,:,14); GRAV => HYD(:,:,15)

! I/O units

  SNP => OPT(:,1); PRF => OPT(:,2); VPL => OPT(:,3); CPL => OPT(:,4); SPR => OPT(:,5); FLX => OPT(:,6)

! Zero variables

  ITR  = 0;   JBTR = 0;   KTTR = 0;   KBTR = 0;   QTR  = 0.0; TTR  = 0.0; CTR  = 0.0; QTRF = 0.0; SNPD  = 0.0; TSRD  = 0.0
  PRFD = 0.0; SPRD = 0.0; CPLD = 0.0; VPLD = 0.0; SCRD = 0.0; FLXD = 0.0; WDOD = 0.0; RSOD = 0.0; ELTRB = 0.0; ELTRT = 0.0

! Input file unit numbers

  NUNIT = 40
  DO JW=1,NWB
    BTH(JW) = NUNIT
    VPR(JW) = NUNIT+1
    LPR(JW) = NUNIT+2
    NUNIT   = NUNIT+3
  END DO
  GRF = NUNIT; NUNIT = NUNIT+1

! Time control cards

  READ (CON,'(//8X,2F8.0,I8)')         TMSTRT,   TMEND,    YEAR
  READ (CON,'(//8X,I8,8F8.0)')         NDLT,     DLTMIN
  READ (CON,'(//(:8X,9F8.0))')        (DLTD(J),            J =1,NDLT)
  READ (CON,'(//(:8X,9F8.0))')        (DLTMAX(J),          J =1,NDLT)
  READ (CON,'(//(:8X,9F8.0))')        (DLTF(J),            J =1,NDLT)
  READ (CON,'(//(8X,2A8))')           (VISC(JW), CELC(JW), JW=1,NWB)

! Grid definition cards

  READ (CON,'(//(8X,7I8,F8.3))')      (US(JB),  DS(JB),     UHS(JB),   DHS(JB), UQB(JB), DQB(JB),  NL(JB), SLOPE(JB), JB=1,NBR)
  READ (CON,'(//(8X,3F8.0,3I8))')     (LAT(JW), LONGIT(JW), ELBOT(JW), BS(JW),  BE(JW),  JBDN(JW),                    JW=1,NWB)

! Initial condition cards

  READ (CON,'(//(8X,2F8.0,2A8))')     (T2I(JW),    ICETHI(JW),  WTYPEC(JW),  GRIDC(JW),                               JW=1,NWB)
  READ (CON,'(//(8X,6A8))')           (VBC(JW),    EBC(JW),     MBC(JW),     PQC(JW),   EVC(JW),   PRC(JW),           JW=1,NWB)
  READ (CON,'(//(8X,4A8))')           (WINDC(JW),  QINC(JW),    QOUTC(JW),   HEATC(JW),                               JW=1,NWB)
  READ (CON,'(//(8X,3A8))')           (QINIC(JB),  DTRIC(JB),   HDIC(JB),                                             JB=1,NBR)
  READ (CON,'(//(8X,5A8,4F8.0))')     (SLHTC(JW),  SROC(JW),    RHEVC(JW),   METIC(JW), FETCHC(JW), AFW(JW),                       &
                                       BFW(JW),    CFW(JW),     WINDH(JW),                                            JW=1,NWB)
  READ (CON,'(//(8X,2A8,6F8.0))')     (ICEC(JW),   SLICEC(JW),  ALBEDO(JW),  HWI(JW),   BETAI(JW),  GAMMAI(JW),                    &
                                       ICEMIN(JW), ICET2(JW),                                                         JW=1,NWB)
  READ (CON,'(//(8X,A8,F8.0))')       (SLTRC(JW),  THETA(JW),                                                         JW=1,NWB)
  READ (CON,'(//(8X,6F8.0,A8))')      (AX(JW),     DXI(JW),     CBHE(JW),    TSED(JW),  FI(JW),     TSEDF(JW),                     &
                                       FRICC(JW),                                                                     JW=1,NWB)
  READ (CON,'(//(8X,2A8,2F8.0))')     (AZC(JW),    AZSLC(JW),   AZMAX(JW),   PHISET,                                  JW=1,NWB)

! Inflow-outflow cards

  READ (CON,'(//(8X,I8))')            (NSTR(JB),      JB=1,NBR)
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9A8)')            (STRIC(JS,JB),  JS=1,NSTR(JB))
  END DO
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9I8)')            (KTSWT(JS,JB), JS=1,NSTR(JB))
  END DO
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9I8)')            (KBSWT(JS,JB), JS=1,NSTR(JB))
  END DO
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9A8)')            (SINKCT(JS,JB),JS=1,NSTR(JB))
  END DO
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9F8.0)')          (ESTRT(JS,JB), JS=1,NSTR(JB))
  END DO
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9F8.0)')          (WSTRT(JS,JB), JS=1,NSTR(JB))
  END DO
  READ (CON,'(//(:8X,2I8,6F8.0,A8))') (IUPI(JP),   IDPI(JP),   EUPI(JP),   EDPI(JP),    WPI(JP),                                   &
                                       DLXPI(JP),  FPI(JP),    FMINPI(JP), LATPIC(JP),              JP=1,NPI)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PUPIC(JP),  ETUPI(JP),  EBUPI(JP),  KTUPI(JP),   KBUPI(JP),  JP=1,NPI)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PDPIC(JP),  ETDPI(JP),  EBDPI(JP),  KTDPI(JP),   KBDPI(JP),  JP=1,NPI)
  READ (CON,'(//(:8X,2I8,5F8.0,A8))') (IUSP(JS),   IDSP(JS),   ESP(JS),    A1SP(JS),    B1SP(JS),                                  &
                                       A2SP(JS),   B2SP(JS),   LATSPC(JS),                          JS=1,NSP)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PUSPC(JS),  ETUSP(JS),  EBUSP(JS),  KTUSP(JS),   KBUSP(JS),  JS=1,NSP)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PDSPC(JS),  ETDSP(JS),  EBDSP(JS),  KTDSP(JS),   KBDSP(JS),  JS=1,NSP)
  READ (CON,'(//(:8X,A8,I8,3F8.0))')  (GASSPC(JS), EQSP(JS),   AGASSP(JS), BGASSP(JS),  CGASSP(JS), JS=1,NSP)
  READ (CON,'(//(:8X,2I8,7F8.0,A8))') (IUGT(JG),   IDGT(JG),   EGT(JG),    A1GT(JG),    B1GT(JG),                                  &
                                       G1GT(JG),   A2GT(JG),   B2GT(JG),   G2GT(JG),    LATGTC(JG), JG=1,NGT)
  READ (CON,'(//(:8X,4F8.0,A8))')     (GTA1(JG),   GTB1(JG),   GTA2(JG),   GTB2(JG),    DYNGTC(JG), JG=1,NGT)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PUGTC(JG),  ETUGT(JG),  EBUGT(JG),  KTUGT(JG),   KBUGT(JG),  JG=1,NGT)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PDGTC(JG),  ETDGT(JG),  EBDGT(JG),  KTDGT(JG),   KBDGT(JG),  JG=1,NGT)
  READ (CON,'(//(:8X,A8,I8,3F8.0))')  (GASGTC(JG), EQGT(JG),   AGASGT(JG), BGASGT(JG),  CGASGT(JG), JG=1,NGT)
  READ (CON,'(//(:8X,2I8,6F8.0,A8))') (IUPU(JP),   IDPU(JP),   EPU(JP),    STRTPU(JP),  ENDPU(JP),                                 &
                                       EONPU(JP),  EOFFPU(JP), QPU(JP),    LATPUC(JP),              JP=1,NPU)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PPUC(JP),   ETPU(JP),   EBPU(JP),   KTPU(JP),    KBPU(JP),   JP=1,NPU)
  READ (CON,'(//(:8X,9I8))')          (IWR(JW),    JW=1,NIW)
  READ (CON,'(//(:8X,9I8))')          (KTWR(JW),   JW=1,NIW)
  READ (CON,'(//(:8X,9I8))')          (KBWR(JW),   JW=1,NIW)
  READ (CON,'(//(:8X,9A8))')          (WDIC(JW),   JW=1,NWD)
  READ (CON,'(//(:8X,9I8))')          (IWD(JW),    JW=1,NWD)
  READ (CON,'(//(:8X,9F8.0))')        (EWD(JW),    JW=1,NWD)
  READ (CON,'(//(:8X,9I8))')          (KTWD(JW),   JW=1,NWD)
  READ (CON,'(//(:8X,9I8))')          (KBWD(JW),   JW=1,NWD)
  READ (CON,'(//(:8X,9A8))')          (TRC(JT),    JT=1,NTR)
  READ (CON,'(//(:8X,9A8))')          (TRIC(JT),   JT=1,NTR)
  READ (CON,'(//(:8X,9I8))')          (ITR(JT),    JT=1,NTR)
  READ (CON,'(//(:8X,9F8.0))')        (ELTRT(JT),  JT=1,NTR)
  READ (CON,'(//(:8X,9F8.0))')        (ELTRB(JT),  JT=1,NTR)
  READ (CON,'(//(8X,A8))')            (DTRC(JB),   JB=1,NBR)

! Output control cards (excluding constituents)

  READ (CON,'(/)')
  DO JH=1,NHY
    READ (CON,'(:8X,9A8)')            (HPRWBC(JH,JW),JW=1,NWB)
  END DO
  READ (CON,'(//(8X,A8,2I8))')        (SNPC(JW), NSNP(JW), NISNP(JW), JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SNPD(J,JW),J=1,NSNP(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SNPF(J,JW),J=1,NSNP(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9I8)')            (ISNP(I,JW),I=1,NISNP(JW))
  END DO
  READ (CON,'(//(8X,A8,I8))')         (SCRC(JW), NSCR(JW), JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SCRD(J,JW),J=1,NSCR(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SCRF(J,JW),J=1,NSCR(JW))
  END DO
  READ (CON,'(//(8X,A8,2I8))')        (PRFC(JW), NPRF(JW), NIPRF(JW), JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (PRFD(J,JW),J=1,NPRF(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (PRFF(J,JW),J=1,NPRF(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9I8)')            (IPRF(J,JW),J=1,NIPRF(JW))
  END DO
  READ (CON,'(//(8X,A8,2I8))')        (SPRC(JW), NSPR(JW), NISPR(JW), JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SPRD(J,JW),J=1,NSPR(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SPRF(J,JW),J=1,NSPR(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9I8)')            (ISPR(J,JW), J=1,NISPR(JW))
  END DO
  READ (CON,'(//(8X,A8,I8))')         (VPLC(JW),  NVPL(JW),  JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (VPLD(J,JW), J=1,NVPL(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (VPLF(J,JW), J=1,NVPL(JW))
  END DO
  READ (CON,'(//(8X,A8,I8))')         (CPLC(JW),   NCPL(JW), JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (CPLD(J,JW), J=1,NCPL(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (CPLF(J,JW), J=1,NCPL(JW))
  END DO
  READ (CON,'(//(8X,A8,I8))')         (FLXC(JW),   NFLX(JW), JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (FLXD(J,JW), J=1,NFLX(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (FLXF(J,JW), J=1,NFLX(JW))
  END DO
  READ (CON,'(//8X,A8,2I8)')           TSRC,    NTSR,    NIKTSR; ALLOCATE (ITSR(MAX(1,NIKTSR)), ETSR(MAX(1,NIKTSR)))
  READ (CON,'(//(:8X,9F8.0))')        (TSRD(J), J=1,NTSR)
  READ (CON,'(//(:8X,9F8.0))')        (TSRF(J), J=1,NTSR)
  READ (CON,'(//(:8X,9I8))')          (ITSR(J), J=1,NIKTSR)
  READ (CON,'(//(:8X,9F8.0))')        (ETSR(J), J=1,NIKTSR)
  READ (CON,'(//8X,A8,2I8)')           WDOC,    NWDO,    NIWDO;  ALLOCATE (IWDO(MAX(1,NIWDO)))
  READ (CON,'(//(:8X,9F8.0))')        (WDOD(J), J=1,NWDO)
  READ (CON,'(//(:8X,9F8.0))')        (WDOF(J), J=1,NWDO)
  READ (CON,'(//(8X,9I8))')           (IWDO(J), J=1,NIWDO)
  READ (CON,'(//8X,A8,I8,A8)')         RSOC,    NRSO,    RSIC
  READ (CON,'(//(:8X,9F8.0))')        (RSOD(J), J=1,NRSO)
  READ (CON,'(//(:8X,9F8.0))')        (RSOF(J), J=1,NRSO)

! Constituent control cards

  READ (CON,'(//8X,2A8,I8)')           CCC, LIMC, CUF
  READ (CON,'(//(2A8))')              (CNAME2(JC),  CAC(JC),      JC=1,NCT)
  READ (CON,'(/)')
  DO JD=1,NDC
    READ (CON,'(A8,(:9A8))')           CDNAME2(JD),(CDWBC(JD,JW), JW=1,NWB)
  END DO
  READ (CON,'(/)')
!  DO JF=1,NFL
  do jf=1,73   ! Fix this later
    READ (CON,'(:8X,9A8)')            (KFWBC(JF,JW),  JW=1,NWB)
  END DO
  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(:8X,9F8.0)')          (C2I(JC,JW),    JW=1,NWB)
  END DO
  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(:8X,9A8)')            (CPRWBC(JC,JW), JW=1,NWB)
  END DO
  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(:8X,9A8)')            (CINBRC(JC,JB), JB=1,NBR)
  END DO
  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(:8X,9A8)')            (CTRTRC(JC,JT), JT=1,NTR)
  END DO
  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(:8X,9A8)')            (CDTBRC(JC,JB), JB=1,NBR)
  END DO
  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(:8X,9A8)')            (CPRBRC(JC,JB), JB=1,NBR)
  END DO

! Kinetics coefficients

  READ (CON,'(//(8X,4F8.0,2A8))')     (EXH2O(JW),  EXSS(JW),   EXOM(JW),   BETA(JW),   EXC(JW),   EXIC(JW),    JW=1,NWB)
  READ (CON,'(//(8X,9F8.0))')         (EXA(JA),                                                                JA=1,NAL)
  READ (CON,'(//(8X,9F8.0))')         (EXZ(Jz),                                                                Jz=1,Nzpt)  !v3.5
  READ (CON,'(//(8X,9F8.0))')         (EXM(Jm),                                                                Jm=1,nmct)  !v3.5
  READ (CON,'(//(8X,4F8.0))')         (CGQ10(JG),  CG0DK(JG),  CG1DK(JG),  CGS(JG),                            JG=1,NGC)
  READ (CON,'(//(8X,F8.0,A,F8.0))')   (SSS(JS),    SEDRC(JS),  TAUCR(JS),                                      JS=1,NSS)
  READ (CON,'(//(8X,9F8.0))')         (AG(JA),     AR(JA),     AE(JA),     AM(JA),     AS(JA),                                     &
                                       AHSP(JA),   AHSN(JA),   AHSSI(JA),  ASAT(JA),                           JA=1,NAL)
  READ (CON,'(//(8X,8F8.0))')         (AT1(JA),    AT2(JA),    AT3(JA),    AT4(JA),    AK1(JA),   AK2(JA),                         &
                                       AK3(JA),    AK4(JA),                                                    JA=1,NAL)
  READ (CON,'(//(8X,6F8.0,I8,F8.0))') (AP(JA),     AN(JA),     AC(JA),     ASI(JA),    ACHLA(JA), APOM(JA),                        &
                                       ANEQN(JA),  ANPR(JA),   JA=1,NAL)      
  READ (CON,'(//(8X,9A8))')           (EPIC(JW,1),                                                             JW=1,NWB)
  DO JE=2,NEPT
    READ (CON,'(8X,9A8)')             (EPIC(JW,JE),                                                            JW=1,NWB)
  END DO
  READ (CON,'(//(8X,9A8))')           (EPIPRC(JW,1),                                                           JW=1,NWB)
  DO JE=2,NEPT
    READ (CON,'(8X,9A8)')             (EPIPRC(JW,JE),                                                          JW=1,NWB)
  END DO
  READ (CON,'(//(8X,9F8.0))')         (EPICI(JW,1),                                                            JW=1,NWB)
  DO JE=2,NEPT
    READ (CON,'(8X,9F8.0)')           (EPICI(JW,JE),                                                           JW=1,NWB)
  END DO
  READ (CON,'(//(8X,8F8.0))')         (EG(JE),     ER(JE),     EE(JE),     EM(JE),     EB(JE),    EHSP(JE),                        &
                                       EHSN(JE),   EHSSI(JE),                                                  JE=1,NEP)
  READ (CON,'(//(8X,2F8.0,I8,F8.0))') (ESAT(JE),   EHS(JE),    ENEQN(JE),  ENPR(JE),                           JE=1,NEP)
  READ (CON,'(//(8X,8F8.0))')         (ET1(JE),    ET2(JE),    ET3(JE),    ET4(JE),    EK1(JE),   EK2(JE),                         &
                                       EK3(JE),    EK4(JE),                                                    JE=1,NEP)
  READ (CON,'(//(8X,6F8.0))')         (EP(JE),     EN(JE),     EC(JE),     ESI(JE),    ECHLA(JE), EPOM(JE),    JE=1,NEP)
! v3.5 start
  READ (CON,'(//(8X,7F8.0))')         (zg(jz),zr(jz),zm(jz),zeff(jz),PREFP(jz),ZOOMIN(jz),ZS2P(jz),            Jz=1,Nzpt)
  
  READ (CON,'(//(8X,8F8.0))')         (PREFA(ja,1),                                                            Ja=1,nal)          ! MM 7/13/06                            
  do jz=2,nzpt
    READ (CON,'((8X,8F8.0))')       (PREFA(ja,jz),                                                           Ja=1,nal)                                       
  end do
  READ (CON,'(//(8X,8F8.0))')       (PREFz(jjz,1),                                                          Jjz=1,nzpt)                                       
  do jz=2,nzpt
    READ (CON,'((8X,8F8.0))')       (PREFz(jjz,jz),                                                          Jjz=1,nzpt)           ! MM 7/13/06                             
  end do
  READ (CON,'(//(8X,8F8.0))')         (zT1(Jz),    zT2(Jz),    zT3(Jz),    zT4(Jz),    zK1(Jz),   zK2(Jz),                         &
                                       zK3(Jz),    zK4(Jz),                                                    Jz=1,Nzpt)
  READ (CON,'(//(8X,3F8.0))')         (zP(Jz),     zN(Jz),     zC(Jz),                                         Jz=1,Nzpt)
  READ (CON,'(//(8X,9A8))')           (macwbC(JW,1),                                                           JW=1,NWB)
  DO Jm=2,nmct
    READ (CON,'(8X,9A8)')             (macwbC(JW,Jm),                                                          JW=1,NWB)
  END DO
  READ (CON,'(//(8X,9A8))')           (mprwbC(JW,1),                                                           JW=1,NWB)
  DO Jm=2,nmct
    READ (CON,'(8X,9A8)')             (mprwbC(JW,Jm),                                                          JW=1,NWB)
  END DO
  READ (CON,'(//(8X,9F8.0))')         (macwbCI(JW,1),                                                          JW=1,NWB)
  DO Jm=2,nmct
    READ (CON,'(8X,9F8.0)')           (macwbcI(JW,Jm),                                                         JW=1,NWB)
  END DO
  READ (CON,'(//(8X,9F8.0))')         (mG(jm), mR(jm), mM(jm), msat(jm),mhsp(jm),mhsn(jm),mhsc(jm),                           &
                                          mpom(jm),lrpmac(jm),     jm=1,nmct)
  READ (CON,'(//(8X,2F8.0))')         (psed(jm), nsed(jm),                                                     jm=1,nmct)
  READ (CON,'(//(8X,2F8.0))')         (mbmp(jm), mmax(jm),                                                     jm=1,nmct)
  READ (CON,'(//(8X,4F8.0))')         (cddrag(jm),dwv(jm),dwsa(jm),anorm(jm),                                 jm=1,nmct)  !cb 6/29/06
  READ (CON,'(//(8X,8F8.0))')         (mT1(Jm),    mT2(Jm),    mT3(Jm),    mT4(Jm),    mK1(Jm),   mK2(Jm),                         &
                                       mK3(Jm),    mK4(Jm),                                                    Jm=1,nmct)
  READ (CON,'(//(8X,3F8.0))')         (mP(Jm),     mN(Jm),     mC(Jm),                                         Jm=1,nmct)
  READ (CON,'(//(8X,3F8.0))')         (LDOMDK(JW), RDOMDK(JW), LRDDK(JW),                                      JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (LPOMDK(JW), RPOMDK(JW), LRPDK(JW),  POMS(JW),                           JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (ORGP(JW),   ORGN(JW),   ORGC(JW),   ORGSI(JW),                          JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (OMT1(JW),   OMT2(JW),   OMK1(JW),   OMK2(JW),                           JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (KBOD(JB),   TBOD(JB),   RBOD(JB), cbods(jb),                           JB=1,NBOD)  !v3.5
  READ (CON,'(//(8X,3F8.0))')         (BODP(JB),   BODN(JB),   BODC(JB),                                       JB=1,NBOD)
  READ (CON,'(//(8X,2F8.0))')         (PO4R(JW),   PARTP(JW),                                                  JW=1,NWB)
  READ (CON,'(//(8X,2F8.0))')         (NH4R(JW),   NH4DK(JW),                                                  JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (NH4T1(JW),  NH4T2(JW),  NH4K1(JW),  NH4K2(JW),                          JW=1,NWB)
  READ (CON,'(//(8X,2F8.0))')         (NO3DK(JW),  NO3S(JW),                                                   JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (NO3T1(JW),  NO3T2(JW),  NO3K1(JW),  NO3K2(JW),                          JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (DSIR(JW),   PSIS(JW),   PSIDK(JW),  PARTSI(JW),                         JW=1,NWB)
  READ (CON,'(//(8X,2F8.0))')         (FER(JW),    FES(JW),                                                    JW=1,NWB)
  READ (CON,'(//(8X,F8.0))')          (CO2R(JW),                                                               JW=1,NWB)
  READ (CON,'(//(8X,2F8.0))')         (O2NH4(JW),  O2OM(JW),                                                   JW=1,NWB)
  READ (CON,'(//(8X,2F8.0))')         (O2AR(JA),   O2AG(JA),                                                   JA=1,NAL)
  READ (CON,'(//(8X,2F8.0))')         (O2ER(JE),   O2EG(JE),                                                   JE=1,NEPT)
  READ (CON,'(//(8X,F8.0))')          (O2zR(Jz),                                                               Jz=1,Nzpt)
  READ (CON,'(//(8X,2F8.0))')         (O2mR(Jm),   O2mG(jm),                                                   Jm=1,nmct)
  READ (CON,'(//(8X,F8.0))')           KDO
  READ (CON,'(//(8X,2A8,6F8.0))')     (SEDCc(JW),   SEDPRC(JW), SEDCI(JW),  SDK(JW), seds(jw),   FSOD(JW),   FSED(JW),  &
                                       sedb(jw),   JW=1,NWB)  ! cb 11/28/06  
  READ (CON,'(//(8X,4F8.0))')         (SODT1(JW),  SODT2(JW),  SODK1(JW),  SODK2(JW),                          JW=1,NWB)
  READ (CON,'(//(8X,9F8.0))')         (SOD(I),                                                                  I=1,IMX)
  READ (CON,'(//(8X,A8,I8,4F8.2))')   (REAERC(JW), NEQN(JW),   RCOEF1(JW), RCOEF2(JW), RCOEF3(JW), RCOEF4(JW), JW=1,NWB)

! Input filenames

  READ (CON,'(//(8X,A72))')  RSIFN
  READ (CON,'(//(8X,A72))')  QWDFN
  READ (CON,'(//(8X,A72))')  QGTFN
  READ (CON,'(//(8X,A72))')  WSCFN
  READ (CON,'(//(8X,A72))')  SHDFN  
  READ (CON,'(//(8X,A72))') (BTHFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (METFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (EXTFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (VPRFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (LPRFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (QINFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (TINFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (CINFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (QOTFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (QTRFN(JT), JT=1,NTR)
  READ (CON,'(//(8X,A72))') (TTRFN(JT), JT=1,NTR)
  READ (CON,'(//(8X,A72))') (CTRFN(JT), JT=1,NTR)
  READ (CON,'(//(8X,A72))') (QDTFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (TDTFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (CDTFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (PREFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (TPRFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (CPRFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (EUHFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (TUHFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (CUHFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (EDHFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (TDHFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (CDHFN(JB), JB=1,NBR)
 
! Output filenames
 
  READ (CON,'(//(8X,A72))') (SNPFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (PRFFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (VPLFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (CPLFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (SPRFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (FLXFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))')  TSRFN
  READ (CON,'(//(8X,A72))')  WDOFN
  CLOSE (CON)

! Bathymetry file
 
  DO JW=1,NWB
    OPEN (BTH(JW),FILE=BTHFN(JW),STATUS='OLD')
    READ (BTH(JW),*)
    READ (BTH(JW),'(//(10F8.0))') (DLX(I),  I=US(BS(JW))-1,DS(BE(JW))+1)
    READ (BTH(JW),'(//(10F8.0))') (ELWS(I), I=US(BS(JW))-1,DS(BE(JW))+1)
    READ (BTH(JW),'(//(10F8.0))') (PHI0(I), I=US(BS(JW))-1,DS(BE(JW))+1)
    READ (BTH(JW),'(//(10F8.0))') (FRIC(I), I=US(BS(JW))-1,DS(BE(JW))+1)
    READ (BTH(JW),'(//(10F8.0))') (H(K,JW), K=1,KMX)
    DO I=US(BS(JW))-1,DS(BE(JW))+1
      READ (BTH(JW),'(//(10F8.0))') (B(K,I), K=1,KMX)
      H2(:,I) = H(:,JW)
    END DO
    CLOSE (BTH(JW))
  END DO
  H1 = H2
  BI = B
  
! Output file unit numbers

  ALLOCATE (TSR(NIKTSR))
  ALLOCATE (WDO(NIWDO,4))
  DO J=1,7
    DO JW=1,NWB
      OPT(JW,J) = NUNIT; NUNIT = NUNIT+1
    END DO
  END DO
  DO J=1,NIKTSR
    TSR(J) = NUNIT; NUNIT = NUNIT+1
  END DO
  DO JW=1,NIWDO
    WDO(JW,1) = NUNIT; NUNIT = NUNIT+1
    WDO(JW,2) = NUNIT; NUNIT = NUNIT+1
    WDO(JW,3) = NUNIT; NUNIT = NUNIT+1
    WDO(JW,4) = NUNIT; NUNIT = NUNIT+1
  END DO

! Variable names, formats, multipliers, and Compaq Visual FORTRAN array viewer controls

  OPEN (GRF,FILE='graph.npt',STATUS='OLD')
  READ (GRF,'(///(A43,1X,A9,3F8.0,A8))') (HNAME(J),  FMTH(J),  HMULT(J),  HYMIN(J), HYMAX(J), HPLTC(J), J=1,NHY)
  READ (GRF,'(// (A43,1X,A9,3F8.0,A8))') (CNAME(J),  FMTC(J),  CMULT(J),  CMIN(J),  CMAX(J),  CPLTC(J), J=1,NCT)
  READ (GRF,'(// (A43,1X,A9,3F8.0,A8))') (CDNAME(J), FMTCD(J), CDMULT(J), CDMIN(J), CDMAX(J), CDPLTC(J),J=1,NDC)
  CLOSE (GRF)
  DO JC=1,NCT
    L3         = 1
    L1         = SCAN (CNAME(JC),',')+2
    L2         = SCAN (CNAME(JC)(L1:43),'  ')+L1
    CUNIT(JC)  = CNAME(JC)(L1:L2)
    CNAME1(JC) = CNAME(JC)(1:L1-3)
    CNAME3(JC) = CNAME1(JC)
    DO WHILE (L3 < L1-3)
      IF (CNAME(JC)(L3:L3) == ' ') CNAME3(JC)(L3:L3) = '_'
      L3 = L3+1
    END DO
    CUNIT1(JC) = CUNIT(JC)(1:1)
    CUNIT2(JC) = CUNIT(JC)
    IF (CUNIT(JC)(1:2) == 'mg') THEN
      CUNIT1(JC) = 'g'
      CUNIT2(JC) = 'g/m^3'
    END IF
    IF (CUNIT(JC)(1:2) /= 'g/' .AND. CUNIT(JC)(1:2) /= 'mg') CUNIT1(JC) = '  '
  END DO
  DO JC=1,NDC
    L1          = 1
    L2          = MAX(4,SCAN (CDNAME(JC),',')-1)
    CDNAME3(JC) = CDNAME(JC)(1:L2)
    DO WHILE (L1 < L2)
      IF (CDNAME(JC)(L1:L1) == ' ') CDNAME3(JC)(L1:L1) = '_'
      L1 = L1+1
    END DO
  END DO
  FMTH(1:NHY) = ADJUSTL (FMTH(1:NHY))

! Initialize logical control variables

  VERT_PROFILE = .FALSE.
  LONG_PROFILE = .FALSE.
  CONSTITUENTS =  CCC  == '      ON'
  RESTART_IN   =  RSIC == '      ON'
  DO JW=1,NWB
    ISO_TEMP(JW)         = T2I(JW)     >=  0
    VERT_TEMP(JW)        = T2I(JW)     == -1
    LONG_TEMP(JW)        = T2I(JW)     <  -1
! v3.5 start
    ISO_SEDIMENT(JW)     = SEDCI(JW)   >=  0   .AND. SEDCc(JW)   == '      ON'
    VERT_SEDIMENT(JW)    = SEDCI(JW)   == -1.0 .AND. SEDCc(JW)   == '      ON'
    LONG_SEDIMENT(JW)    = SEDCI(JW)   <  -1.0 .AND. SEDCc(JW)   == '      ON'
! v3.5 end
    ISO_EPIPHYTON(JW,:)  = EPICI(JW,:) >=  0   .AND. EPIC(JW,:) == '      ON'
    VERT_EPIPHYTON(JW,:) = EPICI(JW,:) == -1.0 .AND. EPIC(JW,:) == '      ON'
    LONG_EPIPHYTON(JW,:) = EPICI(JW,:) <  -1.0 .AND. EPIC(JW,:) == '      ON'
    DO JC=1,NCT
      ISO_CONC(JC,JW)  = C2I(JC,JW) >=  0.0
      VERT_CONC(JC,JW) = C2I(JC,JW) == -1.0 .AND. CAC(JC) == '      ON'
      LONG_CONC(JC,JW) = C2I(JC,JW) <  -1.0 .AND. CAC(JC) == '      ON'
      IF (VERT_CONC(JC,JW)) VERT_PROFILE(JW) = .TRUE.
      IF (LONG_CONC(JC,JW)) LONG_PROFILE(JW) = .TRUE.
    END DO
    IF (VERT_TEMP(JW))             VERT_PROFILE(JW) = .TRUE.
    IF (LONG_TEMP(JW))             LONG_PROFILE(JW) = .TRUE.
    IF (VERT_SEDIMENT(JW))         VERT_PROFILE(JW) = .TRUE.
    IF (LONG_SEDIMENT(JW))         LONG_PROFILE(JW) = .TRUE.
    IF (ANY(VERT_EPIPHYTON(JW,:))) VERT_PROFILE(JW) = .TRUE.
    IF (ANY(LONG_EPIPHYTON(JW,:))) LONG_PROFILE(JW) = .TRUE.
    do m=1,nmc
      macrophyte_calc(jw,m) = constituents.and.macwbc(jw,m).eq.' ON'
      print_macrophyte(jw,m) = macrophyte_calc(jw,m).AND.mprwbc(jw,m).EQ.' ON'
      if(macrophyte_calc(jw,m))macrophyte_on=.true.
    end do
  END DO

! Restart data

  JDAY = TMSTRT
  IF (RESTART_IN) THEN
    VERT_PROFILE = .FALSE.      
    LONG_PROFILE = .FALSE.      
    OPEN  (RSI,FILE=RSIFN,FORM='UNFORMATTED',STATUS='OLD')
    READ  (RSI) NIT,    NV,     KMIN,   IMIN,   NSPRF,  CMBRT,  ZMIN,   IZMIN,  START,  CURRENT
    READ  (RSI) DLTDP,  SNPDP,  TSRDP,  VPLDP,  PRFDP,  CPLDP,  SPRDP,  RSODP,  SCRDP,  FLXDP,  WDODP
    READ  (RSI) JDAY,   YEAR,   ELTM,   ELTMF,  DLT,    DLTAV,  DLTS,   MINDLT, JDMIN,  CURMAX
    READ  (RSI) NXTMSN, NXTMTS, NXTMPR, NXTMCP, NXTMVP, NXTMRS, NXTMSC, NXTMSP, NXTMFL, NXTMWD
    READ  (RSI) VOLIN,  VOLOUT, VOLUH,  VOLDH,  VOLPR,  VOLTRB, VOLDT,  VOLWD,  VOLEV,  VOLSBR, VOLTR, VOLSR
    READ  (RSI) TSSEV,  TSSPR,  TSSTR,  TSSDT,  TSSWD,  TSSIN,  TSSOUT, TSSS,   TSSB,   TSSICE
    READ  (RSI) TSSUH,  TSSDH,  TSSUH2, TSSDH2, CSSUH2, CSSDH2, VOLUH2, VOLDH2, QUH1
    READ  (RSI) ESBR,   ETBR,   EBRI
    READ  (RSI) Z,      SZ,     ELWS,   SAVH2,  SAVHR,  H2
    READ  (RSI) KTWB,   KTI,    SKTI,   SBKT   
    READ  (RSI) ICE,    ICETH,  CUF,    QSUM
    READ  (RSI) U,      W,      SU,     SW,     AZ,     SAZ,    DLTLIM
    READ  (RSI) T1,     T2,     C1,     C2,     C1S,    EPD,    SED,    KFS,    CSSK
    READ  (RSI) SEDC, SEDN, SEDP, ZOO, CD  ! mlm 10/06
    READ  (RSI) sdkv                       ! cb 11/30/06
    READ  (RSI) TKE                        ! sw 10/5/07
    CLOSE (RSI)
  END IF

! Open error and warning files

  OPEN (W2ERR,FILE='w2.err',STATUS='UNKNOWN'); OPEN (WRN,FILE='w2.wrn',STATUS='UNKNOWN')

!***********************************************************************************************************************************
!**                                             Task 1.1: Variable Initialization                                                 **
!***********************************************************************************************************************************
 
!***********************************************************************************************************************************
!**                                                 Task 1.1.1: Zero Variables                                                    **
!***********************************************************************************************************************************

  KB     = 0;   KBR    = 0;   NAC    = 0;   NTAC   = 0;   NACD   = 0;   NACIN  = 0;   NACTR  = 0;   NACDT  = 0;   NACPR  = 0
  NDSP   = 0;   HMAX   = 0;   KBMAX  = 0;   DLXMAX = 0;   KBQIN  = 0;   KTQIN  = 0
  NAF    = 0;   TISS   = 0.0; VSS    = 0.0; VS     = 0.0; YS     = 0.0; YSS    = 0.0; YST    = 0.0; YSTS   = 0.0; VST    = 0.0
  VSTS   = 0.0; DTP    = 0.0; DTPS   = 0.0; QOLD   = 0.0; QOLDS  = 0.0; CSHE   = 0.0; CIN    = 0.0; TIN    = 0.0; EV     = 0.0
  QDTR   = 0.0; DZ     = 0.0; ET     = 0.0; CSHE   = 0.0; A      = 0.0; F      = 0.0; D      = 0.0; C      = 0.0; ELTMF  = 0.0
  EL     = 0.0; DX     = 0.0; ST     = 0.0; SB     = 0.0; VS     = 0.0; YS     = 0.0; DZQ    = 0.0; TSS    = 0.0; TTR    = 0.0
  HSEG   = 0.0; CTR    = 0.0; QSS    = 0.0; HPG    = 0.0; HDG    = 0.0; VSH    = 0.0; QDH1   = 0.0; ADMX   = 0.0; DECAY  = 0.0
  ADMZ   = 0.0; UYBR   = 0.0; GRAV   = 0.0; FETCH  = 0.0; FETCHU = 0.0; FETCHD = 0.0; DLTTVD = 0.0; ICETHU = 0.0; ICETH1 = 0.0
  ICETH2 = 0.0; P      = 0.0; CELRTY = 0.0; TAU1   = 0.0; TAU2   = 0.0; VOLSR  = 0.0; VOLTR  = 0.0; AF     = 0.0; EF     = 0.0
  ELTMS  = 0.0; TDTR   = 0.0; DM     = 0.0; QIN    = 0.0; REAER  = 0.0; ST     = 0.0; SB     = 0.0; ADMX   = 0.0; QWD    = 0.0
  ADMZ   = 0.0; HPG    = 0.0; HDG    = 0.0; RHO    = 0.0; JDAYTS = 0.0; JDAY1  = 0.0; DEPTHB = 0.0; DEPTHM = 0.0; UXBR   = 0.0
  BHRHO  = 0.0; DLMR   = 0.0; SRON   = 0.0; CSSB   = 0.0; Q      = 0.0; BH1    = 0.0; BH2    = 0.0; BHR1   = 0.0; BHR2   = 0.0
  AVHR   = 0.0; GRAV   = 0.0; kbp = 0; DZT=0.0; AZT=0.0   ! SW 1/23/06 10/20/07
  IF (.NOT. RESTART_IN) THEN
    NSPRF  = 0;   IZMIN  = 0;   KTWB   = 2;   KMIN   = 1;   IMIN   = 1
    T1     = 0.0; T2     = 0.0; C1     = 0.0; C2     = 0.0; CD     = 0.0; CIN    = 0.0; C1S    = 0.0; KF     = 0.0; CMBRT  = 0.0
    KFS    = 0.0; U      = 0.0; W      = 0.0; SU     = 0.0; SW     = 0.0; SAZ    = 0.0; AZ     = 0.0; ESBR   = 0.0; EPD    = 0.0
    ETBR   = 0.0; EBRI   = 0.0; DLTLIM = 0.0; VOLEV  = 0.0; VOLPR  = 0.0; VOLDT  = 0.0; VOLWD  = 0.0; CURRENT= 0.0; EPI    = 0.0
    VOLUH  = 0.0; VOLDH  = 0.0; VOLIN  = 0.0; VOLOUT = 0.0; VOLSBR = 0.0; VOLTRB = 0.0; TSSS   = 0.0; TSSB   = 0.0; EF     = 0.0
    TSSEV  = 0.0; TSSPR  = 0.0; TSSTR  = 0.0; TSSDT  = 0.0; TSSWD  = 0.0; TSSUH  = 0.0; TSSDH  = 0.0; TSSIN  = 0.0; CSSK   = 0.0
    TSSOUT = 0.0; TSSICE = 0.0; TSSUH1 = 0.0; TSSUH2 = 0.0; CSSUH1 = 0.0; CSSUH2 = 0.0; TSSDH1 = 0.0; TSSDH2 = 0.0; QC     = 0.0
    CSSDH1 = 0.0; CSSDH2 = 0.0; QIND   = 0.0; TIND   = 0.0; CIND   = 0.0; SAVH2  = 0.0; SAVHR  = 0.0; VOLUH2 = 0.0
    AVH1   = 0.0; AVH2   = 0.0; VOLDH2 = 0.0; ADMZ   = 0.0; Z      = 0.0; QUH1   = 0.0; SED    = 0.0
    ZMIN   = -1000.0
    TKE=0.0                                                     ! SG 10/4/07
! v3.5 start
    sedp=0.0;sedc=0.0;sedn=0.0
    macmbrt=0.0; macrc=0.0;smacrc=0.0;mac=0.0;smac=0.0;macrm=0.0
    kticol=.false.
! v3.5 end
  END IF
  ANLIM = 1.0; APLIM = 1; ASLIM = 1.0; ALLIM = 1.0; ENLIM = 1.0; EPLIM = 1; ESLIM = 1.0; ELLIM = 1.0; KLOC = 1; ILOC = 1
! v3.5 start
  mNLIM = 1.0; mPLIM = 1; mcLIM = 1.0; mLLIM = 1.0
! v3.5 end
  ICESW      = 1.0
  HMIN       = 1.0E10
  DLXMIN     = 1.0E10
  LFPR       = BLANK
  CONV       = BLANK
  CONV1      = BLANK1
  CNAME2     = ADJUSTR(CNAME2)
  CDNAME2    = ADJUSTR(CDNAME2)
  TITLE(11)  = ' '
  TEXT       = ' '
  IF (.NOT. CONSTITUENTS) THEN
    NAL = 0; NEP = 0; NSS = 0; NBOD = 0;
  END IF
  DO JW=1,NWB
    GAMMA(:,US(BS(JW)):DS(BE(JW))) = EXH2O(JW)
  END DO

!***********************************************************************************************************************************
!**                                            Task 1.1.2: Miscellaneous Variables                                                **
!***********************************************************************************************************************************

! Logical controls

  NEW_PAGE              = .TRUE.;  VOLUME_WARNING = .TRUE.;  INITIALIZE_GRAPH = .TRUE.;  UPDATE_GRAPH    = .TRUE.
  ICE                   = .FALSE.; FLUX           = .FALSE.; PUMPON           = .FALSE.; WINTER          = .FALSE.
  TDG_GATE              = .FALSE.; TDG_SPILLWAY   = .FALSE.; INTERNAL_WEIR    = .FALSE.; SURFACE_WARNING = .FALSE.
  WARNING_OPEN          = .FALSE.; PRINT_CONST    = .FALSE.; PRINT_DERIVED    = .FALSE.; ERROR_OPEN      = .FALSE.
  LIMITING_FACTOR       = .FALSE.
  HEAD_BOUNDARY         = .FALSE.; PRINT_HYDRO    = .FALSE.; ONE_LAYER        = .FALSE.; ZERO_SLOPE      = .TRUE.
  INTERNAL_FLOW         = .FALSE.; DAM_INFLOW     = .FALSE.; DAM_OUTFLOW      = .FALSE.; HEAD_FLOW       = .FALSE.     !TC 08/03/04
  UPDATE_RATES          = .FALSE.                                                                                       !TC 08/03/04
  WEIR_CALC             =  NIW > 0; GATES    = NGT > 0; PIPES       = NPI > 0
  PUMPS                 =  NPU > 0; SPILLWAY = NSP > 0; TRIBUTARIES = NTR > 0
  WITHDRAWALS           =  NWD > 0
  VOLUME_BALANCE        = VBC         == '      ON'
  PLACE_QIN             = PQC         == '      ON'; EVAPORATION        = EVC    == '      ON'
  ENERGY_BALANCE        = EBC         == '      ON'; RH_EVAP            = RHEVC  == '      ON'
  PRECIPITATION         = PRC         == '      ON'; RESTART_OUT        = RSOC   == '      ON'
  INTERP_TRIBS          = TRIC        == '      ON'; INTERP_DTRIBS      = DTRIC  == '      ON'
  INTERP_HEAD           = HDIC        == '      ON'; INTERP_INFLOW      = QINIC  == '      ON'
  INTERP_OUTFLOW        = STRIC       == '      ON'; INTERP_WITHDRAWAL  = WDIC   == '      ON'
  INTERP_METEOROLOGY    = METIC       == '      ON'; DOWNSTREAM_OUTFLOW = WDOC   == '      ON'
  CELERITY_LIMIT        = CELC        == '      ON'; VISCOSITY_LIMIT    = VISC   == '      ON'
  HYDRO_PLOT            = HPLTC       == '      ON'; PRINT_HYDRO        = HPRWBC == '      ON'
  LIMITING_DLT          = HPRWBC(1,:) == '      ON'; FETCH_CALC         = FETCHC == '      ON'
  SCREEN_OUTPUT         = SCRC        == '      ON'; SNAPSHOT           = SNPC   == '      ON'
  CONTOUR               = CPLC        == '      ON'; VECTOR             = VPLC   == '      ON'
  PROFILE               = PRFC        == '      ON'; SPREADSHEET        = SPRC   == '      ON'
  TIME_SERIES           = TSRC        == '      ON'; READ_RADIATION     = SROC   == '      ON'
  ICE_CALC              = ICEC        == '      ON'; DIST_TRIBS         = DTRC   == '      ON'
  INTERP_EXTINCTION     = EXIC        == '      ON'; READ_EXTINCTION    = EXC    == '      ON'
  NO_INFLOW             = QINC        == '     OFF'; NO_OUTFLOW         = QOUTC  == '     OFF'
  NO_HEAT               = HEATC       == '     OFF'; NO_WIND            = WINDC  == '     OFF'
  SPECIFY_QTR        = TRC    == ' SPECIFY'
  IMPLICIT_VISC         = AZSLC       == '     IMP'; UPWIND             = SLTRC  == '  UPWIND'
  ULTIMATE              = SLTRC       == 'ULTIMATE'; TERM_BY_TERM       = SLHTC  == '    TERM'
  MANNINGS_N            = FRICC       == '    MANN'; PLACE_QTR          = TRC    == ' DENSITY'
  LATERAL_SPILLWAY      = LATSPC      /= '    DOWN'; LATERAL_PUMP       = LATPUC /= '    DOWN'
  LATERAL_GATE          = LATGTC      /= '    DOWN'; LATERAL_PIPE       = LATPIC /= '    DOWN'
  TRAPEZOIDAL           = GRIDC       == '    TRAP'                                                                    !SW 07/16/04
  EPIPHYTON_CALC        = CONSTITUENTS .AND. EPIC        == '      ON'
  MASS_BALANCE          = CONSTITUENTS .AND. MBC         == '      ON'
  SUSP_SOLIDS           = CONSTITUENTS .AND. CAC(NSSS)   == '      ON'
  OXYGEN_DEMAND         = CONSTITUENTS .AND. CAC(NDO)    == '      ON'
! v3.5
  SEDIMENT_CALC         = CONSTITUENTS .AND. SEDCc        == '      ON'
  zooplankton_CALC      = CONSTITUENTS .AND. cac(nzooS) == '      ON'
! v3.5 end
  SEDIMENT_RESUSPENSION = CONSTITUENTS .AND. SEDRC       == '      ON'
  DERIVED_PLOT          = CONSTITUENTS .AND. CDPLTC      == '      ON'
  DERIVED_CALC          = CONSTITUENTS .AND. ANY(CDWBC   == '      ON')
  PH_CALC               = CONSTITUENTS .AND. CDWBC(20,:) == '      ON'
  PRINT_EPIPHYTON       = CONSTITUENTS .AND. EPIPRC      == '      ON' .AND. EPIPHYTON_CALC
  PRINT_SEDIMENT        = CONSTITUENTS .AND. SEDPRC      == '      ON' .AND. SEDIMENT_CALC
  FRESH_WATER           = CONSTITUENTS .AND. WTYPEC      == '   FRESH' .AND. CAC(NTDS) == '      ON'
  SALT_WATER            = CONSTITUENTS .AND. WTYPEC      == '    SALT' .AND. CAC(NTDS) == '      ON'
  CONSTITUENT_PLOT      = CONSTITUENTS .AND. CPLTC       == '      ON' .AND. CAC       == '      ON'
  DETAILED_ICE          = ICE_CALC     .AND. SLICEC      == '  DETAIL'
  LEAP_YEAR             = MOD(YEAR,4) == 0
  ICE_COMPUTATION       = ANY(ICE_CALC)
  END_RUN               = JDAY > TMEND
  UPDATE_KINETICS       = CONSTITUENTS
  IF (WEIR_CALC) THEN
    DO JWR=1,NIW
      DO K=2,KMX-1
        IF ((K >= KTWR(JWR) .AND. K <= KBWR(JWR))) INTERNAL_WEIR(K,IWR(JWR)) = .TRUE.
      END DO
    END DO
  END IF
  IF (RESTART_IN) THEN
    IF (JDAY > 300.0 .OR. JDAY < 40.0)     WINTER = .TRUE.
  ELSE
    IF (TMSTRT > 300.0 .OR. TMSTRT < 40.0) WINTER = .TRUE.
  END IF
  WHERE (READ_EXTINCTION)      
    EXOM = 0.0
    EXSS = 0.0
  ENDWHERE
  IF (CONSTITUENTS) THEN
    SUSP_SOLIDS   = .FALSE.
    FLUX          =  FLXC   == '      ON'
    PRINT_CONST   =  CPRWBC == '      ON'
    PRINT_DERIVED =  CDWBC  == '      ON'
    IF (ANY(CAC(NSSS:NSSE)  == '      ON')) SUSP_SOLIDS  = .TRUE.
    IF (ANY(CAC(NSSS:NCT)   == '      ON')) UPDATE_RATES = .TRUE.      
    DO JA=1,NAL
      LIMITING_FACTOR(JA) = CONSTITUENTS .AND. CAC(NAS-1+JA) == '      ON' .AND. LIMC == '      ON'      
    END DO
  END IF
  JBDAM = 0
  CDHS  = DHS
  DO JB=1,NBR
    UP_FLOW(JB)     = UHS(JB) == 0
    DN_FLOW(JB)     = DHS(JB) == 0
    UP_HEAD(JB)     = UHS(JB) /= 0
    UH_INTERNAL(JB) = UHS(JB) >  0
    IF (UP_HEAD(JB)) THEN
      DO JJB=1,NBR
        IF (ABS(UHS(JB)) >= US(JJB) .AND. ABS(UHS(JB)) <= DS(JJB)) THEN
          IF (ABS(UHS(JB)) == DS(JJB)) THEN
            IF (DHS(JJB) == US(JB)) THEN
              UP_FLOW(JB)       = .TRUE.
              HEAD_FLOW(JB)     = .TRUE.
              INTERNAL_FLOW(JB) = .TRUE.
              UP_HEAD(JB)       = .FALSE.
              UH_INTERNAL(JB)   = .FALSE.
            END IF
            IF (UHS(JB) < 0) THEN
              DO JJJB=1,NBR
                IF (ABS(UHS(JB)) == DS(JJJB)) EXIT                                                                     ! CB 1/2/05
              END DO
              UP_FLOW(JB)       = .TRUE.
              DAM_INFLOW(JB)    = .TRUE.                                                                               !TC 08/03/04
              DAM_OUTFLOW(JJJB) = .TRUE.                                                                               !TC 08/03/04
              INTERNAL_FLOW(JB) = .TRUE.
              UP_HEAD(JB)       = .FALSE.
              UHS(JB)           =  ABS(UHS(JB))
              JBDAM(JJJB)       =  JB
            END IF
          END IF
          EXIT
        END IF
      END DO
    END IF
    DH_INTERNAL(JB) = DHS(JB)  >   0; DN_HEAD(JB)     = DHS(JB)  /=  0; UH_EXTERNAL(JB) = UHS(JB)  == -1
    DH_EXTERNAL(JB) = DHS(JB)  == -1; UQ_EXTERNAL(JB) = UHS(JB)  ==  0; DQ_EXTERNAL(JB) = DHS(JB)  ==  0
    DQ_INTERNAL(JB) = DQB(JB)  >   0; UQ_INTERNAL(JB) = UQB(JB)  >   0 .AND. .NOT. DAM_INFLOW(JB)                      !TC 08/03/04
  END DO
  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      IF (UH_EXTERNAL(JB) .OR. DH_EXTERNAL(JB)) HEAD_BOUNDARY(JW) = .TRUE.
      IF (SLOPE(JB) /= 0.0)                     ZERO_SLOPE(JW)    = .FALSE.
    END DO
  END DO
  WHERE (CAC == '     OFF') CPLTC = '     OFF'

! Kinetic flux variables

  KFNAME(1)  = 'TISS settling in - source, kg/day            '; KFNAME(2)  = 'TISS settling out - sink, kg/day             '
  KFNAME(3)  = 'PO4 algal respiration - source, kg/day       '; KFNAME(4)  = 'PO4 algal growth - sink, kg/day              '
  KFNAME(5)  = 'PO4 algal net- source/sink, kg/day           '; KFNAME(6)  = 'PO4 epiphyton respiration - source, kg/day   '
  KFNAME(7)  = 'PO4 epiphyton growth - sink, kg/day          '; KFNAME(8)  = 'PO4 epiphyton net- source/sink, kg/day       '
  KFNAME(9)  = 'PO4 POM decay - source, kg/day               '; KFNAME(10) = 'PO4 DOM decay - source, kg/day               '
  KFNAME(11) = 'PO4 OM decay - source, kg/day                '; KFNAME(12) = 'PO4 sediment decay - source, kg/day          '
  KFNAME(13) = 'PO4 SOD release - source, kg/day             '; KFNAME(14) = 'PO4 net settling  - source/sink, kg/day      '
  KFNAME(15) = 'NH4 nitrification - sink, kg/day             '; KFNAME(16) = 'NH4 algal respiration - source, kg/day       '
  KFNAME(17) = 'NH4 algal growth - sink, kg/day              '; KFNAME(18) = 'NH4 algal net - source/sink, kg/day          '
  KFNAME(19) = 'NH4 epiphyton respiration - source, kg/day   '; KFNAME(20) = 'NH4 epiphyton growth - sink, kg/day          '
  KFNAME(21) = 'NH4 epiphyton net - source/sink, kg/day      '; KFNAME(22) = 'NH4 POM decay - source, kg/day               '
  KFNAME(23) = 'NH4 DOM decay  - source, kg/day              '; KFNAME(24) = 'NH4 OM decay - source, kg/day                '
  KFNAME(25) = 'NH4 sediment decay - source, kg/day          '; KFNAME(26) = 'NH4 SOD release - source, kg/day             '
  KFNAME(27) = 'NO3 denitrification - sink, kg/day           '; KFNAME(28) = 'NO3 algal growth - sink, kg/day              '
  KFNAME(29) = 'NO3 epiphyton growth - sink, kg/day          '; KFNAME(30) = 'NO3 sediment uptake - sink, kg/day           '
  KFNAME(31) = 'DSi algal growth - sink, kg/day              '; KFNAME(32) = 'DSi epiphyton growth - sink, kg/day          '
  KFNAME(33) = 'DSi PBSi decay - source, kg/day              '; KFNAME(34) = 'DSi sediment decay - source, kg/day          '
  KFNAME(35) = 'DSi SOD release  - source, kg/day            '; KFNAME(36) = 'DSi net settling - source/sink, kg/day       '
  KFNAME(37) = 'PBSi algal mortality  - source, kg/day       '; KFNAME(38) = 'PBSi net settling - source/sink, kg/day      '
  KFNAME(39) = 'PBSi decay - sink, kg/day                    '; KFNAME(40) = 'Fe net settling - source/sink, kg/day        '
  KFNAME(41) = 'Fe sediment release - source, kg/day         '; KFNAME(42) = 'LDOM decay - sink, kg/day                    '
  KFNAME(43) = 'LDOM decay to RDOM - sink, kg/day            '; KFNAME(44) = 'RDOM decay - sink, kg/day                    '
  KFNAME(45) = 'LDOM algal mortality - source, kg/day        '; KFNAME(46) = 'LDOM epiphyton mortality - source, kg/day    '
  KFNAME(47) = 'LPOM decay - sink, kg/day                    '; KFNAME(48) = 'LPOM decay to RPOM - sink, kg/day            '
  KFNAME(49) = 'RPOM decay - sink, kg/day                    '; KFNAME(50) = 'LPOM algal production - source, kg/day       '
  KFNAME(51) = 'LPOM epiphyton production - source, kg/day   '; KFNAME(52) = 'LPOM net settling - source/sink, kg/day      '
  KFNAME(53) = 'RPOM net settling - source/sink, kg/day      '; KFNAME(54) = 'CBOD decay - sink, kg/day                    '
  KFNAME(55) = 'DO algal production  - source, kg/day        '; KFNAME(56) = 'DO epiphyton production  - source, kg/day    '
  KFNAME(57) = 'DO algal respiration - sink, kg/day          '; KFNAME(58) = 'DO epiphyton respiration - sink, kg/day      '
  KFNAME(59) = 'DO POM decay - sink, kg/day                  '; KFNAME(60) = 'DO DOM decay - sink, kg/day                  '
  KFNAME(61) = 'DO OM decay - sink, kg/day                   '; KFNAME(62) = 'DO nitrification - sink, kg/day              '
  KFNAME(63) = 'DO CBOD uptake - sink, kg/day                '; KFNAME(64) = 'DO rearation - source, kg/day                '
  KFNAME(65) = 'DO sediment uptake - sink, kg/day            '; KFNAME(66) = 'DO SOD uptake - sink, kg/day                 '
  KFNAME(67) = 'TIC algal uptake - sink, kg/day              '; KFNAME(68) = 'TIC epiphyton uptake - sink, kg/day          '
  KFNAME(69) = 'Sediment decay - sink, kg/day                '; KFNAME(70) = 'Sediment algal settling - sink, kg/day       '
  KFNAME(71) = 'Sediment LPOM settling - source,kg/day       '; KFNAME(72) = 'Sediment net settling - source/sink, kg/day  '
  KFNAME(73) = 'SOD decay - sink, kg/day                     '

! v3.5 start
  KFNAME(74) = 'LDOM P algal mortality - source, kg/day      '; KFNAME(75) = 'LDOM P epiphyton mortality - source, kg/day  '
  KFNAME(76) = 'LPOM P algal production- source, kg/day      '; KFNAME(77) = 'LPOM P net settling - source/sink, kg/day    '
  KFNAME(78) = 'RPOM P net settling - source/sink, kg/day    '
  KFNAME(79) = 'LDOM P algal mortality - source, kg/day      '; KFNAME(80) = 'LDOM P epiphyton mortality - source, kg/day  '
  KFNAME(81) = 'LPOM P algal production- source, kg/day      '; KFNAME(82) = 'LPOM P net settling - source/sink, kg/day    '
  KFNAME(83) = 'RPOM P net settling - source/sink, kg/day    '
  KFNAME(84) = 'Sediment P decay - sink, kg/day              '; KFNAME(85) = 'Sediment algal P settling - source, kg/day   '
  KFNAME(86) = 'Sediment P LPOM settling - source,kg/day     '; KFNAME(87) = 'Sediment net P settling - source/sink, kg/day'
  KFNAME(88) = 'Sediment epiphyton P settling - source,kg/day'
  KFNAME(89) = 'Sediment N decay - sink, kg/day              '; KFNAME(90) = 'Sediment algal N settling - source, kg/day   '
  KFNAME(91) = 'Sediment N LPOM settling - source,kg/day     '; KFNAME(92) = 'Sediment net N settling - source/sink, kg/day'
  KFNAME(93) = 'Sediment epiphyton N settling - source,kg/day'
  KFNAME(94) = 'Sediment C decay - sink, kg/day              '; KFNAME(95) = 'Sediment algal C settling - source, kg/day   '
  KFNAME(96) = 'Sediment C LPOM settling - source,kg/day     '; KFNAME(97) = 'Sediment net C settling - source/sink, kg/day'
  KFNAME(98) = 'Sediment epiphyton C settling - source,kg/day'
  KFNAME(99) = 'Sediment N denitrification - source, kg/day  ';
  KFNAME(100) = 'PO4 macrophyte resp - source, kg/day         '
  KFNAME(101) = 'PO4 macrophyte growth - sink, kg/day         '
  KFNAME(102) = 'NH4 macrophyte resp - source, kg/day         '
  KFNAME(103) = 'NH4 macrophyte growth - sink, kg/day         '
  KFNAME(104) = 'LDOM macrophyte mort  - source, kg/day       '
  KFNAME(105) = 'LPOM macrophyte mort  - source, kg/day       '
  KFNAME(106) = 'RPOM macrophyte mort  - source, kg/day       '
  KFNAME(107) = 'DO  macrophyte production  - source, kg/day  '
  KFNAME(108) = 'DO  macrophyte respiration - sink, kg/day    '
  KFNAME(109) = 'TIC macrophyte growth/resp  - S/S, kg/day    '
  KFNAME(110) = 'CBOD settling - sink, kg/day                 '
  KFNAME(111) = 'Sediment CBOD settling - source, kg/day      '
  KFNAME(112) = 'Sediment CBOD P settling - source, kg/day    '
  KFNAME(113) = 'Sediment CBOD N settling - source, kg/day    '
  KFNAME(114) = 'Sediment CBOD C settling - source, kg/day    '
! v3.5 end

! Convert rates from per-day to per-second
 
  IF (CONSTITUENTS) THEN
    AE     = AE    /DAY; AM     = AM    /DAY; AR     = AR    /DAY; AG    = AG    /DAY; AS     = AS    /DAY
    EE     = EE    /DAY; EM     = EM    /DAY; ER     = ER    /DAY; EG    = EG    /DAY; EB     = EB    /DAY
    CGS    = CGS   /DAY; CG0DK  = CG0DK /DAY; CG1DK  = CG1DK /DAY; SSS   = SSS   /DAY; FES    = FES   /DAY
    PSIS   = PSIS  /DAY; POMS   = POMS  /DAY; SDK    = SDK   /DAY; NH4DK = NH4DK /DAY; NO3DK  = NO3DK /DAY
    NO3S   = NO3S  /DAY; PSIDK  = PSIDK /DAY; LRDDK  = LRDDK /DAY; LRPDK = LRPDK /DAY; LDOMDK = LDOMDK/DAY
    LPOMDK = LPOMDK/DAY; RDOMDK = RDOMDK/DAY; RPOMDK = RPOMDK/DAY; KBOD  = KBOD  /DAY; seds   = seds/day   !v3.5
    sedb=sedb/day   !cb 11/27/06
    cbods  = cbods/day   !cb 7/23/07
    DO JW=1,NWB
      SOD(US(BS(JW))-1:DS(BE(JW))+1) = (SOD(US(BS(JW))-1:DS(BE(JW))+1)/DAY)*FSOD(JW)
    END DO
    DO J=1,NEP
      EBR(:,:,J) = EB(J)
    END DO
!  v3.5 start    
    mg=mg/day
    mr=mr/day
    mm=mm/day    
    zg=zg/day
    zr=zr/day
    zm=zm/day
! v3.5 end

  END IF

! Convert slope to angle alpha in radians

  ALPHA = ATAN(SLOPE)
  SINA  = SIN(ALPHA)
  COSA  = COS(ALPHA)

! Time and printout control variables

  IF (.NOT. RESTART_IN) THEN
    JDAY   = TMSTRT
    ELTM   = TMSTRT*DAY
    DLT    = DLTMAX(1)
    DLTS   = DLT
    MINDLT = DLT
    NIT    = 0
    NV     = 0
    DLTDP  = 1; RSODP = 1; TSRDP = 1; SNPDP = 1; VPLDP = 1; PRFDP = 1
    SPRDP  = 1; CPLDP = 1; SCRDP = 1; FLXDP = 1; WDODP = 1
    DO JW=1,NWB
      DO J=1,NOD
        IF (TMSTRT > SNPD(J,JW)) SNPD(J,JW) = TMSTRT; IF (TMSTRT > PRFD(J,JW)) PRFD(J,JW) = TMSTRT
        IF (TMSTRT > SPRD(J,JW)) SPRD(J,JW) = TMSTRT; IF (TMSTRT > CPLD(J,JW)) CPLD(J,JW) = TMSTRT
        IF (TMSTRT > VPLD(J,JW)) VPLD(J,JW) = TMSTRT; IF (TMSTRT > SCRD(J,JW)) SCRD(J,JW) = TMSTRT
        IF (TMSTRT > FLXD(J,JW)) FLXD(J,JW) = TMSTRT
      END DO
      NXTMSN(JW) = SNPD(SNPDP(JW),JW); NXTMPR(JW) = PRFD(PRFDP(JW),JW); NXTMSP(JW) = SPRD(SPRDP(JW),JW)
      NXTMCP(JW) = CPLD(CPLDP(JW),JW); NXTMVP(JW) = VPLD(VPLDP(JW),JW); NXTMSC(JW) = SCRD(SCRDP(JW),JW)
      NXTMFL(JW) = FLXD(FLXDP(JW),JW)
    END DO
    DO J=1,NOD
      IF (TMSTRT > TSRD(J)) TSRD(J) = TMSTRT; IF (TMSTRT > WDOD(J)) WDOD(J) = TMSTRT
      IF (TMSTRT > RSOD(J)) RSOD(J) = TMSTRT; IF (TMSTRT > DLTD(J)) DLTD(J) = TMSTRT
    END DO
    NXTMTS = TSRD(TSRDP); NXTMWD = WDOD(WDODP); NXTMRS = RSOD(RSODP)
  END IF
  TSRD(NTSR+1:NOD) = TMEND+1.0; WDOD(NWDO+1:NOD) = TMEND+1.0; RSOD(NRSO+1:NOD) = TMEND+1.0; DLTD(NDLT+1:NOD) = TMEND+1.0
  DO JW=1,NWB
    SNPD(NSNP(JW)+1:NOD,JW) = TMEND+1.0; PRFD(NPRF(JW)+1:NOD,JW) = TMEND+1.0; SPRD(NSPR(JW)+1:NOD,JW) = TMEND+1.0
    VPLD(NVPL(JW)+1:NOD,JW) = TMEND+1.0; CPLD(NCPL(JW)+1:NOD,JW) = TMEND+1.0; SCRD(NSCR(JW)+1:NOD,JW) = TMEND+1.0
    FLXD(NFLX(JW)+1:NOD,JW) = TMEND+1.0      
  END DO
  JDAYG  = JDAY
  JDAYNX = JDAYG+1
  NXTVD  = JDAY
  CURMAX = DLTMAX(DLTDP)/DLTF(DLTDP)
 
! Hydraulic structures

  IF (SPILLWAY) THEN
    DO JS=1,NSP
      IF (LATERAL_SPILLWAY(JS)) THEN
        IF (IDSP(JS) /= 0) THEN
          TRIBUTARIES = .TRUE.
          WITHDRAWALS = .TRUE.
        ELSE
          WITHDRAWALS = .TRUE.
        END IF
      END IF     
      DO JB=1,NBR
        IF (IUSP(JS) >= US(JB) .AND. IUSP(JS) <= DS(JB)) EXIT
      END DO
      JBUSP(JS) = JB
      IF (IUSP(JS) == DS(JBUSP(JS)) .AND. .NOT. LATERAL_SPILLWAY(JS)) NST = NST+1
      DO JW=1,NWB
        IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
      END DO
      JWUSP(JS) = JW
      IF (IDSP(JS) > 0) THEN
        DO JB=1,NBR
          IF (IDSP(JS) >= US(JB) .AND. IDSP(JS) <= DS(JB)) EXIT
        END DO
        JBDSP(JS) = JB
        DO JW=1,NWB
          IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
        END DO
        JWDSP(JS) = JW
      END IF
    END DO
  END IF
  IF (PIPES) THEN
    DO JP=1,NPI
      IF (LATERAL_PIPE(JP)) THEN
        IF (IDPI(JP) /= 0) THEN
          TRIBUTARIES = .TRUE.
          WITHDRAWALS = .TRUE.
        ELSE
          WITHDRAWALS = .TRUE.
        END IF
      END IF
      DO JB=1,NBR
        IF (IUPI(JP) >= US(JB) .AND. IUPI(JP) <= DS(JB)) EXIT
      END DO
      JBUPI(JP) = JB
      IF (IUPI(JP) == DS(JBUPI(JP)) .AND. .NOT. LATERAL_PIPE(JP)) NST = NST+1
      DO JW=1,NWB
        IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
      END DO
      JWUPI(JP) = JW
      IF (IDPI(JP) > 0) THEN
        DO JB=1,NBR
          IF (IDPI(JP) >= US(JB) .AND. IDPI(JP) <= DS(JB)) EXIT
        END DO
        JBDPI(JP) = JB
        DO JW=1,NWB
          IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
        END DO
        JWDPI(JP) = JW
      END IF
    END DO
  END IF
  IF (GATES) THEN
    DO JG=1,NGT
      IF (LATERAL_GATE(JG)) THEN
        IF (IDGT(JG) /= 0) THEN
          TRIBUTARIES = .TRUE.
          WITHDRAWALS = .TRUE.
        ELSE
          WITHDRAWALS = .TRUE.
        END IF
      END IF     
      DO JB=1,NBR
        IF (IUGT(JG) >= US(JB) .AND. IUGT(JG) <= DS(JB)) EXIT
      END DO
      JBUGT(JG) = JB
      IF (IUGT(JG) == DS(JBUGT(JG)) .AND. .NOT. LATERAL_GATE(JG)) NST = NST+1
      DO JW=1,NWB
        IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
      END DO
      JWUGT(JG) = JW
      IF (IDGT(JG) > 0) THEN
        DO JB=1,NBR
          IF (IDGT(JG) >= US(JB) .AND. IDGT(JG) <= DS(JB)) EXIT
        END DO
        JBDGT(JG) = JB
        DO JW=1,NWB
          IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
        END DO
        JWDGT(JG) = JW
      END IF
    END DO
  END IF
  IF (PUMPS) THEN
    DO JP=1,NPU
      IF (LATERAL_PUMP(JP)) THEN
        IF (IDPU(JP) /= 0) THEN
          TRIBUTARIES = .TRUE.
          WITHDRAWALS = .TRUE.
        ELSE
          WITHDRAWALS = .TRUE.
        END IF
      END IF     
      DO JB=1,NBR
        IF (IUPU(JP) >= US(JB) .AND. IUPU(JP) <= DS(JB)) EXIT
      END DO
      JBUPU(JP) = JB
      IF (IUPU(JP) ==  DS(JBUPU(JP)) .AND. .NOT. LATERAL_PUMP(JP)) NST = NST+1
      DO JW=1,NWB
        IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
      END DO
      JWUPU(JP) = JW
      IF (IDPU(JP) > 0) THEN
        DO JB=1,NBR
          IF (IDPU(JP) >= US(JB) .AND. IDPU(JP) <= DS(JB)) EXIT
        END DO
        JBDPU(JP) = JB
        DO JW=1,NWB
          IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
        END DO
        JWDPU(JP) = JW
      END IF
    END DO
  END IF

  ALLOCATE (ESTR(NST,NBR),WSTR(NST,NBR),QSTR(NST,NBR),KTSW(NST,NBR),KBSW(NST,NBR),SINKC(NST,NBR),POINT_SINK(NST,NBR),QNEW(KMX))
  QSTR = 0.0
  DO JB=1,NBR
    ESTR(1:NSTR(JB),JB)  = ESTRT(1:NSTR(JB),JB)
    KTSW(1:NSTR(JB),JB)  = KTSWT(1:NSTR(JB),JB)
    KBSW(1:NSTR(JB),JB)  = KBSWT(1:NSTR(JB),JB)
    WSTR(1:NSTR(JB),JB)  = WSTRT(1:NSTR(JB),JB)
    SINKC(1:NSTR(JB),JB) = SINKCT(1:NSTR(JB),JB)
  END DO
  POINT_SINK = SINKC == '   POINT'
  DEALLOCATE (ESTRT,KBSWT,KTSWT,WSTRT,SINKCT)

! Active constituents, derived constituents, and fluxes

  IF (CONSTITUENTS) THEN
    DO JC=1,NCT
      IF (CAC(JC) == '      ON') THEN
        NAC     = NAC+1
        CN(NAC) = JC
      END IF
      DO JB=1,NBR
        IF (CINBRC(JC,JB) == '      ON') THEN
          NACIN(JB)       = NACIN(JB)+1
          INCN(NACIN(JB),JB) = JC
        END IF
        IF (CDTBRC(JC,JB) == '      ON') THEN
          NACDT(JB)       = NACDT(JB)+1
          DTCN(NACDT(JB),JB) = JC
        END IF
        IF (CPRBRC(JC,JB) == '      ON') THEN
          NACPR(JB)       = NACPR(JB)+1
          PRCN(NACPR(JB),JB) = JC
        END IF
      END DO
      DO JT=1,NTR
        IF (CTRTRC(JC,JT) == '      ON') THEN
          NACTR(JT)          = NACTR(JT)+1
          TRCN(NACTR(JT),JT) = JC
        END IF
      END DO
    END DO
    DO JW=1,NWB
      DO JD=1,NDC
        IF (CDWBC(JD,JW) == '      ON') THEN
          NACD(JW)         = NACD(JW)+1
          CDN(NACD(JW),JW) = JD
        END IF
      END DO
      DO JF=1,NFL
        IF (KFWBC(JF,JW) == '      ON') THEN
          NAF(JW)          = NAF(JW)+1
          KFCN(NAF(JW),JW) = JF
        END IF
      END DO
    END DO
  END IF

! Starting time

  DEG = CHAR(248)//'C'
  ESC = CHAR(027)
  CALL DATE_AND_TIME (CDATE,CTIME)
  DO JW=1,NWB
    TITLE(11) = 'Model run at '//CTIME(1:2)//':'//CTIME(3:4)//':'//CTIME(5:6)//' on '//CDATE(5:6)//'/'//CDATE(7:8)//'/'//CDATE(3:4)
    IF (RESTART_IN) TITLE(11) = 'Model restarted at '//CTIME(1:2)//':'//CTIME(3:4)//':'//CTIME(5:6)//' on '//CDATE(5:6)//'/'//     &
                                 CDATE(7:8)//'/'//CDATE(3:4)
  END DO

!***********************************************************************************************************************************
!**                                                    Task 1.1.3: Geometry                                                       **
!***********************************************************************************************************************************

! Layer elevations

  NPOINT = 0
  DO JW=1,NWB
    IF (ZERO_SLOPE(JW)) THEN
      DO I=US(BS(JW))-1,DS(BE(JW))+1
        EL(KMX,I) = ELBOT(JW)
        DO K=KMX-1,1,-1
          EL(K,I) = EL(K+1,I)+H(K,JW)
        END DO
      END DO
    ELSE
      EL(KMX,DS(JBDN(JW))+1) = ELBOT(JW)
      JB                     = JBDN(JW)
      NPOINT(JB)             = 1
      NNBP                   = 1
      NCBP                   = 0
      NINTERNAL              = 0
      NUP                    = 0
      DO WHILE (NNBP <= (BE(JW)-BS(JW)+1))
        NCBP = NCBP+1
        IF (NINTERNAL == 0) THEN
          IF (NUP == 0) THEN
            DO I=DS(JB),US(JB),-1
              IF (I /= DS(JB)) THEN
                EL(KMX,I) = EL(KMX,I+1)+SINA(JB)*(DLX(I)+DLX(I+1))*0.5
              ELSE
                EL(KMX,I) = EL(KMX,I+1)
              END IF
              DO K=KMX-1,1,-1
                EL(K,I) = EL(K+1,I)+H(K,JW)*COSA(JB)
              END DO
            END DO
          ELSE
            DO I=US(JB),DS(JB)
              IF (I /= US(JB)) THEN
                EL(KMX,I) = EL(KMX,I-1)-SINA(JB)*(DLX(I)+DLX(I-1))*0.5
              ELSE
                EL(KMX,I) = EL(KMX,I-1)
              END IF
              DO K=KMX-1,1,-1
                EL(K,I) = EL(K+1,I)+H(K,JW)*COSA(JB)
              END DO
            END DO
            NUP = 0
          END IF
          DO K=KMX,1,-1
            IF (UP_HEAD(JB)) THEN
              EL(K,US(JB)-1) = EL(K,US(JB))+SINA(JB)*DLX(US(JB))
            ELSE
              EL(K,US(JB)-1) = EL(K,US(JB))
            END IF
            IF (DN_HEAD(JB)) THEN
              EL(K,DS(JB)+1) = EL(K,DS(JB))-SINA(JB)*DLX(DS(JB))
            ELSE
              EL(K,DS(JB)+1) = EL(K,DS(JB))
            END IF
          END DO
        ELSE
          DO K=KMX-1,1,-1
            EL(K,UHS(JJB)) = EL(K+1,UHS(JJB))+H(K,JW)*COSA(JB)
          END DO
          DO I=UHS(JJB)+1,DS(JB)
            EL(KMX,I) = EL(KMX,I-1)-SINA(JB)*(DLX(I)+DLX(I-1))*0.5
            DO K=KMX-1,1,-1
              EL(K,I) = EL(K+1,I)+H(K,JW)*COSA(JB)
            END DO
          END DO
          DO I=UHS(JJB)-1,US(JB),-1
            EL(KMX,I) = EL(KMX,I+1)+SINA(JB)*(DLX(I)+DLX(I+1))*0.5
            DO K=KMX-1,1,-1
              EL(K,I) = EL(K+1,I)+H(K,JW)*COSA(JB)
            END DO
          END DO
          NINTERNAL = 0
        END IF
        IF (NNBP == (BE(JW)-BS(JW)+1)) EXIT

!****** Find next branch connected to furthest downstream branch

        DO JB=BS(JW),BE(JW)
          IF (NPOINT(JB) /= 1) THEN
            DO JJB=BS(JW),BE(JW)
              IF (DHS(JB) >= US(JJB) .AND. DHS(JB) <=DS (JJB) .AND. NPOINT(JJB) == 1) THEN
                NPOINT(JB)       = 1
                EL(KMX,DS(JB)+1) = EL(KMX,DHS(JB))+SINA(JB)*(DLX(DS(JB))+DLX(DHS(JB)))*0.5
                NNBP             = NNBP+1; EXIT
              END IF
              IF (UHS(JJB) == DS(JB) .AND. NPOINT(JJB) == 1) THEN
                NPOINT(JB)       = 1      
                EL(KMX,DS(JB)+1) = EL(KMX,US(JJB))+(SINA(JJB)*DLX(US(JJB))+SINA(JB)*DLX(DS(JB)))*0.5
                NNBP             = NNBP+1; EXIT
              END IF
              IF (UHS(JJB) >= US(JB) .AND. UHS(JJB) <= DS(JB) .AND. NPOINT(JJB)==1) THEN
                NPOINT(JB)       = 1
                EL(KMX,UHS(JJB)) = EL(KMX,US(JJB))+SINA(JJB)*DLX(US(JJB))*0.5
                NNBP             = NNBP+1
                NINTERNAL        = 1; EXIT
              END IF
              IF (UHS(JB) >= US(JJB) .AND. UHS(JB) <= DS(JJB) .AND. NPOINT(JJB) == 1) THEN
                NPOINT(JB)       = 1
                EL(KMX,US(JB)-1) = EL(KMX,UHS(JB))-SINA(JB)*DLX(US(JB))*0.5
                NNBP             = NNBP+1
                NUP              = 1; EXIT
              END IF
            END DO
            IF (NPOINT(JB)==1) EXIT
          END IF
        END DO
      END DO
    END IF
  END DO

! Minimum/maximum layer heights

  DO JW=1,NWB
    DO K=KMX-1,1,-1
      HMIN = MIN(H(K,JW),HMIN)
      HMAX = MAX(H(K,JW),HMAX)
    END DO
  END DO
  HMAX2 = HMAX**2

! Water surface and bottom layers

  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      DO I=US(JB)-1,DS(JB)+1
        IF (.NOT. RESTART_IN) THEN
          KTI(I) = 2
          DO WHILE (EL(KTI(I),I) > ELWS(I))
            KTI(I) = KTI(I)+1
          END DO
          Z(I)     = (EL(KTI(I),I)-ELWS(I))/COSA(JB)
          ZMIN(JW) =  MAX(ZMIN(JW),Z(I))
          KTMAX    =  MAX(2,KTI(I))
          KTWB(JW) =  MAX(KTMAX,KTWB(JW))
          KTI(I)   =  MAX(KTI(I)-1,2)
          IF (Z(I) > ZMIN(JW)) IZMIN(JW) = I
        END IF
        K = 2
        DO WHILE (B(K,I) > 0.0)
          KB(I) = K
          K     = K+1
        END DO
        KBMAX(JW) = MAX(KBMAX(JW),KB(I))
      END DO
      KB(US(JB)-1) = KB(US(JB))
      KB(DS(JB)+1) = KB(DS(JB))
    END DO
    KBI = KB

!** Correct for water surface going over several layers

    IF (.NOT. RESTART_IN) THEN
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        DO I=US(JB)-1,DS(JB)+1
          H2(KT,I) = H(KT,JW)-Z(I)
          K        = KTI(I)+1
          DO WHILE (KT > K)
            Z(I)     = Z(I)-H(K,JW)
            H2(KT,I) = H(KT,JW)-Z(I)
            K        = K+1
          END DO
        END DO
      END DO
    END IF
    ELKT(JW) = EL(KTWB(JW),DS(BE(JW)))-Z(DS(BE(JW)))*COSA(BE(JW))
  END DO
  DO JW=1,NWB
    KT = KTWB(JW)
    DO JB=BS(JW),BE(JW)
      IU = US(JB)
      ID = DS(JB)

!**** Boundary bottom layers

      IF (UH_EXTERNAL(JB)) KB(IU-1) = KB(IU)                                        ! CB 6/12/07
      IF (DH_EXTERNAL(JB)) KB(ID+1)  = KB(ID)

!**** Branch numbers corresponding to tributaries, withdrawals, and head

      IF (TRIBUTARIES) THEN
        DO JT=1,NTR
          IF (ITR(JT) >= US(JB) .AND. ITR(JT) <= DS(JB)) JBTR(JT) = JB
        END DO
      END IF
      IF (WITHDRAWALS) THEN
        DO JWD=1,NWD
          IF (IWD(JWD) >= US(JB) .AND. IWD(JWD) <= DS(JB)) JBWD(JWD) = JB
        END DO
      END IF
      IF (UH_INTERNAL(JB)) THEN
        DO JJB=1,NBR
          JBUH(JB) = JJB
          IF (UHS(JB) >= US(JJB) .AND. UHS(JB) <= DS(JJB)) EXIT
        END DO
        DO JJW=1,NWB
          JWUH(JB) = JJW
          IF (JBUH(JB) >= BS(JJW) .AND. JBUH(JB) <= BE(JJW)) EXIT
        END DO
      END IF
      IF (INTERNAL_FLOW(JB)) THEN
        DO JJB=1,NBR
          JBUH(JB) = JJB
          IF (UHS(JB) >= US(JJB) .AND. UHS(JB) <= DS(JJB)) EXIT
        END DO
        DO JJW=1,NWB
          JWUH(JB) = JJW
          IF (JBUH(JB) >= BS(JJW) .AND. JBUH(JB) <= BE(JJW)) EXIT
        END DO
      END IF
      IF (DH_INTERNAL(JB)) THEN
        DO JJB=1,NBR
          JBDH(JB) = JJB
          IF (DHS(JB) >= US(JJB) .AND. DHS(JB) <= DS(JJB)) EXIT
        END DO
        DO JJW=1,NWB
          JWDH(JB) = JJW
          IF (JBDH(JB) >= BS(JJW) .AND. JBDH(JB) <= BE(JJW)) EXIT
        END DO
      END IF

!**** Bottom boundary cells

      IF (UH_INTERNAL(JB)) THEN
        IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
          KB(IU-1) = MIN(KB(UHS(JB)),KB(IU))                                 ! CB 6/12/07
        ELSE
          IF (EL(KB(IU),IU) >= EL(KB(UHS(JB)),UHS(JB))) THEN                 ! CB 6/12/07
            KB(IU-1) = KB(IU)                                                ! CB 6/12/07
          ELSE
            DO K=KT,KB(IU)                                                   ! CB 6/12/07
              IF (EL(KB(UHS(JB)),UHS(JB)) >= EL(K,IU)) THEN                  ! CB 6/12/07
                KB(IU-1) = K; EXIT                                           ! CB 6/12/07
              END IF
            END DO
          END IF
        END IF
      END IF
      IF (DH_INTERNAL(JB)) THEN
        IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
          KB(ID+1) = MIN(KB(DHS(JB)),KB(ID))
        ELSE
          IF (EL(KB(ID),ID) >= EL(KB(DHS(JB)),DHS(JB))) THEN
            KB(ID+1) = KB(ID)
          ELSE
            DO K=KT,KB(ID)
              IF (EL(KB(DHS(JB)),DHS(JB)) >= EL(K,ID)) THEN
                KB(ID+1) = K; EXIT
              END IF
            END DO
          END IF
        END IF
      END IF

!**** Boundary segment lengths

      DLX(IU-1) = DLX(IU)
      DLX(ID+1) = DLX(ID)

!**** Minimum bottom layers and average segment lengths

      DO I=IU-1,ID
        KBMIN(I) =  MIN(KB(I),KB(I+1))
        DLXR(I)  = (DLX(I)+DLX(I+1))*0.5
      END DO
      KBMIN(ID+1) = KBMIN(ID)
      DLXR(ID+1)  = DLX(ID)

!**** Minimum/maximum segment lengths

      DO I=IU,ID
        DLXMIN = MIN(DLXMIN,DLX(I))
        DLXMAX = MAX(DLXMAX,DLX(I))
      END DO
    END DO
  END DO

! Boundary widths

  DO JW=1,NWB
    KT = KTWB(JW)
    DO JB=BS(JW),BE(JW)
      IU = US(JB)
      ID = DS(JB)
      DO I=IU-1,ID+1
        B(1,I) = B(2,I)
        DO K=KB(I)+1,KMX
          B(K,I) = B(KB(I),I)
        END DO
      END DO
    END DO
  END DO

  DO JW=1,NWB
    KT = KTWB(JW)
    DO JB=BS(JW),BE(JW)
      IU    = US(JB)
      ID    = DS(JB)
      IEXIT = 0
      DO K=1,KMX-1
        B(K,IU-1) = B(K,IU)
        IF (UH_INTERNAL(JB) .OR. HEAD_FLOW(JB)) THEN   
          IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
            B(K,IU-1) = B(K,UHS(JB))
          ELSE
            ELR = EL(K,IU)+SINA(JB)*DLX(IU)*0.5
            ELL = EL(2,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
            IF (ELR >= ELL) THEN
              B(K,IU-1) = B(2,UHS(JB))
            ELSE
              DO KUP=2,KMX-1
                ELL1 = EL(KUP,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                ELL2 = EL(KUP+1,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                IF (ELL1 > ELR .AND. ELL2 <= ELR) THEN
                  IF (KUP > KB(UHS(JB))) THEN
                    KB(IU-1)    = K-1
                    KBMIN(IU-1) = MIN(KB(IU),KB(IU-1))
                    IEXIT       = 1
                    EXIT
                  END IF
                  ELR2 = EL(K+1,IU)+SINA(JB)*DLX(IU)*0.5
                  IF (ELR2 >= ELL2) THEN
                    B(K,IU-1) = B(KUP,UHS(JB)); EXIT
                  ELSE
                    K1 = KUP+1
                    IF (K1 > KMX) EXIT
                    B11 = 0.0
                    EL1 = ELR
                    EL2 = EL(K1,UHS(JB))-SINA(JBUH(JB))*DLX(IU)*0.5
                    DO WHILE (ELR2 <= EL2)
                      B11 = B11+(EL1-EL2)*B(K1-1,UHS(JB))
                      EL1 = EL2
                      K1  = K1+1
                      IF (K1 >= KMX+1 .OR. EL2 == ELR2) EXIT
                      EL2 = EL(K1,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                      IF (EL2 <= ELR2) EL2 = ELR2   
                    END DO
                    B(K,IU-1) = B11/H(K,JW); EXIT
                  END IF
                END IF
              END DO
              IF (EL(KMX,UHS(JB)) > EL(K,IU)) B(K,IU-1) = B(K-1,IU-1)
              IF (B(K,IU-1) == 0.0) B(K,IU-1) = B(K-1,IU-1)
              IF (IEXIT == 1) EXIT
            END IF
          END IF
        END IF
      END DO
      IEXIT = 0 
      DO K=1,KMX-1
        B(K,ID+1) = B(K,ID)
        IF (DH_INTERNAL(JB)) THEN
          IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
            B(K,ID+1) = B(K,DHS(JB))
          ELSE
            ELL = EL(K,ID)-SINA(JB)*DLX(ID)*0.5
            ELR = EL(2,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
            IF (ELL >= ELR) THEN
              B(K,ID+1) = B(2,DHS(JB))
            ELSE
              DO KDN=2,KMX-1
                ERR1 = EL(KDN,DHS(JB))  +SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                ERR2 = EL(KDN+1,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                IF (ERR1 >= ELL .AND. ERR2 < ELL) THEN
                  IF (KDN > KB(DHS(JB))) THEN
                    KB(ID+1)  = K-1
                    KBMIN(ID) = MIN(KB(ID),KB(ID+1))
                    IEXIT     = 1
                    EXIT
                  END IF
                  ELL2 = EL(K+1,ID)-SINA(JB)*DLX(ID)*0.5
                  IF (ELL2 >= ERR2) THEN
                    B(K,ID+1) = B(KDN,DHS(JB)); EXIT
                  ELSE
                    K1  = KDN+1
                    IF (K1 > KMX) EXIT
                    B11 = 0.0
                    EL2 = ELL
                    EL1 = EL(K1,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                    DO WHILE (ELL2 <= EL1)      
                      B11 = B11+(EL2-EL1)*B(K1-1,DHS(JB))
                      EL2 = EL1
                      K1  = K1+1
                      IF (K1 >= KMX+1 .OR. EL1 == ELL2) EXIT
                      EL1 = EL(K1,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                      IF (EL1 <= ELL2) EL1 = ELL2   
                    END DO
                    B(K,ID+1) = B11/H(K,JW); EXIT
                  END IF
                END IF
              END DO
              IF (EL(KMX,DHS(JB)) > EL(K,ID)) B(K,ID+1) = B(K-1,ID+1)
              IF (B(K,ID+1) == 0.0) B(K,ID+1) = B(K-1,ID+1)
              IF (IEXIT == 1) EXIT
            END IF
          END IF
        END IF
      END DO
    END DO   ! SW 1/23/06
   END DO    ! SW 1/23/06
   bnew=b    ! SW 1/23/06

!**** Upstream active segment and single layer  ! 1/23/06 entire section moved SW
  DO JW=1,NWB
    KT = KTWB(JW)
    DO JB=BS(JW),BE(JW)
      IU    = US(JB)
      ID    = DS(JB)
      IEXIT = 0
      IF (SLOPE(JB) /= 0.0) THEN
        DO I=US(JB)-1,DS(JB)+1
          IF (KBi(I) < KT ) THEN     ! SW 1/23/06                                              
            do k=kbi(i)+1,kt
            Bnew(K,I) = 0.000001      ! SW 1/23/06
            end do
            KB(I)   = KT
          END IF
        END DO
      END IF
      IUT = IU
      DO I=IU,ID
        IF (KB(I)-KT < NL(JB)-1) IUT = I+1
        ONE_LAYER(I) = KT == KB(I)
      END DO
      CUS(JB) = IUT
      ! reinitialize KBMIN
      DO I=IU-1,ID
        KBMIN(I) =  MIN(KB(I),KB(I+1))
      END DO
      KBMIN(ID+1) = KBMIN(ID)

!**** Areas and bottom widths

      IF (.NOT. TRAPEZOIDAL(JW)) THEN                                                                                  !SW 07/16/04
        DO I=IU-1,ID+1
          DO K=1,KMX-1
            BH2(K,I) = B(K,I)*H(K,JW)
            BH(K,I)  = B(K,I)*H(K,JW)
            BB(K,I)  = B(K,I)-(B(K,I)-B(K+1,I))/(0.5*(H(K,JW)+H(K+1,JW)))*H(K,JW)*0.5                                  !SW 08/02/04
          END DO
          BH(KMX,I) = BH(KMX-1,I)
        END DO

! v3.5 start
! column widths
        DO I=IU-1,ID+1
          cw(kb(i),i)=b(kb(i),i)
          do k=1,kb(i)-1
              cw(k,i)=b(k,i)-b(k+1,i)            
          end do
        end do
! v3.5 end

!****** Derived geometry

        DO I=IU-1,ID+1
          BH2(KT,I) = B(KTI(I),I)*(EL(KT,I)-EL(KTI(I)+1,I)-Z(I)*COSA(JB))/COSA(JB)
          IF (KT == KTI(I)) BH2(KT,I) = H2(KT,I)*B(KT,I)
          DO K=KTI(I)+1,KT
            BH2(KT,I) = BH2(KT,I)+Bnew(K,I)*H(K,JW)                                      ! sw 1/23/06    BH(K,I)
          END DO
          BKT(I)   = BH2(KT,I)/H2(KT,I)
          BI(KT,I) = B(KTI(I),I)
        END DO
      ELSE                                                                                                             !SW 07/16/04
        DO I=IU-1,ID+1
          DO K=1,KMX-1
            BB(K,I)  = B(K,I)-(B(K,I)-B(K+1,I))/(0.5*(H(K,JW)+H(K+1,JW)))*H(K,JW)*0.5 
          END DO
          BB(KB(I),I) = B(KB(I),I)*0.5
          BH2(1,I)    = B(1,I)*H(1,JW)  
          BH(1,I)     = BH2(1,I)   
          DO K=2,KMX-1
            BH2(K,I) = 0.25*H(K,JW)*(BB(K-1,I)+2.*B(K,I)+BB(K,I))
            BH(K,I)  = BH2(K,I)
          END DO
          BH(KMX,I) = BH(KMX-1,I)
        END DO
        DO I=IU-1,ID+1
          ELT = EL(KT,I)-Z(I)
          CALL GRID_AREA1(ELT,EL(KT+1,I),BH2(KT,I),BI(KT,I))
          BKT(I) = BH2(KT,I)/H2(KT,I)
        END DO
      END IF
      DO I=IU-1,ID
        DO K=1,KMX-1
          AVH2(K,I) = (H2(K,I) +H2(K+1,I)) *0.5
          AVHR(K,I) =  H2(K,I)+(H2(K,I+1)-H2(K,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                  !SW 07/29/04
        END DO
        AVH2(KMX,I) = H2(KMX,I)
        DO K=1,KMX
          BR(K,I)   = B(K,I)  +(B(K,I+1)  -B(K,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                !SW 07/29/04
          BHR(K,I)  = BH(K,I) +(BH(K,I+1) -BH(K,I)) /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                !SW 07/29/04
          BHR2(K,I) = BH2(K,I)+(BH2(K,I+1)-BH2(K,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                !SW 07/29/04
        END DO
      END DO
      DO K=1,KMX-1  
        AVH2(K,ID+1) = (H2(K,ID+1)+H2(K+1,ID+1))*0.5
        BR(K,ID+1)   =   B(K,ID+1)
        BHR(K,ID+1)  =   BH(K,ID+1)
      END DO
      AVH2(KMX,ID+1) = H2(KMX,ID+1)
      AVHR(KT,ID+1)  = H2(KT,ID+1)
      BHR2(KT,ID+1)  = BH2(KT,ID+1)
      IUT = IU
      IF (UP_HEAD(JB)) IUT = IU-1
      DO I=IUT,ID
        DO K=1,KMX-1
          VOL(K,I) = B(K,I)*H2(K,I)*DLX(I)
        END DO
        VOL(KT,I)    = BH2(KT,I)*DLX(I)
        DEPTHB(KT,I) = H2(KT,I)
        DEPTHM(KT,I) = H2(KT,I)*0.5
        DO K=KT+1,KMX
          DEPTHB(K,I) = DEPTHB(K-1,I)+ H2(K,I)
          DEPTHM(K,I) = DEPTHM(K-1,I)+(H2(K-1,I)+H2(K,I))*0.5
        END DO
      END DO
    END DO
  END DO
  H1   = H2
  BH1  = BH2
  BHR1 = BHR2
  AVH1 = AVH2

! Temporary downstream head segment

  DO JB=1,NBR
   IF (DHS(JB).GT.0) THEN
     DO JJB=1,NBR
       IF (DHS(JB) >= US(JJB) .AND. DHS(JB) <= DS(JJB)) EXIT
     END DO
     IF (CUS(JJB) > DHS(JB)) CDHS(JB) = CUS(JJB)
   END IF
  END DO

! Total active cells

  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      DO I=CUS(JB),DS(JB)
        DO K=KTWB(JW),KB(I)
          NTAC = NTAC+1
        END DO
      END DO
      NTACMX = NTAC
      NTACMN = NTAC

!**** Wind fetch lengths

      DO I=US(JB),DS(JB)
        FETCHD(I,JB) = FETCHD(I-1,JB)+DLX(I)
      END DO
      DO I=DS(JB),US(JB),-1
        FETCHU(I,JB) = FETCHU(I+1,JB)+DLX(I)
      END DO
    END DO
  END DO

! Segment heights

  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      DO I=US(JB)-1,DS(JB)+1
        DO K=MIN(KMX-1,KB(I)),2,-1
          HSEG(K,I) = HSEG(K+1,I)+H2(K,I)
        END DO
      END DO
    END DO
  END DO

! Beginning and ending segments/layers for snapshots

  DO JW=1,NWB
    DO I=1,NISNP(JW)
      KBR(JW) = MAX(KB(ISNP(I,JW)),KBR(JW))
    END DO
  END DO

! Density related derived constants

  RHOWCP  = RHOW*CP
  RHOICP  = RHOI*CP
  RHOIRL1 = RHOI*RL1
  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      DLXRHO(US(JB):DS(JB)) = 0.5/(DLXR(US(JB):DS(JB))*RHOW)
      IF (UP_HEAD(JB)) DLXRHO(US(JB)-1) = 0.5/(DLXR(US(JB))*RHOW)
    END DO
  END DO

! Transport interpolation multipliers

  DO JW=1,NWB
    CALL INTERPOLATION_MULTIPLIERS
  END DO

!***********************************************************************************************************************************
!**                                              Task 1.4.4: Initial conditions                                                   **
!***********************************************************************************************************************************

! macrophytes start
  do i=1,imx
    do k=1,kmx
      bic(k,i)=b(k,i)
    end do
  end do
! macrophytes sned

  DO JW=1,NWB
    KT = KTWB(JW)
    IF (VERT_PROFILE(JW)) THEN

!**** Temperature and water quality

      OPEN (VPR(JW),FILE=VPRFN(JW),STATUS='OLD')
      READ (VPR(JW),*)
      IF (VERT_TEMP(JW)) READ (VPR(JW),'(//(8X,9F8.0))') (TVP(K,JW),K=KT,KBMAX(JW))
      IF (CONSTITUENTS) THEN
        DO JC=1,NCT
          IF (VERT_CONC(JC,JW))      READ (VPR(JW),'(//(8X,9F8.0))') (CVP(K,JC,JW),  K=KT,KBMAX(JW))
        END DO
        DO JE=1,NEP
          IF (VERT_EPIPHYTON(JW,JE)) READ (VPR(JW),'(//(8X,9F8.0))') (EPIVP(K,JW,JE),K=KT,KBMAX(JW))
        END DO
        IF (VERT_SEDIMENT(JW))       READ (VPR(JW),'(//(8X,9F8.0))') (SEDVP(K,JW),   K=KT,KBMAX(JW))
      END IF
    END IF

!** Longitudinal/vertical initial profiles

    IF (LONG_PROFILE(JW)) THEN
      OPEN (LPR(JW),FILE=LPRFN(JW),STATUS='OLD')
      READ (LPR(JW),*)
    END IF

!** Branch related variables

    IF (.NOT. RESTART_IN) THEN
      DO JB=BS(JW),BE(JW)

!****** Temperature

        DO I=CUS(JB),DS(JB)
          IF (LONG_TEMP(JW)) READ (LPR(JW),'(//(8X,9F8.0))') (T1(K,I),K=KT,KB(I))
          DO K=KT,KB(I)
            IF (ISO_TEMP(JW))  T1(K,I) = T2I(JW)
            IF (VERT_TEMP(JW)) T1(K,I) = TVP(K,JW)
            T2(K,I) = T1(K,I)
          END DO
        END DO
      END DO

!**** Constituents

      DO JC=1,NAC
        DO JB=BS(JW),BE(JW)
          DO I=CUS(JB),DS(JB)
            JAC = CN(JC)
            IF (LONG_CONC(JAC,JW)) READ (LPR(JW),'(//(8X,9F8.0))') (C2(K,I,JAC),K=KT,KB(I))
            DO K=KT,KB(I)
              IF (ISO_CONC(JAC,JW))  C2(K,I,JAC) = C2I(JAC,JW)
              IF (VERT_CONC(JAC,JW)) C2(K,I,JAC) = CVP(K,JAC,JW)
              C1(K,I,JAC)  = C2(K,I,JAC)
              C1S(K,I,JAC) = C1(K,I,JAC)
            END DO
          END DO
        END DO
      END DO

!**** Epiphyton

      DO JB=BS(JW),BE(JW)
        DO JE=1,NEP
          IF (EPIPHYTON_CALC(JW,JE)) THEN
            DO I=CUS(JB),DS(JB)
              IF (LONG_EPIPHYTON(JW,JE)) READ (LPR(JW),'(//(8X,9F8.0))') (EPD(K,I,JE),K=KT,KB(I))
              IF (ISO_EPIPHYTON(JW,JE))  EPD(:,I,JE) = EPICI(JW,JE)
              IF (VERT_EPIPHYTON(JW,JE)) EPD(:,I,JE) = EPIVP(K,JW,JE)
            END DO
          END IF
        END DO
      END DO

!**** Sediments

      DO JB=BS(JW),BE(JW)
        IF (SEDIMENT_CALC(JW)) THEN
          DO I=CUS(JB),DS(JB)
            IF (LONG_SEDIMENT(JW)) READ (LPR(JW),'(//(8X,9F8.0))') (SED(K,I),K=KT,KB(I))
            DO K=KT,KB(I)
              IF (ISO_SEDIMENT(JW))  SED(K,I) = SEDCI(JW)
              IF (VERT_SEDIMENT(JW)) SED(K,I) = SEDVP(K,JW)
!v3.5 start
            END DO
            SED(KT,I)         = SED(KT,I)/H2(KT,I)
            SED(KT+1:KB(I),I) = SED(KT+1:KB(I),I)/H2(KT+1:KB(I),I)
            do k=kt,kb(i)                         !cb 11/30/06
              sdkv(k,i)=sdk(jw)                   !cb 11/30/06
              seddktot(k,i)=sed(k,i)*sdkv(k,i)    !cb 11/30/06
            end do                                !cb 11/30/06
          END DO
        END IF
      END DO
        DO JB=BS(JW),BE(JW)
          IF (SEDIMENT_CALC(JW)) THEN
            DO I=CUS(JB),DS(JB)
              DO K=KT,KB(I)
                IF (ISO_SEDIMENT(JW))SEDp(K,I) = orgp(JW)*sedci(jw)              
                IF (VERT_SEDIMENT(JW))SEDp(K,I) = SEDVP(K,JW)*orgp(jw)
                IF (LONG_SEDIMENT(JW)) sedp(k,i)=orgp(jw)*sed(k,i)
              END DO                       
              SEDp(KT,I)         = SEDp(KT,I)/H2(KT,I)
              SEDp(KT+1:KB(I),I) = SEDp(KT+1:KB(I),I)/H2(KT+1:KB(I),I)                          
            END DO
          END IF
        END DO
        DO JB=BS(JW),BE(JW)
          IF (SEDIMENT_CALC(JW)) THEN
            DO I=CUS(JB),DS(JB)
              DO K=KT,KB(I)
                IF (ISO_SEDIMENT(JW))SEDn(K,I) = orgn(JW)*sedci(jw)              
                IF (VERT_SEDIMENT(JW))SEDn(K,I) = SEDVP(K,JW)*orgn(jw)
                IF (LONG_SEDIMENT(JW)) sedn(k,i)=orgn(jw)*sed(k,i)
              END DO                       
              SEDn(KT,I)         = SEDn(KT,I)/H2(KT,I)
              SEDn(KT+1:KB(I),I) = SEDn(KT+1:KB(I),I)/H2(KT+1:KB(I),I)                          
          END DO
        END IF
      END DO
        DO JB=BS(JW),BE(JW)
          IF (SEDIMENT_CALC(JW)) THEN
            DO I=CUS(JB),DS(JB)
              DO K=KT,KB(I)
                IF (ISO_SEDIMENT(JW))SEDc(K,I) = SEDCI(JW)*orgc(jw)              
                IF (VERT_SEDIMENT(JW))SEDc(K,I) = SEDVP(K,JW)*orgc(jw)
                IF (LONG_SEDIMENT(JW)) sedc(k,i)=orgc(jw)*sed(k,i)
              END DO                       
              SEDc(KT,I)         = SEDc(KT,I)/H2(KT,I)
              SEDc(KT+1:KB(I),I) = SEDc(KT+1:KB(I),I)/H2(KT+1:KB(I),I)                          
            END DO
          END IF
        END DO      
! v3.5 end

      SED(:,US(BS(JW)):DS(BE(JW))) = SED(:,US(BS(JW)):DS(BE(JW)))*FSED(JW)
      SEDp(:,US(BS(JW)):DS(BE(JW))) = SEDp(:,US(BS(JW)):DS(BE(JW)))*FSED(JW)
      SEDn(:,US(BS(JW)):DS(BE(JW))) = SEDn(:,US(BS(JW)):DS(BE(JW)))*FSED(JW)
      SEDc(:,US(BS(JW)):DS(BE(JW))) = SEDc(:,US(BS(JW)):DS(BE(JW)))*FSED(JW)
      
      DO JB=BS(JW),BE(JW)  
        do m=1,nmc
          IF (macrophyte_CALC(jw,m)) THEN

!c distributing initial macrophyte conc to bottom column cells; macwbci = g/m^3              
            DO I=cus(jb),ds(jb)

              wlelev=(el(kt,i)-z(i)*cos(alpha(jb)))
              depkti=wlelev-el(kti(i)+1,i)

              if(depkti.ge.thrkti)then
                kticol(i)=.true.
                jt=kti(i)
              else
                kticol(i)=.false.
                jt=kti(i)+1
              end if
              
              je=kb(i)
              DO j=jt,je
                if(j.le.kt)then
                  k=kt
                else
                  k=j
                end if
                macrc(j,K,I,m) = macwbci(jw,m)
                smacrc(j,K,I,m) = macwbci(jw,m)
              END DO
            END DO

            DO I=cus(jb),ds(jb)
              tmac=0.0
              xsar=0.0
              do k=kti(i),kt
                jt=k
                je=kb(i)
                wlelev=(el(kt,i)-z(i)*cosa(jb))
                colb=el(k+1,i)
                coldep=wlelev-colb
                do j=jt,je                                
                  tmac=tmac+macrc(j,kt,i,m)*cw(j,i)*coldep
                  xsar=xsar+cw(j,i)*coldep
                end do
              end do
              mac(kt,i,m)=tmac/xsar
              smac(kt,i,m)=mac(kt,i,m)
                
              DO K=KT+1,KB(I)                  
                jt=k                    
                je=kb(i)
                tmac=0.0
                do j=jt,je
                  tmac=tmac+macrc(j,k,i,m)*cw(j,i)                      
                end do
                mac(k,i,m)=tmac/b(k,i)
                smac(k,i,m)=mac(k,i,m)
              end do
            end do

            do i=cus(jb),ds(jb)
              jt=kti(i)
              je=kb(i)
              wlelev=(el(kt,i)-z(i)*cosa(jb))                    
              do j=jt,je
                if(j.lt.kt)then               
                  colb=el(j+1,i)
                else
                  colb=el(kt+1,i)
                end if
                coldep=wlelev-colb               
                macrm(j,kt,i,m)=macrc(j,kt,i,m)*coldep*cw(j,i)*dlx(i)
                smacrm(j,kt,i,m)=macrm(j,kt,i,m)
              end do
                                 
              do K=KT+1,KB(I)
                
                jt=k                    
                je=kb(i)
                
                do j=jt,je
                   
                  macrm(j,k,i,m)=macrc(j,k,i,m)*h2(k,i)*cw(j,i)*dlx(i)     
                  smacrm(j,k,i,m)=macrm(j,k,i,m)
                end do                
                
              END DO
            end do

          END IF
        end do
      end do
! v3.5 end

!**** Energy

      DO JB=BS(JW),BE(JW)
        DO I=CUS(JB),DS(JB)
          IF (ENERGY_BALANCE(JW)) THEN
            DO K=KT,KB(I)
              EBRI(JB) = EBRI(JB)+T2(K,I)*DLX(I)*BH2(K,I)
            END DO
          END IF
          DO K=KT,KB(I)
            CMBRT(CN(1:NAC),JB) = CMBRT(CN(1:NAC),JB)+C2(K,I,CN(1:NAC))*DLX(I)*BH2(K,I)
          END DO
        END DO

! v3.5 start
!c   initializing macrophyte temporal mass balance term....
          do m=1,nmc
            if(macrophyte_calc(jw,m))then                
              DO I=CUS(JB),DS(JB)
                if(kticol(i))then               
                  jt=kti(i)
                else                
                  jt=kti(i)+1
                end if                                                                                              
                je=kb(i)                    
                do j=jt,je               
                  maCMBRT(JB,m) = maCMBRT(JB,m)+macrm(j,Kt,I,m)     
                end do                               
                DO K=KT+1,KB(I)               
                  jt=k                    
                  je=kb(i)              
                  do j=jt,je                    
                    maCMBRT(JB,m) = maCMBRT(JB,m)+macrm(j,K,I,m)     
                  end do                                  
                END DO
              END DO
            end if
          end do
! v3.5 end

!****** Ice cover

        IF (ICE_CALC(JW)) THEN
          ICETH(US(JB):DS(JB)) = ICETHI(JW)
          ICE(US(JB):DS(JB))   = ICETH(US(JB):DS(JB)) > 0.0
        END IF

!****** Vertical eddy viscosity

        IUT = CUS(JB)
        IDT = DS(JB)-1
        IF (UP_HEAD(JB)) IUT = IU-1
        IF (DN_HEAD(JB)) IDT = ID
        DO I=IUT,IDT
!         DO K=KT,KBMIN(I)-1 
          DO K=KT,KB(I)-1       !SG 10/4/07
            AZ(K,I)    = AZMIN
            TKE(K,I,1) = 1.25E-7
            TKE(K,I,2) = 1.0E-9
          END DO
        END DO
        DO JWR=1,NIW
          IF (WEIR_CALC) AZ(MAX(KT,KTWR(JWR)-1):KBWR(JWR),IWR(JWR)) = 0.0
        END DO
      END DO
    END IF

!** Horizontal diffusivities

    DO JB=BS(JW),BE(JW)
      DO I=CUS(JB),DS(JB)-1
        DO K=KT,KBMIN(I)
          DX(K,I) = DXI(JW)
          IF (INTERNAL_WEIR(K,I)) DX(K,I) = 0.0
        END DO
      END DO
    END DO
    IF (VERT_PROFILE(JW)) CLOSE (VPR(JW))
    IF (LONG_PROFILE(JW)) CLOSE (LPR(JW))
  END DO

! Atmospheric pressure

  IF (CONSTITUENTS) PALT(:) = (1.0-EL(KT,:)/1000.0/44.3)**5.25

! Saved variables for autostepping

  IF (.NOT. RESTART_IN) THEN
    SZ    = Z
    SU    = U
    SW    = W
    SAZ   = AZ
    SKTI  = KTI
    SBKT  = BKT
    SAVH2 = AVH2
    SAVHR = AVHR
  END IF
  CALL GREGORIAN_DATE
  CALL TIME_VARYING_DATA
  CALL READ_INPUT_DATA (NXTVD)
  IF (CONSTITUENTS) THEN
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IU = US(JB)
        ID = DS(JB)
        CALL TEMPERATURE_RATES
        CALL KINETIC_RATES
      END DO
    END DO
  END IF

!***********************************************************************************************************************************
!*                                                           Task 1.5: Outputs                                                    **
!***********************************************************************************************************************************

! Open output files

  IF (RESTART_IN) THEN
    DO JW=1,NWB
      IF (SNAPSHOT(JW))    OPEN (SNP(JW),FILE=SNPFN(JW),POSITION='APPEND')
      IF (VECTOR(JW))      OPEN (VPL(JW),FILE=VPLFN(JW),POSITION='APPEND')
      IF (SPREADSHEET(JW)) OPEN (SPR(JW),FILE=SPRFN(JW),POSITION='APPEND')
      IF (CONTOUR(JW))     OPEN (CPL(JW),FILE=CPLFN(JW),POSITION='APPEND')
      IF (PROFILE(JW))     OPEN (PRF(JW),FILE=PRFFN(JW),POSITION='APPEND')
      IF (FLUX(JW))        OPEN (FLX(JW),FILE=FLXFN(JW),POSITION='APPEND')
      IF (SNAPSHOT(JW)) THEN
        REWIND (SNP(JW))
        DO WHILE (.TRUE.)
          READ (SNP(JW),'(A72)',END=100) LINE
          IF (LINE(26:28) == 'NIT') THEN
            BACKSPACE SNP(JW)
            READ (SNP(JW),'(31X,I10)',END=100) NIT1
            IF (NIT1 > NIT) THEN
              DO J=1,24
                BACKSPACE (SNP(JW))
              END DO
              EXIT
            END IF
          END IF
        END DO
      END IF
100   CONTINUE
      IF (SPREADSHEET(JW)) THEN
        REWIND (SPR(JW))
        READ (SPR(JW),*)
        DO WHILE (JDAY1 < JDAY)
          READ (SPR(JW),'(A,F10.0)',END=101) LINE(1:38),JDAY1
        END DO
        BACKSPACE (SPR(JW))
101     CONTINUE
        JDAY1 = 0.0
      END IF
      IF (PROFILE(JW)) THEN
        REWIND (PRF(JW))
        READ   (PRF(JW),'(A)')        (LINE,J=1,11)
        READ   (PRF(JW),'(8I8)')       I
        READ   (PRF(JW),'(10I8)')     (I,J=1,NIPRF(JW))
        READ   (PRF(JW),'(20(1X,A))')  LINE (1:8), (LINE (1:3), JC=1,NCT),(LINE (1:3), JD=1,NDC)
        READ   (PRF(JW),'(2A)')        LINE (1:26),(LINE (1:26),JC=1,NCT),(LINE (1:43),JD=1,NDC)
        DO WHILE (JDAY1 < JDAY)
          READ (PRF(JW),'(A72)',END=102) LINE
          L1 = 0      
          L1 = SCAN(LINE,',')
          IF (L1 /= 0) THEN
            BACKSPACE (PRF(JW))
            READ (PRF(JW),'(F8.0)',END=102) JDAY1
          END IF
        END DO
        JDAY1 = 0.0
      END IF
102   CONTINUE
      JDAY1 = 0.0
      IF (CONTOUR(JW)) THEN
        REWIND (CPL(JW))
        DO WHILE (JDAY1 < JDAY)
          READ (CPL(JW),'(A72)',END=103) LINE
          IF (LINE(1:8) == 'New date') THEN
            BACKSPACE (CPL(JW))
            READ (CPL(JW),'(A,F12.4)',END=103) LINE(1:9),JDAY1
          END IF
        END DO
        BACKSPACE (CPL(JW))
        JDAY1 = 0.0
      END IF
103   CONTINUE
      IF (VECTOR(JW)) THEN
        REWIND (VPL(JW))
        DO WHILE (JDAY1 < JDAY)
          READ (VPL(JW),'(A72)',END=104) LINE
          IF (LINE(1:8) == 'New date') THEN
            BACKSPACE (VPL(JW))
            READ (VPL(JW),*,END=104) LINE(1:9),JDAY1
          END IF
        END DO
        BACKSPACE (VPL(JW))
        JDAY1 = 0.0
      END IF
104   CONTINUE
      IF (FLUX(JW)) THEN
        REWIND (FLX(JW))
        DO WHILE (JDAY1 < JDAY)
          READ (FLX(JW),'(A72)',END=105) LINE
          IF (LINE(1:8) == 'New date') THEN
            BACKSPACE (FLX(JW))
            READ (FLX(JW),'(8X,F10.0)',END=105) JDAY1
          END IF
        END DO
        BACKSPACE (FLX(JW))
      END IF
105   CONTINUE
    END DO     
    IF (DOWNSTREAM_OUTFLOW) THEN
      DO JWD=1,NIWDO
        WRITE (SEGNUM,'(I0)') IWDO(JWD)
        SEGNUM = ADJUSTL(SEGNUM)
        L      = LEN_TRIM(SEGNUM)
        OPEN   (WDO(JWD,1),FILE='qwo_'//SEGNUM(1:L)//'.opt',POSITION='APPEND')
        REWIND (WDO(JWD,1))
        READ   (WDO(JWD,1),'(//)')
        DO WHILE (JDAY1 < JDAY)
          READ (WDO(JWD,1),'(F8.0)',END=106) JDAY1
        END DO
        BACKSPACE (WDO(JWD,1))
106     CONTINUE
        OPEN   (WDO(JWD,2),FILE='two_'//SEGNUM(1:L)//'.opt',POSITION='APPEND')
        REWIND (WDO(JWD,2))
        READ   (WDO(JWD,2),'(//)')
        DO WHILE (JDAY1 < JDAY)
          READ (WDO(JWD,2),'(F8.0)',END=107) JDAY1
        END DO
        BACKSPACE (WDO(JWD,2))
107     CONTINUE
        IF (CONSTITUENTS) THEN
          OPEN   (WDO(JWD,3),FILE='cwo_'//SEGNUM(1:L)//'.opt',POSITION='APPEND')
          REWIND (WDO(JWD,3))
          READ   (WDO(JWD,3),'(//)')
          DO WHILE (JDAY1 < JDAY)
            READ (WDO(JWD,3),'(F8.0)',END=108) JDAY1
          END DO
          BACKSPACE (WDO(JWD,3))
108       CONTINUE
        END IF
        IF (DERIVED_CALC) THEN
          OPEN   (WDO(JWD,4),FILE='dwo_'//SEGNUM(1:L)//'.opt',POSITION='APPEND')
          REWIND (WDO(JWD,4))
          READ   (WDO(JWD,4),'(//)')
          DO WHILE (JDAY1 < JDAY)
            READ (WDO(JWD,4),'(F8.0)',END=109) JDAY1
          END DO
          BACKSPACE (WDO(JWD,4))
109       CONTINUE
        END IF
      END DO
    END IF
    IF (TIME_SERIES) THEN
      L1 = SCAN(TSRFN,'.')
      DO J=1,NIKTSR
        WRITE (SEGNUM,'(I0)') J
        SEGNUM = ADJUSTL(SEGNUM)
        L      = LEN_TRIM(SEGNUM)
        TSRFN  = TSRFN(1:L1-1)//'_'//SEGNUM(1:L)//'.opt'
        OPEN   (TSR(J),FILE=TSRFN,POSITION='APPEND')
        REWIND (TSR(J))
        READ   (TSR(J),'(A72)')   (LINE,I=1,11)
        READ   (TSR(J),'(/F10.3)',END=110) JDAYTS
        DO WHILE (JDAYTS < JDAY)
          READ (TSR(J),'(F10.0)',END=110) JDAYTS
        END DO
        BACKSPACE (TSR(J))
110     CONTINUE
      END DO
    END IF
  ELSE
    DO JW=1,NWB
      IF (SNAPSHOT(JW))    OPEN (SNP(JW),FILE=SNPFN(JW),STATUS='UNKNOWN')
      IF (VECTOR(JW))      OPEN (VPL(JW),FILE=VPLFN(JW),STATUS='UNKNOWN')
      IF (PROFILE(JW))     OPEN (PRF(JW),FILE=PRFFN(JW),STATUS='UNKNOWN')
      IF (SPREADSHEET(JW)) OPEN (SPR(JW),FILE=SPRFN(JW),STATUS='UNKNOWN')
      IF (CONTOUR(JW))     OPEN (CPL(JW),FILE=CPLFN(JW),STATUS='UNKNOWN')
      IF (FLUX(JW))        OPEN (FLX(JW),FILE=FLXFN(JW),STATUS='UNKNOWN')

!**** Output files

      IF (PROFILE(JW)) THEN
        TIME = TMSTRT
        DO WHILE (TIME <= TMEND)
          NDSP = NDSP+1
          TIME = TIME+PRFF(PRFDP(JW),JW)
          IF (TIME >= PRFD(PRFDP(JW)+1,JW)) PRFDP(JW) = PRFDP(JW)+1
        END DO
        PRFDP(JW) = 1
        WRITE (PRF(JW),'(A)')         TITLE
        WRITE (PRF(JW),'(8I8,L2)')    KMX,NIPRF(JW),NDSP,NCT,NDC,NAC+NACD(JW)+1,PRFDP(JW),KTWB(JW),CONSTITUENTS
        WRITE (PRF(JW),'(10I8)')      IPRF(1:NIPRF(JW),JW)
        WRITE (PRF(JW),'(20(1X,A))') ' ON',CPRWBC(:,JW)(6:8),CDWBC(:,JW)(6:8)
        WRITE (PRF(JW),'(2A)')       'Temperature, \F8C                            ',ADJUSTL(CNAME),ADJUSTL(CDNAME)
        WRITE (PRF(JW),'(20I4)')      1,CN(1:NAC)+1,CDN(1:NACD(JW),JW)+NCT+1
        WRITE (PRF(JW),'(10F8.0)')    1.0,CMULT,CDMULT
        WRITE (PRF(JW),'(20I4)')      KB(IPRF(1:NIPRF(JW),JW))
        WRITE (PRF(JW),'(10F8.2)')    H
        DO JP=1,NIPRF(JW)
          NRS = KB(IPRF(JP,JW))-KTWB(JW)+1
          WRITE (PRF(JW),'(A8,I4/(8(F10.2)))') 'TEMP    ',NRS,(T2(K,IPRF(JP,JW)),K=KTWB(JW),KB(IPRF(JP,JW)))
        END DO
        DO JC=1,NAC
          IF (PRINT_CONST(CN(JC),JW)) THEN
            DO JP=1,NIPRF(JW)
              NRS = KB(IPRF(JP,JW))-KTWB(JW)+1
              WRITE (PRF(JW),'(A,I4/(8(E13.6,1x)))') ADJUSTL(CNAME2(CN(JC))),NRS,(C2(K,IPRF(JP,JW),CN(JC))*CMULT(CN(JC)),&
                                                     K=KTWB(JW),KB(IPRF(JP,JW)))    
            END DO
          END IF
        END DO
        DO JD=1,NACD(JW)
          DO JP=1,NIPRF(JW)
            NRS = KB(IPRF(JP,JW))-KTWB(JW)+1
            WRITE (PRF(JW),'(A,I4/(8(E13.6,1x)))') ADJUSTL(CDNAME2(CDN(JD,JW))),NRS,&
                                                   (CD(K,IPRF(JP,JW),CDN(JD,JW))*CDMULT(CDN(JD,JW)),&
                                                    K=KTWB(JW),KB(IPRF(JP,JW)))    
          END DO
        END DO
      END IF
      IF (SPREADSHEET(JW)) THEN
        DO J=1,NISPR(JW)
          WRITE (SEGNUM,'(I0)') ISPR(J,JW)
          SEGNUM = ADJUSTL(SEGNUM)
          L      = LEN_TRIM(SEGNUM)
          SEG(J) = 'Seg_'//SEGNUM(1:L)
        END DO
        WRITE (SPR(JW),'(A,27X,A,5X,A,1000(1X,"Elevation",2X,A7))') 'Constituent','Julian_day','Depth',(SEG(J),J=1,NISPR(JW))
      END IF
      IF (CONTOUR(JW)) THEN
        WRITE (CPL(JW),'(A)')           TITLE
        WRITE (CPL(JW),'(8(I8,2X))')    NBR
        WRITE (CPL(JW),'(8(I8,2X))')    IMX,KMX
        DO JB=BS(JW),BE(JW)
          WRITE (CPL(JW),'(9(I8,2X))')  US(JB),DS(JB)
          WRITE (CPL(JW),'(9(I8,2X))')  KB(US(JB):DS(JB))
        END DO
        WRITE (CPL(JW),'(8(E13.6,2X))') DLX
        WRITE (CPL(JW),'(8(E13.6,2X))') H
        WRITE (CPL(JW),'(8(I8,2X))')    NAC
        WRITE (CPL(JW),'(A)')           CNAME1(CN(1:NAC))
      END IF
      IF (VECTOR(JW)) THEN
        WRITE (VPL(JW),*)  TITLE
        WRITE (VPL(JW),*)  H,KB,US,DS,DLX
      END IF
    END DO
    IF (TIME_SERIES) THEN
      L1 = SCAN(TSRFN,'.')
      DO J=1,NIKTSR
        WRITE (SEGNUM,'(I0)') J
        SEGNUM = ADJUSTL(SEGNUM)
        L      = LEN_TRIM(SEGNUM)
        TSRFN  = TSRFN(1:L1-1)//'_'//SEGNUM(1:L)//'.opt'
        OPEN  (TSR(J),FILE=TSRFN,STATUS='UNKNOWN')
        WRITE (TSR(J),'(A)') (TITLE(I),I=1,11)
        I = ITSR(J)                                                                                            ! SR 5/10/05 
        DO JW=1,NWB
          IF (I >= US(BS(JW)) .AND. I <= DS(BE(JW))) EXIT
        END DO
! v3.5 start
        IF (ICE_COMPUTATION) THEN
          if(sediment_calc(jw))then
          WRITE (TSR(J),'(1000(2X,A))') '    JDAY','     DLT','    ELWS','      T2','       U','       Q','    SRON','      ET',   &
                                        '   DEPTH','   WIDTH','   SHADE','   ICETH',(CNAME2(CN(JC)),JC=1,NAC),                     &
                                        ('     EPI',JE=1,NEP),('     MAC',Jm=1,nmc),'     SED','    SEDP','    SEDN','    SEDC',   &
                                        (CDNAME2(CDN(JD,JW)),JD=1,NACD(JW))
          else
            WRITE (TSR(J),'(1000(2X,A))') '    JDAY','     DLT','    ELWS','      T2','       U','       Q','    SRON','      ET', &
                                        '   DEPTH','   WIDTH','   SHADE','   ICETH',(CNAME2(CN(JC)),JC=1,NAC),                     &
                                        ('     EPI',JE=1,NEP),('     MAC',Jm=1,nmc),(CDNAME2(CDN(JD,JW)),JD=1,NACD(JW))     
          end if
        ELSE
          if(sediment_calc(jw))then                                                                                                    !mlm 7/25/06
          WRITE (TSR(J),'(1000(2X,A))') '    JDAY','     DLT','    ELWS','      T2','       U','       Q','    SRON','      ET',   &
                                        '   DEPTH','   WIDTH','   SHADE',(CNAME2(CN(JC)),JC=1,NAC),                     &
                                        ('     EPI',JE=1,NEP),('     MAC',Jm=1,nmc),'     SED','    SEDP','    SEDN','    SEDC',   &
                                        (CDNAME2(CDN(JD,JW)),JD=1,NACD(JW))
          else
            WRITE (TSR(J),'(1000(2X,A))') '    JDAY','     DLT','    ELWS','      T2','       U','       Q','    SRON','      ET', &
                                        '   DEPTH','   WIDTH','   SHADE',(CNAME2(CN(JC)),JC=1,NAC),                     &
                                        ('     EPI',JE=1,NEP),('     MAC',Jm=1,nmc),(CDNAME2(CDN(JD,JW)),JD=1,NACD(JW))    
          end if
        END IF
! v3.5 end

      END DO
    END IF
    IF (DOWNSTREAM_OUTFLOW) THEN
      DO JWD=1,NIWDO
        WRITE (SEGNUM,'(I0)') IWDO(JWD)
        SEGNUM = ADJUSTL(SEGNUM)
        L      = LEN_TRIM(SEGNUM)
        OPEN  (WDO(JWD,1),FILE='qwo_'//SEGNUM(1:L)//'.opt',STATUS='UNKNOWN')
        OPEN  (WDO(JWD,2),FILE='two_'//SEGNUM(1:L)//'.opt',STATUS='UNKNOWN')
        WRITE (WDO(JWD,1),'(A,I0//A)') 'Flow file for segment ',       IWDO(JWD),'    JDAY     QWD'
        WRITE (WDO(JWD,2),'(A,I0//A)') 'Temperature file for segment ',IWDO(JWD),'    JDAY       T'
        DO JW=1,NWB
          IF (IWDO(JWD) >= US(BS(JW)) .AND. IWDO(JWD) <= DS(BE(JW))) EXIT
        END DO
        IF (CONSTITUENTS) THEN
          OPEN  (WDO(JWD,3),FILE='cwo_'//SEGNUM(1:L)//'.opt',STATUS='UNKNOWN')
          WRITE (WDO(JWD,3),'(A,I0//,(100A))') 'Concentration file for segment ',     IWDO(JWD),'    JDAY',CNAME2(CN(1:NAC))
        END IF
        IF (DERIVED_CALC) THEN
          OPEN  (WDO(JWD,4),FILE='dwo_'//SEGNUM(1:L)//'.opt',STATUS='UNKNOWN')
          WRITE (WDO(JWD,4),'(A,I0//(100A))') 'Derived constituent file for segment ',IWDO(JWD),'    JDAY',                        &
                                               CDNAME2(CDN(1:NACD(JW),JW))
        END IF
      END DO
    END IF
  END IF
  IF (.NOT. RESTART_IN) CALL CPU_TIME (START)

  if (macrophyte_on.and.constituents) call porosity  !v3.5

!***********************************************************************************************************************************
!**                                                   Task 2: Calculations                                                        **
!***********************************************************************************************************************************

  DO WHILE (.NOT. END_RUN)
    IF (JDAY >= NXTVD) CALL READ_INPUT_DATA (NXTVD)
    CALL INTERPOLATE_INPUTS
    DLTTVD = (NXTVD-JDAY)*DAY
    DLT    =  MIN(DLT,DLTTVD+1.0)
    DLTS1  =  DLT
    IF (DLT <= DLTTVD+0.999) THEN
      DLTS = DLT
    ELSE
      KLOC = 1
      ILOC = 1
    END IF

!***********************************************************************************************************************************
!**                                            Task 2.1: Hydrodynamic sources/sinks                                               **
!***********************************************************************************************************************************

!** Timestep violation entry point

210 CONTINUE
    QINSUM = 0.0; TINSUM = 0.0; CINSUM = 0.0; UXBR = 0.0; UYBR = 0.0
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IU   = CUS(JB)
        ID   = DS(JB)
        TSUM = 0.0; CSUM = 0.0; QSUM(JB) = 0.0; QOUT(:,JB) = 0.0

!****** Densities

        DO I=IU-1,ID+1
          DO K=KT,KB(I)
            TISS(K,I) = 0.0
            DO JS=1,NSS
              TISS(K,I) = TISS(K,I)+SS(K,I,JS)
            END DO
            RHO(K,I) = DENSITY(T2(K,I),MAX(TDS(K,I),0.0),MAX(TISS(K,I),0.0))
          END DO
        END DO
! v3.5 deleted pumpback code from v3.2
        DO JS=1,NSTR(JB)
          IF (QSTR(JS,JB) /= 0.0) THEN
            CALL DOWNSTREAM_WITHDRAWAL (JS)
          END IF
        END DO
        DO K=KT,KB(ID)
          QSUM(JB)        = QSUM(JB)       +QOUT(K,JB)
          TSUM            = TSUM           +QOUT(K,JB)*T2(K,ID)
          CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QOUT(K,JB)*C2(K,ID,CN(1:NAC))
        END DO
        IF (QSUM(JB) /= 0.0) THEN
          TOUT(JB)           = TSUM           /QSUM(JB)
          COUT(CN(1:NAC),JB) = CSUM(CN(1:NAC))/QSUM(JB)
        END IF
        IF (QSUM(JB) /= 0.0 .AND. DAM_OUTFLOW(JB)) THEN                                                                !TC 08/03/04
          TINSUM(JBDAM(JB))           = (TSUM           +QINSUM(JBDAM(JB))*TINSUM(JBDAM(JB)))          /(QSUM(JB)                  &
                                        +QINSUM(JBDAM(JB)))
          CINSUM(CN(1:NAC),JBDAM(JB)) = (CSUM(CN(1:NAC))+QINSUM(JBDAM(JB))*CINSUM(CN(1:NAC),JBDAM(JB)))/(QSUM(JB)                  &
                                        +QINSUM(JBDAM(JB)))
          QINSUM(JBDAM(JB))           =  QINSUM(JBDAM(JB))+QSUM(JB)
        END IF     
      END DO
    END DO
    ILAT = 0
    JWW  = NWD
    JTT  = NTR
    JSS  = NSTR
    IF (SPILLWAY) THEN
      CALL SPILLWAY_FLOW
      DO JS=1,NSP

!****** Positive flows

        JLAT = 0
        JBU  = JBUSP(JS)
        JBD  = JBDSP(JS)
        IF (QSP(JS) > 0.0) THEN
          IF (LATERAL_SPILLWAY(JS)) THEN
            JLAT      = 1
            JWW       = JWW+1
            IWD(JWW)  = IUSP(JS)
            QWD(JWW)  = QSP(JS)
            KTWD(JWW) = KTUSP(JS)
            KBWD(JWW) = KBUSP(JS)
            EWD(JWW)  = ESP(JS)
            JBWD(JWW) = JBU
            I         = MAX(CUS(JBWD(JWW)),IWD(JWW))
            JB        = JBWD(JWW)
            JW        = JWUSP(JS)
            KT        = KTWB(JW)
            CALL LATERAL_WITHDRAWAL (JWW)
            DO K=KTW(JWW),KBW(JWW)
              QSS(K,I) = QSS(K,I)-QSW(K,JWW)
            END DO
          ELSE
            JSS(JBU)                 =  JSS(JBU)+1
            KTSW(JSS(JBU),JBU)       =  KTUSP(JS)
            KBSW(JSS(JBU),JBU)       =  KBUSP(JS)
            JB                       =  JBU
            POINT_SINK(JSS(JBU),JBU) = .TRUE.
            ID                       =  IUSP(JS)
            QSTR(JSS(JBU),JBU)       =  QSP(JS)
            ESTR(JSS(JBU),JBU)       =  ESP(JS)
            KT                       =  KTWB(JWUSP(JS))
            JW                       =  JWUSP(JS)
            CALL DOWNSTREAM_WITHDRAWAL(JSS(JBU))
            IF (IDSP(JS) /= 0 .AND. US(JBD) == IDSP(JS)) THEN
              QSUMM = 0.0; TSUM  = 0.0;  CSUM  = 0.0
              DO K=KT,KB(ID)
                QSUMM = QSUMM+QNEW(K)
                TSUM  = TSUM+QNEW(K)*T2(K,ID)
                DO JC=1,NAC
                  IF (CAC(NDO) == '      ON' .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN
                    CALL TOTAL_DISSOLVED_GAS (PALT(ID),0,JS,T2(K,ID),CGAS)
                    CSUM(CN(JC)) = CSUM(CN(JC))+QNEW(K)*CGAS
                  ELSE
                    CSUM(CN(JC)) = CSUM(CN(JC))+QNEW(K)*C2(K,ID,CN(JC))
                  END IF
                END DO
              END DO
              IF (QSUMM /= 0.0) THEN
                TINSUM(JBD)           = (TSUM           +QINSUM(JBD)*TINSUM(JBD))          /(QSUMM+QINSUM(JBD))
                CINSUM(CN(1:NAC),JBD) = (CSUM(CN(1:NAC))+QINSUM(JBD)*CINSUM(CN(1:NAC),JBD))/(QSUMM+QINSUM(JBD))
                QINSUM(JBD)           =  QINSUM(JBD)+QSUMM
              END IF   
            END IF  
            QSUM(JB) = 0.0; TSUM = 0.0; CSUM = 0.0
            DO K=KT,KB(ID)
              QSUM(JB) = QSUM(JB)+QOUT(K,JB)
              TSUM     = TSUM+QOUT(K,JB)*T2(K,ID)
              DO JC=1,NAC
                IF (CAC(NDO) == '      ON' .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN
                  CALL TOTAL_DISSOLVED_GAS (PALT(ID),0,JS,T2(K,ID),CGAS)
                  CSUM(CN(JC)) = CSUM(CN(JC))+QOUT(K,JB)*CGAS
                ELSE
                  CSUM(CN(JC)) = CSUM(CN(JC))+QOUT(K,JB)*C2(K,ID,CN(JC))
                END IF
              END DO
            END DO
            IF (QSUM(JB) /= 0.0) THEN
              TOUT(JB)           = TSUM           /QSUM(JB)
              COUT(CN(1:NAC),JB) = CSUM(CN(1:NAC))/QSUM(JB)
            END IF   
          END IF   
          IF (IDSP(JS) /= 0) THEN
            IF (US(JBD) /= IDSP(JS) .OR. HEAD_FLOW(JBD) .OR. UP_HEAD(JBD)) THEN
              JTT              = JTT+1
              QTR(JTT)         = QSP(JS)
              ITR(JTT)         = IDSP(JS)
              PLACE_QTR(JTT)   = PDSPC(JS) == ' DENSITY'
              SPECIFY_QTR(JTT) = PDSPC(JS) == ' SPECIFY'
              IF (SPECIFY_QTR(JTT)) THEN
                ELTRT(JTT) = ETDSP(JS)
                ELTRB(JTT) = EBDSP(JS)
              END IF
              JBTR(JTT) = JBD
              IF (JLAT == 1) THEN   
                CSUM(CN(1:NAC)) =  0.0; TSUM = 0.0; QSUMM = 0.0
                DO K=KTW(JWW),KBW(JWW)
                  QSUMM           = QSUMM          +QSW(K,JWW)
                  TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
                  CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
                END DO
                TTR(JTT) = TSUM/QSUMM
                DO JC=1,NAC
                  CTR(CN(JC),JTT) = CSUM(CN(JC))/QSUMM
                  IF (CN(JC) == NDO .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN
                    TDG_SPILLWAY(JWW,JS) = .TRUE.      
                    CALL TOTAL_DISSOLVED_GAS (PALT(ITR(JTT)),0,JS,TTR(JTT),CTR(CN(JC),JTT))
                  END IF
                END DO
              ELSE   
                TTR(JTT) =  TOUT(JB)
                DO JC=1,NAC
                  CTR(CN(JC),JTT) = COUT(CN(JC),JB)
                  IF (CN(JC) == NDO .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN
                    CALL TOTAL_DISSOLVED_GAS (PALT(ITR(JTT)),0,JS,TTR(JTT),CTR(CN(JC),JTT))
                  END IF
                END DO
              END IF     
            ELSE IF (LATERAL_SPILLWAY(JS)) THEN
              TSUM      = 0.0; QSUMM = 0.0; CSUM = 0.0
              ILAT(JWW) = 1
              DO K=KTW(JWW),KBW(JWW)
                QSUMM           = QSUMM          +QSW(K,JWW)
                TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
                CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
              END DO
              TINSUM(JBD)           = (TINSUM(JBD)          *QINSUM(JBD)+TSUM)           /(QSUMM+QINSUM(JBD))
              CINSUM(CN(1:NAC),JBD) = (CINSUM(CN(1:NAC),JBD)*QINSUM(JBD)+CSUM(CN(1:NAC)))/(QSUMM+QINSUM(JBD))
              QINSUM(JBD)           =  QSUMM                +QINSUM(JBD)
            ELSE IF (CAC(NDO) == '      ON' .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN
              TDG_SPILLWAY(JWW,JS) = .TRUE.
            END IF
          END IF
        ELSE IF (QSP(JS) < 0.0) THEN
          JTT              =  JTT+1
          JWW              =  JWW+1
          IWD(JWW)         =  IDSP(JS)
          ITR(JTT)         =  IUSP(JS)
          QTR(JTT)         = -QSP(JS)
          QWD(JWW)         = -QSP(JS)
          KTWD(JWW)        =  KTDSP(JS)
          KBWD(JWW)        =  KBDSP(JS)
          EWD(JWW)         =  ESP(JS)
          PLACE_QTR(JTT)   =  PUSPC(JS) == ' DENSITY'
          SPECIFY_QTR(JTT) =  PUSPC(JS) == ' SPECIFY'
          IF (SPECIFY_QTR(JTT)) THEN
            ELTRT(JTT) = ETUSP(JS)
            ELTRB(JTT) = EBUSP(JS)
          END IF
          JBTR(JTT) = JBU
          JBWD(JWW) = JBD
          I         = MAX(CUS(JBWD(JWW)),IWD(JWW))
          JB        = JBWD(JWW)
          JW        = JWDSP(JS)
          KT        = KTWB(JW)
          CALL LATERAL_WITHDRAWAL (JWW)
          DO K=KTW(JWW),KBW(JWW)
            QSS(K,I) = QSS(K,I)-QSW(K,JWW)
          END DO
          IF (IDSP(JS) /= 0) THEN
            TSUM  =  0.0; QSUMM = 0.0; CSUM = 0.0
            DO K=KTW(JWW),KBW(JWW)
              QSUMM           = QSUMM          +QSW(K,JWW)
              TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
              CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
            END DO
            TTR(JTT) = TSUM/QSUMM
            DO JC=1,NAC
              CTR(CN(JC),JTT) = CSUM(CN(JC))/QSUMM
              IF (CN(JC) == NDO .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN
                TDG_SPILLWAY(JWW,JS) = .TRUE.      
                CALL TOTAL_DISSOLVED_GAS (PALT(I),0,JS,TTR(JTT),CTR(CN(JC),JTT))
              END IF
            END DO
          ELSE IF (CAC(NDO) == '      ON' .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN
            TDG_SPILLWAY(JWW,JS) = .TRUE.
          END IF
        END IF
      END DO
    END IF
    IF (PUMPS) THEN
      DO JP=1,NPU
        JLAT = 0
        JWU  = JWUPU(JP)
        JBU  = JBUPU(JP)
        JBD  = JBDPU(JP)
        IF (JDAY >= ENDPU(JP)) PUMPON(JP) = .FALSE.                                                        !  CB 1/13/06
        IF (JDAY >= STRTPU(JP) .AND. JDAY < ENDPU(JP)) THEN
          IF (LATERAL_PUMP(JP)) THEN
            ELW = EL(KTWB(JWU),IUPU(JP))-Z(IUPU(JP))*COSA(JBU)
          ELSE
            ELW = EL(KTWB(JWU),IUPU(JP))-Z(IUPU(JP))*COSA(JBU)-SINA(JBU)*DLX(IUPU(JP))*0.5
          END IF
          IF (ELW <= EOFFPU(JP)) PUMPON(JP) = .FALSE.                                                       ! CB 1/13/06
          IF (ELW > EOFFPU(JP) .AND. QPU(JP) > 0.0) THEN
            IF (ELW >= EONPU(JP)) PUMPON(JP) = .TRUE.
            IF (PUMPON(JP)) THEN
              IF (LATERAL_PUMP(JP)) THEN
                JLAT      = 1      
                JWW       = JWW+1
                JBWD(JWW) = JBU
                IWD(JWW)  = IUPU(JP)
                QWD(JWW)  = QPU(JP)
                KTWD(JWW) = KTPU(JP)
                KBWD(JWW) = KBPU(JP)
                EWD(JWW)  = EPU(JP)
                I         = MAX(CUS(JBWD(JWW)),IWD(JWW))
                JB        = JBWD(JWW)
                JW        = JWU
                KT        = KTWB(JW)
                CALL LATERAL_WITHDRAWAL (JWW)
                DO K=KTW(JWW),KBW(JWW)
                  QSS(K,I) = QSS(K,I)-QSW(K,JWW)
                END DO
              ELSE
                JSS(JBU)                 =  JSS(JBU)+1
                KTSW(JSS(JBU),JBU)       =  KTPU(JP)
                KBSW(JSS(JBU),JBU)       =  KBPU(JP)
                JB                       =  JBU
                POINT_SINK(JSS(JBU),JBU) = .TRUE.
                ID                       =  IUPU(JP)
                QSTR(JSS(JBU),JBU)       =  QPU(JP)
                ESTR(JSS(JBU),JBU)       =  EPU(JP)
                KT                       =  KTWB(JWU)
                JW                       =  JWU
                CALL DOWNSTREAM_WITHDRAWAL (JSS(JBU))
                IF (IDPU(JP) /= 0 .AND. US(JBD) == IDPU(JP)) THEN
                  QSUMM = 0.0; TSUM  = 0.0; CSUM  = 0.0
                  DO K=KT,KB(ID)      
                    QSUMM           = QSUMM          +QNEW(K)
                    TSUM            = TSUM           +QNEW(K)*T2(K,ID)
                    CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QNEW(K)*C2(K,ID,CN(1:NAC))
                  END DO     
                  IF (QSUMM /= 0.0) THEN
                    TINSUM(JBD)           = (TSUM           +TINSUM(JBD)          *QINSUM(JBD))/(QSUMM+QINSUM(JBD))
                    CINSUM(CN(1:NAC),JBD) = (CSUM(CN(1:NAC))+CINSUM(CN(1:NAC),JBD)*QINSUM(JBD))/(QSUMM+QINSUM(JBD))
                    QINSUM(JBD)           =  QINSUM(JBD)    +QSUMM  
                  END IF     
                END IF     
                QSUM(JB) = 0.0; TSUM = 0.0; CSUM = 0.0
                DO K=KT,KB(ID)      
                  QSUM(JB)        = QSUM(JB)       +QOUT(K,JB)
                  TSUM            = TSUM           +QOUT(K,JB)*T2(K,ID)
                  CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QOUT(K,JB)*C2(K,ID,CN(1:NAC))
                END DO     
                IF (QSUM(JB) /= 0.0) THEN
                  TOUT(JB)           = TSUM           /QSUM(JB)
                  COUT(CN(1:NAC),JB) = CSUM(CN(1:NAC))/QSUM(JB)
                END IF
              END IF
              IF (IDPU(JP) /= 0) THEN
                IF (US(JBD) /= IDPU(JP) .OR. HEAD_FLOW(JBD) .OR. UP_HEAD(JBD)) THEN
                  JTT              = JTT+1
                  QTR(JTT)         = QPU(JP)
                  ITR(JTT)         = IDPU(JP)
                  PLACE_QTR(JTT)   = PPUC(JP) == ' DENSITY'
                  SPECIFY_QTR(JTT) = PPUC(JP) == ' SPECIFY'
                  IF (SPECIFY_QTR(JTT)) THEN
                    ELTRT(JTT) = ETPU(JP)
                    ELTRB(JTT) = EBPU(JP)
                  END IF
                  JBTR(JTT) = JBD
                  IF (JLAT == 1) THEN
                    TSUM = 0.0; QSUMM = 0.0; CSUM(CN(1:NAC)) = 0.0
                    DO K=KTW(JWW),KBW(JWW)
                      QSUMM           = QSUMM          +QSW(K,JWW)
                      TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
                      CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
                    END DO
                    TTR(JTT)           = TSUM           /QSUMM
                    CTR(CN(1:NAC),JTT) = CSUM(CN(1:NAC))/QSUMM
                  ELSE   
                    TTR(JTT)          = TOUT(JB)
                    CTR(CN(1:NAC),JTT)= COUT(CN(1:NAC),JB)
                  END IF     
                ELSE IF (LATERAL_PUMP(JP)) THEN
                  TSUM      = 0.0; QSUMM = 0.0; CSUM = 0.0
                  ILAT(JWW) = 1
                  DO K=KTW(JWW),KBW(JWW)
                    QSUMM           = QSUMM          +QSW(K,JWW)
                    TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
                    CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
                  END DO
                  TINSUM(JBD)           = (TSUM           +TINSUM(JBD)          *QINSUM(JBD))/(QSUMM+QINSUM(JBD))
                  CINSUM(CN(1:NAC),JBD) = (CSUM(CN(1:NAC))+CINSUM(CN(1:NAC),JBD)*QINSUM(JBD))/(QSUMM+QINSUM(JBD))
                  QINSUM(JBD)           =  QINSUM(JBD)+QSUMM
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
    END IF
    IF (PIPES) THEN
      YSS   = YS
      VSS   = VS
      VSTS  = VST
      YSTS  = YST
      DTPS  = DTP
      QOLDS = QOLD
      CALL PIPE_FLOW (NIT,JDAY)
      DO JP=1,NPI

!****** Positive flows

        JLAT = 0
        JBU  = JBUPI(JP)
        JBD  = JBDPI(JP)
        IF (QPI(JP) >= 0.0) THEN
          IF (LATERAL_PIPE(JP)) THEN
            JLAT      = 1      
            JWW       = JWW+1
            IWD(JWW)  = IUPI(JP)
            QWD(JWW)  = QPI(JP)
            KTWD(JWW) = KTUPI(JP)
            KBWD(JWW) = KBUPI(JP)
            EWD(JWW)  = EUPI(JP)
            JBWD(JWW) = JBU
            I         = MAX(CUS(JBWD(JWW)),IWD(JWW))
            JB        = JBWD(JWW)
            JW        = JWUPI(JP)
            KT        = KTWB(JW)
            CALL LATERAL_WITHDRAWAL (JWW)
            DO K=KTW(JWW),KBW(JWW)
              QSS(K,I) = QSS(K,I)-QSW(K,JWW)
            END DO
          ELSE
            JSS(JBU)                 =  JSS(JBU)+1
            KTSW(JSS(JBU),JBU)       =  KTDPI(JP)
            KBSW(JSS(JBU),JBU)       =  KBDPI(JP)
            JB                       =  JBU
            POINT_SINK(JSS(JBU),JBU) = .TRUE.
            ID                       =  IUPI(JP)
            QSTR(JSS(JBU),JBU)       =  QPI(JP)
            ESTR(JSS(JBU),JBU)       =  EUPI(JP)
            KT                       =  KTWB(JWUPI(JP))
            JW                       =  JWUPI(JP)
            CALL DOWNSTREAM_WITHDRAWAL(JSS(JBU))
            IF (IDPI(JP) /= 0 .AND. US(JBD) == IDPI(JP)) THEN
              QSUMM = 0.0; TSUM  = 0.0; CSUM  = 0.0
              DO K=KT,KB(ID)
                QSUMM           = QSUMM          +QNEW(K)
                TSUM            = TSUM           +QNEW(K)*T2(K,ID)
                CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QNEW(K)*C2(K,ID,CN(1:NAC))
              END DO
              IF (QSUMM /= 0.0) THEN
                TINSUM(JBD)           = (TSUM           +QINSUM(JBD)*TINSUM(JBD))          /(QSUMM+QINSUM(JBD))
                CINSUM(CN(1:NAC),JBD) = (CSUM(CN(1:NAC))+QINSUM(JBD)*CINSUM(CN(1:NAC),JBD))/(QSUMM+QINSUM(JBD))
                QINSUM(JBD)           =  QINSUM(JBD)    +QSUMM  
              END IF 
            END IF     
            QSUM(JB) = 0.0; TSUM = 0.0; CSUM = 0.0
            DO K=KT,KB(ID)
              QSUM(JB)        = QSUM(JB)       +QOUT(K,JB)
              TSUM            = TSUM           +QOUT(K,JB)*T2(K,ID)
              CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QOUT(K,JB)*C2(K,ID,CN(1:NAC))
            END DO
            IF (QSUM(JB) /= 0.0) THEN
              TOUT(JB)           = TSUM           /QSUM(JB)
              COUT(CN(1:NAC),JB) = CSUM(CN(1:NAC))/QSUM(JB)
            END IF   
          END IF
          IF (IDPI(JP) /= 0) THEN
            IF (US(JBD) /= IDPI(JP) .OR. HEAD_FLOW(JBD) .OR. UP_HEAD(JBD)) THEN
              JTT              = JTT+1
              QTR(JTT)         = QPI(JP)
              ITR(JTT)         = IDPI(JP)
              PLACE_QTR(JTT)   = PDPIC(JP) == ' DENSITY'
              SPECIFY_QTR(JTT) = PDPIC(JP) == ' SPECIFY'
              IF (SPECIFY_QTR(JTT)) THEN
                ELTRT(JTT) = ETDPI(JP)
                ELTRB(JTT) = EBDPI(JP)
              END IF
              JBTR(JTT) = JBD
              IF (JLAT == 1) THEN   
                TSUM = 0.0; QSUMM = 0.0; CSUM(CN(1:NAC)) = 0.0
                DO K=KTW(JWW),KBW(JWW)
                  QSUMM           = QSUMM          +QSW(K,JWW)
                  TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
                  CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
                END DO
                TTR(JTT)           = TSUM           /QSUMM
                CTR(CN(1:NAC),JTT) = CSUM(CN(1:NAC))/QSUMM
              ELSE   
                TTR(JTT)           = TOUT(JB)
                CTR(CN(1:NAC),JTT) = COUT(CN(1:NAC),JB)
              END IF     
            ELSE 
              IF (LATERAL_PIPE(JP)) THEN
                TSUM      = 0.0; QSUMM = 0.0; CSUM = 0.0
                ILAT(JWW) = 1
                JB        = JBD
                DO K=KTW(JWW),KBW(JWW)
                  QSUMM           = QSUMM          +QSW(K,JWW)
                  TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
                  CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
                END DO
                TINSUM(JB)           = (TINSUM(JB)          *QINSUM(JB)+TSUM)           /(QSUMM+QINSUM(JB))
                CINSUM(CN(1:NAC),JB) = (CINSUM(CN(1:NAC),JB)*QINSUM(JB)+CSUM(CN(1:NAC)))/(QSUMM+QINSUM(JB))
                QINSUM(JB)           =  QSUMM               +QINSUM(JB)
              END IF
            END IF
          END IF
        ELSE
          JTT              =  JTT+1
          JWW              =  JWW+1
          IWD(JWW)         =  IDPI(JP)
          ITR(JTT)         =  IUPI(JP)
          QTR(JTT)         = -QPI(JP)
          QWD(JWW)         = -QPI(JP)
          KTWD(JWW)        =  KTDPI(JP)
          KBWD(JWW)        =  KBDPI(JP)
          EWD(JWW)         =  EDPI(JP)
          PLACE_QTR(JTT)   =  PUPIC(JP) == ' DENSITY'
          SPECIFY_QTR(JTT) =  PUPIC(JP) == ' SPECIFY'
          IF (SPECIFY_QTR(JTT)) THEN
            ELTRT(JTT) = ETUPI(JP)
            ELTRB(JTT) = EBUPI(JP)
          END IF
          JBTR(JTT) = JBU
          JBWD(JWW) = JBD
          I         = MAX(CUS(JBWD(JWW)),IWD(JWW))
          JB        = JBWD(JWW)
          JW        = JWDPI(JP)
          KT        = KTWB(JW)
          CALL LATERAL_WITHDRAWAL (JWW)
          DO K=KTW(JWW),KBW(JWW)
            QSS(K,I) = QSS(K,I)-QSW(K,JWW)
          END DO
          IF (IDPI(JP) /= 0) THEN
            TSUM  = 0.0; QSUMM = 0.0; CSUM  = 0.0
            DO K=KTW(JWW),KBW(JWW)
              QSUMM           = QSUMM          +QSW(K,JWW)
              TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
              CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
            END DO
            TTR(JTT)           = TSUM           /QSUMM
            CTR(CN(1:NAC),JTT) = CSUM(CN(1:NAC))/QSUMM
          END IF
        END IF
      END DO
    END IF
    IF (GATES) THEN
      CALL GATE_FLOW
      DO JG=1,NGT

!****** Positive flows

        JLAT = 0
        JBU  = JBUGT(JG)
        JBD  = JBDGT(JG)
        IF (QGT(JG) > 0.0) THEN
          IF (LATERAL_GATE(JG)) THEN
            JLAT      = 1
            JWW       = JWW+1
            IWD(JWW)  = IUGT(JG)
            QWD(JWW)  = QGT(JG)
            KTWD(JWW) = KTUGT(JG)
            KBWD(JWW) = KBUGT(JG)
            EWD(JWW)  = EGT(JG)
            JBWD(JWW) = JBU
            I         = MAX(CUS(JBWD(JWW)),IWD(JWW))
            JW        = JWUGT(JG)
            JB        = JBWD(JWW)
            KT        = KTWB(JW)
            CALL LATERAL_WITHDRAWAL (JWW)
            DO K=KTW(JWW),KBW(JWW)
              QSS(K,I) = QSS(K,I)-QSW(K,JWW)
            END DO
            IF (IDGT(JG) /= 0) THEN
              CSUM(CN(1:NAC)) = 0.0; TSUM = 0.0; QSUMM = 0.0
              DO K=KTW(JWW),KBW(JWW)
                QSUMM           = QSUMM          +QSW(K,JWW)
                TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
                CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
              END DO
              TTR(JTT) = TSUM/QSUMM
              DO JC=1,NAC
                CTR(CN(JC),JTT) = CSUM(CN(JC))/QSUMM
                IF (CN(JC) == NDO .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN
                  TDG_GATE(JWW,JG) = .TRUE.
                  CALL TOTAL_DISSOLVED_GAS (PALT(ID),1,JG,TTR(JTT),CTR(CN(JC),JTT))
                END IF
              END DO
            ELSE IF (CAC(NDO) == '      ON' .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN
              TDG_GATE(JWW,JG) = .TRUE.
            END IF
          ELSE
            JSS(JBU)                 =  JSS(JBU)+1
            KTSW(JSS(JBU),JBU)       =  KTUGT(JG)
            KBSW(JSS(JBU),JBU)       =  KBUGT(JG)
            JB                       =  JBU
            POINT_SINK(JSS(JBU),JBU) = .TRUE.
            ID                       =  IUGT(JG)
            ESTR(JSS(JBU),JBU)       =  EGT(JG)
            QSTR(JSS(JBU),JBU)       =  QGT(JG)
            KT                       =  KTWB(JWUGT(JG))
            JW                       =  JWUGT(JG)
            CALL DOWNSTREAM_WITHDRAWAL (JSS(JBU))
            IF (IDGT(JG) /= 0 .AND. US(JBD) == IDGT(JG)) THEN
              QSUMM = 0.0
              TSUM  = 0.0
              CSUM  = 0.0
              DO K=KT,KB(ID)
                QSUMM = QSUMM+QNEW(K)
                TSUM  = TSUM+QNEW(K)*T2(K,ID)
                DO JC=1,NAC
                  IF (CAC(NDO) == '      ON' .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN
                    CALL TOTAL_DISSOLVED_GAS(PALT(ID),1,JG,T2(K,ID),CGAS)
                    CSUM(CN(JC)) = CSUM(CN(JC))+QNEW(K)*CGAS
                  ELSE
                    CSUM(CN(JC)) = CSUM(CN(JC))+QNEW(K)*C2(K,ID,CN(JC))
                  END IF
                END DO
              END DO
              IF (QSUMM /= 0.0) THEN
                TINSUM(JBD)           = (TSUM           +QINSUM(JBD)*TINSUM(JBD))          /(QSUMM+QINSUM(JBD))
                CINSUM(CN(1:NAC),JBD) = (CSUM(CN(1:NAC))+QINSUM(JBD)*CINSUM(CN(1:NAC),JBD))/(QSUMM+QINSUM(JBD))
                QINSUM(JBD)           =  QINSUM(JBD)    +QSUMM
              END IF
            END IF
            QSUM(JB) = 0.0; TSUM = 0.0; CSUM = 0.0
            DO K=KT,KB(ID)
              QSUM(JB) = QSUM(JB)+QOUT(K,JB)
              TSUM     = TSUM+QOUT(K,JB)*T2(K,ID)
              DO JC=1,NAC
                IF (CAC(NDO) == '      ON' .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN
                  CALL TOTAL_DISSOLVED_GAS (PALT(ID),1,JG,T2(K,ID),CGAS)
                  CSUM(CN(JC)) = CSUM(CN(JC))+QOUT(K,JB)*CGAS
                ELSE
                  CSUM(CN(JC)) = CSUM(CN(JC))+QOUT(K,JB)*C2(K,ID,CN(JC))
                END IF
              END DO
            END DO
            IF (QSUM(JB) /= 0.0) THEN
              TOUT(JB)           = TSUM           /QSUM(JB)
              COUT(CN(1:NAC),JB) = CSUM(CN(1:NAC))/QSUM(JB)
            END IF   
          END IF     
          IF (IDGT(JG) /= 0) THEN
            IF (US(JBD) /= IDGT(JG) .OR. HEAD_FLOW(JBD) .OR. UP_HEAD(JBD)) THEN
              JTT              = JTT+1
              QTR(JTT)         = QGT(JG)
              ITR(JTT)         = IDGT(JG)
              PLACE_QTR(JTT)   = PDGTC(JG) == ' DENSITY'
              SPECIFY_QTR(JTT) = PDGTC(JG) == ' SPECIFY'
              IF (SPECIFY_QTR(JTT)) THEN
                ELTRT(JTT) = ETDGT(JG)
                ELTRB(JTT) = EBDGT(JG)
              END IF
              JBTR(JTT) = JBD
              IF (JLAT == 1) THEN
                CSUM(CN(1:NAC)) =  0.0; TSUM = 0.0; QSUMM = 0.0
                DO K=KTW(JWW),KBW(JWW)
                  QSUMM           = QSUMM          +QSW(K,JWW)
                  TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
                  CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
                END DO
                TTR(JTT) = TSUM/QSUMM
                DO JC=1,NAC
                  CTR(CN(JC),JTT) = CSUM(CN(JC))/QSUMM
                  IF (CN(JC) == NDO .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN
                    TDG_GATE(JWW,JG) = .TRUE.
                    CALL TOTAL_DISSOLVED_GAS (PALT(ID),1,JG,TTR(JTT),CTR(CN(JC),JTT))
                  END IF
                END DO
              ELSE   
                TTR(JTT) =  TOUT(JB)
                DO JC=1,NAC
                  CTR(CN(JC),JTT) = COUT(CN(JC),JB)
                  IF (CN(JC) == NDO .AND. GASGTC(JS) == '      ON' .AND. QGT(JG) > 0.0) THEN
                    CALL TOTAL_DISSOLVED_GAS (PALT(ID),0,JS,TTR(JTT),CTR(CN(JC),JTT))
                  END IF
                END DO
              END IF     
            ELSE IF (CAC(NDO) == '      ON' .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN
              TDG_GATE(JWW,JG) = .TRUE.      
            ELSE IF (LATERAL_GATE(JG)) THEN
              TSUM      = 0.0; QSUMM = 0.0; CSUM = 0.0
              ILAT(JWW) = 1
              DO K=KTW(JWW),KBW(JWW)
                QSUMM           = QSUMM          +QSW(K,JWW)
                TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
                CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
              END DO
              JB                   =  JBD
              TINSUM(JB)           = (TINSUM(JB)          *QINSUM(JB)+TSUM)           /(QSUMM+QINSUM(JB))
              CINSUM(CN(1:NAC),JB) = (CINSUM(CN(1:NAC),JB)*QINSUM(JB)+CSUM(CN(1:NAC)))/(QSUMM+QINSUM(JB))
              QINSUM(JB)           =  QSUMM+QINSUM(JB)
            END IF   
          END IF
        ELSE IF (QGT(JG) < 0.0) THEN
          JTT              =  JTT+1
          JWW              =  JWW+1
          IWD(JWW)         =  IDGT(JG)
          ITR(JTT)         =  IUGT(JG)
          QTR(JTT)         = -QGT(JG)
          QWD(JWW)         = -QGT(JG)
          KTWD(JWW)        =  KTDGT(JG)
          KBWD(JWW)        =  KBDGT(JG)
          EWD(JWW)         =  EGT(JG)
          PLACE_QTR(JTT)   =  PUGTC(JG) == ' DENSITY'
          SPECIFY_QTR(JTT) =  PUGTC(JG) == ' SPECIFY'
          IF (SPECIFY_QTR(JTT)) THEN
            ELTRT(JTT) = ETUGT(JG)
            ELTRB(JTT) = EBUGT(JG)
          END IF
          JBTR(JTT) = JBU
          JBWD(JWW) = JBD
          I         = MAX(CUS(JBWD(JWW)),IWD(JWW))
          JW        = JWDGT(JG)
          JB        = JBWD(JWW)
          KT        = KTWB(JW)
          CALL LATERAL_WITHDRAWAL (JWW)
          DO K=KTW(JWW),KBW(JWW)
            QSS(K,I) = QSS(K,I)-QSW(K,JWW)
          END DO
          IF (IDGT(JG) /= 0) THEN
            CSUM(CN(1:NAC)) = 0.0; TSUM = 0.0; QSUMM = 0.0
            DO K=KTW(JWW),KBW(JWW)
              QSUMM           = QSUMM          +QSW(K,JWW)
              TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
              CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
            END DO
            TTR(JTT) = TSUM/QSUMM
            DO JC=1,NAC
              CTR(CN(JC),JTT) = CSUM(CN(JC))/QSUMM
              IF (CN(JC) == NDO .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN
                TDG_GATE(JWW,JG) = .TRUE.
                CALL TOTAL_DISSOLVED_GAS (PALT(ID),1,JG,TTR(JTT),CTR(CN(JC),JTT))
              END IF
            END DO
          ELSE IF (CAC(NDO) == '      ON' .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN
            TDG_GATE(JWW,JG) = .TRUE.
          END IF
        END IF
      END DO
    END IF
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IU = CUS(JB)
        ID = DS(JB)
        IF (EVAPORATION(JW)) THEN
          EVBR(JB) = 0.0
          DO I=IU,ID
            FW = AFW(JW)+BFW(JW)*WIND2(I)**CFW(JW)
            IF (RH_EVAP(JW)) THEN
              EA = EXP(2.3026*(7.5*TDEW(JW)/(TDEW(JW)+237.3)+0.6609))
              ES = EXP(2.3026*(7.5*T2(KT,I)/(T2(KT,I)+237.3)+0.6609))
              IF (TDEW(JW) < 0.0) EA = EXP(2.3026*(9.5*TDEW(JW)/(TDEW(JW)+265.5)+0.6609))
              IF (T2(KT,I) < 0.0) ES = EXP(2.3026*(9.5*T2(KT,I)/(T2(KT,I)+265.5)+0.6609))
              TAIRV = (TAIR(JW)+273.0)/(1.0-0.378*EA/760.0)
              DTV   = (T2(KT,I)+273.0)/(1.0-0.378*ES/760.0)-TAIRV
              DTVL  =  0.0084*WIND2(I)**3      
              IF (DTV < DTVL) DTV = DTVL
              FW = (3.59*DTV**0.3333333+4.26*WIND2(I))
            END IF
            TM    = (T2(KT,I)+TDEW(JW))*0.5
            VPTG  =  0.35+0.015*TM+0.0012*TM*TM
            EV(I) =  VPTG*(T2(KT,I)-TDEW(JW))*FW*BI(KT,I)*DLX(I)/2.45E9
            IF (EV(I) < 0.0 .OR. ICE(I)) EV(I) = 0.0
            QSS(KT,I) = QSS(KT,I)-EV(I)
            EVBR(JB)  = EVBR(JB)+EV(I)
          END DO
        END IF
        IF (PRECIPITATION(JW)) THEN
          QPRBR(JB) = 0.0
          DO I=IU,ID
            QPR(I)    = PR(JB)*BI(KT,I)*DLX(I)
            QPRBR(JB) = QPRBR(JB)+QPR(I)
            QSS(KT,I) = QSS(KT,I)+QPR(I)
          END DO
        END IF
        IF (TRIBUTARIES) THEN
          DO JT=1,JTT

!********** Inflow fractions

            IF (JB == JBTR(JT)) THEN
              I = MAX(ITR(JT),IU)
              QTRF(KT:KB(I),JT) = 0.0
              IF (PLACE_QTR(JT)) THEN

!************** Inflow layer

                SSTOT = 0.0
                DO J=NSSS,NSSE
                  SSTOT = SSTOT+CTR(J,JT)
                END DO
                RHOTR = DENSITY(TTR(JT),CTR(NTDS,JT),SSTOT)
                K     = KT
                DO WHILE (RHOTR > RHO(K,I) .AND. K < KB(I))
                  K = K+1
                END DO
                KTTR(JT) = K
                KBTR(JT) = K

!************** Layer inflows

                VQTR  =  QTR(JT)*DLT
                VQTRI =  VQTR
                QTRFR =  1.0
                INCR  = -1
                DO WHILE (QTRFR > 0.0)
                  IF (K <= KB(I)) THEN
                    V1 = VOL(K,I)
                    IF (VQTR > 0.5*V1) THEN
                      QTRF(K,JT) = 0.5*V1/VQTRI
                      QTRFR      = QTRFR-QTRF(K,JT)
                      VQTR       = VQTR-QTRF(K,JT)*VQTRI
                      IF (K == KT) THEN
                        K    = KBTR(JT)
                        INCR = 1
                      END IF
                    ELSE
                      QTRF(K,JT) = QTRFR
                      QTRFR      = 0.0
                    END IF
                    IF (INCR < 0) KTTR(JT) = K
                    IF (INCR > 0) KBTR(JT) = MIN(KB(I),K)
                    K = K+INCR
                  ELSE
                    QTRF(KT,JT) = QTRF(KT,JT)+QTRFR
                    QTRFR       = 0.0
                  END IF
                END DO
              ELSE
                IF (SPECIFY_QTR(JT)) THEN
                  KTTR(JT) = 2
                  DO WHILE (EL(KTTR(JT),I) > ELTRT(JT))
                    KTTR(JT) = KTTR(JT)+1
                  END DO
                  KBTR(JT) = KMX-1
                  DO WHILE (EL(KBTR(JT),I) < ELTRB(JT))
                    KBTR(JT) = KBTR(JT)-1
                  END DO
                ELSE
                  KTTR(JT) = KT
                  KBTR(JT) = KB(I)
                END IF
                KTTR(JT) = MAX(KT,KTTR(JT))
                KBTR(JT) = MIN(KB(I),KBTR(JT))
                IF (KBTR(JT) < KTTR(JT)) KBTR(JT) = KTTR(JT)
                BHSUM = 0.0
                DO K=KTTR(JT),KBTR(JT)
                  BHSUM = BHSUM+BH2(K,I)
                END DO
                DO K=KTTR(JT),KBTR(JT)
                  QTRF(K,JT) = BH2(K,I)/BHSUM
                END DO
              END IF
              DO K=KTTR(JT),KBTR(JT)
                QSS(K,I) = QSS(K,I)+QTR(JT)*QTRF(K,JT)
              END DO
            END IF
          END DO
        END IF
        IF (DIST_TRIBS(JB)) THEN
          AKBR = 0.0
          DO I=IU,ID
            AKBR = AKBR+BI(KT,I)*DLX(I)
          END DO
          DO I=IU,ID
            QDT(I)    = QDTR(JB)*BI(KT,I)*DLX(I)/AKBR
            QSS(KT,I) = QSS(KT,I)+QDT(I)
          END DO
        END IF
        IF (WITHDRAWALS) THEN
          DO JWD=1,NWD
            IF (JB == JBWD(JWD)) THEN
              I = MAX(CUS(JBWD(JWD)),IWD(JWD))
              CALL LATERAL_WITHDRAWAL (JWD)
              DO K=KTW(JWD),KBW(JWD)
                QSS(K,I) = QSS(K,I)-QSW(K,JWD)
              END DO
            END IF
          END DO
        END IF
        IF (UH_INTERNAL(JB)) THEN
          IF (UHS(JB) /= DS(JBUH(JB)) .OR. DHS(JBUH(JB)) /= US(JB)) THEN
            IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
              DO K=KT,KB(IU-1)
                QSS(K,UHS(JB)) = QSS(K,UHS(JB))-VOLUH2(K,JB)/DLT
              END DO
            ELSE
              CALL UPSTREAM_FLOW
            END IF
          END IF
        END IF
        IF (DH_INTERNAL(JB)) THEN
          IF (DHS(JB) /= US(JBDH(JB)) .OR. UHS(JBDH(JB)) /= DS(JB)) THEN
            IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
              DO K=KT,KB(ID+1)
                QSS(K,CDHS(JB)) = QSS(K,CDHS(JB))+VOLDH2(K,JB)/DLT
              END DO
            ELSE
              CALL DOWNSTREAM_FLOW
            END IF
          END IF
        END IF
      END DO
    END DO

!** Compute tributary contribution to cross-shear

    IF (TRIBUTARIES) THEN
      DO JW=1,NWB
        DO JB=BS(JW),BE(JW)
          DO JT=1,JTT
            IF (JB == JBTR(JT)) THEN
              I = MAX(CUS(JB),ITR(JT))
              DO K=KTWB(JW),KBMIN(I)
                UYBR(K,I) = UYBR(K,I)+ABS(QTR(JT))*QTRF(K,JT)
              END DO
            END IF
          END DO
        END DO
      END DO
    END IF

!***********************************************************************************************************************************
!**                                           Task 2.2: Hydrodynamic calculations                                                 **
!***********************************************************************************************************************************

    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IU = CUS(JB)
        ID = DS(JB)

!***********************************************************************************************************************************
!**                                Task 2.2.1: Boundary concentrations, temperatures, and densities                               **
!***********************************************************************************************************************************

        IUT = IU
        IDT = ID
        IF (UP_FLOW(JB)) THEN
          IF (.NOT. INTERNAL_FLOW(JB)) THEN
            DO K=KT,KB(IU)
                IF (QIND(JB)+QINSUM(JB).GT.0.0) THEN
                  TIN(JB)               = (TINSUM(JB)               *QINSUM(JB)+TIND(JB)          *QIND(JB))/(QIND(JB)+QINSUM(JB))
                  CIN(CN(1:NAC),JB)     =  MAX((CINSUM(CN(1:NAC),JB)*QINSUM(JB)+CIND(CN(1:NAC),JB)*QIND(JB))/(QIND(JB)+QINSUM(JB)),&
                                                0.0)
                  T1(K,IU-1)            =  TIN(JB)
                  T2(K,IU-1)            =  TIN(JB)
                  C1S(K,IU-1,CN(1:NAC)) =  CIN(CN(1:NAC),JB)
                  QIN(JB)               =  QIND(JB)+QINSUM(JB)
                ELSE      
                  QIN(JB)               =  0.0
                  TIN(JB)               =  TIND(JB)
                  T1(K,IU-1)            =  TIND(JB)
                  T2(K,IU-1)            =  TIND(JB)
                  C1S(K,IU-1,CN(1:NAC)) =  CIND(CN(1:NAC),JB)
                END IF
            END DO
          ELSE IF (.NOT. DAM_INFLOW(JB)) THEN                                                                          !TC 08/03/04
            IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
              TIN(JB)           = T1(KT,UHS(JB))
              CIN(CN(1:NAC),JB) = MAX(C1S(KT,UHS(JB),CN(1:NAC)),0.0)
              DO K=KT,KB(IU)
                T1(K,IU-1)            = T1(K,UHS(JB))
                T2(K,IU-1)            = T1(K,UHS(JB))
                C1S(K,IU-1,CN(1:NAC)) = C1S(K,UHS(JB),CN(1:NAC))
                C1(K,IU-1,CN(1:NAC))  = C1S(K,UHS(JB),CN(1:NAC))
                C2(K,IU-1,CN(1:NAC))  = C1S(K,UHS(JB),CN(1:NAC))
              END DO
            ELSE
              CALL UPSTREAM_WATERBODY
              TIN(JB)           = T1(KT,IU-1)
              CIN(CN(1:NAC),JB) = MAX(C1(KT,IU-1,CN(1:NAC)),0.0)
            END IF
          ELSE
            TIN(JB)           = TINSUM(JB)
            QIN(JB)           = QINSUM(JB)
            CIN(CN(1:NAC),JB) = MAX(CINSUM(CN(1:NAC),JB),0.0)
            DO K=KT,KB(ID)
              T1(K,IU-1)            = TIN(JB)
              T2(K,IU-1)            = TIN(JB)
              C1S(K,IU-1,CN(1:NAC)) = CIN(CN(1:NAC),JB)
            END DO
           END IF
        END IF
        IF (DN_FLOW(JB)) THEN
          DO K=KT,KB(ID)
              T1(K,ID+1)            = T2(K,ID)
              T2(K,ID+1)            = T2(K,ID)
              C1S(K,ID+1,CN(1:NAC)) = C1S(K,ID,CN(1:NAC))
          END DO
        END IF
        IF (UP_HEAD(JB)) THEN
          IUT = IU-1
          IF (UH_INTERNAL(JB)) THEN
            IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
              DO K=KT,KB(IUT)
                RHO(K,IUT)           = RHO(K,UHS(JB))
                T1(K,IUT)            = T2(K,UHS(JB))
                T2(K,IUT)            = T2(K,UHS(JB))
                C1S(K,IUT,CN(1:NAC)) = C1S(K,UHS(JB),CN(1:NAC))
                C1(K,IUT,CN(1:NAC))  = C1S(K,UHS(JB),CN(1:NAC))
                C2(K,IUT,CN(1:NAC))  = C1S(K,UHS(JB),CN(1:NAC))
              END DO
            ELSE
              CALL UPSTREAM_WATERBODY
            END IF
            DO K=KT,KB(IUT)      
              RHO(K,IUT) = DENSITY(T2(K,IUT),MAX(TDS(K,IUT),0.0),MAX(TISS(K,IUT),0.0))
            END DO
          ELSE IF (UH_EXTERNAL(JB)) THEN
            DO K=KT,KB(IUT)
              RHO(K,IUT)           = DENSITY(TUH(K,JB),MAX(TDS(K,IUT),0.0),MAX(TISS(K,IUT),0.0))
              T1(K,IUT)            = TUH(K,JB)
              T2(K,IUT)            = TUH(K,JB)
              C1S(K,IUT,CN(1:NAC)) = CUH(K,CN(1:NAC),JB)
              C1(K,IUT,CN(1:NAC))  = CUH(K,CN(1:NAC),JB)
              C2(K,IUT,CN(1:NAC))  = CUH(K,CN(1:NAC),JB)
            END DO
          END IF
        END IF
        IF (DN_HEAD(JB)) THEN
          IDT = ID+1
          IF (DH_INTERNAL(JB)) THEN
            IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
              DO K=KT,KB(IDT)
                RHO(K,IDT)           = RHO(K,DHS(JB))
                T1(K,IDT)            = T2(K,DHS(JB))
                T2(K,IDT)            = T2(K,DHS(JB))
                C1S(K,IDT,CN(1:NAC)) = C1S(K,DHS(JB),CN(1:NAC))
                C1(K,IDT,CN(1:NAC))  = C1S(K,DHS(JB),CN(1:NAC))
                C2(K,IDT,CN(1:NAC))  = C1S(K,DHS(JB),CN(1:NAC))
              END DO
            ELSE
              CALL DOWNSTREAM_WATERBODY
            END IF
            DO K=KT,KB(ID)      
              RHO(K,IDT) = DENSITY(T2(K,IDT),MAX(TDS(K,IDT),0.0),MAX(TISS(K,IDT),0.0))
            END DO
          ELSE IF (DH_EXTERNAL(JB)) THEN
            DO K=KT,KB(IDT)
              RHO(K,IDT)           = DENSITY(TDH(K,JB),MAX(TDS(K,IDT),0.0),MAX(TISS(K,IDT),0.0))
              T1(K,IDT)            = TDH(K,JB)
              T2(K,IDT)            = TDH(K,JB)
              C1S(K,IDT,CN(1:NAC)) = CDH(K,CN(1:NAC),JB)
              C1(K,IDT,CN(1:NAC))  = CDH(K,CN(1:NAC),JB)
              C2(K,IDT,CN(1:NAC))  = CDH(K,CN(1:NAC),JB)
            END DO
          END IF
        END IF

!***********************************************************************************************************************************
!**                                                 Task 2.2.2: Momentum terms                                                    **
!***********************************************************************************************************************************

!****** Density pressures

        DO I=IUT,IDT
          DO K=KT,KB(I)
            P(K,I) = P(K-1,I)+RHO(K,I)*G*H(K,JW)*COSA(JB)
          END DO
        END DO

!****** Horizontal density gradients

        DO I=IUT,IDT-1
          HDG(KT,I) = DLXRHO(I)*(BKT(I)+BKT(I+1))*0.5*H(KT,JW)*(P(KT,I+1)-P(KT,I))
          DO K=KT+1,KBMIN(I)
            HDG(K,I) = DLXRHO(I)*BHR2(K,I)*((P(K-1,I+1)-P(K-1,I))+(P(K,I+1)-P(K,I)))
          END DO
        END DO

!****** Adjusted wind speed and surface wind shear drag coefficient
        
        DO I=IU-1,ID+1
          WIND10(I) = WIND(JW)*WSC(I)*LOG(10.0/0.01)/LOG(WINDH(JW)/0.01)
          FETCH(I)  = FETCHD(I,JB)
          IF (COS(PHI(JW)-PHI0(I)) < 0.0) FETCH(I) = FETCHU(I,JB)
          IF (FETCH(I) <= 0.0) FETCH(I) = DLX(I)
          IF (FETCH_CALC(JW)) THEN
            ZB        = 0.8*LOG(FETCH(I)*0.5)-1.0718
            WIND10(I) = WIND10(I)*(5.0*ZB+4.6052)/(3.0*ZB+9.2103)
          END IF
          CZ(I) = 0.0
          IF (WIND10(I) >= 1.0)  CZ(I) = 0.0005*SQRT(WIND10(I))
          IF (WIND10(I) >= 15.0) CZ(I) = 0.0026
        END DO     

!****** Longitudinal and lateral surface wind shear and exponential decay

        DO I=IUT,IDT-1
          WSHX(I) = CZ(I)*WIND10(I)**2*RHOA/RHOW*    COS(PHI(JW)-PHI0(I))* ICESW(I)
          WSHY(I) = CZ(I)*WIND10(I)**2*RHOA/RHOW*ABS(SIN(PHI(JW)-PHI0(I)))*ICESW(I)
          WWT     = 0.0
          IF (WIND10(I) /= 0.0) WWT = 6.95E-2*(FETCH(I)**0.233)*WIND10(I)**0.534
          DFC = -8.0*PI*PI/(G*WWT*WWT+NONZERO)
          DO K=KT,KBMIN(I)
            DECAY(K,I) = EXP(MAX(DFC*DEPTHB(K,I),-30.0))
          END DO

!******** Branch inflow lateral shear and friction

          DO JJB=1,NBR
            IF (I == UHS(JJB) .AND. .NOT. INTERNAL_FLOW(JJB)) THEN
              BETABR = (PHI0(I)-PHI0(US(JJB)))
              IF (JJB >= BS(JW) .AND. JJB <= BE(JW)) THEN
                DO K=KT,KBMIN(I)
                  IF (U(K,US(JJB)) < 0.0) THEN
                    UXBR(K,I) = UXBR(K,I)+ABS(U(K,US(JJB)))*COS(BETABR)     *VOLUH2(K,JJB)/(DLT*DLX(I))
                    UYBR(K,I) = UYBR(K,I)              +ABS(SIN(BETABR))*ABS(VOLUH2(K,JJB))/DLT
                  END IF
                END DO
              ELSE
                CALL UPSTREAM_BRANCH
              END IF
            END IF
            IF (I == DHS(JJB)) THEN
              BETABR = (PHI0(I)-PHI0(DS(JJB)))
              IF (I == US(JB) .AND. UHS(JB) /= DS(JJB)) THEN
                IF (JJB >= BS(JW) .AND. JJB <= BE(JW)) THEN
                  DO K=KT,KBMIN(I)
                    IF (U(K,DS(JJB)) >= 0.0) THEN
                      UXBR(K,I) = UXBR(K,I)+U(K,DS(JJB))*   COS(BETABR) *VOLDH2(K,JJB)/(DLT*DLX(I))
                      UYBR(K,I) = UYBR(K,I)            +ABS(SIN(BETABR))*VOLDH2(K,JJB)/DLT
                    END IF
                  END DO
                ELSE
                  CALL DOWNSTREAM_BRANCH
                END IF
              ELSE IF (I /= US(JB)) THEN
                IF (JJB >= BS(JW) .AND. JJB <= BE(JW)) THEN
                  DO K=KT,KBMIN(I)
                    IF (U(K,DS(JJB)) >= 0.0) THEN
                      UXBR(K,I) = UXBR(K,I)+U(K,DS(JJB))*   COS(BETABR) *VOLDH2(K,JJB)/(DLT*DLX(I))     
                      UYBR(K,I) = UYBR(K,I)            +ABS(SIN(BETABR))*VOLDH2(K,JJB)/DLT
                    END IF
                  END DO
                ELSE
                  CALL DOWNSTREAM_BRANCH
                END IF
              END IF
            END IF
          END DO
          DO K=KT,KBMIN(I)
            FRICBR(K,I) = (FI(JW)/8.0)*RHO(K,I)*(UYBR(K,I)/(DLX(I)*H2(K,I)))**2
          END DO
        END DO

!****** Vertical eddy viscosities/diffusivities

        DO I=IUT,IDT-1
          CALL CALCULATE_AZ
          IF (KBMIN(I) <= KT+1 .AND. KB(I) > KBMIN(I)) THEN
            AZ(KBMIN(I),I) = AZMIN
            DZ(KBMIN(I),I) = DZMIN
          END IF
        END DO
        IF (AZC(JW) == '     TKE') THEN
          DO I=IUT,IDT-1
            DO K=KT,KBMIN(I)
              AZT(K,I)  = 0.5*(AZ(K,I)+AZ(K,I+1))    ! SW 10/18/07
            END DO
          AZ(KBMIN(I),I) = AZMIN              !SG 10/4/07 SW 10/4/07
          END DO
          AZ(KT:KMX-1,IUT:IDT-1)=AZT(KT:KMX-1,IUT:IDT-1)                  ! SW 10/18/07
        END IF
        DO JWR=1,NIW
        IF (WEIR_CALC) AZ(KTWR(JWR)-1:KBWR(JWR),IWR(1:NIW)) = 0.0
        END DO

!****** Average eddy diffusivities

        IF(AZC(JW) == '     TKE')THEN 
        DZ(KT:KB(IDT)-1,IDT) = DZ(KT:KB(IDT)-1,IDT-1)
        ELSE
        DZ(KT:KB(IDT)-1,IDT) = DZT(KT:KB(IDT)-1,IDT-1)
        ENDIF
        DO I=IUT,IDT-1
          DO K=KT,KB(I)-1
            IF (K >= KBMIN(I)) THEN
              IF (KB(I-1) >= KB(I) .AND. I /= IUT) THEN
                DZ(K,I) = DZ(K,I-1)
              ELSE
                DZ(K,I) = DZMIN
              END IF
            ELSEIF(AZC(JW) /= '     TKE')THEN              ! SW 10/20/07 Average for non-TKE algorithm, for TKE DZ is at seg center-vertical average
              if(i == iut)then                             ! SW 10/20/07
              dz(k,i)=dzt(k,i)
              else
              DZ(K,I) = (DZT(K,I)+DZT(K,I-1))*0.5          ! SW 10/20/07  (DZT(K,I)+DZT(K+1,I))*0.5
              endif
            END IF
          END DO
        END DO

!****** Density inversions

        DO I=IUT,IDT
          DO K=KT,KB(I)-1
            DZQ(K,I) = MIN(1.0E-4,DZ(K,I))
            IF (RHO(K,I) > RHO(K+1,I)) DZ(K,I) = DZMAX
          END DO
        END DO

!****** Wind, velocity, and bottom shear stresses @ top and bottom of cell

        SB(:,IUT:IDT-1) = 0.0
        DO I=IUT,IDT-1
          ST(KT,I) = WSHX(I)*BR(KTI(I),I)
          DO K=KT+1,KBMIN(I)
            ST(K,I) = WSHX(I)*DECAY(K-1,I)*BR(K,I)
            IF (.NOT. IMPLICIT_VISC(JW)) ST(K,I) = ST(K,I)+AZ(K-1,I)*(BR(K-1,I)+BR(K,I))*0.5*(U(K-1,I)-U(K,I))/((AVH2(K-1,I)       &
                                                   +AVH2(K-1,I+1))*0.5)
          END DO
          GC2 = 0.0
          IF (FRIC(I) /= 0.0) GC2 = G/(FRIC(I)*FRIC(I))
! v3.5 start
          hrad=BHR2(KT,I)/(BR(KTI(I),I)-BR(KT+1,I)+2.0*AVHR(KT,I))
          if(macrophyte_on.and.mannings_n(jw))then
            call macrophyte_friction(hrad,fric(i),effric,kt,i)
            gc2=g*effric*effric/hrad**0.33333333
          else if(.not.macrophyte_on.and.mannings_n(jw))then
            gc2=g*fric(i)*fric(i)/hrad**0.33333333
          end if
          IF (ONE_LAYER(I)) THEN
            SB(KT,I) = ST(KT+1,I)+GC2*(BR(KTI(I),I)+2.0*AVHR(KT,I))*U(KT,I)*ABS(U(KT,I))
          ELSE
            SB(KT,I) = GC2*(BR(KTI(I),I)-BR(KT+1,I)+2.0*AVHR(KT,I))*U(KT,I)*ABS(U(KT,I))
            DO K=KT+1,KBMIN(I)-1
              hrad=(BHR2(K,I)/(BR(K,I)-BR(K+1,I)+2.0*H(K,JW)))
              if(macrophyte_on.and.mannings_n(jw))then
                call macrophyte_friction(hrad,fric(i),effric,k,i)
                gc2=g*effric*effric/hrad**0.33333333
              else if(.not.macrophyte_on.and.mannings_n(jw))then
                gc2=g*fric(i)*fric(i)/hrad**0.33333333
              end if
              SB(K,I) = GC2*(BR(K,I)-BR(K+1,I)+2.0*H(K,JW))*U(K,I)*ABS(U(K,I))
            END DO
            IF (KT /= KBMIN(I)) THEN
              hrad=(BHR2(KBMIN(I),I)/(BR(KBMIN(I),I)+2.0*H(KBMIN(I),JW)))
              if(macrophyte_on.and.mannings_n(jw))then
                call macrophyte_friction(hrad,fric(i),effric,kbmin(i),i)
                gc2=g*effric*effric/hrad**0.33333333
              else if(.not.macrophyte_on.and.mannings_n(jw))then
                gc2=g*fric(i)*fric(i)/hrad**0.33333333
              end if
! v3.5 end
              IF (KBMIN(I) /= KB(I)) THEN
                SB(KBMIN(I),I) = GC2*(BR(KBMIN(I),I)-BR(KBMIN(I)+1,I)+2.0*H2(K,I))*U(KBMIN(I),I)*ABS(U(KBMIN(I),I))
              ELSE
                SB(KBMIN(I),I) = GC2*(BR(KBMIN(I),I)+2.0*H2(K,I))*U(KBMIN(I),I)*ABS(U(KBMIN(I),I))
              END IF
            END IF
          END IF
          DO K=KT,KBMIN(I)-1
            SB(K,I) = SB(K,I)+ST(K+1,I)
          END DO
          SB(KBMIN(I),I) = SB(KBMIN(I),I)+WSHX(I)*DECAY(KBMIN(I),I)*(BR(KBMIN(I)-1,I)+BR(KBMIN(I),I))*0.5
        END DO

!****** Horizontal advection of momentum

        DO I=IU,ID-1
          DO K=KT,KBMIN(I)
            UDR       = (1.0+SIGN(1.0,(U(K,I)+U(K,I+1))*0.5))*0.5
            UDL       = (1.0+SIGN(1.0,(U(K,I)+U(K,I-1))*0.5))*0.5
            ADMX(K,I) = (BH2(K,I+1)*(U(K,I+1)+U(K,I))*0.5*(UDR*U(K,I)+(1.0-UDR)*U(K,I+1))-BH2(K,I)*(U(K,I)+U(K,I-1))               &
                        *0.5*(UDL*U(K,I-1)+(1.0-UDL)*U(K,I)))/DLXR(I)
          END DO
        END DO

!****** Horizontal dispersion of momentum

        DO I=IU,ID-1
          DO K=KT,KBMIN(I)
            DM(K,I) = AX(JW)*(BH2(K,I+1)*(U(K,I+1)-U(K,I))/DLX(I+1)-BH2(K,I)*(U(K,I)-U(K,I-1))/DLX(I))/DLXR(I)
          END DO
        END DO

!****** Vertical advection of momentum

        DO I=IU,ID-1
          DO K=KT,KB(I)-1
            AB        = (1.0+SIGN(1.0,(W(K,I+1)+W(K,I))*0.5))*0.5
            ADMZ(K,I) = (BR(K,I)+BR(K+1,I))*0.5*(W(K,I+1)+W(K,I))*0.5*(AB*U(K,I)+(1.0-AB)*U(K+1,I))
          END DO
        END DO

!****** Gravity force due to channel slope

        DO I=IU-1,ID
          GRAV(KT,I) = AVHR(KT,I)*(BKT(I)+BKT(I+1))*0.5*G*SINA(JB)                                                     !SW 09/09/04
          DO K=KT+1,KB(I)                                                                                              !SW 09/09/04
            GRAV(K,I) = BHR2(K,I)*G*SINA(JB)
          END DO
        END DO

!***********************************************************************************************************************************
!**                                            Task 2.2.3: Water surface elevation                                                **
!***********************************************************************************************************************************

!****** Tridiagonal coefficients

        BHRHO(IU-1:ID+1) = 0.0; D(IU-1:ID+1) = 0.0; F(IU-1:ID+1) = 0.0
        DO I=IU,ID-1
          DO K=KT,KBMIN(I)
            BHRHO(I) = BHRHO(I)+(BH2(K,I+1)/RHO(K,I+1)+BH2(K,I)/RHO(K,I))
          END DO
          DO K=KT,KB(I)
            D(I) = D(I)+(U(K,I)*BHR2(K,I)-U(K,I-1)*BHR2(K,I-1)-QSS(K,I)+(UXBR(K,I)-UXBR(K,I-1))*DLT)
            F(I) = F(I)+(-SB(K,I)+ST(K,I)-ADMX(K,I)+DM(K,I)-HDG(K,I)+GRAV(K,I))
          END DO
        END DO

!****** Boundary tridiagonal coefficients

        D(IU) = 0.0
        DO K=KT,KB(IU)
          D(IU) = D(IU)+(U(K,IU)*BHR2(K,IU)-QSS(K,IU))+UXBR(K,IU)*DLT
        END DO
        IF (DN_FLOW(JB)) THEN
          DO K=KT,KB(ID)
            D(ID) = D(ID)-U(K,ID-1)*BHR2(K,ID-1)-QSS(K,ID)+(UXBR(K,ID)-UXBR(K,ID-1))*DLT+QOUT(K,JB)
          END DO
        END IF
        IF (UP_HEAD(JB)) THEN
          DO K=KT,KBMIN(IU-1)
            BHRHO(IU-1) = BHRHO(IU-1)+(BH2(K,IU)/RHO(K,IU)+BH2(K,IU-1)/RHO(K,IU-1))
          END DO
          DO K=KT,KB(IU)
            D(IU)   = D(IU)-U(K,IU-1)*BHR2(K,IU-1)
            F(IU-1) = F(IU-1)-(SB(K,IU-1)-ST(K,IU-1)+HDG(K,IU-1)-GRAV(K,IU-1))
          END DO
        END IF
        IF (DN_HEAD(JB)) THEN
          DO K=KT,KBMIN(ID)
            BHRHO(ID) = BHRHO(ID)+(BH2(K,ID+1)/RHO(K,ID+1)+BH2(K,ID)/RHO(K,ID))
          END DO
          DO K=KT,KB(ID)
            D(ID) = D(ID)+(U(K,ID)*BHR2(K,ID)-U(K,ID-1)*BHR2(K,ID-1)-QSS(K,ID))+(UXBR(K,ID)-UXBR(K,ID-1))*DLT
            F(ID) = F(ID)+(-SB(K,ID)+ST(K,ID)-HDG(K,ID)+GRAV(K,ID))
          END DO
        END IF
      END DO
    END DO
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IU = CUS(JB)
        ID = DS(JB)
        IF (INTERNAL_FLOW(JB) .AND. .NOT. DAM_INFLOW(JB)) THEN                                                         !TC 08/03/04
          QIN(JB) = 0.0  
          DO K=KTWB(JWUH(JB)),KB(UHS(JB))
            QIN(JB) = QIN(JB)+U(K,UHS(JB))*BHR2(K,UHS(JB))
          END DO
        END IF
        IF (UP_FLOW(JB)) D(IU) = D(IU)-QIN(JB)

!****** Boundary surface elevations

        IF (UH_INTERNAL(JB)) THEN
          Z(IU-1)    = ((-EL(KTWB(JWUH(JB)),UHS(JB))+Z(UHS(JB))*COSA(JBUH(JB)))+EL(KT,IU-1)+SINA(JB)*DLXR(IU-1))/COSA(JB)
          ELWS(IU-1) = EL(KT,IU-1)-Z(IU-1)*COSA(JB)
          KTI(IU-1)  = 2
          DO WHILE (EL(KTI(IU-1),IU-1) > ELWS(IU-1))
            KTI(IU-1) = KTI(IU-1)+1
          END DO
          KTI(IU-1) = MAX(KTI(IU-1)-1,2)
        END IF
        IF (UH_EXTERNAL(JB)) Z(IU-1) = (EL(KT,IU-1)-(ELUH(JB)+SINA(JB)*DLX(IU)*0.5))/COSA(JB)
        IF (DH_INTERNAL(JB)) THEN
          Z(ID+1)    = ((-EL(KTWB(JWDH(JB)),DHS(JB))+Z(DHS(JB))*COSA(JBDH(JB)))+EL(KT,ID+1))/COSA(JB)
          ELWS(ID+1) = EL(KT,ID+1)-Z(ID+1)*COSA(JB)
          KTI(ID+1)  = 2      
          DO WHILE (EL(KTI(ID+1),ID+1) > ELWS(ID+1))
            KTI(ID+1) = KTI(ID+1)+1
          END DO     
          KTI(ID+1) = MAX(KTI(ID+1)-1,2)
          IF (KTI(ID+1) >= KB(ID)) THEN
            Z(ID+1)    = Z(ID)-SLOPE(JB)*DLX(ID)/2.0
            ELWS(ID+1) = EL(KT,ID+1)-Z(ID+1)*COSA(JB)
            KTI(ID+1)  = 2
            DO WHILE (EL(KTI(ID+1),ID+1) > ELWS(ID+1))
              KTI(ID+1) = KTI(ID+1)+1
            END DO
            KTI(ID+1) = MAX(KTI(ID+1)-1,2)
          END IF
        END IF
        IF (DH_EXTERNAL(JB)) Z(ID+1) = (EL(KT,ID+1)-(ELDH(JB)-SINA(JB)*DLX(ID)*0.5))/COSA(JB)

!****** Implicit water surface elevation solution

        DO I=IU,ID
          A(I) = -RHO(KT,I-1)*G*COSA(JB)*DLT*DLT* BHRHO(I-1)*0.5/DLXR(I-1)
          C(I) = -RHO(KT,I+1)*G*COSA(JB)*DLT*DLT* BHRHO(I)  *0.5/DLXR(I)
          V(I) =  RHO(KT,I)  *G*COSA(JB)*DLT*DLT*(BHRHO(I)  *0.5/DLXR(I)+BHRHO(I-1)*0.5/DLXR(I-1))+DLX(I)*BI(KT,I)
          D(I) =  DLT*(D(I)+DLT*(F(I)-F(I-1)))+DLX(I)*BI(KT,I)*Z(I)
        END DO
        IF (UP_HEAD(JB)) D(IU) = D(IU)-A(IU)*Z(IU-1)
        IF (DN_HEAD(JB)) D(ID) = D(ID)-C(ID)*Z(ID+1)
        BTA(IU) = V(IU)
        GMA(IU) = D(IU)
        DO I=IU+1,ID
          BTA(I) = V(I)-A(I)/BTA(I-1)*C(I-1)
          GMA(I) = D(I)-A(I)/BTA(I-1)*GMA(I-1)
        END DO
        Z(ID) = GMA(ID)/BTA(ID)
        DO I=ID-1,IU,-1
          Z(I) = (GMA(I)-C(I)*Z(I+1))/BTA(I)
        END DO

!****** Boundary water surface elevations

        IF (UP_FLOW(JB) .AND. .NOT. HEAD_FLOW(JB)) Z(IU-1) = Z(IU)      
        IF (UP_FLOW(JB) .AND.       HEAD_FLOW(JB)) Z(IU-1) = (-EL(KTWB(JWUH(JB)),UHS(JB))+Z(UHS(JB))*COSA(JBUH(JB))+EL(KT,IU-1)  &
                                                             +SINA(JBUH(JB))*DLXR(IU-1))/COSA(JBUH(JB))
        IF (DN_FLOW(JB))                           Z(ID+1) = Z(ID)

!****** Updated surface layer and geometry

        IF (.NOT. TRAPEZOIDAL(JW)) THEN                                                                                !SW 07/16/04
          DO I=IU-1,ID+1
            IF (EL(KT,I)-Z(I)*COSA(JB) > EL(KTI(I),I)) THEN
              DO WHILE (EL(KT,I)-Z(I)*COSA(JB) > EL(KTI(I),I) .AND. KTI(I) /= 2)
                Z(I)   = (EL(KT,I)-EL(KTI(I),I)-(EL(KT,I)-EL(KTI(I),I)-Z(I)*COSA(JB))*(B(KTI(I),I)/B(KTI(I)-1,I)))/COSA(JB)
! v3.5 start
                if(macrophyte_on)then
                  ktip=kti(i)
!c  keeping track if column kti has macrophytes
                  if(ktip.gt.2)kticol(i)=.false.
                end if
! v3.5 end
                KTI(I) =  MAX(KTI(I)-1,2)
              END DO
            ELSE IF (EL(KT,I)-Z(I)*COSA(JB) < EL(KTI(I)+1,I)) THEN
              DO WHILE (EL(KT,I)-Z(I)*COSA(JB) < EL(KTI(I)+1,I))
                Z(I)   = (EL(KT,I)-EL(KTI(I)+1,I)-(EL(KT,I)-EL(KTI(I)+1,I)-Z(I)*COSA(JB))*(B(KTI(I),I)/B(KTI(I)+1,I)))/COSA(JB)
                KTI(I) =  KTI(I)+1
                if(macrophyte_on)kticol(i)=.true.   !v3.5
                IF (KTI(I) >= KB(I)) EXIT
              END DO
            END IF
            BI(KT:KB(I),I) =  B(KT:KB(I),I)
            BI(KT,I)       =  B(KTI(I),I)
            H1(KT,I)       =  H(KT,JW)-Z(I)
            AVH1(KT,I)     = (H1(KT,I)+H1(KT+1,I))*0.5
            IF (KT == KTI(I) .OR. KTI(I) >= KB(I)) THEN
              BH1(KT,I) = B(KT,I)*H1(KT,I)
            ELSE
              BH1(KT,I) = BI(KT,I)*(EL(KT,I)-Z(I)*COSA(JB)-EL(KTI(I)+1,I))/COSA(JB)
            END IF
            DO K=KTI(I)+1,KT
              BH1(KT,I) = BH1(KT,I)+Bnew(K,I)*H(K,JW) !Bnew(K,I)*H(K,JW)   ! sw 1/23/06
            END DO
            BKT(I)    = BH1(KT,I)/H1(KT,I)
            if(kbi(i) < kb(i))bkt(i)=bh1(kt,i)/(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))    ! SW 1/23/06
            VOL(KT,I) = BH1(KT,I)*DLX(I)
          END DO
          DO I=IU-1,ID
            AVHR(KT,I) = H1(KT,I)  +(H1(KT,I+1) -H1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                          !SW 07/29/04
            if(kbi(i) < kb(i))avhr(kt,i)=(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))  &
               +(H1(KT,I+1)-(el(kbi(i)+1,i+1)-el(kb(i)+1,i+1)) -H1(KT,I)+(el(kbi(i)+1,i)&
               -el(kb(i)+1,i)))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)        ! SW 1/23/06
            BHR1(KT,I) = BH1(KT,I)+(BH1(KT,I+1)-BH1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                          !SW 07/29/04
          END DO
          AVHR(KT,ID+1) = H1(KT,ID+1)
          BHR1(KT,ID+1) = BH1(KT,ID+1)
          DLVOL(JB)        = 0.0
        ELSE                                                                                                           !SW 07/16/04
          DO I=IU-1,ID+1
            BI(KT:KB(I),I) =  B(KT:KB(I),I)
            CALL GRID_AREA2
            H1(KT,I)   =  H(KT,JW)-Z(I)
            AVH1(KT,I) = (H1(KT,I)+H1(KT+1,I))*0.5
            ELT        =  EL(KT,I)-Z(I)
            CALL GRID_AREA1 (ELT,EL(KT+1,I),BH1(KT,I),BI(KT,I))
            BKT(I)    = BH1(KT,I)/H1(KT,I)
            if(kbi(i) < kb(i))bkt(i)=bh1(kt,i)/(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))    ! SW 1/23/06
            VOL(KT,I) = BH1(KT,I)*DLX(I)
          END DO
          DO I=IU-1,ID
            AVHR(KT,I) = H1(KT,I)  +(H1(KT,I+1) -H1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                          !SW 07/29/04
            if(kbi(i) < kb(i))avhr(kt,i)=(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i))) &
               +(H1(KT,I+1)-(el(kbi(i)+1,i+1)-el(kb(i)+1,i+1)) -H1(KT,I)+(el(kbi(i)+1,i)&
               -el(kb(i)+1,i)))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                                     ! SW 1/23/06
            BHR1(KT,I) = BH1(KT,I)+(BH1(KT,I+1)-BH1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                          !SW 07/29/04
          END DO
          AVHR(KT,ID+1) = H1(KT,ID+1)
          BHR1(KT,ID+1) = BH1(KT,ID+1)
          DLVOL(JB)     = 0.0
        END IF
        DO I=IU,ID
          DLVOL(JB) = DLVOL(JB)+(BH1(KT,I)-BH2(KT,I))*DLX(I)
          IF (KT == 2 .AND. H1(KT,I) > H(2,JW) .AND. .NOT. SURFACE_WARNING) THEN
            WRITE (WRN,'(A,I0,A,F0.3)') 'Water surface is above the top of layer 2 in segment ',I,' at day ',JDAY
            WARNING_OPEN    = .TRUE.
            SURFACE_WARNING = .TRUE.
          END IF
        END DO

! v3.5
        if(macrophyte_on)then
!c  if depth in kti layer becomes greater than threshold, setting
!c      macrophyte conc. in kti column to initial conc.
          do i=iu,id
            wlelev=(el(kt,i)-z(i)*cos(alpha(jb)))
            depkti=wlelev-el(kti(i)+1,i)

!******* macrophytes, setting conc. of macrophytes in new columns to
!********* initial concentration if column depth is greater than 'thrkti'
            if(.not.kticol(i).and.depkti.ge.thrkti)then
              kticol(i)=.true.
              jt=kti(i)
              mact(jt,kt,i)=0.0
              do m=1,nmc
                macrc(jt,kt,i,m)=macwbci(jw,m)
                colb=el(kti(i)+1,i)
                coldep=wlelev-colb
                macrm(jt,kt,i,m)=macwbci(jw,m)*coldep*cw(jt,i)*dlx(i)
                mact(jt,kt,i)=mact(jt,kt,i)+macwbci(jw,m)
                maCMBRT(JB,m) = maCMBRT(JB,m)+macrm(jt,kt,i,m)                   
              end do
            end if

!****** macrophytes, when column depth is less than 'thrkti', zeroing out conc.
            if(kticol(i).and.depkti.lt.thrkti)then
              kticol(i)=.false.
              jt=kti(i)
              mact(jt,kt,i)=0.0
              do m=1,nmc
                maCMBRT(JB,m) = maCMBRT(JB,m)-macrm(jt,kt,i,m) 
                macrc(jt,kt,i,m)=0.0
                macrm(jt,kt,i,m)=0.0
              end do        
            end if
          end do
        end if
! v3.5 end

!***********************************************************************************************************************************
!**                                             Task 2.2.4: Longitudinal velocities                                               **
!***********************************************************************************************************************************

        IUT = IU
        IDT = ID
        IF (UP_HEAD(JB)) IUT = IU-1
        IF (DN_HEAD(JB)) IDT = ID+1

!****** Pressures

        DO I=IUT,IDT
          DO K=KT,KB(I)
            P(K,I) = P(K-1,I)+RHO(K,I)*G*H1(K,I)*COSA(JB)
          END DO
        END DO

!****** Horizontal pressure gradients

        DO I=IUT,IDT-1
          HPG(KT,I) = DLXRHO(I)*(BKT(I)+BKT(I+1))*0.5*(H1(KT,I+1)*P(KT,I+1)-H1(KT,I)*P(KT,I))
          DO K=KT+1,KBMIN(I)
            HPG(K,I) = DLXRHO(I)*BHR2(K,I)*((P(K-1,I+1)-P(K-1,I))+(P(K,I+1)-P(K,I)))
          END DO
        END DO

!****** Boundary horizontal velocities

        IF (UP_FLOW(JB)) THEN
          IF (.NOT. HEAD_FLOW(JB)) THEN
            QINF(:,JB) = 0.0
            IF (PLACE_QIN(JW)) THEN

!************ Inflow layer

              K     = KT
              SSTOT = 0.0
              DO JC=NSSS,NSSE   
                SSTOT = SSTOT+CIN(JC,JB)
              END DO
              RHOIN = DENSITY(TIN(JB),MAX(CIN(1,JB),0.0),MAX(SSTOT,0.0))
              DO WHILE (RHOIN > RHO(K,IU) .AND. K < KB(IU))
                K = K+1
              END DO
              KTQIN(JB) = K
              KBQIN(JB) = K

!************ Layer inflows

              VQIN  =  QIN(JB)*DLT
              VQINI =  VQIN
              QINFR =  1.0
              INCR  = -1
              DO WHILE (QINFR > 0.0)
                V1 = VOL(K,IU)
                IF (K <= KB(IU)) THEN
                  IF (VQIN > 0.5*V1) THEN
                    QINF(K,JB) = 0.5*V1/VQINI
                    QINFR      = QINFR-QINF(K,JB)
                    VQIN       = VQIN-QINF(K,JB)*VQINI
                    IF (K == KT) THEN
                      K    = KBQIN(JB)
                      INCR = 1
                    END IF
                  ELSE
                    QINF(K,JB) = QINFR
                    QINFR      = 0.0
                  END IF
                  IF (INCR < 0) KTQIN(JB) = K
                  IF (INCR > 0) KBQIN(JB) = MIN(KB(IU),K)
                  K = K+INCR
                ELSE
                  QINF(KT,JB) = QINF(KT,JB)+QINFR
                  QINFR       = 0.0
                END IF
              END DO
            ELSE
              KTQIN(JB) = KT
              KBQIN(JB) = KB(IU)
              BHSUM     = 0.0
              DO K=KT,KB(IU)
                BHSUM = BHSUM+BH1(K,IU)
              END DO
              DO K=KT,KB(IU)
                QINF(K,JB) = BH1(K,IU)/BHSUM
              END DO
            END IF
            DO K=KT,KB(IU)
              U(K,IU-1) = QINF(K,JB)*QIN(JB)/BHR1(K,IU-1)
            END DO
          ELSE
            KTQIN(JB) = KT     
            KBQIN(JB) = KB(IU)      
            IF (JBUH(JB) <= BE(JW) .AND. JBUH(JB) >= BS(JW)) THEN
              DO K=KT,KB(IU)
                U(K,IU-1) = U(K,UHS(JB))*BHR1(K,UHS(JB))/BHR1(K,IU-1)
              END DO
            ELSE
              CALL UPSTREAM_VELOCITY
            END IF
          END IF
        END IF
        IF (DN_FLOW(JB)) THEN
          DO K=KT,KB(ID)
            U(K,ID) = QOUT(K,JB)/BHR1(K,ID)
          END DO
        END IF
        IF (UP_HEAD(JB)) THEN
          DO K=KT,KB(IU-1)
            U(K,IU-1) = (BHR2(K,IU-1)*U(K,IU-1)+DLT*(-SB(K,IU-1)+ST(K,IU-1)-HPG(K,IU-1)+GRAV(K,IU-1)))/BHR1(K,IU-1)
          END DO
        END IF
        IF (DN_HEAD(JB)) THEN
          DO K=KT,KB(ID+1)
            U(K,ID) = (BHR2(K,ID)*U(K,ID)+DLT*(-SB(K,ID)+ST(K,ID)-HPG(K,ID)+GRAV(K,ID)))/BHR1(K,ID)
          END DO
        END IF

!****** Horizontal velocities

        DO I=IU,ID-1
          DO K=KT,KBMIN(I)
            U(K,I) = (BHR2(K,I)*U(K,I))/BHR1(K,I)+(DLT*(-SB(K,I)+ST(K,I)-ADMZ(K,I)+ADMZ(K-1,I)-ADMX(K,I)+DM(K,I)-HPG(K,I)+GRAV(K,I)&
                     +UXBR(K,I)/H2(K,I)))/BHR1(K,I)
            IF (INTERNAL_WEIR(K,I)) U(K,I) = 0.0
          END DO
        END DO

!****** Implicit vertical eddy viscosity

        IF (IMPLICIT_VISC(JW)) THEN
          AT = 0.0; CT = 0.0; VT = 0.0; DT = 0.0
          DO I=IUT,IDT-1
            DO K=KT,KB(I)
              AT(K,I) = -DLT/BHR1(K,I)*AZ(K-1,I)*(BHR1(K-1,I)/AVHR(K-1,I)+BR(K,I))  /(AVH1(K-1,I)+AVH1(K-1,I+1))
              CT(K,I) = -DLT/BHR1(K,I)*AZ(K,I)  *(BHR1(K,I)  /AVHR(K,I)  +BR(K+1,I))/(AVH1(K,I)  +AVH1(K,I+1))
              VT(K,I) =  1.0-AT(K,I)-CT(K,I)
              DT(K,I) =  U(K,I)
            END DO
            CALL TRIDIAG(AT(:,I),VT(:,I),CT(:,I),DT(:,I),KT,KBMIN(I),KMX,U(:,I))
          END DO
        END IF

!****** Corrected horizontal velocities

        IF (UP_HEAD(JB)) THEN
          IS    =  ID
          IE    =  IU-1
          INCR  = -1
          Q(IS) =  0.0
          DO K=KT,KB(ID)
            Q(IS) = Q(IS)+U(K,IS)*BHR1(K,IS)
          END DO
          QSSUM(IS) = 0.0
          DO K=KT,KB(IS)
            QSSUM(IS) = QSSUM(IS)+QSS(K,IS)
          END DO
        ELSE
          IS   = IU-1
          IE   = ID
          INCR = 1
          IF (DN_FLOW(JB)) IE = ID-1
          Q(IS) = 0.0
          DO K=KT,KB(IU)
            Q(IS) = Q(IS)+U(K,IS)*BHR1(K,IS)
          END DO
        END IF
        QC(IS) = Q(IS)
        DO I=IS+INCR,IE,INCR
          QSSUM(I) = 0.0
          DO K=KT,KB(I)
            QSSUM(I) = QSSUM(I)+QSS(K,I)
          END DO
          BHRSUM = 0.0
          Q(I)   = 0.0
          DO K=KT,KBMIN(I)
            IF (.NOT. INTERNAL_WEIR(K,I)) THEN
              BHRSUM = BHRSUM+BHR1(K,I)
              Q(I)   = Q(I)+U(K,I)*BHR1(K,I)
            END IF
          END DO
          IF (UP_HEAD(JB)) THEN
            QC(I) = QC(I+1)+(BH1(KT,I+1)-BH2(KT,I+1))*DLX(I+1)/DLT-QSSUM(I+1)
          ELSE
            QC(I) = QC(I-1)-(BH1(KT,I)  -BH2(KT,I))  *DLX(I)  /DLT+QSSUM(I)
          END IF
          DO K=KT,KBMIN(I)
            IF (INTERNAL_WEIR(K,I)) THEN
              U(K,I) = 0.0
            ELSE
              U(K,I) =  U(K,I)+(QC(I)-Q(I))/BHRSUM
              IF (Q(I) /= 0.0) QERR(I) = (Q(I)-QC(I))/Q(I)*100.0
            END IF
          END DO
        END DO

!****** Head boundary flows

        IF (UP_HEAD(JB)) QUH1(KT:KB(IU-1),JB) = U(KT:KB(IU-1),IU-1)*BHR1(KT:KB(IU-1),IU-1)
        IF (DN_HEAD(JB)) QDH1(KT:KB(ID+1),JB) = U(KT:KB(ID+1),ID)  *BHR1(KT:KB(ID+1),ID)      

!***********************************************************************************************************************************
!**                                              Task 2.2.5: Vertical velocities                                                  **
!***********************************************************************************************************************************

        DO I=IU,ID
          DO K=KB(I)-1,KT,-1
            WT1    =  W(K+1,I)*BB(K+1,I)
            WT2    = (BHR(K+1,I)*U(K+1,I)-BHR(K+1,I-1)*U(K+1,I-1)-QSS(K+1,I))/DLX(I)
            W(K,I) = (WT1+WT2)/BB(K,I)
          END DO
        END DO
      END DO
    END DO

!***********************************************************************************************************************************
!**                                                  Task 2.2.6: Autostepping                                                     **
!***********************************************************************************************************************************

    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        DO I=CUS(JB),DS(JB)
          IF (H1(KT,I) < 0.0) THEN
            WRITE (WRN,'(A,F0.3,A,I0/4(A,F0.3))') 'Computational warning at Julian day = ',JDAY,' at segment ',I,'timestep = ',DLT,&
                                                  ' water surface deviation [Z] = ',Z(I),' m  layer thickness = ',H1(KT,I),' m'
            WARNING_OPEN = .TRUE.
            IF (DLT > DLTMIN) THEN
              WRITE (WRN,'(A,I0/2(A,F0.3),A,I0)') 'Negative surface layer thickness in segment ',I,'  time step reduced 90% to ',  &
                                                   DLT*0.1,' s on day ',JDAY,' at iteration ',NIT
              WARNING_OPEN = .TRUE.
              CURMAX       =  DLTMIN
              GO TO 220
            ELSE
              WRITE (W2ERR,'(A,F0.3/A,I0)') 'Unstable water surface elevation on day ',JDAY,'negative surface layer thickness '//  &
                                            'using minimum timestep at iteration ',NIT 
              WRITE (W2ERR,'(A)') 'Segment, Surface layer thickness, m'
              DO II=MAX(CUS(JB),I-2),MIN(DS(JB),I+2)
                WRITE (W2ERR,'(T6,I3,T21,F10.2)') II,H1(KT,I)
              END DO
              TEXT = 'Runtime error - see w2.err'
              ERROR_OPEN = .TRUE.
              GO TO 230
            END IF
          END IF
        END DO
        DO I=CUS(JB),DS(JB)
          IF (VISCOSITY_LIMIT(JW)) TAU1   = 2.0*MAX(AX(JW),DXI(JW))/(DLX(I)*DLX(I))
          IF (CELERITY_LIMIT(JW))  CELRTY = SQRT((ABS(RHO(KB(I),I)-RHO(KT,I)))/1000.0*G*DEPTHB(KBI(I),I)*0.5)               ! SW 1/23/06
          DO K=KT,KB(I)
            IF (VISCOSITY_LIMIT(JW) .AND. .NOT. IMPLICIT_VISC(JW)) TAU2 = 2.0*AZ(K,I)/(H1(K,JW)*H1(K,JW))
            QTOT(K,I) = (ABS(U(K,I))*BHR1(K,I)+ABS(U(K,I-1))*BHR1(K,I-1)+(ABS(W(K,I))*BB(K,I)+ABS(W(K-1,I))*BB(K-1,I))*DLX(I)      &
                        +DLX(I)*ABS(BH2(K,I)-BH1(K,I))/DLT+ABS(QSS(K,I)))*0.5
            DLTCAL    = 1.0/((QTOT(K,I)/BH1(K,I)+CELRTY)/DLX(I)+TAU1+TAU2+NONZERO)
            IF (DLTCAL < CURMAX) THEN
              KLOC   = K
              ILOC   = I
              CURMAX = DLTCAL 
              IF (DLTF(DLTDP)*CURMAX < MINDLT) THEN
                KMIN = K
                IMIN = I
              END IF
            END IF
          END DO
        END DO
      END DO
    END DO

!** Restore timestep dependent variables and restart calculations

220 CONTINUE
    IF (CURMAX < DLT .AND. DLT > DLTMIN) THEN
      DLT = DLTF(DLTDP)*CURMAX
      IF (DLT <= DLTMIN) THEN
        WRITE (WRN,'(A,F0.3/A,F0.3,A)') 'Computational warning at Julian day = ',JDAY,' timestep = ',DLT,' sec'
        WARNING_OPEN = .TRUE.
        DLT          =  DLTMIN
      END IF
      NV        = NV+1
      Z         = SZ
      U         = SU
      W         = SW
      AZ        = SAZ
      AVH2      = SAVH2
      AVHR      = SAVHR
      KTI       = SKTI
      BKT       = SBKT
      QSS       = 0.0
      SB        = 0.0
      DLTS      = DLT

        do jw=1,nwb                                                                            ! SW 8/25/05
        do jb=bs(jw),be(jw)
        do i=us(jb)-1,ds(jb)+1
            VOL(KTWB(JW),I) = BH2(KTWB(JW),I)*DLX(I)
            BI(KTWB(JW),I) = B(KTI(I),I)
        end do
        end do
        end do
      
      
      CURMAX    = DLTMAX(DLTDP)/DLTF(DLTDP)
      IF (PIPES) THEN   
        YS   = YSS    
        VS   = VSS    
        VST  = VSTS   
        YST  = YSTS   
        DTP  = DTPS   
        QOLD = QOLDS
      END IF     
! v3.5 start
!********** Macrophytes
      do jw=1,nwb
        do m=1,nmc
          IF (macrophyte_CALC(jw,m)) THEN          
            KT = KTWB(JW)
              DO JB=BS(Jw),BE(Jw)
                DO I=CUS(JB),DS(JB)              
                  DO K=kt,kb(i)
                    mac(k,i,m)=smac(k,i,m)
                    if(kticol(i))then
                      jt=kti(i)
                    else
                      jt=kti(i)+1
                    end if
                    je=kb(i)                
                    do j=jt,je                              
                      macrc(j,K,I,m)=smacrc(j,K,I,m)
                      macrm(j,K,I,m)=smacrm(j,K,I,m)
                    end do
                  END DO
                END DO
             end do          
          end if
        end do
      end do
! v3.5 end
      GO TO 210
    END IF
    DLTLIM(KMIN,IMIN) = DLTLIM(KMIN,IMIN)+1.0

!** Layer bottom and middle depths

    DO JW=1,NWB
      DO JB=BS(JW),BE(JW)
        DO I=CUS(JB)-1,DS(JB)
          DEPTHB(KTWB(JW),I) = H1(KTWB(JW),I)
             if(kbi(i) < kb(i))depthb(ktwb(jw),i)=(h1(ktwb(jw),i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))    ! SW 1/23/06
          DEPTHM(KTWB(JW),I) = H1(KTWB(JW),I)*0.5
             if(kbi(i) < kb(i))depthm(ktwb(jw),i)=(h1(ktwb(jw),i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))*0.5    ! SW 1/23/06
          DO K=KTWB(JW)+1,KMX    
            DEPTHB(K,I) = DEPTHB(K-1,I)+ H1(K,I)
            DEPTHM(K,I) = DEPTHM(K-1,I)+(H1(K-1,I)+H1(K,I))*0.5
          END DO
        END DO
      END DO
    END DO

!***********************************************************************************************************************************
!**                                      Task 2.3: Temporal balance terms and temperatures                                        **
!***********************************************************************************************************************************

    DO JW=1,NWB
      KT = KTWB(JW)
      IF (.NOT. NO_HEAT(JW)) THEN
        IF (.NOT. READ_RADIATION(JW)) CALL SHORT_WAVE_RADIATION (JDAY)
        IF (TERM_BY_TERM(JW))then                                      ! SW 1/25/05
           if(tair(jw).ge.5.0)then
           RAN(JW) = 5.31E-13*(273.15+TAIR(JW))**6*(1.0+0.0017*CLOUD(JW)**2)*0.97
           else
           RAN(JW) = 5.62E-8*(273.15+TAIR(JW))**4*(1.-0.261*exp(-7.77E-4*TAIR(JW)**2))*(1.0+0.0017*CLOUD(JW)**2)*0.97
           endif
        ENDIF
      END IF
      DO JB=BS(JW),BE(JW)
        IU = CUS(JB)
        ID = DS(JB)

!****** Heat exchange

        IF (.NOT. NO_HEAT(JW)) THEN
          DO I=IU,ID
            IF (DYNAMIC_SHADE(I)) CALL SHADING

!********** Surface

            IF (.NOT. ICE(I)) THEN
              IF (TERM_BY_TERM(JW)) THEN
                CALL SURFACE_TERMS (T2(KT,I))
                RS(I)     = SRON(JW)*SHADE(I)
                RN(I)     = RS(I)+RAN(JW)-RB(I)-RE(I)-RC(I)
                HEATEX    = RN(I)/RHOWCP*BI(KT,I)*DLX(I)
              ELSE
                CALL EQUILIBRIUM_TEMPERATURE
                HEATEX = (ET(I)-T2(KT,I))*CSHE(I)*BI(KT,I)*DLX(I)
              END IF
              TSS(KT,I) =  TSS(KT,I)+HEATEX
              TSSS(JB)  =  TSSS(JB) +HEATEX*DLT
              SROOUT    = (1.0-BETA(JW))*(SRON(JW)*SHADE(I)/RHOWCP)*BI(KT,I)*DLX(I)*EXP(-GAMMA(KT,I)*DEPTHB(KT,I))
              TSS(KT,I) =  TSS(KT,I)-SROOUT
              TSSS(JB)  =  TSSS(JB) -SROOUT*DLT
              if(kt == kb(i))then    ! SW 4/18/07
              SROSED    =  SROOUT*TSEDF(JW)
              else
              SROSED    =  SROOUT*(1.0-BI(KT+1,I)/BI(KT,I))*TSEDF(JW)
              endif
              TSS(KT,I) =  TSS(KT,I)+SROSED
              TSSS(JB)  =  TSSS(JB) +SROSED*DLT
              SROIN     =  SROOUT*B(KT+1,I)/BI(KT,I)
              DO K=KT+1,KB(I)
                SROOUT   = SROIN*EXP(-GAMMA(K,I)*(H1(K,I)))
                SRONET   = SROIN-SROOUT
                SROSED   = SROOUT*(1.0-BI(K+1,I)/BI(K,I))*TSEDF(JW)
                TSS(K,I) = TSS(K,I)+ SRONET+SROSED
                TSSS(JB) = TSSS(JB)+(SRONET+SROSED)*DLT
                SROIN    = SROOUT*B(K+1,I)/B(K,I)
              END DO
            END IF

!********** Sediment/water

            DO K=KT,KB(I)
              if(k==kb(i))then                ! SW 4/18/07
              TFLUX    = CBHE(JW)/RHOWCP*(TSED(JW)-T2(K,I))*BI(K,I)*DLX(I)
              else
              TFLUX    = CBHE(JW)/RHOWCP*(TSED(JW)-T2(K,I))*(BI(K,I)-BI(K+1,I))*DLX(I)
              endif
              TSS(K,I) = TSS(K,I)+TFLUX
              TSSB(JB) = TSSB(JB)+TFLUX*DLT
            END DO
          END DO

!******** Ice cover

          IF (ICE_CALC(JW)) THEN
            HIA = 0.2367*CSHE(I)/5.65E-8
            DO I=IU,ID
              ALLOW_ICE(I) = .TRUE.
              DO K=KT,KB(I)
                IF (T2(K,I) > ICET2(JW)) ALLOW_ICE(I) = .FALSE.
              END DO
            END DO
            ICE_IN(JB) = .TRUE.
            DO I=IU,ID
              IF (ICETH(I) < ICEMIN(JW)) ICE_IN(JB) = .FALSE.
            END DO
            DO I=IU,ID
              IF (DETAILED_ICE(JW)) THEN
                IF (T2(KT,I) < 0.0) THEN
                  IF (.NOT. ICE(I)) THEN
                    ICETH2 = -T2(KT,I)*RHO(KT,I)*CP*H2(KT,I)/RHOIRL1
                    IF (ICETH2 < ICE_TOL) THEN
                      ICETH2 = 0.0
                    ELSE
                      TFLUX      = T2(KT,I)*RHO(KT,I)*CP*H2(KT,I)*BI(KT,I)/(RHOWCP*DLT)*DLX(I)
                      TSS(KT,I)  = TSS(KT,I) -TFLUX
                      TSSICE(JB) = TSSICE(JB)-TFLUX*DLT
                    END IF
                  END IF
                END IF

!************** Ice balance

                IF (ICE(I)) THEN
                  TICE = TAIR(JW)
                  DEL  = 2.0
                  J    = 1
                  DO WHILE (DEL > 1.0 .AND. J < 500)
                    CALL SURFACE_TERMS (TICE)
                    RN(I) = SRON(JW)/(REFL*RHOWCP)*SHADE(I)*(1.0-ALBEDO(JW))*BETAI(JW)+RAN(JW)-RB(I)-RE(JW)-RC(I)
                    DEL   = RN(I)+RK1*(RIMT-TICE)/ICETH(I)
                    IF (ABS(DEL) > 1.0) TICE = TICE+DEL/500.0
                    J = J+1
                  END DO

!**************** Solar radiation attenuation

                  TFLUX      = DLX(I)*SRON(JW)/(REFL*RHOWCP)*SHADE(I)*(1.0-ALBEDO(JW))*(1.0-BETAI(JW))                             &
                               *EXP(-GAMMAI(JW)*ICETH(I))*BI(KT,I)
                  TSS(KT,I)  = TSS(KT,I) +TFLUX
                  TSSICE(JB) = TSSICE(JB)+TFLUX*DLT
                  IF (TICE > 0.0) THEN
                    HICE   =  RHOICP*0.5*TICE*0.5*ICETH(I)*BI(KT,I)/(RHOWCP*DLT)
                    ICETHU = -DLT*HICE/B(KTI(I),I)*RHOWCP/RHOIRL1
                    TICE   =  0.0
                  END IF

!**************** Ice growth

                  IF (TICE < 0.0) ICETH1 = DLT*(RK1*(RIMT-TICE)/ICETH(I))/RHOIRL1

!**************** Ice melt

                  IF (T2(KT,I) > 0.0) THEN
                    ICETH2     = -DLT*HWI(JW)*(T2(KT,I)-RIMT)/RHOIRL1
                    TFLUX      =  2.392E-7*HWI(JW)*(RIMT-T2(KT,I))*BI(KT,I)*DLX(I)
                    TSS(KT,I)  =  TSS(KT,I) +TFLUX
                    TSSICE(JB) =  TSSICE(JB)+TFLUX*DLT
                  END IF
                END IF

!************** Ice thickness

                ICETH(I) = ICETH(I)+ICETHU+ICETH1+ICETH2
                IF (ICETH(I) < ICE_TOL) ICETH(I) = 0.0
                IF (WINTER .AND. (.NOT. ICE_IN(JB))) THEN
                  IF (.NOT. ALLOW_ICE(I)) ICETH(I) = 0.0
                END IF
                ICE(I)   = ICETH(I) > 0.0
                IF (ICE(I)) ICESW(I) = 0.0
                ICETHU   = 0.0
                ICETH1   = 0.0
                ICETH2   = 0.0
                ICESW(I) = 1.0
                IF (ICETH(I) < ICE_TOL .AND. ICETH(I) > 0.0) ICETH(I) = ICE_TOL
              ELSE
                HIA      = 0.2367*CSHE(I)/5.65E-8
                ICETH(I) = MAX(0.0,ICETH(I)+DLT*((RIMT-ET(I))/(ICETH(I)/RK1+1.0/HIA)-(T2(KT,I)-RIMT))/RHOIRL1)
                ICE(I)   = ICETH(I) > 0.0
                ICESW(I) = 1.0
                IF (ICE(I)) THEN
                  TFLUX      = 2.392E-7*(RIMT-T2(KT,I))*BI(KT,I)*DLX(I)
                  TSS(KT,I)  = TSS(KT,I) +TFLUX
                  TSSICE(JB) = TSSICE(JB)+TFLUX*DLT
                END IF
              END IF
            END DO
          END IF
        END IF

!****** Heat sources/sinks and total inflow/outflow

        IF (EVAPORATION(JW)) THEN
          DO I=IU,ID
            TSS(KT,I) = TSS(KT,I)-EV(I)*T2(KT,I)
            TSSEV(JB) = TSSEV(JB)-EV(I)*T2(KT,I)*DLT
            VOLEV(JB) = VOLEV(JB)-EV(I)         *DLT
          END DO
        END IF
        IF (PRECIPITATION(JW)) THEN
          DO I=IU,ID
            TSS(KT,I) = TSS(KT,I)+QPR(I)*TPR(JB)
            TSSPR(JB) = TSSPR(JB)+QPR(I)*TPR(JB)*DLT
            VOLPR(JB) = VOLPR(JB)+QPR(I)        *DLT
          END DO
        END IF
        IF (TRIBUTARIES) THEN
          DO JT=1,JTT
            IF (JB == JBTR(JT)) THEN
              I = ITR(JT)
              IF (I < CUS(JB)) I = CUS(JB)
              DO K=KTTR(JT),KBTR(JT)
                IF (QTR(JT) < 0) THEN
                  TSS(K,I)  = TSS(K,I) +T2(K,I)*QTR(JT)*QTRF(K,JT)
                  TSSTR(JB) = TSSTR(JB)+T2(K,I)*QTR(JT)*QTRF(K,JT)*DLT
                ELSE
                  TSS(K,I)  = TSS(K,I) +TTR(JT)*QTR(JT)*QTRF(K,JT)
                  TSSTR(JB) = TSSTR(JB)+TTR(JT)*QTR(JT)*QTRF(K,JT)*DLT
                END IF
              END DO
              VOLTRB(JB) = VOLTRB(JB)+QTR(JT)*DLT
            END IF
          END DO
        END IF
        IF (DIST_TRIBS(JB)) THEN
          DO I=IU,ID
            IF (QDT(I) < 0) THEN
              TSS(KT,I) = TSS(KT,I)+T2(KT,I)*QDT(I)
              TSSDT(JB) = TSSDT(JB)+T2(KT,I)*QDT(I)*DLT
            ELSE
              TSS(KT,I) = TSS(KT,I)+TDTR(JB)*QDT(I)
              TSSDT(JB) = TSSDT(JB)+TDTR(JB)*QDT(I)*DLT
            END IF
            VOLDT(JB) = VOLDT(JB)+QDT(I)*DLT
          END DO
        END IF
        IF (WITHDRAWALS) THEN
          DO JWD=1,JWW
            IF (QWD(JWD) /= 0.0) THEN
              IF (JB == JBWD(JWD)) THEN
                I = MAX(CUS(JBWD(JWD)),IWD(JWD))
                DO K=KTW(JWD),KBW(JWD)
                  TSS(K,I)  = TSS(K,I) -T2(K,I)*QSW(K,JWD)
                  TSSWD(JB) = TSSWD(JB)-T2(K,I)*QSW(K,JWD)*DLT
                END DO
                VOLWD(JB) = VOLWD(JB)-QWD(JWD)*DLT
              END IF
            END IF
          END DO
        END IF
        IF (UP_FLOW(JB)) THEN
            DO K=KT,KB(IU)
              IF (.NOT. HEAD_FLOW(JB)) THEN
                TSS(K,IU) = TSS(K,IU)+QINF(K,JB)*QIN(JB)*TIN(JB)
                TSSIN(JB) = TSSIN(JB)+QINF(K,JB)*QIN(JB)*TIN(JB)*DLT
              ELSE
                IF (U(K,IU-1) >= 0.0) THEN   
                  TSS(K,IU) = TSS(K,IU)+U(K,IU-1)*BHR1(K,IU-1)*T1(K,IU-1)
                  TSSIN(JB) = TSSIN(JB)+U(K,IU-1)*BHR1(K,IU-1)*T1(K,IU-1)*DLT
                ELSE
                  TSS(K,IU) = TSS(K,IU)+U(K,IU-1)*BHR1(K,IU-1)*T1(K,IU)
                  TSSIN(JB) = TSSIN(JB)+U(K,IU-1)*BHR1(K,IU-1)*T1(K,IU)*DLT
                END IF
              END IF
            END DO
          VOLIN(JB) = VOLIN(JB)+QIN(JB)*DLT
        END IF
        IF (DN_FLOW(JB)) THEN
          DO K=KT,KB(ID)
            TSS(K,ID)  = TSS(K,ID) -QOUT(K,JB)*T2(K,ID+1)
            TSSOUT(JB) = TSSOUT(JB)-QOUT(K,JB)*T2(K,ID+1)*DLT
            VOLOUT(JB) = VOLOUT(JB)-QOUT(K,JB)           *DLT
          END DO
        END IF
        IF (UP_HEAD(JB)) THEN
          DO K=KT,KB(IU)
            IUT = IU
            IF (QUH1(K,JB) >= 0.0) IUT = IU-1
            TSSUH1(K,JB) = T2(K,IUT)*QUH1(K,JB)
            TSS(K,IU)    = TSS(K,IU)+TSSUH1(K,JB)
            TSSUH(JB)    = TSSUH(JB)+TSSUH1(K,JB)*DLT
            VOLUH(JB)    = VOLUH(JB)+QUH1(K,JB)  *DLT
          END DO
        END IF
        IF (UH_INTERNAL(JB)) THEN
          IF (UHS(JB) /= DS(JBUH(JB)) .OR. DHS(JBUH(JB)) /= US(JB)) THEN
            IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
              DO K=KT,KB(IU-1)
                TSS(K,UHS(JB))  = TSS(K,UHS(JB)) -TSSUH2(K,JB)/DLT
                TSSUH(JBUH(JB)) = TSSUH(JBUH(JB))-TSSUH2(K,JB)
                VOLUH(JBUH(JB)) = VOLUH(JBUH(JB))-VOLUH2(K,JB)
              END DO
            ELSE
              CALL UPSTREAM_CONSTITUENT(T2,TSS)
              DO K=KT,KB(IU-1)
                TSSUH(JBUH(JB)) = TSSUH(JBUH(JB))-TSSUH2(K,JB)
                VOLUH(JBUH(JB)) = VOLUH(JBUH(JB))-VOLUH2(K,JB)
              END DO
            END IF
          END IF
        END IF
        IF (DN_HEAD(JB)) THEN
          DO K=KT,KB(ID+1)
            IDT = ID+1
            IF (QDH1(K,JB) >= 0.0) IDT = ID
            TSSDH1(K,JB) = T2(K,IDT)*QDH1(K,JB)
            TSS(K,ID)    = TSS(K,ID)-TSSDH1(K,JB)
            TSSDH(JB)    = TSSDH(JB)-TSSDH1(K,JB)*DLT
            VOLDH(JB)    = VOLDH(JB)-QDH1(K,JB)  *DLT
          END DO
        END IF
        IF (DH_INTERNAL(JB)) THEN
          IF (DHS(JB) /= US(JBDH(JB)) .OR. UHS(JBDH(JB)) /= DS(JB)) THEN
            IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
              DO K=KT,KB(ID+1)
                TSS(K,CDHS(JB)) = TSS(K,CDHS(JB))+TSSDH2(K,JB)/DLT
                TSSDH(JBDH(JB)) = TSSDH(JBDH(JB))+TSSDH2(K,JB)
                VOLDH(JBDH(JB)) = VOLDH(JBDH(JB))+VOLDH2(K,JB)
              END DO
            ELSE
              CALL DOWNSTREAM_CONSTITUENT(T2,TSS)
              DO K=KT,KB(ID+1)
                TSSDH(JBDH(JB)) = TSSDH(JBDH(JB))+TSSDH2(K,JB)
                VOLDH(JBDH(JB)) = VOLDH(JBDH(JB))+VOLDH2(K,JB)
              END DO
            END IF
          END IF
        END IF
      END DO
    END DO

!** Temperature transport

    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IU =  CUS(JB)
        ID =  DS(JB)
        COLD => HYD(:,:,4)
        CALL HORIZONTAL_MULTIPLIERS
        CALL VERTICAL_MULTIPLIERS
        CNEW => T1(:,:)
        SSB  => TSS(:,:)
        SSK  => CSSB(:,:,1)
        CALL HORIZONTAL_TRANSPORT
        DO I=IU,ID
          AT(:,I) = 0.0; CT(:,I) = 0.0; VT(:,I) = 0.0; DT(:,I) = 0.0
          DO K=KT,KB(I)
            AT(K,I) = -DLT/BH1(K,I)*(BB(K-1,I)*(DZ(K-1,I)/AVH1(K-1,I)+THETA(JW)*0.5*W(K-1,I)))
            CT(K,I) =  DLT/BH1(K,I)*(BB(K,I)*(THETA(JW)*0.5*W(K,I)-DZ(K,I)/AVH1(K,I)))
            VT(K,I) =  1.0+DLT/BH1(K,I)*(BB(K,I)*(DZ(K,I)/AVH1(K,I)+THETA(JW)*0.5*W(K,I))+BB(K-1,I)*(DZ(K-1,I)/AVH1(K-1,I)         &
                      -THETA(JW)*0.5*W(K-1,I)))
            DT(K,I) =  CNEW(K,I)
          END DO
          CALL TRIDIAG(AT(:,I),VT(:,I),CT(:,I),DT(:,I),KT,KB(I),KMX,CNEW(:,I))
        END DO
      END DO
    END DO

!***********************************************************************************************************************************
!**                                                    Task 2.4:  Constituents                                                    **
!***********************************************************************************************************************************

    IF (CONSTITUENTS) THEN
      if(macrophyte_on.and.update_kinetics)call porosity !v3.5
      DO JW=1,NWB
        KT = KTWB(JW)
        DO JB=BS(JW),BE(JW)
          IU = CUS(JB)
          ID = DS(JB)

!******** Kinetic sources/sinks

! v3.5 start          
          IF (SEDIMENT_CALC(JW))then
            CALL SEDIMENT
            CALL SEDIMENTp
            CALL SEDIMENTn
            CALL SEDIMENTc
            CALL SEDIMENT_decay_rate
          end if
          do m=1,nmc
            if (macrophyte_calc(jw,m))then
              call macrophyte(m)                              
            end if
          end do
! v3.5 end              
          IF (UPDATE_KINETICS) THEN
            IF (UPDATE_RATES) THEN
              CALL TEMPERATURE_RATES
              CALL KINETIC_RATES
            END IF
            DO JAC=1,NAC
              JC = CN(JAC)
              IF (JC == NPO4)                    CALL PHOSPHORUS
              IF (JC == NNH4)                    CALL AMMONIUM
              IF (JC == NNO3)                    CALL NITRATE
              IF (JC == NDSI)                    CALL DISSOLVED_SILICA
              IF (JC == NPSI)                    CALL PARTICULATE_SILICA
              IF (JC == NFE)                     CALL IRON
              IF (JC == NLDOM)                   CALL LABILE_DOM
              IF (JC == NRDOM)                   CALL REFRACTORY_DOM
              IF (JC == NLPOM)                   CALL LABILE_POM
              IF (JC == NRPOM)                   CALL REFRACTORY_POM
              IF (JC == NDO)                     CALL DISSOLVED_OXYGEN
              IF (JC >= NGCS  .AND. JC <= NGCE)  CALL GENERIC_CONST(JC-NGCS+1)
              IF (JC >= NSSS  .AND. JC <= NSSE)  CALL SUSPENDED_SOLIDS(JC-NSSS+1)
              IF (JC >= NAS   .AND. JC <= NAE)   CALL ALGAE(JC-NAS+1)
              IF (JC >= NBODS .AND. JC <= NBODE) CALL BIOCHEMICAL_O2_DEMAND(JC-NBODS+1)
! v3.5 start             
			  IF (JC >= NZOOS  .AND. JC <= NZOOE)  CALL zooplankton  !v3.5			               
              IF (JC == NLDOMP)                CALL LABILE_DOM_p
              IF (JC == NRDOMP)                CALL refractory_DOM_p
              IF (JC == NLPOMP)                CALL labile_POM_p
              IF (JC == NRPOMP)                CALL refractory_POM_p
              IF (JC == NLDOMN)                CALL LABILE_DOM_n
              IF (JC == NRDOMN)                CALL refractory_DOM_n
              IF (JC == NLPOMN)                CALL labile_POM_n
              IF (JC == NRPOMN)                CALL refractory_POM_n              
            END DO
! v3.5 end
            IF (PH_CALC(JW)) CALL INORGANIC_CARBON
            IF (PH_CALC(JW)) CALL PH_CO2
          END IF
          DO JE=1,NEP   ! sw 5/16/06
            IF (EPIPHYTON_CALC(JW,JE)) CALL EPIPHYTON(JE)
          END DO

!******** External sources/sinks

          DO JAC=1,NAC
            JC = CN(JAC)
            IF (TRIBUTARIES) THEN
              DO JT=1,JTT
                IF (JB == JBTR(JT)) THEN
                  I = ITR(JT)
                  IF (I < CUS(JB)) I = CUS(JB)
                  DO K=KTTR(JT),KBTR(JT)
                    IF (QTR(JT) < 0.0) THEN
                      CSSB(K,I,JC) = CSSB(K,I,JC)+C1(K,I,JC)*QTR(JT)*QTRF(K,JT)
                    ELSE
                      CSSB(K,I,JC) = CSSB(K,I,JC)+CTR(JC,JT)*QTR(JT)*QTRF(K,JT)
                    END IF
                  END DO
                END IF
              END DO
            END IF
            IF (DIST_TRIBS(JB)) THEN
              DO I=IU,ID
                IF (QDT(I) < 0.0) THEN
                  CSSB(KT,I,JC) = CSSB(KT,I,JC)+C1(KT,I,JC)*QDT(I)
                ELSE
                  CSSB(KT,I,JC) = CSSB(KT,I,JC)+CDTR(JC,JB)*QDT(I)
                END IF
              END DO
            END IF
            IF (WITHDRAWALS) THEN
              DO JWD=1,JWW
                IF (QWD(JWD) /= 0.0) THEN
                  IF (JB == JBWD(JWD)) THEN
                    I = MAX(CUS(JBWD(JWD)),IWD(JWD))
                    DO K=KTW(JWD),KBW(JWD)
                      CSSB(K,I,JC) = CSSB(K,I,JC)-C1S(K,I,JC)*QSW(K,JWD)
                    END DO
                  END IF
                END IF
              END DO
            END IF
            IF (PRECIPITATION(JW)) THEN
              DO I=IU,ID
                CSSB(KT,I,JC) = CSSB(KT,I,JC)+CPR(JC,JB)*QPR(I)
              END DO
            END IF
            IF (UP_FLOW(JB)) THEN
              DO K=KT,KB(IU)
                IF (.NOT. HEAD_FLOW(JB)) THEN
                  CSSB(K,IU,JC) = CSSB(K,IU,JC)+QINF(K,JB)*QIN(JB)*CIN(JC,JB)
                ELSE
                  IF (U(K,IU-1) >= 0.0) THEN   
                    CSSB(K,IU,JC) = CSSB(K,IU,JC)+U(K,IU-1)*BHR1(K,IU-1)*C1S(K,IU-1,JC)
                  ELSE
                    CSSB(K,IU,JC) = CSSB(K,IU,JC)+U(K,IU-1)*BHR1(K,IU-1)*C1S(K,IU,JC)
                  END IF
                END IF
              END DO
            END IF
            IF (DN_FLOW(JB)) CSSB(KT:KB(ID),ID,JC) = CSSB(KT:KB(ID),ID,JC)-QOUT(KT:KB(ID),JB)*C1S(KT:KB(ID),ID,JC)
            IF (UP_HEAD(JB)) THEN
              DO K=KT,KB(IU)
                IUT = IU
                IF (QUH1(K,JB) >= 0.0) IUT = IU-1
                CSSUH1(K,JC,JB) = C1S(K,IUT,JC)*QUH1(K,JB)
                CSSB(K,IU,JC)   = CSSB(K,IU,JC)+CSSUH1(K,JC,JB)
              END DO
              IF (UH_INTERNAL(JB)) THEN
                IF (UHS(JB) /= DS(JBUH(JB)) .OR. DHS(JBUH(JB)) /= US(JB)) THEN
                  IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
                    I = UHS(JB)
                    DO K=KT,KB(IU)
                      CSSB(K,I,JC) = CSSB(K,I,JC)-CSSUH2(K,JC,JB)/DLT
                    END DO
                  ELSE
                    CALL UPSTREAM_CONSTITUENT(C2(:,:,JC),CSSB(:,:,JC))
                  END IF
                END IF
              END IF
            END IF
            IF (DN_HEAD(JB)) THEN
              DO K=KT,KB(ID+1)
                IDT = ID+1
                IF (QDH1(K,JB) >= 0.0) IDT = ID     
                CSSDH1(K,JC,JB) = C1S(K,IDT,JC)*QDH1(K,JB)
                CSSB(K,ID,JC)   = CSSB(K,ID,JC)-CSSDH1(K,JC,JB)
              END DO
              IF (DH_INTERNAL(JB)) THEN
                IF (DHS(JB) /= US(JBDH(JB)) .OR. UHS(JBDH(JB)) /= DS(JB)) THEN
                  IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
                    I = DHS(JB)
                    DO K=KT,KB(ID+1)
                      CSSB(K,I,JC) = CSSB(K,I,JC)+CSSDH2(K,JC,JB)/DLT
                    END DO
                  ELSE
                    CALL DOWNSTREAM_CONSTITUENT(C2(:,:,JC),CSSB(:,:,JC))
                  END IF
                END IF
              END IF
            END IF
          END DO
        END DO
      END DO

!**** Kinetic fluxes

      DO JW=1,NWB
        IF (FLUX(JW)) CALL KINETIC_FLUXES
      END DO

!**** Constituent transport

      DO JW=1,NWB
        KT = KTWB(JW)
        DO JB=BS(JW),BE(JW)
          IU = CUS(JB)
          ID = DS(JB)
          DO JAC=1,NAC
            JC   =  CN(JAC)
            COLD => C1S(:,:,JC)
            CALL HORIZONTAL_MULTIPLIERS
            CALL VERTICAL_MULTIPLIERS
            CNEW => C1(:,:,JC)
            SSB  => CSSB(:,:,JC)
            SSK  => CSSK(:,:,JC)
            CALL HORIZONTAL_TRANSPORT
            DT = CNEW
            DO I=IU,ID
              CALL TRIDIAG(AT(:,I),VT(:,I),CT(:,I),DT(:,I),KT,KB(I),KMX,CNEW(:,I))
            END DO
          END DO
        END DO
      END DO
      IF (DERIVED_CALC) CALL DERIVED_CONSTITUENTS
    END IF

!***********************************************************************************************************************************
!**                                       Task 2.5: Layer - Segment Additions and Subtractions                                    **
!***********************************************************************************************************************************

!** Water surface minimum thickness

    DO JW=1,NWB
      KT       =  KTWB(JW)
      ZMIN(JW) = -1000.0
      KTMAX    =  2                                                                                                 ! SR 10/17/05
      DO JB=BS(JW),BE(JW)
        DO I=CUS(JB),DS(JB)
          IF(KB(I) > KTMAX) KTMAX = KB(I)                                                                           ! SR 10/17/05
          IF (Z(I) > ZMIN(JW)) THEN
            IZMIN(JW) = I
            JBIZ      = JB
          END IF
          ZMIN(JW) = MAX(ZMIN(JW),Z(I))
        END DO
      END DO
      ADD_LAYER = ZMIN(JW) < -0.85*H(KT-1,JW) .AND. KT /= 2
      SUB_LAYER = ZMIN(JW) >  0.60*H(KT,JW)   .AND. KT < KTMAX                                                       ! SR 10/17/05
      IF (KTWB(JW) == KMX-1 .AND. SLOPE(JBIZ) > 0.0 .AND. SUB_LAYER .AND. ONE_LAYER(IZMIN(JW))) THEN
        IF (ZMIN(JW) > 0.99*H(KT,JW)) WRITE (WRN,'(A,I0,2(A,F0.3))') 'Low water in segment ',IZMIN(JW),&
                                              ' water surface deviation'//' = ',ZMIN(JW),' at day ',JDAY
        WARNING_OPEN = .TRUE.
        SUB_LAYER    = .FALSE.
      END IF

!**** Add layers

      DO WHILE (ADD_LAYER)
        IF (SNAPSHOT(JW)) WRITE (SNP(JW),'(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,13("*"))') '   Add layer ',KT-1,&
                                                        ' at Julian day = ',JDAY,'   NIT = ',NIT,' IZMIN =',IZMIN(JW)   ! SW 1/23/06

!****** Variable initialization

        KTWB(JW) = KTWB(JW)-1
        KT       = KTWB(JW)
        DO JB=BS(JW),BE(JW)
          IU = CUS(JB)
          ID = DS(JB)
          DO I=IU-1,ID+1
            Z(I)          =  H(KT,JW)+Z(I)
            H1(KT,I)      =  H(KT,JW)-Z(I)
            H1(KT+1,I)    =  H(KT+1,JW)
            H2(KT+1,I)    =  H(KT+1,JW)
            AVH1(KT,I)    = (H1(KT,I)  +H1(KT+1,I))*0.5
            AVH1(KT+1,I)  = (H1(KT+1,I)+H1(KT+2,I))*0.5
            IF (.NOT. TRAPEZOIDAL(JW)) THEN
              BH1(KT,I)   = BH1(KT+1,I)-Bnew(KT+1,I)*H1(KT+1,I)                              ! SW 1/23/06
              BH1(KT+1,I) = Bnew(KT+1,I)*H1(KT+1,I)                                          ! SW 1/23/06
            ELSE
              ELT = EL(KT,I)-Z(I)
              CALL GRID_AREA1(ELT,EL(KT+1,I),BH1(KT,I),DUMMY)                                                          !SW 08/03/04
              BH1(KT+1,I) = 0.25*H1(KT+1,JW)*(BB(KT,I)+2.*B(KT+1,I)+BB(KT+1,I))
            ENDIF
            VOL(KT,I)     = BH1(KT,I)  *DLX(I)
            VOL(KT+1,I)   = BH1(KT+1,I)*DLX(I)
            BKT(I)        = BH1(KT,I)/H1(KT,I)
            DEPTHB(KT,I)  = H1(KT,I)
            DEPTHM(KT,I)  = H1(KT,I)*0.5
            BI(KT:KB(I),I) =  B(KT:KB(I),I)   ! SW 8/26/05
            BI(KT,I)       =  B(KTI(I),I)
            T1(KT,I)      = T1(KT+1,I)
!            RHO(KT,I)     = DENSITY(T1(KT,I),MAX(TDS(KT,I),0.0),MAX(TISS(KT,I),0.0))    ! SR 5/15/06
            DO K=KT+1,KMX
              DEPTHB(K,I) = DEPTHB(K-1,I)+ H1(K,I)
              DEPTHM(K,I) = DEPTHM(K-1,I)+(H1(K-1,I)+H1(K,I))*0.5
            END DO
            C1(KT,I,CN(1:NAC))             = C1(KT+1,I,CN(1:NAC))
            CSSK(KT,I,CN(1:NAC))           = CSSK(KT+1,I,CN(1:NAC))
            KF(KT,I,KFCN(1:NAF(JW),JW))    = KF(KT+1,I,KFCN(1:NAF(JW),JW))
            KFS(KT,I,KFCN(1:NAF(JW),JW))   = KF(KT+1,I,KFCN(1:NAF(JW),JW))
            KF(KT+1,I,KFCN(1:NAF(JW),JW))  = 0.0
            KFS(KT+1,I,KFCN(1:NAF(JW),JW)) = 0.0
            if(kt >= kbi(i))then         ! CB 5/24/06
            adx(kt+1,i)=0.0             ! CB 5/15/06
            c1(kt+1,i,cn(1:nac))=0.0    ! CB 5/15/06
            cssk(kt+1,i,cn(1:nac))=0.0  ! CB 5/15/06
            endif                       ! CB 5/15/06
            DO JE=1,NEP
              if(kt < kbi(i))then    ! CB 4/28/06
              EPD(KT,I,JE)   = EPD(KT+1,I,JE)
              EPM(KT,I,JE)   = EPD(KT,I,JE)*(Bi(KT,I)-B(KT+1,I)+2.0*H1(KT,I))*DLX(I)             ! SR 5/15/06
              EPM(KT+1,I,JE) = EPM(KT+1,I,JE)-EPM(KT,I,JE)
              EPC(KT,I,JE)   = EPM(KT,I,JE)/VOL(KT,I)
              EPC(KT+1,I,JE) = EPM(KT+1,I,JE)/VOL(KT+1,I)
              else
              EPD(KT,I,JE)   = EPD(KT+1,I,JE) ! SW 5/15/06
              EPM(KT,I,JE)   = EPM(KT+1,I,JE)
              EPC(KT,I,JE) =   EPC(KT+1,I,JE)
              EPD(KT+1,I,JE) = 0.0
              EPM(KT+1,I,JE) = 0.0
              EPC(KT+1,I,JE) = 0.0
              end if                ! CB 4/28/06
            END DO
          END DO
! v3.5 start
!********macrophytes...
          do i=iu,id
            jt=kti(i)
            je=kb(i)
            wlelev=(el(kt,i)-z(i)*cosa(jb))           
            do j=jt,je
              if(j.lt.kt)then                 
                colb=el(j+1,i)
              else
                colb=el(kt+1,i)
              end if
              coldep=wlelev-colb
              mact(j,kt,i)=mact(j,kt+1,i)
              do m=1,nmc
                if(macrophyte_calc(jw,m))then
                  macrc(j,kt,i,m)=macrc(j,kt+1,i,m)
                  macrm(j,kt,i,m)=macrc(j,kt,i,m)*cw(j,i)*coldep*dlx(i)
                end if
              end do
            end do

            jt=kt+1
            je=kb(i)              
            do j=jt,je
              do m=1,nmc
                if(macrophyte_calc(jw,m))then
                  macrm(j,kt+1,i,m)=macrc(j,kt+1,i,m)*cw(j,i)*h(kt+1,jw)*dlx(i)
                end if
              end do
            end do
          end do
! v3.5 end          
          DO I=IU-1,ID
            AVHR(KT+1,I) = H1(KT+1,I) +(H1(KT+1,I+1) -H1(KT+1,I)) /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
            AVHR(KT,I)   = H1(KT,I)   +(H1(KT,I+1)   -H1(KT,I))   /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
            BHR1(KT,I)   = BH1(KT,I)  +(BH1(KT,I+1)  -BH1(KT,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
            BHR1(KT+1,I) = BH1(KT+1,I)+(BH1(KT+1,I+1)-BH1(KT+1,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
            U(KT,I)      = U(KT+1,I)
          END DO
          DO I=IU,ID
            IF (ONE_LAYER(I)) THEN
              W(KT,I) = 0.0
            ELSE
              W1      =  W(KT+1,I)*BB(KT+1,I)
              W2      = (BHR1(KT+1,I)*U(KT+1,I)-BHR1(KT+1,I-1)*U(KT+1,I-1))/DLX(I)
              W3      = (-QSS(KT+1,I)*BH1(KT+1,I)/(BH1(KT+1,I)+BH1(KT,I)))/DLX(I)
              W(KT,I) = (W1+W2+W3)/BB(KT,I)
            END IF
          END DO
          IF (UP_HEAD(JB)) THEN
            BHSUM                     = BHR1(KT,IU-1)            +BHR1(KT+1,IU-1)
            QUH1(KT,JB)               = QUH1(KT+1,JB)            *BHR1(KT,IU-1)  /BHSUM
            QUH1(KT+1,JB)             = QUH1(KT+1,JB)            *BHR1(KT+1,IU-1)/BHSUM
            TSSUH1(KT,JB)             = TSSUH1(KT+1,JB)          *BHR1(KT,IU-1)  /BHSUM
            TSSUH1(KT+1,JB)           = TSSUH1(KT+1,JB)          *BHR1(KT+1,IU-1)/BHSUM
            CSSUH1(KT,CN(1:NAC),JB)   = CSSUH1(KT+1,CN(1:NAC),JB)*BHR1(KT,IU-1)  /BHSUM
            CSSUH1(KT+1,CN(1:NAC),JB) = CSSUH1(KT+1,CN(1:NAC),JB)*BHR1(KT+1,IU-1)/BHSUM
          END IF
          IF (DN_HEAD(JB)) THEN
            BHSUM                     = BHR1(KT,ID)              +BHR1(KT+1,ID)
            QDH1(KT,JB)               = QDH1(KT+1,JB)            *BHR1(KT,ID)    /BHSUM
            QDH1(KT+1,JB)             = QDH1(KT+1,JB)            *BHR1(KT+1,ID)  /BHSUM
            TSSDH1(KT,JB)             = TSSDH1(KT+1,JB)          *BHR1(KT,ID)    /BHSUM
            TSSDH1(KT+1,JB)           = TSSDH1(KT+1,JB)          *BHR1(KT+1,ID)  /BHSUM
            CSSDH1(KT,CN(1:NAC),JB)   = CSSDH1(KT+1,CN(1:NAC),JB)*BHR1(KT,ID)    /BHSUM
            CSSDH1(KT+1,CN(1:NAC),JB) = CSSDH1(KT+1,CN(1:NAC),JB)*BHR1(KT+1,ID)  /BHSUM
          END IF
          DO I=IU,ID-1
            DX(KT,I) = DXI(JW)
            IF (INTERNAL_WEIR(KT,I)) DX(KT,I) = 0.0
          END DO
          IUT = IU
          IDT = ID-1
          IF (UP_HEAD(JB)) IUT = IU-1
          IF (DN_HEAD(JB)) IDT = ID
          DO I=IUT,IDT
            AZ(KT,I)  = AZMIN
            TKE(KT,I,1) = 1.25E-7 !sg 10/4/07
            TKE(KT,I,2) = 1.0E-9  !sg 10/4/07
            SAZ(KT,I) = AZMIN
            IF (INTERNAL_WEIR(KT,I)) THEN
              AZ(KT,I)  = 0.0
              TKE(KT,I,1) = 0.0 !sg  10/4/07
              TKE(KT,I,2) = 0.0 !sg  10/4/07
              SAZ(KT,I) = 0.0
            END IF
          END DO
          IF (CONSTITUENTS) THEN
            CALL TEMPERATURE_RATES
            CALL KINETIC_RATES
          END IF

!******** Upstream active segment

          IUT = US(JB)
          IF (SLOPE(JB) == 0.0) THEN
            DO I=US(JB),DS(JB)
              IF (KB(I)-KT < NL(JB)-1) IUT = I+1
            END DO
          ELSE
          ilayer=0
            DO I=US(JB)-1,DS(JB)+1
              IF (KB(I) > KBI(I)) THEN
                Bnew(KB(I),I)  = b(kb(i),i)                                                    ! SW 1/23/06                                           ! SW 3/2/05  
                DX(KB(I),I) = 0.0
                KB(I)       = KB(I)-1
                ilayer(i) = 1
                u(kb(i)+1,i)=0.0                                                               ! SW 1/23/06
                WRITE (WRN,'(2(A,I8),A,F0.3)') 'Raising bottom layer at segment ',I,' at iteration ',NIT,' at Julian day ',JDAY
              END IF
            END DO                   ! SW 1/23/06
            DO I=US(JB)-1,DS(JB)+1   ! SW 1/23/06
!              IF (KB(I)-KT < NL(JB)-1) IUT = I+1    ! SW 1/23/06
!                IF (I /= DS(JB)+1) KBMIN(I)   = MIN(KB(I),KB(I+1))                    ! SW 1/23/06                                 ! SW 3/2/05
                IF (I /= US(JB)-1) KBMIN(I-1) = MIN(KB(I-1),KB(I))                     ! SW 1/23/06  
                if(kbi(i) < kb(i))then 
                bkt(i)=bh1(kt,i)/(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))    ! SW 1/23/06
                depthb(ktwb(jw),i)=(h1(ktwb(jw),i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))    ! SW 1/23/06
                depthm(ktwb(jw),i)=(h1(ktwb(jw),i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))*0.5    ! SW 1/23/06
                avhr(kt,i)=(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))   +(H1(KT,I+1)-(el(kbi(i)+1,i+1)-el(kb(i)+1,i+1))&
                 -H1(KT,I)+(el(kbi(i)+1,i)-el(kb(i)+1,i)))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                                                           ! SW 1/23/06
                end if
            ENDDO
        DO I=US(JB)-1,DS(JB)+1   ! SW 1/23/06
         
         do k=kbmin(i)+1,kb(i)
         u(k,i)=0.0
         end do

         if(ilayer(i).eq.1.and.ilayer(i+1).eq.0)then
         bhrsum=0.0
         q(i)=0.0
          DO K=KT,KBMIN(I)
            IF (.NOT. INTERNAL_WEIR(K,I)) THEN
              BHRSUM = BHRSUM+BHR1(K,I)
              Q(I)   = Q(I)+U(K,I)*BHR1(K,I)
            END IF
          END DO
          DO K=KT,KBMIN(I)
            IF (INTERNAL_WEIR(K,I)) THEN
              U(K,I) = 0.0
            ELSE
              U(K,I) =  U(K,I)+(QC(I)-Q(I))/BHRSUM
            END IF
          END DO
          elseif(ilayer(i).eq.1.and.ilayer(i-1).eq.0)then
          bhrsum=0.0
          q(i-1)=0.0
          DO K=KT,KBMIN(I-1)
            IF (.NOT. INTERNAL_WEIR(K,I-1)) THEN
              BHRSUM = BHRSUM+BHR1(K,I-1)
              Q(I-1)   = Q(I-1)+U(K,I-1)*BHR1(K,I-1)
            END IF
          END DO
          DO K=KT,KBMIN(I-1)
            IF (INTERNAL_WEIR(K,I-1)) THEN
              U(K,I-1) = 0.0
            ELSE
              U(K,I-1) =  U(K,I-1)+(QC(I-1)-Q(I-1))/BHRSUM
            END IF
          END DO
          endif

            END DO
          END IF

!******** Segment addition

          IF (IUT /= IU) THEN
            IF (SNAPSHOT(JW)) WRITE (SNP(JW),'(/17X,2(A,I0))') ' Add segments ',IUT,' through ',IU-1
            DO I=IUT-1,IU-1
              Z(I)         =  Z(IU)
              KTI(I)       =  KTI(IU)
              H1(KT+1,I)   =  H(KT+1,JW)
              AVH1(KT+1,I) = (H1(KT+1,I)+H1(KT+2,I))*0.5
              AVHR(KT+1,I) =  H1(KT+1,I)+(H1(KT+1,I+1)-H1(KT+1,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
              IF (.NOT. TRAPEZOIDAL(JW)) THEN
                BH1(KT+1,I) =  B(KT+1,I)*H(KT+1,JW)
                H1(KT,I)    =  H(KT,JW)-Z(I)
! Check                BI(KT:KB(I),I) =  B(KT:KB(I),I)                                                                        ! SW 4/18/07
                BI(KT,I)    =  B(KTI(I),I)
                AVH1(KT,I)  = (H1(KT,I)+H1(KT+1,I))*0.5
                AVHR(KT,I)  =  H1(KT,I)+(H1(KT,I+1)-H1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                       !SW 07/29/04
                BH1(KT,I)   =  BI(KT,I)*(EL(KT,I)-Z(I)*COSA(JB)-EL(KTI(I)+1,I))/COSA(JB)
                IF (KTI(I) >= KB(I)) BH1(KT,I) = B(KT,I)*H1(KT,I)
                DO K=KTI(I)+1,KT
                  BH1(KT,I) = BH1(KT,I)+BH(K,I)
                END DO
              ELSE
                ELT = EL(KT,I)-Z(I)
                CALL GRID_AREA1 (ELT,EL(KT+1,I),BH1(KT,I),BI(KT,I))                                                    !SW 08/03/04
                BH1(KT+1,I) =  0.25*H(KT+1,JW)*(BB(KT,I)+2.*B(KT+1,I)+BB(KT+1,I))
                H1(KT,I)    =  H(KT,JW)-Z(I)
                AVH1(KT,I)  = (H1(KT,I)+H1(KT+1,I))*0.5
                AVHR(KT,I)  =  H1(KT,I)+(H1(KT,I+1)-H1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                       !SW 07/29/04
              END IF
              BKT(I)      = BH1(KT,I)/H1(KT,I)
              DEPTHB(K,I) = H1(KT,I)
              DEPTHM(K,I) = H1(KT,I)*0.5
              DO K=KT+1,KB(I)    
                DEPTHB(K,I) = DEPTHB(K-1,I)+ H1(K,I)
                DEPTHM(K,I) = DEPTHM(K-1,I)+(H1(K-1,I)+H1(K,I))*0.5
              END DO
            END DO
            DO I=IUT-1,IU-1
              BHR1(KT+1,I) = BH1(KT+1,I)+(BH1(KT+1,I+1)-BH1(KT+1,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                !SW 07/29/04
              BHR1(KT,I)   = BH1(KT,I)  +(BH1(KT,I+1)  -BH1(KT,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                !SW 07/29/04
            END DO
            DO I=IUT,IU-1
              ICE(I)   = ICE(IU)
              ICETH(I) = ICETH(IU)
              WIND2(I) = WIND2(IU)
              IF (DYNAMIC_SHADE(I)) CALL SHADING
              DO K=KT,KB(I)
                DX(K,I) = DXI(JW)
                U(K,I)  = U(K,IU)
                IF (INTERNAL_WEIR(K,I)) THEN
                  DX(K,I) = 0.0
                  U(K,I)  = 0.0
                END IF
                T1(K,I)           = T1(K,IU)
                T2(K,I)           = T1(K,IU)
                SU(K,I)           = U(K,IU)
                C1(K,I,CN(1:NAC)) = C1(K,IU,CN(1:NAC))
                C2(K,I,CN(1:NAC)) = C1(K,IU,CN(1:NAC))
                DO JE=1,NEP
                  EPD(K,I,JE) = 0.01
                  EPC(K,I,JE) = 0.01/H1(K,I)
                END DO
                CMBRT(CN(1:NAC),JB) = CMBRT(CN(1:NAC),JB)+C1(K,IU,CN(1:NAC))*DLX(I)*BH1(K,I)
                EBRI(JB)            = EBRI(JB)           +T1(K,IU)          *DLX(I)*BH1(K,I)
              END DO
              DO K=KT,KB(I)-1
                AZ(K,I)  = AZ(K,IU)
		        TKE(K,I,1) = TKE(K,IU,1) !sg 10/4/07
                TKE(K,I,2) = TKE(K,IU,2) !sg 10/4/07
                SAZ(K,I) = AZ(K,IU)
                IF (INTERNAL_WEIR(K,I)) THEN
                  AZ(K,I)  = 0.0
	              TKE(K,I,1) = 0.0 !sg 10/4/07
                  TKE(K,I,2) = 0.0 !sg 10/4/07
                  SAZ(K,I) = 0.0
                END IF     
              END DO
            END DO
! v3.5 start
!*********macrophytes
              do m=1,nmc
                if(macrophyte_calc(jw,m))then                                    
                  jt=kti(i)
                  je=kb(i)
                  wlelev=(el(kt,i)-z(i)*cosa(jb))                                     
                  do j=jt,je
                    if(j.lt.kt)then               
                      colb=el(j+1,i)
                    else
                      colb=el(kt+1,i)
                    end if
                    coldep=wlelev-colb               
                    macrc(j,kt,i,m)=macwbci(jw,m)
                    macrm(j,kt,i,m)=macrc(j,kt,i,m)*cw(j,i)*coldep*dlx(i)     
                    maCMBRT(JB,m) = maCMBRT(JB,m)+macrm(j,kt,i,m)                   
                  end do                                  
                  DO K=KT+1,KB(I)                 
                    jt=k                    
                    je=kb(i)                  
                    do j=jt,je
                      macrc(j,k,i,m)=macwbci(jw,m)
                      macrm(j,k,i,m)=macrc(j,k,i,m)*cw(j,i)*h2(k,i)*dlx(i)
                      maCMBRT(JB,m) = maCMBRT(JB,m)+macrm(j,k,i,m)
                    end do                                
                  END DO
                end if
              end do                 
! v3.5 end      
            U(KB(IUT):KB(IU),IU-1)  = 0.0
            SU(KB(IUT):KB(IU),IU-1) = 0.0
            ADX(KB(IUT):KB(IU),IU)  = 0.0
            IU                      = IUT
            CUS(JB)                 = IU
            IF (UH_EXTERNAL(JB)) KB(IU-1) = KB(IU)
            IF (UH_INTERNAL(JB)) THEN
              IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
                KB(IU-1) = MIN(KB(UHS(JB)),KB(IU))
              ELSE
                DO KKB=KT,KMX
                  IF (EL(KKB,IU) <= EL(KB(UHS(JB)),UHS(JB))) EXIT
                END DO
                KB(IU-1) = MIN(KKB,KB(IU))
              END IF
            END IF
            IF (UP_HEAD(JB)) THEN
              AZ(KT:KB(IU-1)-1,IU-1)  = AZMIN
              TKE(KT:KB(IU-1)-1,IU-1,1) = 1.25E-7 !sg 10/4/07
              TKE(KT:KB(IU-1)-1,IU-1,2) = 1.0E-9  !sg 10/4/07
              SAZ(KT:KB(IU-1)-1,IU-1) = AZMIN
            END IF
          END IF
          IF (CONSTITUENTS) THEN   
            CALL TEMPERATURE_RATES
            CALL KINETIC_RATES
          END IF     

!******** Total active cells and single layers

          DO I=IU,ID
            NTAC         = NTAC+1
            ONE_LAYER(I) = KTWB(JW) == KB(I)
          END DO
          NTACMX = MAX(NTAC,NTACMX)
        END DO
        CALL INTERPOLATION_MULTIPLIERS

!****** Additional layers

        ZMIN(JW) = -1000.0
        DO JB=BS(JW),BE(JW)
          DO I=CUS(JB),DS(JB)
            ZMIN(JW) = MAX(ZMIN(JW),Z(I))
          END DO
        END DO
        ADD_LAYER = ZMIN(JW) < -0.80*H(KT-1,JW) .AND. KT /= 2
      END DO

!**** Subtract layers

      DO WHILE (SUB_LAYER)
        IF (SNAPSHOT(JW)) WRITE (SNP(JW),'(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,1x,13("*"))') 'Subtract layer ',&
        KT,' at Julian day = ', JDAY,' NIT = ',NIT,' IZMIN =',IZMIN(JW)                                                                    ! SW 1/23/06

!****** Variable initialization

        KTWB(JW) = KTWB(JW)+1
        KT       = KTWB(JW)
        DO JB=BS(JW),BE(JW)
          IU = CUS(JB)
          ID = DS(JB)
          IF (CONSTITUENTS) DO1(KT-1,IU-1:ID+1) = 0.0
          DO I=IU-1,ID+1
            Z(I)         =  Z(I)-H(KT-1,JW)
            H1(KT-1,I)   =  H(KT-1,JW)
            H1(KT,I)     =  H(KT,JW)-Z(I)
            BI(KT,I)     =  B(KTI(I),I)
            BI(KT-1,I)   =  B(KT-1,I)
            AVH1(KT-1,I) = (H(KT-1,JW)+H(KT,JW))*0.5
            AVH1(KT,I)   = (H1(KT,I)+H1(KT+1,I))*0.5
            IF(.NOT. TRAPEZOIDAL(JW))THEN
              if(kb(i) >= kt)then                            ! sw 1/23/06
              BH1(KT,I)   =  BH1(KT-1,I)+BH1(KT,I)           ! sw 1/23/06
              BH1(KT-1,I) =  B(KT-1,I)*H(KT-1,JW)
              else  ! sw 1/23/06
              BH1(KT,I)   =  BH1(KT-1,I)         ! sw 1/23/06
              end if  ! sw 1/23/06
            ELSE
              ELT = EL(KT,I)-Z(I)
              CALL GRID_AREA1 (ELT,EL(KT+1,I),BH1(KT,I),DUMMY)                                                         !SW 08/03/04
              BH1(KT-1,I) = 0.25*H1(KT-1,JW)*(BB(KT-2,I)+2.*B(KT-1,I)+BB(KT-1,I))
            ENDIF
            VOL(KT,I)   =  BH1(KT,I)*DLX(I)
            VOL(KT-1,I) =  BH1(KT-1,I)*DLX(I)
            BKT(I)      =  BH1(KT,I)/H1(KT,I)
            if(kb(i) >= kt)then                             ! SW 1/23/06 
              U(KT,I)  = (U(KT-1,I)*BHR1(KT-1,I)+U(KT,I)*BHR(KT,I))/(BHR1(KT-1,I)+BHR(KT,I))
              T1(KT,I) = (T1(KT-1,I)*(BH1(KT,I)-BH(KT,I))+T1(KT,I)*BH(KT,I))/BH1(KT,I)
            ELSE
!              EBRI(JB) = EBRI(JB)-T1(KT,I)*VOL(KT,I)   1/23/06
              u(kt,i) = u(kt-1,i)    ! SW 1/23/06
              t1(kt,i) = t1(kt-1,i)  ! SW 1/23/06
            END IF
            if(kb(i) >= kt)then       ! sw 1/23/06
            C1(KT,I,CN(1:NAC))             = (C1(KT-1,I,CN(1:NAC))  *(BH1(KT,I)-BH(KT,I))+C1(KT,I,CN(1:NAC))  *BH(KT,I))/BH1(KT,I)
            CSSK(KT,I,CN(1:NAC))           = (CSSK(KT-1,I,CN(1:NAC))*(BH1(KT,I)-BH(KT,I))+CSSK(KT,I,CN(1:NAC))*BH(KT,I))/BH1(KT,I)
            else
            C1(KT,I,CN(1:NAC))             = C1(KT-1,I,CN(1:NAC))    ! SW 1/23/06
            CSSK(KT,I,CN(1:NAC))           = CSSK(KT-1,I,CN(1:NAC))   ! sw 1/23/06
            endif   ! sw 1/23/06
            CSSB(KT,I,CN(1:NAC))           =  CSSB(KT-1,I,CN(1:NAC))+CSSB(KT,I,CN(1:NAC))
            KF(KT,I,KFCN(1:NAF(JW),JW))    =  KF(KT-1,I,KFCN(1:NAF(JW),JW))
            KFS(KT,I,KFCN(1:NAF(JW),JW))   =  KFS(KT-1,I,KFCN(1:NAF(JW),JW))
            C1(KT-1,I,CN(1:NAC))           =  0.0
            C2(KT-1,I,CN(1:NAC))           =  0.0
            CSSB(KT-1,I,CN(1:NAC))         =  0.0
            CSSK(KT-1,I,CN(1:NAC))         =  0.0
            KF(KT-1,I,KFCN(1:NAF(JW),JW))  =  0.0
            KFS(KT-1,I,KFCN(1:NAF(JW),JW)) =  0.0
            DO JE=1,NEP
              if(kt <= kbi(i))then    ! CB 4/28/06
              EPM(KT,I,JE)   = EPM(KT-1,I,JE)+EPM(KT,I,JE)
              EPD(KT,I,JE)   = EPM(KT,I,JE)/((BI(KT,I)-BI(KT+1,I)+2.0*H1(KT,I))*DLX(I))
              EPC(KT,I,JE)   = EPM(KT,I,JE)/VOL(KT,I)
              EPM(KT-1,I,JE) = 0.0
              EPD(KT-1,I,JE) = 0.0
              EPC(KT-1,I,JE) = 0.0
              else                   ! SW 5/15/06
              EPM(KT,I,JE)   = EPM(KT-1,I,JE)
              EPD(KT,I,JE)   = EPD(KT-1,I,JE)
              EPC(KT,I,JE)   = EPC(KT-1,I,JE)
              EPM(KT-1,I,JE) = 0.0
              EPD(KT-1,I,JE) = 0.0
              EPC(KT-1,I,JE) = 0.0
              endif                   ! CB 4/28/06
            END DO
          END DO
! v3.5 start
          do i=iu,id
            do m=1,nmc
              if(macrophyte_calc(jw,m))then
                if(kticol(i))then
                  jt=kti(i)
                else
                  jt=kti(i)+1
                end if
                je=kb(i)
                wlelev=(el(kt,i)-z(i)*cosa(jb))           
                do j=jt,je
                  if(j.lt.kt)then                 
                    colb=el(j+1,i)
                  else
                    colb=el(kt+1,i)
                  end if
                  coldep=wlelev-colb
                  if(j.lt.kt)then
                    macrm(j,kt,i,m)=macrm(j,kt-1,i,m)
                  else
                    macrm(j,kt,i,m)=macrm(j,kt-1,i,m)+macrm(j,kt,i,m)
                  end if
                  if(cw(j,i).gt.0.0)then                      
                    macrc(j,kt,i,m)=macrm(j,kt,i,m)/(cw(j,i)*coldep*dlx(i))     
                  else
                    macrc(j,kt,i,m)=0.0
                  end if
                  macrm(j,kt-1,i,m)=0.0
                  macrc(j,kt-1,i,m)=0.0
                  
                end do
              end if
            end do                                
            jt=kti(i)
            je=kb(i)
            do j=jt,je
              mact(j,kt,i)=0.0
              mact(j,kt-1,i)=0.0
            end do            
            do m=1,nmc                  
              if(macrophyte_CALC(jw,m))then
                 do j=jt,je
                   mact(j,kt,i)=macrc(j,Kt,I,m)+mact(j,kt,i)
                 end do
              end if
            end do                
            do m=1,nmc
              tmac=0.0
              if(macrophyte_CALC(jw,m))then
                jt=kti(i)                  
                je=kb(i)                  
                do j=jt,je                                
                  tmac=tmac+macrm(j,kt,i,m)
                end do                
              end if
              mac(kt,i,m)=tmac/(bh1(kt,i)*dlx(i))
            end do                                        
            do m=1,nmc
              if(macrophyte_CALC(jw,m))then
                mac(kt-1,i,m)=0.0
              end if
            end do
          end do
! v3.5 end
          DO I=IU-1,ID
            AVHR(KT-1,I) =  H1(KT-1,I)+( H1(KT-1,I+1)- H1(KT-1,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
            AVHR(KT,I)   =  H1(KT,I)  +( H1(KT,I+1)  - H1(KT,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
            BHR1(KT,I)   = BH1(KT,I)  +(BH1(KT,I+1)  -BH1(KT,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
            BHR1(KT-1,I) = BH1(KT-1,I)+(BH1(KT-1,I+1)-BH1(KT-1,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
          END DO
          U(KT-1,IU-1:ID+1)     = 0.0
          W(KT-1,IU-1:ID+1)     = 0.0
          SW(KT-1,IU-1:ID+1)    = 0.0                                                                                  !TC 3/9/05
          P(KT-1,IU-1:ID+1)     = 0.0
          AZ(KT-1,IU-1:ID+1)    = 0.0
	      TKE(KT-1,IU-1:ID+1,1) = 0.0 !sg 10/4/07
	      TKE(KT-1,IU-1:ID+1,2) = 0.0 !sg 10/4/07
          DZ(KT-1,IU-1:ID+1)    = 0.0
          ADMZ(KT-1,IU-1:ID+1)  = 0.0
          ADZ(KT-1,IU-1:ID+1)   = 0.0
          DECAY(KT-1,IU-1:ID+1) = 0.0
          IF (UP_HEAD(JB)) THEN
            QUH1(KT,JB)             = QUH1(KT,JB)+QUH1(KT-1,JB)
            TSSUH1(KT,JB)           = TSSUH1(KT-1,JB)          +TSSUH1(KT,JB)
            CSSUH1(KT,CN(1:NAC),JB) = CSSUH1(KT-1,CN(1:NAC),JB)+CSSUH1(KT,CN(1:NAC),JB)
          END IF
          IF (DN_HEAD(JB)) THEN
            QDH1(KT,JB)             = QDH1(KT,JB)+QDH1(KT-1,JB)      
            TSSDH1(KT,JB)           = TSSDH1(KT-1,JB)          +TSSDH1(KT,JB)
            CSSDH1(KT,CN(1:NAC),JB) = CSSDH1(KT-1,CN(1:NAC),JB)+CSSDH1(KT,CN(1:NAC),JB)
          END IF

!******** Upstream active segment

          IUT = US(JB)
          IF (SLOPE(JB) /= 0.0) THEN
          ilayer=0   ! SW 1/23/06
            DO I=US(JB)-1,DS(JB)+1
              IF (KB(I) < KT ) THEN                                                                                  ! SR 10/17/05
                KB(I)                 = KT
                Bnew(KB(I),I)         = 0.000001   ! sw 1/23/06
                DX(KB(I),I)           = DXI(JW)
                ilayer(i)=1
                T1(KB(I),I)           = T1(KT,I)                   !    SW 5/15/06    T1(KB(I)-1,I) 
                C1(KB(I),I,CN(1:NAC)) = C1(KT,I,CN(1:NAC))         !    SW 5/15/06    C1(KB(I)-1,I,CN(1:NAC))
                WRITE (WRN,'(2(A,I8),A,F0.3,A,F0.3)') 'Lowering bottom segment ',I,' at iteration ',NIT,' at Julian day ',&
                                                       JDAY,' Z(I)=',Z(I)
              END IF
            ENDDO                    ! SW 1/23/06
            DO I=US(JB)-1,DS(JB)+1   ! SW 1/23/06
!               IF (I /= DS(JB)+1) KBMIN(I)   = MIN(KB(I),KB(I+1))  ! SW 1/23/06
                IF (I /= US(JB)-1) KBMIN(I-1) = MIN(KB(I-1),KB(I))  ! SW 1/23/06
                if(kbi(i) < kb(i))then 
                bkt(i)=bh1(kt,i)/(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))    ! SW 1/23/06
                depthb(ktwb(jw),i)=(h1(ktwb(jw),i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))    ! SW 1/23/06
                depthm(ktwb(jw),i)=(h1(ktwb(jw),i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))*0.5    ! SW 1/23/06
                avhr(kt,i)=(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))   +(H1(KT,I+1)-(el(kbi(i)+1,i+1)-el(kb(i)+1,i+1))&
                                           -H1(KT,I)+(el(kbi(i)+1,i)-el(kb(i)+1,i)))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                                                           ! SW 1/23/06
                end if
            ENDDO

        DO I=US(JB)-1,DS(JB)+1   ! SW 1/23/06
         
         if(ilayer(i).eq.1.and.ilayer(i+1).eq.0)then  ! SW 1/23/06
         bhrsum=0.0
         q(i)=0.0
          DO K=KT,KBMIN(I)
            IF (.NOT. INTERNAL_WEIR(K,I)) THEN
              BHRSUM = BHRSUM+BHR1(K,I)
              Q(I)   = Q(I)+U(K,I)*BHR1(K,I)
            END IF
          END DO
          DO K=KT,KBMIN(I)
            IF (INTERNAL_WEIR(K,I)) THEN
              U(K,I) = 0.0
            ELSE
              U(K,I) =  U(K,I)+(QC(I)-Q(I))/BHRSUM
            END IF
          END DO
          elseif(ilayer(i).eq.1.and.ilayer(i-1).eq.0)then
          bhrsum=0.0
          q(i-1)=0.0
          DO K=KT,KBMIN(I-1)
            IF (.NOT. INTERNAL_WEIR(K,I-1)) THEN
              BHRSUM = BHRSUM+BHR1(K,I-1)
              Q(I-1)   = Q(I-1)+U(K,I-1)*BHR1(K,I-1)
            END IF
          END DO
          DO K=KT,KBMIN(I-1)
            IF (INTERNAL_WEIR(K,I-1)) THEN
              U(K,I-1) = 0.0
            ELSE
              U(K,I-1) =  U(K,I-1)+(QC(I-1)-Q(I-1))/BHRSUM
            END IF
          END DO
          endif    ! SW 1/23/06
 
            END DO
          END IF
          DO I=US(JB),DS(JB)
            IF (KB(I)-KT < NL(JB)-1) IUT = I+1
            ONE_LAYER(I) = KTWB(JW) == KB(I)
          END DO
          IF (IUT > DS(JB)) THEN
            WRITE (W2ERR,'(A,I0/A,F0.2,2(A,I0))') 'Fatal error - insufficient segments in branch ',JB,'Julian day = ',JDAY,      &
                                                  ' at iteration ',NIT,' with water surface layer = ',KT
            WRITE (W2ERR,'(2(A,I0))')             'Minimum water surface located at segment ',IZMIN(JW),' with bottom layer at ',&
                                                   KB(IZMIN(JW))     
            TEXT = 'Runtime error - see w2.err'
            ERROR_OPEN = .TRUE.
            GO TO 230
          END IF

!******** Segment subtraction

          IF (IUT /= IU) THEN
            IF (SNAPSHOT(JW)) WRITE (SNP(JW),'(/17X,A,I0,A,I0)') ' Subtract segments ',IU,' through ',IUT-1
            DO I=IU,IUT-1
              DO K=KT,KB(I)
                EBRI(JB)            = EBRI(JB)-T1(K,I)*VOL(K,I)
                CMBRT(CN(1:NAC),JB) = CMBRT(CN(1:NAC),JB)-C1(K,I,CN(1:NAC))*VOL(K,I)+(CSSB(K,I,CN(1:NAC))+CSSK(K,I,CN(1:NAC))   &
                                      *VOL(K,I))*DLT
              END DO
            END DO
! v3.5 start
            do i=iu,iut-1
              do m=1,nmc
                if(macrophyte_calc(jw,m))then                               
                  jt=kti(i)
                  je=kb(i)
                  wlelev=(el(kt,i)-z(i)*cosa(jb))                 
                  do j=jt,je
                    if(j.lt.kt)then               
                      colb=el(j+1,i)
                    else
                      colb=el(kt+1,i)
                    end if                  
                    coldep=wlelev-colb               
!                    maCMBRT(JB,m) = maCMBRT(JB,m)-macrm(j,kt,i,m)+(macSS(j,Kt,I,m)*coldep*cw(j,i)*DLX(I))*DLT                               
                    maCMBRT(JB,m) = maCMBRT(JB,m)-macrm(j,kt,i,m)
                  end do                                 
                  DO K=KT+1,KB(I)                 
                    jt=k                    
                    je=kb(i)                
                    do j=jt,je                  
!                      maCMBRT(JB,m) = maCMBRT(JB,m)-macrm(j,k,i,m)+(macSS(j,K,I,m)*H2(k,i)*cw(j,i)*DmX(I))*DLT          
                      maCMBRT(JB,m) = maCMBRT(JB,m)-macrm(j,k,i,m)
                    end do                                
                  END DO                
                end if
              end do
            end do
! v3.5 end
            F(IU-1:IUT-1)     =  0.0
            Z(IU-1:IUT-1)     =  0.0
            ICETH(IU-1:IUT-1) =  0.0
            BHRHO(IU-1:IUT-1) =  0.0
            ICE(IU-1:IUT-1)   = .FALSE.
            DO K=KT,KB(IUT)
              ADX(K,IU-1:IUT-1)            = 0.0
              DX(K,IU-1:IUT-1)             = 0.0
              AZ(K,IU-1:IUT-1)             = 0.0
	          TKE(K,IU-1:IUT-1,1)          = 0.0 !SG  10/4/07
	          TKE(K,IU-1:IUT-1,2)          = 0.0 !SG  10/4/07
              SAZ(K,IU-1:IUT-1)            = 0.0
              U(K,IU-1:IUT-1)              = 0.0
              SU(K,IU-1:IUT-1)             = 0.0
              T1(K,IU-1:IUT-1)             = 0.0
              TSS(K,IU-1:IUT-1)            = 0.0
              QSS(K,IU-1:IUT-1)            = 0.0
              C1(K,IU-1:IUT-1,CN(1:NAC))   = 0.0
              C2(K,IU-1:IUT-1,CN(1:NAC))   = 0.0
              C1S(K,IU-1:IUT-1,CN(1:NAC))  = 0.0
              CSSB(K,IU-1:IUT-1,CN(1:NAC)) = 0.0
              CSSK(K,IU-1:IUT-1,CN(1:NAC)) = 0.0
            END DO
! v3.5 start
            do m=1,nmc
              IF (macrophyte_CALC(jw,m)) THEN              
                mac(k,i,m)=0.0
                mact(j,k,i)=0.0
              end if
            end do
            jt=kti(i)
            je=kb(i)
            do j=jt,je
              do m=1,nmc
                if(macrophyte_calc(jw,m))then
                  macrc(j,K,I,m) = 0.0                      
                end if
              end do
            end do
! v3.5 end

            IU           =  IUT
            CUS(JB)      =  IU
            Z(IU-1)      = (EL(KT,IU-1)-(EL(KT,IU)-Z(IU)*COSA(JB)))/COSA(JB)
            SZ(IU-1)     =  Z(IU)
            KTI(IU-1)    =  KTI(IU)
            IF (.NOT. TRAPEZOIDAL(JW)) THEN
              BI(KT,IU-1)  = B(KTI(IU-1),I)
              H1(KT,IU-1)  = H(KT,JW)-Z(IU-1)
              BH1(KT,IU-1) = Bnew(KTI(IU-1),IU-1)*(EL(KT,IU-1)-EL(KTI(IU-1)+1,IU-1)-Z(IU-1)*COSA(JB))/COSA(JB)   ! sw 1/23/06  Bnew(KTI(IU-1),IU-1)*(EL(KT,IU-1)-EL(KTI(IU-1)+1,IU-1)-Z(IU-1)*COSA(JB))/COSA(JB)     ! SR 10/17/05
              IF (KT >= KB(IU-1)) BH1(KT,IU-1) = Bnew(KT,IU-1)*H1(KT,IU-1)   ! sw 1/23/06 
              DO K=KTI(IU-1)+1,KT
                BH1(KT,IU-1) = BH1(KT,IU-1)+BH1(K,IU-1)
              END DO
            ELSE
              ELT = EL(KT,I)-Z(I)
              CALL GRID_AREA1 (ELT,EL(KT+1,IU-1),BH1(KT,IU-1),BI(KT,IU-1))                                             !SW 08/03/04
              BH1(KT,I) = 0.25*H(KT,JW)*(BB(KT-1,I)+2.*B(KT,I)+BB(KT,I))
              H1(KT,I)  = H(KT,JW)-Z(I)
            END IF
            BKT(IU-1)     =  BH1(KT,IU-1)/H1(KT,IU-1)
            BHR1(KT,IU-1) =  BH1(KT,IU-1)+(BH1(KT,IU)-BH1(KT,IU-1))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                 !SW 07/29/04
            IF (UH_EXTERNAL(JB)) KB(IU-1) = KB(IU)
            IF (UH_INTERNAL(JB)) THEN
              IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
                KB(IU-1) = MIN(KB(UHS(JB)),KB(IU))
              ELSE
                DO KKB = KT, KMX
                  IF (EL(KKB,IU) <= EL(KB(UHS(JB)),UHS(JB))) EXIT
                END DO
                KB(IU-1) = MIN(KKB,KB(IU))
              END IF
            END IF
          END IF
          IF (CONSTITUENTS) THEN      ! SW 5/15/06 
            CALL TEMPERATURE_RATES
            CALL KINETIC_RATES
          END IF     


!******** Total active cells

          DO I=IU,ID
            NTAC = NTAC-1
          END DO
          NTACMN = MIN(NTAC,NTACMN)
        END DO
        CALL INTERPOLATION_MULTIPLIERS

!****** Additional layer subtractions

        ZMIN(JW) = -1000.0
        DO JB=BS(JW),BE(JW)
          DO I=CUS(JB),DS(JB)
            ZMIN(JW) = MAX(ZMIN(JW),Z(I))
          END DO
        END DO
        SUB_LAYER = ZMIN(JW) > 0.60*H(KT,JW) .AND. KT < KTMAX                                                         ! SR 10/17/05
      END DO
    END DO

!** Temporary downstream head segment

    DO JB=1,NBR
      IF (DHS(JB) > 0) THEN
        DO JJB=1,NBR
          IF (DHS(JB) >= US(JJB) .AND. DHS(JB) <= DS(JJB)) EXIT
        END DO
        IF (CUS(JJB) > DHS(JB)) CDHS(JB) = CUS(JJB)
      END IF
    END DO

!***********************************************************************************************************************************
!*                                                    Task 2.6: Balances                                                          **
!***********************************************************************************************************************************

    QINT  = 0.0
    QOUTT = 0.0
    VOLSR = 0.0
    VOLTR = 0.0
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IF (VOLUME_BALANCE(JW)) THEN
          VOLSBR(JB) = VOLSBR(JB)+DLVOL(JB)
          VOLTBR(JB) = VOLEV(JB)+VOLPR(JB)+VOLTRB(JB)+VOLDT(JB)+VOLWD(JB)+VOLUH(JB)+VOLDH(JB)+VOLIN(JB)+VOLOUT(JB)
          VOLSR(JW)  = VOLSR(JW)+VOLSBR(JB)
          VOLTR(JW)  = VOLTR(JW)+VOLTBR(JB)
          QINT(JW)   = QINT(JW) +VOLIN(JB)+VOLTRB(JB)+VOLDT(JB)+VOLPR(JB)
          QOUTT(JW)  = QOUTT(JW)-VOLEV(JB)-VOLWD(JB) -VOLOUT(JB)
          IF (ABS(VOLSBR(JB)-VOLTBR(JB)) > VTOL .AND. VOLTBR(JB) > 100.0*VTOL) THEN
            IF (VOLUME_WARNING) THEN
              WRITE (WRN,'(A,F0.3,3(:/A,E15.8,A))') 'Computational warning at Julian day = ',JDAY,'spatial change  =', VOLSBR(JB), &
                                                    ' m^3','temporal change =',VOLTBR(JB),' m^3','volume error    =',              &
                                                     VOLSBR(JB)-VOLTBR(JB),' m^3'
              WARNING_OPEN   = .TRUE.
              VOLUME_WARNING = .FALSE.
            END IF
          END IF
        END IF
        IF (VOLSR(JW) /= 0.0) DLVR(JW) = (VOLTR(JW)-VOLSR(JW))/VOLSR(JW)*100.0
      END DO
      IF (ENERGY_BALANCE(JW)) THEN
        ESR(JW) = 0.0
        ETR(JW) = 0.0
        DO JB=BS(JW),BE(JW)
          ETBR(JB) = EBRI(JB)+TSSEV(JB)+TSSPR(JB)+TSSTR(JB)+TSSDT(JB)+TSSWD(JB)+TSSUH(JB)+TSSDH(JB)+TSSIN(JB)+TSSOUT(JB)+TSSS(JB)  &
                     +TSSB(JB)+TSSICE(JB)
          ESBR(JB) = 0.0
          DO I=CUS(JB),DS(JB)
            DO K=KT,KB(I)
              ESBR(JB) = ESBR(JB)+T1(K,I)*DLX(I)*BH1(K,I)
            END DO
          END DO
          ETR(JW) = ETR(JW)+ETBR(JB)
          ESR(JW) = ESR(JW)+ESBR(JB)
        END DO
      END IF
      IF (MASS_BALANCE(JW)) THEN
        DO JB=BS(JW),BE(JW)
          DO JC=1,NAC
            CMBRS(CN(JC),JB) = 0.0
            DO I=CUS(JB),DS(JB)
              DO K=KT,KB(I)
                CMBRS(CN(JC),JB) = CMBRS(CN(JC),JB)+C1(K,I,CN(JC))*DLX(I)*BH1(K,I)
                CMBRT(CN(JC),JB) = CMBRT(CN(JC),JB)+(CSSB(K,I,CN(JC))+CSSK(K,I,CN(JC))*BH1(K,I)*DLX(I))*DLT
              END DO
            END DO
          END DO
! v3.5 start          
          do m=1,nmc
            if(macrophyte_calc(jw,m))then
              maCMBRS(JB,m) = 0.0
              DO I=CUS(JB),DS(JB)
                if(kticol(i))then
                  jt=kti(i)
                else
                  jt=kti(i)+1
                end if
                je=kb(i)
                wlelev=(el(kt,i)-z(i)*cosa(jb))               
                do j=jt,je
                  if(j.lt.kt)then                 
                    colb=el(j+1,i)
                  else
                     colb=el(kt+1,i)
                  end if
                  coldep=wlelev-colb  
                  maCMBRS(JB,m) = maCMBRS(jb,m)+macrm(j,Kt,I,m)
                  maCMBRT(JB,m) = maCMBRT(JB,m)+(macSS(j,Kt,I,m)*coldep*cw(j,i)*DLX(I))*DLT       
                end do                          
                DO K=KT+1,KB(I)               
                  jt=k                    
                  je=kb(i)              
                  do j=jt,je
                    maCMBRS(JB,m) =maCMBRS(JB,m)+macrm(j,K,I,m)
!                    maCMBRT(JB,m) = maCMBRT(JB,m)+(macSS(j,K,I,m)*H2(k,i)*cw(j,i)*DLX(I))*DLT     
                    maCMBRT(JB,m) = maCMBRT(JB,m)+(macSS(j,K,I,m)*(cw(j,i)/b(k,i))*bh1(k,i)*DLX(I))*DLT     
                  end do                                  
                END DO
              END DO
            end if
          end do
! v3.5 end
        END DO
      END IF
    END DO

!***********************************************************************************************************************************
!*                                       Task 2.7: Variable updates for next timestep                                             **
!***********************************************************************************************************************************

    SZ     = Z
    SKTI   = KTI
    SBKT   = BKT
    VOLUH2 = QUH1  *DLT
    VOLDH2 = QDH1  *DLT
    TSSUH2 = TSSUH1*DLT
    TSSDH2 = TSSDH1*DLT
    DO JW=1,NWB
     KT = KTWB(JW)
      ELKT(JW) = EL(KT,DS(BS(JW)))-Z(DS(BS(JW)))*COSA(BS(JW))
      DO JB=BS(JW),BE(JW)
        ELWS(CUS(JB):DS(JB)+1) = EL(KT,CUS(JB):DS(JB)+1)-Z(CUS(JB):DS(JB)+1)*COSA(JB)
        DO I=US(JB)-1,DS(JB)
          AVHR(KT,I) = H1(KT,I)+(H1(KT,I+1)-H1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                               !SW 07/29/04
        END DO
        AVHR(KT,DS(JB)+1)=H1(KT,DS(JB)+1)                                                                              !SW 03/08/05
        DO I=CUS(JB)-1,DS(JB)+1
          DO K=KTWB(JW),KB(I)
            QSS(K,I)   = 0.0
            TSS(K,I)   = 0.0
            SU(K,I)    = U(K,I)
            SW(K,I)    = W(K,I)
            T2(K,I)    = T1(K,I)
            SAZ(K,I)   = AZ(K,I)
            H2(K,I)    = H1(K,I)
            BH2(K,I)   = BH1(K,I)
            BHR2(K,I)  = BHR1(K,I)
            AVH2(K,I)  = AVH1(K,I)
            SAVH2(K,I) = AVH2(K,I)
            SAVHR(K,I) = AVHR(K,I)
          END DO
        END DO
      END DO
    END DO
    IF (CONSTITUENTS) THEN
      CSSUH2 = CSSUH1*DLT
      CSSDH2 = CSSDH1*DLT
      DO JW=1,NWB
        KT = KTWB(JW)
        DO JB=BS(JW),BE(JW)
          DO JC=1,NAC
            DO I=US(JB)-1,DS(JB)+1
              DO K=KT,KB(I)
                DO JE=1,NEP
                  IF (EPIPHYTON_CALC(JW,JE)) EPD(K,I,JE) = MAX(EPD(K,I,JE),0.0)
                END DO
! v3.5 start
                IF (SEDIMENT_CALC(JW))then
                  SED(K,I) = MAX(SED(K,I),0.0)                  
                  SEDp(K,I) = MAX(SEDp(K,I),0.0)
                  SEDn(K,I) = MAX(SEDn(K,I),0.0)
                  SEDc(K,I) = MAX(SEDc(K,I),0.0)                  
                end if
! v3.5 end
                CSSB(K,I,CN(JC)) = 0.0
                C1S(K,I,CN(JC))  = C1(K,I,CN(JC))
                C2(K,I,CN(JC))   = MAX(C1(K,I,CN(JC)),0.0)
              END DO
            END DO
          END DO
        END DO
      END DO
    END IF
! v3.5 start
    do jw = 1,nwb
      KT = KTWB(JW)
      do m=1,nmc
        IF (macrophyte_CALC(jw,m)) THEN
          do jb=bs(jw),be(jw)
            DO I=US(JB),DS(JB)
              DO K=kt,kb(i)
                smac(k,i,m)=mac(k,i,m)             
                do j=1,kmx
                  smacrc(j,K,I,m)=macrc(j,K,I,m)
                  smacrm(j,K,I,m)=macrm(j,K,I,m)
                end do
              END DO
            END DO          
          end do
        end if
      end do
    end do
! v3.5 end
  do jw = 1,nwb !mlm
    IF (ULTIMATE(JW)) THEN   ! SR 5/15/06
      DO K=2,KMX
      RATZ(K,JW)  =  AVH2(K-1,DS(BE(JW)))/AVH2(K,DS(BE(JW)))                                         ! SW 5/20/05
      CURZ1(K,JW) =  2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))+AVH2(K,DS(BE(JW))))/AVH2(K-1,DS(BE(JW)))   ! SW 5/20/05
      CURZ2(K,JW) = -2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))*AVH2(K,DS(BE(JW))))                        ! SW 5/20/05
      CURZ3(K,JW) =  2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))+AVH2(K,DS(BE(JW))))/AVH2(K,DS(BE(JW)))     ! SW 5/20/05
      END DO
    END IF
    END DO
    NIT     =  NIT+1
    ELTM    =  ELTM+DLT
    ELTMS   =  ELTMS+DLT
    ELTMF   =  ELTMF+DLT
    JDAY    =  ELTM/DAY
    ELTMJD  =  JDAY-TMSTRT
    END_RUN =  JDAY >= TMEND
    DLT     =  MAX(DLTMIN,DLTF(DLTDP)*CURMAX)
    DLT     =  MIN(DLT,1.1*DLTS)
    DLTAV   = (ELTM-TMSTRT*DAY)/NIT
    IF (DLT <  MINDLT) THEN
      MINDLT = DLTS
      JDMIN  = JDAY
    END IF
    IF (JDAY >= DLTD(DLTDP+1)) DLTDP = DLTDP+1
    IF (DLT  >  DLTMAX(DLTDP)) DLT   = DLTMAX(DLTDP)
    CURMAX = DLTMAX(DLTDP)/DLTF(DLTDP)
    IF (INT(JDAY) == JDAYNX) THEN
      JDAYG  = JDAYG+1
      JDAYNX = JDAYNX+1
    END IF
    IF (JDAYG > 300) WINTER = .TRUE.
    IF (JDAYG < 40)  WINTER = .FALSE.
    WRITE (GDCH,'(I3)') GDAY
    CALL GREGORIAN_DATE
    UPDATE_KINETICS = .FALSE.
    IF (MOD(NIT,CUF) == 0) UPDATE_KINETICS = .TRUE.

!***********************************************************************************************************************************
!*                                                    Task 2.8: Output Results                                                    **
!***********************************************************************************************************************************

! v3.5 start  - keep this code for error checking, but delete later
!   if(macrophyte_on)then
!     if(nit.eq.1)open(689,file='totmac.dat',status='unknown')
!     if((nit/200)*200.eq.nit)then        
!       bmass=0.0
!       totar=0.0
!       DO Jw=1,Nwb
!         KT = KTwb(Jw)
!         DO JB=BS(Jw),BE(Jw)
!           IU = CUS(JB)
!           ID = DS(JB)
!           DO I=IU,ID
!             smass=0.0
!             do k=kt,kb(i)
!               do m=1,nmc               
!                 bmass=bmass+mac(k,i,m)*bh2(k,i)*dlx(i)
!                 smass=smass+mac(k,i,m)*bh2(k,i)*dlx(i)                         
!               end do
!             end do
!             totar=totar+dlx(i)*b(kti(i),i) 
!             armac(i)=smass/(dlx(i)*b(kti(i),i))
!           end do
!         end do
!       end do
!       bmassar=bmass/totar
!       write(689,'(f10.4,2f10.3)')jday,bmass/1000.0,bmassar
!     end if
!   end if
! v3.5 end

    DO JW=1,NWB
    
!**** Inactive segments

      JB       = BS(JW)
      NBL(JW)  = 1
      IBPR(JW) = 1
      DO I=1,NISNP(JW)-1
        IF (CUS(JB) > ISNP(I,JW)) THEN
          BL(NBL(JW),JW) = I
          NBL(JW)        = NBL(JW)+1
          IBPR(JW)       = I+1
        END IF
        IF (ISNP(I+1,JW) > DS(JB)) JB = JB+1
      END DO
      NBL(JW) = NBL(JW)-1

!**** Snapshots

      IF (SNAPSHOT(JW)) THEN
        IF (JDAY >= NXTMSN(JW) .OR. JDAY >= SNPD(SNPDP(JW)+1,JW)) THEN
          IF (JDAY >= SNPD(SNPDP(JW)+1,JW)) THEN
            SNPDP(JW)  = SNPDP(JW)+1
            NXTMSN(JW) = SNPD(SNPDP(JW),JW)
          END IF
          NXTMSN(JW) = NXTMSN(JW)+SNPF(SNPDP(JW),JW)
          WRITE (SNP(JW),10490) (TITLE(J),J=1,10)
          WRITE (SNP(JW),10500) 'Time Parameters',MONTH,GDAY,YEAR,INT(JDAY),(JDAY-INT(JDAY))*24.0,INT(ELTMJD),                     &
                                (ELTMJD-INT(ELTMJD))*24.0,INT(DLTS1),KLOC,ILOC,INT(MINDLT),INT(JDMIN),(JDMIN-INT(JDMIN))*24.0,     &
                                 KMIN,IMIN
          IF (LIMITING_DLT(JW))  WRITE (SNP(JW),10510) KMIN,IMIN
          WRITE (SNP(JW),10520)  INT(DLTAV),NIT,NV
          WRITE (SNP(JW),10530) 'Meteorological Parameters'
          WRITE (SNP(JW),10540)  TAIR(JW),DEG,TDEW(JW),DEG,PHI(JW),CLOUD(JW),ET(DS(1)),DEG,CSHE(DS(1)),SRON(JW),DEG
          WRITE (SNP(JW),10550) 'Inflows','Upstream inflows'
          DO JB=BS(JW),BE(JW)
            IF (UP_FLOW(JB)) WRITE (SNP(JW),10560) JB,KTQIN(JB),KBQIN(JB),QIN(JB),TIN(JB),DEG
          END DO
          DO JB=BS(JW),BE(JW)
            IF (DIST_TRIBS(JB)) THEN
              WRITE (SNP(JW),10570)
              WRITE (SNP(JW),10580) JB,QDTR(JB),TDTR(JB),DEG
            END IF
          END DO
          IF (TRIBUTARIES) THEN
            WRITE (SNP(JW),10590) (ITR(JT),          JT=1,JTT)
            WRITE (SNP(JW),10600) (KTTR(JT),KBTR(JT),JT=1,JTT)
            WRITE (SNP(JW),10610) (QTR(JT),          JT=1,JTT)
            WRITE (SNP(JW),10620) (TTR(JT),          JT=1,JTT)
          END IF
          WRITE (SNP(JW),10630)
          DO JB=BS(JW),BE(JW)
            IF (DN_FLOW(JB)) THEN
              WRITE (SNP(JW),10640)  JB,(QSTR(JS,JB),JS=1,JSS(JB))
              WRITE (SNP(JW),10650)  QSUM(JB),(K,K=KTWB(JW),KB(DS(JB)))
              WRITE (SNP(JW),10660) (QOUT(K,JB), K=KTWB(JW),KB(DS(JB)))
            END IF
          END DO
          IF (WITHDRAWALS) THEN
            DO JWD=1,JWW
              WRITE (SNP(JW),10670) MAX(CUS(JBWD(JWD)),IWD(JWD)),QWD(JWD)
              IF (QWD(JWD) /= 0.0) THEN
                WRITE (SNP(JW),10680) (K,         K=KTW(JWD),KBW(JWD))
                WRITE (SNP(JW),10690) (QSW(K,JWD),K=KTW(JWD),KBW(JWD))
              ELSE
                WRITE (SNP(JW),10680)
                WRITE (SNP(JW),10690) QWD(JWD)
              END IF
            END DO
          END IF
          IF (CONSTITUENTS) THEN
            WRITE (SNP(JW),10700) 'Constituent Inflow Concentrations'
            DO JB=BS(JW),BE(JW)
              IF (UP_FLOW(JB) .AND. NACIN(JB) > 0)    WRITE (SNP(JW),10710) JB,(CNAME1(INCN(JC,JB))(1:18),CIN(INCN(JC,JB),JB),     &
                                                                            CUNIT2(INCN(JC,JB)),JC=1,NACIN(JB))
              IF (DIST_TRIBS(JB) .AND. NACDT(JB) > 0) WRITE (SNP(JW),10730) JB,(CNAME1(DTCN(JC,JB))(1:18),CDTR(DTCN(JC,JB),JB),    &
                                                                            CUNIT2(DTCN(JC,JB)),JC=1,NACDT(JB))     
            END DO
            DO JT=1,NTR
              IF (NACTR(JT) > 0) WRITE (SNP(JW),10720) JT,(CNAME1(TRCN(JC,JT))(1:18),CTR(TRCN(JC,JT),JT),CUNIT2(TRCN(JC,JT)),      &
                                                       JC=1,NACTR(JT))
            END DO
          END IF
          IF (EVAPORATION(JW) .OR. PRECIPITATION(JW)) WRITE (SNP(JW),10740)
          IF (EVAPORATION(JW)) THEN
            WRITE (SNP(JW),10750) (JB,EVBR(JB),JB=BS(JW),BE(JW))                                                              ! SW 9/15/05
            WRITE (SNP(JW),10755) (JB,-VOLEV(JB),JB=BS(JW),BE(JW))
          END IF
          IF (PRECIPITATION(JW)) WRITE (SNP(JW),10760) (JB,PR(JB),JB=BS(JW),BE(JW))
          IF (HEAD_BOUNDARY(JW)) THEN
            WRITE (SNP(JW),10770)
            DO JB=BS(JW),BE(JW)
              IF (UH_EXTERNAL(JB)) WRITE (SNP(JW),10780) JB,ELUH(JB)
              IF (DH_EXTERNAL(JB)) WRITE (SNP(JW),10790) JB,ELDH(JB)
            END DO
          END IF
          IF (VOLUME_BALANCE(JW)) THEN
            WRITE (SNP(JW),10800)
            WRITE (SNP(JW),10810) JW,VOLSR(JW),VOLTR(JW),VOLTR(JW)-VOLSR(JW),DLVR(JW)
            DO JB=BS(JW),BE(JW)
              IF (VOLSBR(JB) /= 0.0) DLVBR = (VOLTBR(JB)-VOLSBR(JB))/VOLSBR(JB)
              WRITE (SNP(JW),10820) JB,VOLSBR(JB),VOLTBR(JB),VOLTBR(JB)-VOLSBR(JB),DLVBR*100.0
            END DO
          END IF
          IF (ENERGY_BALANCE(JW)) THEN
            WRITE (SNP(JW),10830)
            IF (ESR(JW) /= 0.0) DLE = (ESR(JW)-ETR(JW))/ESR(JW)
            WRITE (SNP(JW),10840) JW,ESR(JW)*4.184E3,ETR(JW)*4.184E3,(ESR(JW)-ETR(JW))*4.184E3,DLE*100.0
            DO JB=BS(JW),BE(JW)
              WRITE (SNP(JW),10870) JB
              IF (ESBR(JB) /= 0.0) DLE = (ESBR(JB)-ETBR(JB))/ESBR(JB)
              WRITE (SNP(JW),10850) ESBR(JB)*4.184E3,ETBR(JB)*4.1843E3,(ESBR(JB)-ETBR(JB))*4.1843E3,DLE*100.0
            END DO
          END IF
          IF (MASS_BALANCE(JW)) THEN
            WRITE (SNP(JW),10860)
            DO JB=BS(JW),BE(JW)
              WRITE (SNP(JW),10870) JB
              DO JC=1,NAC
                IF (CMBRS(CN(JC),JB) /= 0.0) DLMR = (CMBRT(CN(JC),JB)-CMBRS(CN(JC),JB))/(CMBRS(CN(JC),JB)+NONZERO)*100.0
                WRITE (SNP(JW),10880) CNAME1(CN(JC)),CMBRS(CN(JC),JB),CUNIT1(CN(JC)),CMBRT(CN(JC),JB),CUNIT1(CN(JC)),              &
                                     (CMBRT(CN(JC),JB)-CMBRS(CN(JC),JB)),CUNIT1(CN(JC)),DLMR
              END DO
! v3.5 start
              DO m=1,Nmc
                if(macrophyte_calc(jw,m))then
                  IF (maCMBRS(JB,m).NE.0.0) THEN
                    DLMR = (maCMBRT(JB,m)-maCMBRS(JB,m))/(maCMBRS(JB,m)+NONZERO)                           
                  END IF
                  WRITE (SNP(Jw),3312) m,maCMBRS(JB,m),maCMBRT(JB,m),(maCMBRT(JB,m)-maCMBRS(JB,m)),DLMR*100.0     
3312  FORMAT(5X,'Macrophyte spec ',i2,/7X,'Spatially integrated mass [MACMBRS] = ',1PE15.8E2,1X,'g ',/7X,   &
                'Temporally integrated mass [MACMBRT] = ',1PE15.8E2,1X,'g ',/7X,'Mass error                         = ',  &
                 1PE15.8E2,1X,'g ',/7X,'Percent error                      = ',1PE15.8E2,' %')     
                end if
              END DO
! v3.5 end
            END DO
          END IF
          WRITE (SNP(JW),10890) 'Geometry',KTWB(JW),ELKT(JW)
          WRITE (SNP(JW),10900) (JB,CUS(JB),JB=BS(JW),BE(JW))
          CALL OUTPUT (JDAY,IBPR(JW),NISNP(JW),KBR(JW),ISNP,BL(1,JW),NBL(JW))
        END IF
      END IF

!**** Vertical profiles

      IF (PROFILE(JW)) THEN
        IF (JDAY >= NXTMPR(JW) .OR. JDAY >= PRFD(PRFDP(JW)+1,JW)) THEN
          IF (JDAY >= PRFD(PRFDP(JW)+1,JW)) THEN
            PRFDP(JW)  = PRFDP(JW)+1
            NXTMPR(JW) = PRFD(PRFDP(JW),JW)
          END IF
          NXTMPR(JW) = NXTMPR(JW)+PRFF(PRFDP(JW),JW)
          NSPRF(JW)  = NSPRF(JW)+1
          WRITE (PRF(JW),'(F8.3,1X,A3,I3,A,2I4,F8.4,I8)')JDAY,ADJUSTL(MONTH),GDAY,', ',YEAR,KTWB(JW),SNGL(Z(DS(BS(JW)))),NSPRF(JW)
          DO JP=1,NIPRF(JW)
            NRS = KB(IPRF(JP,JW))-KTWB(JW)+1
            WRITE (PRF(JW),'(A8,I4/(8F10.2))') 'TEMP    ',NRS,(T2(K,IPRF(JP,JW)),K=KTWB(JW),KB(IPRF(JP,JW)))
          END DO
          DO JC=1,NAC
            IF (PRINT_CONST(CN(JC),JW)) THEN
              DO JP=1,NIPRF(JW)
                NRS = KB(IPRF(JP,JW))-KTWB(JW)+1
                WRITE (PRF(JW),'(A8,I4/(8(E13.6,2X)))') ADJUSTL(CNAME2(CN(JC))),NRS,(C2(K,IPRF(JP,JW),CN(JC))*CMULT(CN(JC)),       &     ! CB 1/24/05
                                                            K=KTWB(JW),KB(IPRF(JP,JW)))
              END DO
            END IF
          END DO
          IF (CONSTITUENTS) THEN
            DO JD=1,NACD(JW)
              DO JP=1,NIPRF(JW)
                NRS = KB(IPRF(JP,JW))-KTWB(JW)+1
                WRITE (PRF(JW),'(A8,I4/(8(E13.6,2X)))') ADJUSTL(CDNAME2(CDN(JD,JW))),NRS,(CD(K,IPRF(JP,JW),CDN(JD,JW))             &      ! CB 1/24/05
                                                           *CDMULT(CDN(JD,JW)),K=KTWB(JW), KB(IPRF(JP,JW)))
              END DO
            END DO
          END IF
        END IF
      END IF

!**** Spreadsheet

      IF (SPREADSHEET(JW)) THEN
        IF (JDAY >= NXTMSP(JW) .OR. JDAY >= SPRD(SPRDP(JW)+1,JW)) THEN
          IF (JDAY >= SPRD(SPRDP(JW)+1,JW)) THEN
            SPRDP(JW)  = SPRDP(JW)+1
            NXTMSP(JW) = SPRD(SPRDP(JW),JW)
          END IF
          CONV1      = BLANK1
          NXTMSP(JW) = NXTMSP(JW)+SPRF(SPRDP(JW),JW)
          DO J=1,NISPR(JW)
            KBMAX(JW) = MAX(KB(ISPR(J,JW)),KBMAX(JW))
            DO K=KTWB(JW),KB(ISPR(J,JW))
              WRITE (CONV1(K,J),'(F10.2)') T2(K,ISPR(J,JW))
            END DO
          END DO
          DO K=KTWB(JW),KBMAX(JW)
            WRITE (SPR(JW),'(A18,20X,2F10.3,1000(F10.3,A))') 'Temperature       ',JDAY,-DEPTHM(K,DS(BS(JW))),                      &
                                                             (ELWS(ISPR(J,JW))-DEPTHM(K,ISPR(J,JW)),CONV1(K,J),J=1,NISPR(JW))
          END DO
          DO JC=1,NAC
            IF (PRINT_CONST(CN(JC),JW)) THEN
              DO J=1,NISPR(JW)
                DO K=KTWB(JW),KB(ISPR(J,JW))
                  WRITE (CONV1(K,J),'(F10.3)') C2(K,ISPR(J,JW),CN(JC))*CMULT(CN(JC))                                                ! SW 8/13/06
                END DO
              END DO
              DO K=KTWB(JW),KBMAX(JW)
                WRITE (SPR(JW),'(A38,2F10.3,1000(F10.3,A))') CNAME3(CN(JC)),JDAY,-DEPTHM(K,DS(BS(JW))),                            &
                                                            (ELWS(ISPR(J,JW))-DEPTHM(K,ISPR(J,JW)),CONV1(K,J),J=1,NISPR(JW))
              END DO
            END IF
          END DO
          IF (CONSTITUENTS) THEN
            DO JD=1,NACD(JW)
              IF (PRINT_DERIVED(CDN(JD,JW),JW)) THEN
                DO J=1,NISPR(JW)
                  DO K=KTWB(JW),KB(ISPR(J,JW))
                    WRITE (CONV1(K,J),'(F10.3)') CD(K,ISPR(J,JW),CDN(JD,JW))*CDMULT(CDN(JD,JW))                                     ! SW 8/13/06
                  END DO
                END DO
                DO K=KTWB(JW),KBMAX(JW)
                  WRITE (SPR(JW),'(A38,2F10.3,1000(F10.3,A))') CDNAME3(CDN(JD,JW)),JDAY,-DEPTHM(K,DS(BS(JW))),                     &
                                                              (ELWS(ISPR(J,JW))-DEPTHM(K,ISPR(J,JW)),CONV1(K,J),J=1,NISPR(JW))
                END DO
              END IF
            END DO
          END IF
        END IF
      END IF

!**** Velocity vectors

      IF (VECTOR(JW)) THEN
        IF (JDAY >= NXTMVP(JW) .OR. JDAY >= VPLD(VPLDP(JW)+1,JW)) THEN
          IF (JDAY >= VPLD(VPLDP(JW)+1,JW)) THEN
            VPLDP(JW)  = VPLDP(JW)+1
            NXTMVP(JW) = VPLD(VPLDP(JW),JW)
          END IF
          NXTMVP(JW) = NXTMVP(JW)+VPLF(VPLDP(JW),JW)
          WRITE (VPL(JW),*)  'New date ',JDAY,MONTH//GDCH//',',YEAR,KTWB(JW),(US(JB),JB=BS(JW),BE(JW))
          WRITE (VPL(JW),*) ((Z(I)*COSA(BS(JW))),     I=US(BS(JW)),DS(BE(JW)))
          WRITE (VPL(JW),*) ((EL(K,I),K=KTWB(JW),KMX),I=US(BS(JW)),DS(BE(JW)))
          WRITE (VPL(JW),*) ((U(K,I), K=KTWB(JW),KMX),I=US(BS(JW)),DS(BE(JW)))
          WRITE (VPL(JW),*) ((W(K,I), K=KTWB(JW),KMX),I=US(BS(JW)),DS(BE(JW)))
        END IF
      END IF

!**** Contours

      IF (CONTOUR(JW)) THEN
        IF (JDAY >= NXTMCP(JW) .OR. JDAY >= CPLD(CPLDP(JW)+1,JW)) THEN
          IF (JDAY >= CPLD(CPLDP(JW)+1,JW)) THEN
            CPLDP(JW)  = CPLDP(JW)+1
            NXTMCP(JW) = CPLD(CPLDP(JW),JW)
          END IF
          NXTMCP(JW) = NXTMCP(JW)+CPLF(CPLDP(JW),JW)
          WRITE (CPL(JW),'(A,F12.4,5X,A9,5X,I2,5X,I4)') 'New date ',JDAY,MONTH,GDAY,YEAR
          WRITE (CPL(JW),'(9(I8,2X))')                   KTWB(JW)
          WRITE (CPL(JW),'(9(E13.6,2X))')               (QTR(JT),JT=1,NTR)
          WRITE (CPL(JW),'(9(E13.6,2X))')               (TTR(JT),JT=1,NTR)
          DO JT=1,NTR
            DO JAC=1,NACTR(JT)
              IF (PRINT_CONST(TRCN(JAC,JT),JW)) WRITE (CPL(JW),'(9(E13.6,2X))') CTR(TRCN(JAC,JT),JT)
            END DO
          END DO
          DO JB=BS(JW),BE(JW)
            WRITE (CPL(JW),'(9(I8,2X))')             CUS(JB)
            WRITE (CPL(JW),'(9(E13.6,2X))')          QIN(JB), QSUM(JB)
            DO I=CUS(JB),DS(JB)
              WRITE (CPL(JW),'(A38/(9(E13.6,2X)))') 'BHR', (BHR1(K,I),K=KTWB(JW)+1,KB(I))
            END DO
            DO I=CUS(JB),DS(JB)
              WRITE (CPL(JW),'(A38/(9(E13.6,2X)))') 'U',   (U(K,I),   K=KTWB(JW),KB(I))
            END DO
            WRITE (CPL(JW),'(A38/(9(E13.6,2X)))')   'QC',  (QC(I),    I=CUS(JB),DS(JB))
            WRITE (CPL(JW),'(A38/(9(E13.6,2X)))')   'Z',   (Z(I),     I=CUS(JB),DS(JB))
            WRITE (CPL(JW),'(A38/(9(I8,2X)))')   'KTI',   (kti(I),     I=CUS(JB),DS(JB))  ! v3.5
            DO I=CUS(JB),DS(JB)
              WRITE (CPL(JW),'(A38/(9(E13.6,2X)))') 'Temperature',(T2(K,I),K=KTWB(JW),KB(I))
            END DO
            DO JC=1,NAC
              IF (PRINT_CONST(CN(JC),JW)) THEN
                DO I=CUS(JB),DS(JB)
                  WRITE (CPL(JW),'(A38/(9(E13.6,2X)))') CNAME(CN(JC)),(C2(K,I,CN(JC))*CMULT(CN(JC)),K=KTWB(JW),KB(I))
                END DO
              END IF
            END DO
            DO JE=1,NEP
              DO I=CUS(JB),DS(JB)
                IF (PRINT_EPIPHYTON(JW,JE)) WRITE (CPL(JW),'(A38/(9(E13.6,2X)))') 'Epiphyton',(EPD(K,I,JE),K=KTWB(JW),KB(I))
              END DO
            END DO
! v3.5 start
            if(print_sediment(jw))then
              DO I=CUS(JB),DS(JB)
                WRITE (CPL(Jw),'(A38/(9(E13.6,2X)))')'Sediment',(seD(K,I),K=KTWB(JW),KB(I))     
              end do
              DO I=CUS(JB),DS(JB)
                WRITE (CPL(Jw),'(A38/(9(E13.6,2X)))')'Sediment P',(seDp(K,I),K=KTWB(JW),KB(I))     
              end do
              DO I=CUS(JB),DS(JB)
                WRITE (CPL(Jw),'(A38/(9(E13.6,2X)))')'Sediment N',(seDn(K,I),K=KTWB(JW),KB(I))     
              end do                                
              DO I=CUS(JB),DS(JB)
                WRITE (CPL(Jw),'(A38/(9(E13.6,2X)))')'Sediment C',(seDc(K,I),K=KTWB(JW),KB(I))     
              end do                                
            end if
            do m=1,nmc
              IF (print_macrophyte(jw,m)) THEN                                                  
                DO I=CUS(JB),DS(JB)
                  WRITE (CPL(Jw),'(A38/(9(E13.6,2X)))')'Macrophytes',((macrc(j,K,I,m),j=kti(i),kb(i)), K=KTwb(Jw),KB(I))     
                end do                                
              end if
            end do
! v3.5 end
          IF (CONSTITUENTS) THEN
            DO JD=1,NACD(JW)
              IF (PRINT_DERIVED(CDN(JD,JW),JW)) THEN
                  WRITE (CPL(JW),'(A38/(9(F10.3,2X)))') CDNAME(CDN(JD,JW)),((CD(K,I,CDN(JD,JW))*CDMULT(CDN(JD,JW)),             &        ! SW 8/12/06
                                                        K=KTWB(JW),KB(I)),I=CUS(JB),DS(JB))                                              ! CB 1/03/05
              END IF
            END DO
          END IF
          END DO
        END IF
      END IF

!**** Fluxes

      IF (FLUX(JW)) THEN
        IF (JDAY >= NXTMFL(JW) .OR. JDAY >= FLXD(FLXDP(JW)+1,JW)) THEN
          IF (JDAY >= FLXD(FLXDP(JW)+1,JW)) THEN
            FLXDP(JW)  = FLXDP(JW)+1
            NXTMFL(JW) = FLXD(FLXDP(JW),JW)
          END IF
          NXTMFL(JW) = NXTMFL(JW)+FLXF(FLXDP(JW),JW)
          CONV       = BLANK
          DO JAF=1,NAF(JW)
            DO JB=BS(JW),BE(JW)
              DO I=CUS(JB),DS(JB)
                DO K=KTWB(JW),KB(I)
                  KFS(K,I,KFCN(JAF,JW)) = KFS(K,I,KFCN(JAF,JW))/ELTMF(JW)*DAY
                END DO
              END DO
            END DO
            DO I=1,NISNP(JW)
              DO K=KTWB(JW),KB(ISNP(I,JW))
                WRITE (CONV(K,I),'(E10.3)') KFS(K,ISNP(I,JW),KFCN(JAF,JW))/1000.0
              END DO
            END DO
            IF (NEW_PAGE) THEN
              WRITE (FLX(JW),'(/(A72))') (TITLE(J),J=1,11)
              NLINES   =  KMX-KTWB(JW)+14
              NEW_PAGE = .FALSE.
            END IF
            NLINES   = NLINES+KMX-KTWB(JW)+11
            NEW_PAGE = NLINES > 72
            WRITE (FLX(JW),'(/A,F10.3,X,3(A,I0),A,F0.2,A/)') 'New date ',JDAY,MONTH//' ',GDAY,', ',YEAR,'   Julian Date = ',       &
                                                              INT(JDAY),' days ',(JDAY-INT(JDAY))*24.0,                            &
                                                            ' hours           '//KFNAME(KFCN(JAF,JW))
            WRITE (FLX(JW),'(3X,2000I10)')                  (ISNP(I,JW),I=1,NISNP(JW))
            DO K=KTWB(JW),KBR(JW)
              WRITE (FLX(JW),'(1X,I2,200A)') K,(CONV(K,I),I=1,NISNP(JW))
            END DO
          END DO
          ELTMF(JW)                    = 0.0
          KF(:,KFCN(1:NAF(JW),JW),JW)  = 0.0
          KFS(:,KFCN(1:NAF(JW),JW),JW) = 0.0
        END IF
      END IF

!**** Screen output

      IF (SCREEN_OUTPUT(JW)) THEN
        IF (JDAY >= NXTMSC(JW) .OR. JDAY >= SCRD(SCRDP(JW)+1,JW)) THEN
          IF (JDAY >= SCRD(SCRDP(JW)+1,JW)) THEN
            SCRDP(JW)  = SCRDP(JW)+1
            NXTMSC(JW) = SCRD(SCRDP(JW),JW)
          END IF
          KT         = KTWB(JW)
          NXTMSC(JW) = NXTMSC(JW)+SCRF(SCRDP(JW),JW)
          WRITE (*,'(A,F10.3,2(A,F6.0),A,I0)') 'JDAY = ',JDAY,' DLT = ',DLT,' DLTAV = ',DLTAV,' NV = ',NV
        END IF
      END IF
    END DO

!** Time series

    IF (TIME_SERIES) THEN
      IF (JDAY.GE.NXTMTS.OR.JDAY.GE.TSRD(TSRDP+1)) THEN
        IF (JDAY.GE.TSRD(TSRDP+1)) THEN
          TSRDP  = TSRDP+1
          NXTMTS = TSRD(TSRDP)
        END IF
        NXTMTS = NXTMTS+TSRF(TSRDP)
        DO J=1,NIKTSR
          I = ITSR(J)
          DO JW=1,NWB
            IF (I >= US(BS(JW))-1 .AND. I <= DS(BE(JW))+1) EXIT
          END DO
          IF (ETSR(J) < 0) THEN
            K = INT(ABS(ETSR(J)))
          ELSE
            DO K=KTWB(JW),KB(I)
              IF (DEPTHB(K,I) > ETSR(J)) EXIT
            END DO
            IF (K > KB(I)) CYCLE
          END IF
          DO JAC=1,NAC
            L = LEN_TRIM(FMTC(CN(JAC)))
            WRITE (C2CH(JAC),FMTC(CN(JAC))(1:L)) C2(K,I,CN(JAC))*CMULT(CN(JAC))
          END DO
          DO JAD=1,NACD(JW)
            L = LEN_TRIM(FMTCD(CDN(JAD,JW)))
            WRITE (CDCH(JAD),FMTCD(CDN(JAD,JW))(1:L)) CD(K,I,CDN(JAD,JW))*CDMULT(CDN(JAD,JW))
          END DO
          DO JE=1,NEP
            WRITE (EPCH(JE),'(F10.3)') EPD(K,I,JE)                                                    ! SW 8/13/06
          END DO
          DO Jm=1,Nmc
            WRITE (macCH(Jm),'(F10.3)') mac(K,I,Jm)                                                   ! SW 8/13/06
          END DO
! v3.5 start
          if(sediment_calc(jw))then
            write (sedch,'(F10.3)') sed(K,I)                                                          ! SW 8/13/06
            write (sedpch,'(F10.3)') sedp(K,I)
            write (sednch,'(F10.3)') sedn(K,I)
            write (sedcch,'(F10.3)') sedc(K,I)
          end if 
          IF (ICE_COMPUTATION) THEN
            if(sediment_calc(jw))then
              WRITE (TSR(J),'(f10.3,11F10.2,1000A)') JDAY,DLT,ELWS(I),T1(K,I),U(K,I),QC(I),SRON(JW)*1.06,ET(I),DEPTHB(KB(I),I),    &     ! SW 8/13/06
                                             BI(KTWB(JW),I),SHADE(I),ICETH(I),(ADJUSTR(C2CH(JAC)),JAC=1,NAC),                      &     ! CB 7/26/07
                                            (ADJUSTR(EPCH(JE)),JE=1,NEP),(ADJUSTR(macCH(Jm)),Jm=1,Nmc),sedch,sedpch,sednch,sedcch, &
                                            (ADJUSTR(CDCH(JAD)),JAD=1,NACD(JW))
            else
            WRITE (TSR(J),'(f10.3,11F10.2,1000A)') JDAY,DLT,ELWS(I),T1(K,I),U(K,I),QC(I),SRON(JW)*1.06,ET(I),DEPTHB(KB(I),I),      &     ! SW 8/13/06
                                             BI(KTWB(JW),I),SHADE(I),ICETH(I),(ADJUSTR(C2CH(JAC)),JAC=1,NAC),                      &    ! CB 7/26/07
                                            (ADJUSTR(EPCH(JE)),JE=1,NEP),(ADJUSTR(macCH(Jm)),Jm=1,Nmc),                            &
                                            (ADJUSTR(CDCH(JAD)),JAD=1,NACD(JW))
            end if
          ELSE
            if(sediment_calc(jw))then
              WRITE (TSR(J),'(f10.3,10F10.2,1000A)') JDAY,DLT,ELWS(I),T1(K,I),U(K,I),QC(I),SRON(JW)*1.06,ET(I),DEPTHB(KB(I),I),    &     ! SW 8/13/06
                                             BI(KTWB(JW),I),SHADE(I),(ADJUSTR(C2CH(JAC)),JAC=1,NAC),(ADJUSTR(EPCH(JE)),            &    ! CB 7/26/07
                                             JE=1,NEP),(ADJUSTR(macCH(Jm)),Jm=1,Nmc),sedch,sedpch,sednch,sedcch,                   &
                                             (ADJUSTR(CDCH(JAD)),JAD=1,NACD(JW))                                          
            else
            WRITE (TSR(J),'(f10.3,10F10.2,1000A)') JDAY,DLT,ELWS(I),T1(K,I),U(K,I),QC(I),SRON(JW)*1.06,ET(I),DEPTHB(KB(I),I),      &      ! SW 8/13/06
                                             BI(KTWB(JW),I),SHADE(I),(ADJUSTR(C2CH(JAC)),JAC=1,NAC),(ADJUSTR(EPCH(JE)),            &    ! CB 7/26/07
                                             JE=1,NEP),(ADJUSTR(macCH(Jm)),Jm=1,Nmc),(ADJUSTR(CDCH(JAD)),JAD=1,NACD(JW))
            end if
          END IF
! v3.5 end
        END DO
      END IF
    END IF

!** Downstream flow, temperature, and constituent files

    IF (DOWNSTREAM_OUTFLOW) THEN
      IF (JDAY >= NXTMWD .OR. JDAY >= WDOD(WDODP+1)) THEN
        IF (JDAY >= WDOD(WDODP+1)) THEN
          WDODP  = WDODP+1
          NXTMWD = WDOD(WDODP)
        END IF
        NXTMWD = NXTMWD+WDOF(WDODP)
        DO J=1,NIWDO
          QWDO(J)    = 0.0
          TWDO(J)    = 0.0
          CWDO(:,J)  = 0.0
          CDWDO(:,J) = 0.0
          CDTOT      = 0.0
          DO JWD=1,JWW
            IF (IWD(JWD) == IWDO(J) .AND. QWD(JWD) /= 0.0 .AND. ILAT(JWD) == 0) THEN
              TSUM  = 0.0
              QSUMM = 0.0
              CSUM  = 0.0
              CDSUM = 0.0
              DO JWWD=1,NWB
                IF (JBWD(JWD) >= BS(JWWD) .AND. JBWD(JWD) <= BE(JWWD)) EXIT
              END DO
              DO K=KTW(JWD),KBW(JWD)
                QSUMM = QSUMM+QSW(K,JWD)
                TSUM  = TSUM+T2(K,IWD(JWD))*QSW(K,JWD)
                DO JC=1,NAC
                  if (cn(jc) .ne. ndo) then                 !cb 09/28/04
!				  IF (CN(JC) == NDO) THEN                  
!                    IFLAG = 0
!                    DO JS=1,NSP
!                      IF (TDG_SPILLWAY(JWD,JS)) THEN
!                        CALL TOTAL_DISSOLVED_GAS (PALT(IWD(JWD)),0,JS,T2(K,IWD(JWD)),CGAS)
!                        CGASD        = (CGAS/EXP(7.7117-1.31403*(LOG(T2(K,IWD(JWD))+45.93)))*PALT(IWD(JWD)))*100.0
!                        CDSUM(NDG)   =  CDSUM(NDG)+CGASD*QSW(K,JWD)
!                        CSUM(CN(JC)) =  CSUM(CN(JC))+CGAS*QSW(K,JWD)
!                        IFLAG        =  1; EXIT
!                      END IF
!                    END DO
!                    IF (IFLAG == 0) THEN
!                      DO JG=1,NGT
!                        IF (TDG_GATE(JWD,JG)) THEN
!                          CALL TOTAL_DISSOLVED_GAS (PALT(IWD(JWD)),1,JG,T2(K,IWD(JWD)),CGAS)
!                          CGASD        = (CGAS/EXP(7.7117-1.31403*(LOG(T2(K,IWD(JWD))+45.93)))*PALT(IWD(JWD)))*100.0
!                          CDSUM(NDG)   =  CDSUM(NDG)+CGASD*QSW(K,JWD)
!                          CSUM(CN(JC)) =  CSUM(CN(JC))+CGAS*QSW(K,JWD)
!                          IFLAG        =  1; EXIT
!                        END IF
!                      END DO
!                    END IF
!                    IF (IFLAG == 0) CSUM(CN(JC)) = CSUM(CN(JC))+C2(K,IWD(JWD),CN(JC))*QSW(K,JWD)
!                  ELSE
                    CSUM(CN(JC)) = CSUM(CN(JC))+C2(K,IWD(JWD),CN(JC))*QSW(K,JWD)
                  END IF
                END DO
                DO JC=1,NACD(JWWD)
                  if (cdn(jc,jwwd) .ne. ndg) then                !cb 09/28/04
!                  IF (CDN(JC,JWWD) == NDG) THEN
!                    IFLAG = 0
!                    DO JG=1,NGT
!                      IF (TDG_GATE(JWD,JG)) THEN
!                        IFLAG = 1; EXIT
!                      END IF
!                    END DO
!                    DO JS=1,NSP
!                      IF (TDG_SPILLWAY(JWD,JS)) THEN
!                        IFLAG = 1; EXIT
!                      END IF
!                    END DO
!                    IF (IFLAG == 0) CDSUM(CDN(JC,JWWD)) = CDSUM(CDN(JC,JWWD))+CD(K,IWD(JWD),CDN(JC,JWWD))*QSW(K,JWD)
!                  ELSE
                    CDSUM(CDN(JC,JWWD)) = CDSUM(CDN(JC,JWWD))+CD(K,IWD(JWD),CDN(JC,JWWD))*QSW(K,JWD)
                  END IF
                END DO
              END DO

              iflag = 0                                                                                          !cb 09/28/04
              do js=1,nsp                                                                                        !cb 09/28/04
                if (tdg_spillway(jwd,js)) then                                                                   !cb 09/28/04
                  tspill=tsum/qsp(js)                                                                            !cb 09/28/04
                  call total_dissolved_gas (PALT(IWD(JWD)),0,js,tspill,cgas)                                     !cb 09/28/04
                  cgasd        = (CGAS/EXP(7.7117-1.31403*(LOG(tspill+45.93)))*PALT(IWD(JWD)))*100.0     !cb 09/28/04
                  cdsum(ndg)   =  cdsum(ndg)+cgasd*qsp(js)                                                       !cb 09/28/04
                  csum(ndo) =  csum(ndo)+cgas*qsp(js)                                                            !cb 09/28/04
                  iflag        =  1; exit                                                                        !cb 09/28/04
                end if                                                                                           !cb 09/28/04
              end do                                                                                             !cb 09/28/04
              if (iflag == 0) then                                                                               !cb 09/28/04
                do jg=1,ngt                                                                                      !cb 09/28/04
                  if (tdg_gate(jwd,jg)) then                                                                     !cb 09/28/04
                    tgate=tsum/qgt(jg)                                                                           !cb 09/28/04
                    call total_dissolved_gas (PALT(IWD(JWD)),1,jg,tgate,cgas)                                    !cb 09/28/04
                    cgasd        = (cgas/exp(7.7117-1.31403*(log(tgate+45.93)))*PALT(IWD(JWD)))*100.0                      !cb 09/28/04
                    cdsum(ndg)   =  cdsum(ndg)+cgasd*qgt(jt)                                                     !cb 09/28/04
                    csum(ndo) =  csum(ndo)+cgas*qgt(jg)                                                          !cb 09/28/04
                    iflag        =  1; exit                                                                      !cb 09/28/04
                  end if                                                                                         !cb 09/28/04
                end do                                                                                           !cb 09/28/04
              end if                                                                                             !cb 09/28/04
              if (iflag == 0 .and. cac(ndo)    == '      ON')then                                                !cb 09/28/04 
                do k=ktw(jwd),kbw(jwd)                                                                           !cb 09/28/04
                  csum(ndo) = csum(ndo)+c2(k,iwd(jwd),ndo)*qsw(k,jwd)                                            !cb 09/28/04
                end do                                                                                           !cb 09/28/04
              end if                                                                                             !cb 09/28/04
			  
              iflag = 0                                                                                          !cb 09/28/04
              do jg=1,ngt                                                                                        !cb 09/28/04
                if (tdg_gate(jwd,jg)) then                                                                       !cb 09/28/04
                  iflag = 1; exit                                                                                !cb 09/28/04
                end if                                                                                           !cb 09/28/04
              end do                                                                                             !cb 09/28/04
              do js=1,nsp                                                                                        !cb 09/28/04
                if (tdg_spillway(jwd,js)) then                                                                   !cb 09/28/04
                  iflag = 1; exit                                                                                !cb 09/28/04
                end if                                                                                           !cb 09/28/04
              end do                                                                                             !cb 09/28/04
              if (iflag == 0 .and. cdwbc(ndg,jwwd) == '      ON' ) then                                     !cb 09/28/04
                do k=ktw(jwd),kbw(jwd)                                                                    !cb 09/28/04
                  cdsum(ndg) = cdsum(ndg)+cd(k,iwd(jwd),ndg)*qsw(k,jwd)                                   !cb 09/28/04 
                end do                                                                                    !cb 09/28/04
              end if                                                                                      !cb 09/28/04

              QWDO(J)                         = QWDO(J)                    +QSUMM
              TWDO(J)                         = TSUM                       +TWDO(J)
              CWDO(CN(1:NAC),J)               = CSUM(CN(1:NAC))            +CWDO(CN(1:NAC),J)
              CDWDO(CDN(1:NACD(JWWD),JWWD),J) = CDSUM(CDN(1:NACD(JWWD),JWWD))+CDWDO(CDN(1:NACD(JWWD),JWWD),J)
            END IF
          END DO
          DO JW=1,NWB
            DO JB=BS(JW),BE(JW)
              IF (DS(JB) == IWDO(J)) THEN
                QWDO(J)           = QWDO(J)            +QSUM(JB)
                TWDO(J)           = (TOUT(JB)          *QSUM(JB))+TWDO(J)
                CWDO(CN(1:NAC),J) = (COUT(CN(1:NAC),JB)*QSUM(JB))+CWDO(CN(1:NAC),J)
                DO K=KTWB(JW),KB(DS(JB))
                  CDTOT(CDN(1:NACD(JW),JW)) = CDTOT(CDN(1:NACD(JW),JW))+CD(K,DS(JB),CDN(1:NACD(JW),JW))*QOUT(K,JB)
                END DO
                CDWDO(CDN(1:NACD(JW),JW),J) = CDWDO(CDN(1:NACD(JW),JW),J)+CDTOT(CDN(1:NACD(JW),JW))
              END IF
            END DO
          END DO
          IF (QWDO(J) /= 0.0) TWDO(J) = TWDO(J)/QWDO(J)
          DO JC=1,NAC
            IF (QWDO(J) /= 0.0) CWDO(CN(JC),J) = CWDO(CN(JC),J)/QWDO(J)
            WRITE (CWDOC(CN(JC)),'(G8.3)') CWDO(CN(JC),J)
            CWDOC(CN(JC)) = ADJUSTR(CWDOC(CN(JC)))
          END DO
          DO JW=1,NWB
            IF (IWDO(J) >= US(BS(JW)) .AND. IWDO(J) <= DS(BE(JW))) EXIT
          END DO
          DO JD=1,NACD(JW)
            IF (QWDO(J) /= 0.0) CDWDO(CDN(JD,JW),J) = CDWDO(CDN(JD,JW),J)/QWDO(J)
            WRITE (CDWDOC(CDN(JD,JW)),'(G8.3)') CDWDO(CDN(JD,JW),J)
            CDWDOC(CDN(JD,JW)) = ADJUSTR(CDWDOC(CDN(JD,JW)))
          END DO
          WRITE (WDO(J,1),'(F8.3,F8.2)') JDAY, QWDO(J)
          WRITE (WDO(J,2),'(F8.3,F8.2)') JDAY, TWDO(J)
          IF (CONSTITUENTS) WRITE (WDO(J,3),'(F8.3,1000A8)') JDAY,(CWDOC(CN(JC)),     JC=1,NAC)
          IF (DERIVED_CALC) WRITE (WDO(J,4),'(F8.3,1000A8)') JDAY,(CDWDOC(CDN(JD,JW)),JD=1,NACD(JW))
        END DO
      END IF
    END IF

!** Restart

    IF (RESTART_OUT) THEN
      IF (JDAY >= NXTMRS .OR. JDAY >= RSOD(RSODP+1)) THEN
        IF (JDAY >= RSOD(RSODP+1)) THEN
          RSODP  = RSODP+1
          NXTMRS = RSOD(RSODP)
        END IF
        NXTMRS = NXTMRS+RSOF(RSODP)
        WRITE (EXT,'(I0)') INT(JDAY)
        EXT   = ADJUSTL(EXT)
        L     = LEN_TRIM(EXT)
        RSOFN = 'rso'//EXT(1:L)//'.opt'
        CALL RESTART_OUTPUT (RSOFN)
      END IF
    END IF
  END DO
230 CONTINUE

!***********************************************************************************************************************************
!*                                                    Task 3: End Simulation                                                      **
!***********************************************************************************************************************************

  IF (.NOT. ERROR_OPEN) TEXT = 'Normal termination at '//CTIME(1:2)//':'//CTIME(3:4)//':'//CTIME(5:6)//' on '//CDATE(5:6)//'/'     &
                                                       //CDATE(7:8)//'/'//CDATE(3:4)
  CALL CPU_TIME (CURRENT)
  DO JW=1,NWB
    IF (SNAPSHOT(JW)) THEN
      WRITE (SNP(JW),'(/A/)')            ADJUSTL(TRIM(TEXT))
      WRITE (SNP(JW),'(A)')             'Runtime statistics'
      WRITE (SNP(JW),'(2(A,I0))')       '  Grid                 = ', IMX,' x ',KMX
      WRITE (SNP(JW),'(A,I0)')          '  Maximum active cells = ', NTACMX,'  Minimum active cells = ',NTACMN
      WRITE (SNP(JW),'(3(A,F0.1))')     '  Segment lengths, m   = ', DLXMIN,'-',DLXMAX
      WRITE (SNP(JW),'(3(A,F0.1))')     '  Layer heights, m     = ', HMIN,  '-',HMAX
      WRITE (SNP(JW),'(A)')             '  Timestep'
      WRITE (SNP(JW),'(A,I0)')          '    Total iterations   = ', NIT
      WRITE (SNP(JW),'(A,I0)')          '    # of violations    = ', NV
      WRITE (SNP(JW),'(A,F0.2)')        '    % violations       = ', FLOAT(NV)/FLOAT(NIT)*100.0
      WRITE (SNP(JW),'(A,I0,A)')        '    Average timestep   = ', INT(DLTAV),' sec'
      WRITE (SNP(JW),'(A,I0,A,F0.2,A)') '  Simulation time      = ', INT(ELTMJD),' days ',(ELTMJD-INT(ELTMJD))*24.0,' hours'
      WRITE (SNP(JW),'(A,F0.2,A)')      '  Total CPU runtime    = ',(CURRENT-START)/60.0,' min'
      CLOSE (SNP(JW))
    END IF
    IF (VECTOR(JW))      CLOSE (VPL(JW))
    IF (PROFILE(JW))     CLOSE (PRF(JW))
    IF (SPREADSHEET(JW)) CLOSE (SPR(JW))
    IF (CONTOUR(JW))     CLOSE (CPL(JW))
  END DO
  IF (TIME_SERIES) THEN
    DO J=1,NIKTSR
      CLOSE (TSR(J))
    END DO
  END IF
  IF (WARNING_OPEN) THEN
    CLOSE (WRN)
  ELSE
    CLOSE (WRN,STATUS='DELETE')
  END IF
  IF (ERROR_OPEN) THEN
    CLOSE (W2ERR)
  ELSE
    CLOSE (W2ERR,STATUS='DELETE')
  END IF
  DO J=40,NOPEN
   CLOSE (J)
  END DO
  DEALLOCATE (TSR,    WDO,    ETSR,   IWDO,   ITSR,   TITLE,  CDAC,   WSC,    ESTR,   WSTR,   QSTR,   KTSW,   KBSW,   SINKC)
  DEALLOCATE (EBC,    MBC,    PQC,    EVC,    PRC,    WINDC,  QINC,   QOUTC,  HEATC,  SLHTC,  QINIC,  DTRIC,  TRIC,   WDIC)
  DEALLOCATE (EXC,    EXIC,   VBC,    METIC,  SLTRC,  THETA,  FRICC,  NAF,    ELTMF,  ZMIN,   IZMIN,  C2CH,   CDCH,   EPCH)
  DEALLOCATE (CPLTC,  HPLTC,  CMIN,   CMAX,   HYMIN,  HYMAX,  CDMIN,  CDMAX,  JBDAM,  ILAT,   CDPLTC, QINSUM, TINSUM, TIND)
  DEALLOCATE (QOLD,   DTP,    DTPS,   QOLDS,  QIND,   JSS,    HDIC,   QNEW,   YSS,    VSS,    YS,     VS,     VSTS,   NSPRF)
  DEALLOCATE (LATGTC, LATSPC, LATPIC, LATPUC, DYNGTC, OPT,    CIND,   CINSUM, CDWBC,  KFWBC,  CPRWBC, CINBRC, CTRTRC, CDTBRC)
  DEALLOCATE (YSTS,   YST,    VST,    ALLIM,  APLIM,  ANLIM,  ASLIM,  ELLIM,  EPLIM,  ENLIM,  ESLIM,  CSSK,   C1,     C2)
  DEALLOCATE (KFS,    AF,     EF,     HYD,    KF,     AZSLC,  STRIC,  CPRBRC, CD,     KBMAX,  ELKT,   WIND2,  VISC,   CELC)
  DEALLOCATE (QOAVR,  QIMXR,  QOMXR,  REAERC, LAT,    LONGIT, ELBOT,  BTH,    VPR,    LPR,    NISNP,  NIPRF,  NISPR,  DECL)
  DEALLOCATE (A00,    HH,     T2I,    KTWB,   KBR,    IBPR,   DLVR,   ESR,    ETR,    NBL,    LPRFN,  EXTFN,  BTHFN,  METFN)
  DEALLOCATE (SNPFN,  PRFFN,  SPRFN,  CPLFN,  VPLFN,  FLXFN,  VPRFN,  AFW,    BFW,    CFW,    WINDH,  RHEVC,  FETCHC, JBDN)
! v3.5 start
  DEALLOCATE (SDK,    FSOD,   FSED,   SEDCI,  SEDCc,   SEDPRC, ICEC,   SLICEC, ICETHI, ALBEDO, HWI,    BETAI,  GAMMAI, ICEMIN)
  DEALLOCATE (seds,   sedb)    !cb 11/28/06
! v3.5 end
  DEALLOCATE (EXH2O,  BETA,   EXOM,   EXSS,   DXI,    CBHE,   TSED,   TSEDF,  FI,     ICET2,  AZC,    AZMAX,  QINT,   QOUTT)
  DEALLOCATE (AX,     WTYPEC, TAIR,   TDEW,   WIND,   PHI,    CLOUD,  CSHE,   SRON,   RAN,    RB,     RC,     RE,     SHADE)
  DEALLOCATE (ET,     RS,     RN,     SNPC,   SCRC,   PRFC,   SPRC,   CPLC,   VPLC,   FLXC,   NXTMCP, NXTMVP, NXTMFL, GAMMA)
  DEALLOCATE (NXTMSN, NXTMSC, NXTMPR, NXTMSP, SNPDP,  SCRDP,  PRFDP,  SPRDP,  CPLDP,  VPLDP,  FLXDP,  NCPL,   NVPL,   NFLX)
  DEALLOCATE (NSNP,   NSCR,   NPRF,   NSPR,   NEQN,   PO4R,   PARTP,  NH4DK,  NH4R,   NO3DK,  NO3S,   FER,    FES,    CDSUM)
  DEALLOCATE (CO2R,   SROC,   O2ER,   O2EG,   CAQ10,  CADK,   CAS,    BODP,   BODN,   BODC,   KBOD,   TBOD,   RBOD,   DTRC)
  DEALLOCATE (LDOMDK, RDOMDK, LRDDK,  OMT1,   OMT2,   OMK1,   OMK2,   LPOMDK, RPOMDK, LRPDK,  POMS,   ORGP,   ORGN,   ORGC)
  DEALLOCATE (RCOEF1, RCOEF2, RCOEF3, RCOEF4, ORGSI,  NH4T1,  NH4T2,  NH4K1,  NH4K2,  NO3T1,  NO3T2,  NO3K1,  NO3K2,  NSTR)
  DEALLOCATE (DSIR,   PSIS,   PSIDK,  PARTSI, SODT1,  SODT2,  SODK1,  SODK2,  O2NH4,  O2OM,   O2AR,   O2AG,   CG1DK,  CGS)
  DEALLOCATE (CGQ10,  CG0DK,  CUNIT,  CUNIT1, CUNIT2, CAC,    INCAC,  TRCAC,  DTCAC,  PRCAC,  CNAME,  CNAME1, CNAME2, CMULT)
  DEALLOCATE (CN,     INCN,   DTCN,   PRCN,   CSUM,   DLTMAX, QWDO,   TWDO,   SSS,    SEDRC,  TAUCR,  XBR)
  DEALLOCATE (QTAVB,  QTMXB,  BS,     BE,     JBUH,   JBDH,   TSSS,   TSSB,   TSSICE, ESBR,   ETBR,   EBRI,   QDTR,   EVBR)
  DEALLOCATE (QIN,    PR,     QPRBR,  TIN,    TOUT,   TPR,    TDTR,   TPB,    NACPR,  NACIN,  NACDT,  NACTR,  NACD,   ELDH)
  DEALLOCATE (QSUM,   NOUT,   KTQIN,  KBQIN,  ELUH,   NL,     NPOINT, SLOPE,  ALPHA,  COSA,   SINA,   TDHFN,  QOTFN,  PREFN)
  DEALLOCATE (CPRFN,  EUHFN,  TUHFN,  CUHFN,  EDHFN,  QINFN,  TINFN,  CINFN,  CDHFN,  QDTFN,  TDTFN,  CDTFN,  TPRFN,  VOLEV)
  DEALLOCATE (VOLWD,  VOLSBR, VOLTBR, DLVOL,  VOLG,   VOLSR,  VOLTR,  VOLB,   VOLPR,  VOLTRB, VOLDT,  VOLUH,  VOLDH,  VOLIN)
  DEALLOCATE (US,     DS,     CUS,    UHS,    DHS,    UQB,    DQB,    CDHS,   VOLOUT, TSSWD,  TSSUH,  TSSDH,  TSSIN,  TSSOUT)
  DEALLOCATE (TSSEV,  TSSPR,  TSSTR,  TSSDT,  SOD,    ELWS,   BKT,    REAER,  ICETH,  ICE,    ICESW,  Q,      QC,     QERR)
  DEALLOCATE (KTI,    SROSH,  SEG,    DLXRHO, QSSUM,  DLX,    DLXR,   QUH1,   QDH1,   BI,     JWUH,   JWDH)
  DEALLOCATE (A,      C,      D,      F,      V,      SKTI,   KBMIN,  EV,     QDT,    QPR,    SBKT,   BHRHO)
  DEALLOCATE (SZ,     WSHX,   WSHY,   WIND10, CZ,     FETCH,  PHI0,   FRIC,   ADZ,    HMULT,  FMTC,   FMTCD,  CNAME3, CDNAME3)
  DEALLOCATE (Z,      KB,     VNORM,  ANPR,   ANEQN,  APOM,   ACHLA,  AHSP,   AHSN,   AHSSI)
  DEALLOCATE (AC,     ASI,    AT1,    AT2,    AT3,    AT4,    AK1,    AK2,    AK3,    AK4,    EXA,    ASAT,   AP,     AN)
  DEALLOCATE (AG,     AR,     AE,     AM,     AS,     ENPR,   ENEQN,  EG,     ER,     EE,     EM,     EB,     ESAT,   EP)
  DEALLOCATE (EC,     ESI,    ECHLA,  EHSP,   EHSN,   EHSSI,  EPOM,   EHS,    EN,     ET4,    EK1,    EK2,    EK3,    EK4)
  DEALLOCATE (ET1,    ET2,    ET3,    HNAME,  FMTH,    KFAC,  KFNAME, KFCN,   C2I,    TRCN,   CDN,    CDNAME, CDNAME2,CDMULT)
  DEALLOCATE (CMBRS,  CMBRT,  FETCHU, FETCHD, IPRF,   ISNP,   ISPR,   BL,     LFPR,   DO3,    SED,    TKE,    PALT, AZT)
  DEALLOCATE (ADX,    DO1,    DO2,    B,      CONV,   CONV1,  EL,     DZ,     DZQ,    DX,     SAZ,    T1,TSS,QSS,BNEW, ILAYER)   ! SW 1/23/06
  DEALLOCATE (P,      SU,     SW,     BB,     BR,     BH,     BHR,    VOL,    HSEG,   DECAY,  FPFE,   FRICBR, UXBR,   UYBR)
  DEALLOCATE (DEPTHB, DEPTHM, FPSS,   TUH,    TDH,    TSSUH1, TSSUH2, TSSDH1, TSSDH2, SEDVP,  H,      EPC)
  DEALLOCATE (TVP,    QINF,   QOUT,   KOUT,   VOLUH2, VOLDH2, CWDO,   CDWDO,  CWDOC,  CDWDOC, CDTOT,  CPR,    CPB,    COUT)
  DEALLOCATE (CIN,    CDTR,   RSOD,   RSOF,   DLTD,   DLTF,   TSRD,   TSRF,   WDOD,   WDOF,   SNPD,   SNPF,   SPRD,   SPRF)
  DEALLOCATE (SCRD,   SCRF,   PRFD,   PRFF,   CPLD,   CPLF,   VPLD,   VPLF,   FLXD,   FLXF,   EPIC,   EPICI,  EPIPRC, EPIVP)
  DEALLOCATE (CUH,    CDH,    EPM,    EPD,    C1S,    CSSB,   CVP,    CSSUH1, CSSUH2, CSSDH2, CSSDH1, LNAME,  IWR,    KTWR)
  DEALLOCATE (JWUSP,  JWDSP,  QSP,    KBWR,   KTWD,   KBWD,   JBWD,   GTA1,   GTB1,   GTA2,   GTB2,   BGT,    IUGT,   IDGT)
  DEALLOCATE (QTR,    TTR,    KTTR,   KBTR,   EGT,    AGASGT, BGASGT, CGASGT, GASGTC, PUGTC,  ETUGT,  EBUGT,  KTUGT,  KBUGT)
  DEALLOCATE (PDGTC,  ETDGT,  EBDGT,  KTDGT,  KBDGT,  A1GT,   B1GT,   G1GT,   A2GT,   B2GT,   G2GT,   JWUGT,  JWDGT,  QGT)
  DEALLOCATE (EQGT,   JBUGT,  JBDGT,  JBUPI,  JBDPI,  JWUPI,  JWDPI,  QPI,    IUPI,   IDPI,   EUPI,   EDPI,   WPI,    DLXPI)
  DEALLOCATE (ETUPI,  EBUPI,  KTUPI,  KBUPI,  PDPIC,  ETDPI,  EBDPI,  KTDPI,  KBDPI,  FPI,    FMINPI, PUPIC,  ETDSP,  EBDSP)
  DEALLOCATE (PUSPC,  ETUSP,  EBUSP,  KTUSP,  KBUSP,  PDSPC,  KTDSP,  KBDSP,  IUSP,   IDSP,   ESP,    A1SP,   B1SP,   A2SP)
  DEALLOCATE (B2SP,   AGASSP, BGASSP, CGASSP, EQSP,   GASSPC, JBUSP,  JBDSP,  STRTPU, ENDPU,  EONPU,  EOFFPU, QPU,    PPUC)
  DEALLOCATE (IUPU,   IDPU,   EPU,    ETPU,   EBPU,   KTPU,   KBPU,   JWUPU,  JWDPU,  JBUPU,  JBDPU,  PUMPON, KTW,    KBW)
  DEALLOCATE (IWD,    KWD,    QWD,    EWD,    ITR,    QTRFN,  TTRFN,  CTRFN,  ELTRT,  ELTRB,  TRC,    JBTR,   QTRF,   CLRB)
  DEALLOCATE (TTLB,   TTRB,   CLLB,   SRLB1,  SRRB1,  SRLB2,  SRRB2,  SRFJD1, SHADEI, SRFJD2, TOPO,   QSW,    CTR)    ! SW 10/17/05
  DEALLOCATE (H1,     H2,     BH1,    BH2,    BHR1,   BHR2,   AVH1,   AVH2,   SAVH2,  AVHR,   SAVHR,  CBODD)
  DEALLOCATE (OPEN_VPR,       OPEN_LPR,             POINT_SINK,         HPRWBC,           READ_EXTINCTION, READ_RADIATION)
  DEALLOCATE (DIST_TRIBS,     UPWIND,               ULTIMATE,           FRESH_WATER,      SALT_WATER,      LIMITING_FACTOR)
  DEALLOCATE (UH_EXTERNAL,    DH_EXTERNAL,          UH_INTERNAL,        DH_INTERNAL,      UQ_INTERNAL,     DQ_INTERNAL)
  DEALLOCATE (UQ_EXTERNAL,    DQ_EXTERNAL,          UP_FLOW,            DN_FLOW,          UP_HEAD,         DN_HEAD)
  DEALLOCATE (INTERNAL_FLOW,  DAM_INFLOW,           DAM_OUTFLOW,        HEAD_FLOW,        HEAD_BOUNDARY)      !TC 08/03/04
  DEALLOCATE (ISO_CONC,             VERT_CONC,          LONG_CONC,        VERT_SEDIMENT,   LONG_SEDIMENT)
  DEALLOCATE (ISO_SEDIMENT,   VISCOSITY_LIMIT,      CELERITY_LIMIT,     IMPLICIT_AZ,      ONE_LAYER,       IMPLICIT_VISC)
  DEALLOCATE (FETCH_CALC,     LIMITING_DLT,         TERM_BY_TERM,       MANNINGS_N,       PLACE_QTR,       SPECIFY_QTR)
  DEALLOCATE (PLACE_QIN,      PRINT_CONST,          PRINT_HYDRO,        PRINT_SEDIMENT,   ENERGY_BALANCE,  MASS_BALANCE)
  DEALLOCATE (VOLUME_BALANCE, DETAILED_ICE,         ICE_CALC,           ICE_IN,           ALLOW_ICE,       PH_CALC)
  DEALLOCATE (EVAPORATION,    PRECIPITATION,        RH_EVAP,            NO_INFLOW,        NO_OUTFLOW,      NO_HEAT)
  DEALLOCATE (ISO_TEMP,       VERT_TEMP,            LONG_TEMP,          VERT_PROFILE,     LONG_PROFILE,    NO_WIND)
  DEALLOCATE (SNAPSHOT,       PROFILE,              VECTOR,             CONTOUR,          SPREADSHEET,     INTERNAL_WEIR)
  DEALLOCATE (SCREEN_OUTPUT,  FLUX,                 DYNAMIC_SHADE,      TRAPEZOIDAL)
  DEALLOCATE (SEDIMENT_CALC,  EPIPHYTON_CALC,       PRINT_DERIVED,      PRINT_EPIPHYTON,  TDG_SPILLWAY,    TDG_GATE)
  DEALLOCATE (ISO_EPIPHYTON,  VERT_EPIPHYTON,       LONG_EPIPHYTON,     LATERAL_SPILLWAY, LATERAL_GATE,    LATERAL_PUMP)
  DEALLOCATE (INTERP_HEAD,    INTERP_WITHDRAWAL,    INTERP_EXTINCTION,  INTERP_DTRIBS,    LATERAL_PIPE,    INTERP_TRIBS)
  DEALLOCATE (INTERP_OUTFLOW, INTERP_INFLOW,        INTERP_METEOROLOGY, CONSTITUENT_PLOT, DERIVED_PLOT,    ZERO_SLOPE)
  DEALLOCATE (HYDRO_PLOT,     SEDIMENT_RESUSPENSION)
! v3.5 start
  deallocate (ORGPLD, ORGPRD, ORGPLP, ORGPRP, ORGNLD, ORGNRD, ORGNLP)
  deallocate (ORGNRP)
  deallocate  (print_macrophyte, macrophyte_calc,macwbc,conv2)
  deallocate  (mac, macrc,mact, macrm, macss)
  deallocate  (mgr,mmr, mrr)
  deallocate  (smacrc, smacrm)
  deallocate  (smact, smac)
  deallocate  (mt1,mt2,mt3,mt4,mk1,mk2,mk3,mk4,mg,mr,mm)
  deallocate  (mp, mn, mc,psed,nsed,mhsp,mhsn,mhsc,msat)
  deallocate  (cddrag,kticol,armac,macwbci, anorm, dwv, dwsa)
  deallocate  (mbmp,mmax,mpom,lrpmac,o2mr,o2mg)
  deallocate  (macmbrs, macmbrt,ssmacmb)
  deallocate  (cw, bic)  
  deallocate  (mactrmr, mactrmf,mactrm)     
  deallocate  (mlfpr)
  deallocate  (mllim, mplim,mclim,mnlim)
  deALLOCATE  (GAMMAj)
  deallocate (por,VOLKTi,VOLi,vstem,vstemkt,sarea)  
  DEALLOCATE (IWIND) ! MLM 08/12/05
  deallocate (zg,zm,zeff,prefp,zr,zoomin,zs2p,exz,zt1,zt2,zt3,zt4,zk1,zk2)
  deallocate (ldompmp,ldomnmp,lpompmp,lpomnmp,rpompmp,rpomnmp,o2zr) ! mlm 06/10/06
  deallocate (mprwbc)                                               ! mlm 06/10/06
  deallocate (exm)                                                  ! mlm 06/10/06
! deallocate (sedbr,sedbrp,sedbrn,sedbrc)   ! SW 6/4/07 No need to deallocate pointers
  deallocate (zk3,zk4,zp,zn,zc,prefa,zmu,tgraze,zrt,zmt,zoorm,zoormr,zoormf) ! POINTERS ,zoo,zooss,
  deallocate (lpzooout,lpzooin,po4zr,nh4zr,dozr,ticzr,agz,agzt)
  deallocate (zgz,PREFZ) !omnivorous zooplankton
!  DEALLOCATE (LDOMPSS,  RDOMPSS,  LPOMPSS,  RPOMPSS,  LDOMNSS,  RDOMNSS)
!  DEALLOCATE (LPOMNSS,  RPOMNSS)  
!  DEALLOCATE (LDOMPAP, LDOMPEP, LPOMPAP, LPOMPNS, RPOMPNS)
!  DEALLOCATE (LDOMnAP, LDOMnEP, LPOMnAP, LPOMnNS, RPOMnNS)
!  DEALLOCATE (ldompmp,ldomnmp,lpompmp,lpomnmp,rpompmp,rpomnmp)
  DEALLOCATE (lpzooinp,lpzooinn,lpzoooutp,lpzoooutn)
!  DEALLOCATE (SEDDp,  SEDASp,  SEDOMSp, SEDNSp, lpomepp)
!  DEALLOCATE (SEDDn,  SEDASn,  SEDOMSn, SEDNSn, lpomepn, sedno3)
!  DEALLOCATE (SEDDc,  SEDASc,  SEDOMSc, SEDNSc, lpomepc)
  DEALLOCATE (sedc, sedn, sedp)
  DEALLOCATE (sedvpc, sedvpp, sedvpn)
  DEALLOCATE (sdkv,seddktot)  
  DEALLOCATE (cbods)
  DEALLOCATE(DZT)
! v3.5 end

  CALL DEALLOCATE_TIME_VARYING_DATA
  CALL DEALLOCATE_TRANSPORT
  CALL DEALLOCATE_KINETICS
  CALL DEALLOCATE_WATERBODY
  CALL DEALLOCATE_PIPE_FLOW
  CALL DEALLOCATE_OPEN_CHANNEL
240 CONTINUE
  WRITE (*,'(A)') TEXT

! Snapshot formats

10490 FORMAT ('CE-QUAL-W2 V3.1'/                                                                                                   &
             (1X,A72))
10500 FORMAT (/1X,A/                                                                                                               &
              3X,'Gregorian date      [GDAY] =',A19,1X,I0,', ',I0/                                                                 &
              3X,'Julian date         [JDAY] =',I10,' days',F6.2,' hours'/                                                         &
              3X,'Elapsed time      [ELTMJD] =',I10,' days',F6.2,' hours'/                                                         &
              3X,'Timestep             [DLT] =',I10,' sec'/                                                                        &
              3X,'  at location  [KLOC,ILOC] = (',I0,',',I0,')'/                                                                   &
              3X,'Minimum timestep  [MINDLT] =',I10,' sec '/                                                                       &
              3X,'  at Julian day    [JDMIN] =',I10,' days',F6.2,' hours'/                                                         &
              3X,'  at location  [KMIN,IMIN] = (',I0,',',I0,')')
10510 FORMAT (3X,'Limiting timestep'/                                                                                              &
              3X,'  at location  [KMIN,IMIN] = (',I0,',',I0,')')
10520 FORMAT (3X,'Average timestep   [DLTAV] =',I10,' sec'/                                                                        &
              3X,'Number of iterations [NIT] =',I10/                                                                               &
              3X,'Number of violations  [NV] =',I10/)
10530 FORMAT (1X,A)
10540 FORMAT (3X,'Input'/                                                                                                          &
              3X,'  Air temperature          [TAIR] =',F9.2,1X,A/                                                                  &
              3X,'  Dewpoint temperature     [TDEW] =',F9.2,1X,A/                                                                  &
              3X,'  Wind direction            [PHI] =',F9.2,' rad'/                                                                &
              3X,'  Cloud cover             [CLOUD] =',F9.2/                                                                       &
              3X,'  Calculated'/                                                                                                   &
              5X,'  Equilibrium temperature    [ET] =',F9.2,1X,A/                                                                  &
              5X,'  Surface heat exchange    [CSHE] =',E9.2,' m/sec'/                                                              &
              5X,'  Net short wave radiation [SRON] =',E9.2,1X,A,' W/m^2'/)
10550 FORMAT (1X,A/                                                                                                                &
              3X,A)
10560 FORMAT (5X,'Branch ',I0/                                                                                                     &
              5X,'  Layer       [KQIN] = ',I0,'-',I0/                                                                              &
              5X,'  Inflow       [QIN] =',F8.2,' m^3/sec'/                                                                         &
              5X,'  Temperature  [TIN] =',F8.2,1X,A)
10570 FORMAT (/3X,'Distributed Tributaries')
10580 FORMAT (5X,'Branch ',I0/                                                                                                     &
              5X,'  Inflow      [QDTR] =',F8.2,' m^3/sec'/                                                                         &
              5X,'  Temperature [TDTR] =',F8.2,1X,A)
10590 FORMAT (:/3X,'Tributaries'/                                                                                                  &
              5X,'Segment     [ITR] =',11I8:/                                                                                      &
             (T25,11I8))
10600 FORMAT (:5X,'Layer      [KTWB] = ',11(I0,'-',I0,2X):/                                                                        &
             (T25,11(I0,'-',I0)))
10610 FORMAT (:5X,'Inflow      [QTR] =',11F8.2:/                                                                                   &
             (T25,11F8.1))
10620 FORMAT (:5X,'Temperature [TTR] =',11F8.2:/                                                                                   &
             (T25,11F8.1))
10630 FORMAT (/1X,'Outflows')
10640 FORMAT (3X,'Structure outflows [QSTR]'/                                                                                      &
              3X,'  Branch ',I0,' = ',11F8.2:/                                                                                     &
             (T16,11F8.2))
10650 FORMAT (:/3X,'Total outflow [QOUT] =',F8.2,' m^3/s'/                                                                         &
              5X,'Outlets'/                                                                                                        &
              5X,'  Layer             [KOUT] =',12I7:/                                                                             &
             (33X,12I7))
10660 FORMAT (:7X,'Outflow (m^3/sec) [QOUT] =',12F7.2:/                                                                            &
             (33X,12F7.2))
10670 FORMAT (:5X,'Withdrawals'/                                                                                                   &
              5X,'  Segment            [IWD] =',I7/                                                                                &
              5X,'  Outflow (m^3/sec)  [QWD] =',F7.2)
10680 FORMAT (5X,'  Layer              [KWD] =',12I7/                                                                              &
             (33X,12I7))
10690 FORMAT (:5X,'  Outflow (m^3/sec)  [QSW] =',12F7.2/                                                                           &
             (33X,12F7.2))
10700 FORMAT (/'1',A)
10710 FORMAT (3X,'Branch ',I0,' [CIN]'/                                                                                            &
             (5X,A,T25,'=',F9.3,1X,A))
10720 FORMAT (3X,'Tributary ',I0,' [CTR]'/                                                                                         &
             (5X,A,T25,'=',F9.3,1X,A))
10730 FORMAT (3X,'Distributed tributary ',I0,' [CDT]'/                                                                             &
             (5X,A,T25,'=',F9.3,1X,A))
10740 FORMAT (/'Surface calculations')
10750 FORMAT (3X,'Evaporation rate [EV]'/                                                                                          &                                                                                       
             (:3X,'  Branch ',I0,' = ',E10.3,' m/s'))                                                                                      ! SW 9/15/05
10755 format (3x,'Cumulative evaporation [VOLEV]'/                                                                                 &                                                                                       
             (:3X,'  Branch ',I0,' = ',F0.1,' m^3'))    
10760 FORMAT (3X,'Precipitation [PR]'/                                                                                             &
             (3X,'  Branch ',I0,' = ',F8.6),' m/s')
10770 FORMAT (/1X,'External head boundary elevations'/)
10780 FORMAT (3X,'Branch ',I0/5X,'Upstream elevation   [ELUH] =',F8.3,' m')
10790 FORMAT (3X,'Branch ',I0/5X,'Downstream elevation [ELDH] =',F8.3,' m')
10800 FORMAT (/'Water Balance')
10810 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              3X,'  Spatial change  [VOLSR]  = ',E15.8,' m^3'/                                                                     &
              3X,'  Temporal change [VOLTR]  = ',E15.8,' m^3'/                                                                     &
              3X,'  Volume error             = ',E15.8,' m^3'/                                                                     &
              3X,'  Percent error            = ',E15.8,' %')
10820 FORMAT (3X,'Branch ',I0/                                                                                                     &
              3X,'  Spatial change  [VOLSBR] = ',E15.8,' m^3'/                                                                     &
              3X,'  Temporal change [VOLTBR] = ',E15.8,' m^3'/                                                                     &
              3X,'  Volume error             = ',E15.8,' m^3'/                                                                     &
              3X,'  Percent error            = ',E15.8,' %')
10830 FORMAT (/1X,'Energy Balance')
10840 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              3X,'  Spatially integrated energy   [ESR] = ',E15.8,' kJ'/                                                           &
              3X,'  Temporally integrated energy  [ETR] = ',E15.8,' kJ'/                                                           &
              3X,'  Energy error                        = ',E15.8,' kJ'/                                                           &
              3X,'  Percent error                       = ',E15.8,' %')
10850 FORMAT (3X,'  Spatially integrated energy  [ESBR] = ',E15.8,' kJ'/                                                           &
              3X,'  Temporally integrated energy [ETBR] = ',E15.8,' kJ'/                                                           &
              3X,'  Energy error                        = ',E15.8,' kJ'/                                                           &
              3X,'  Percent error                       = ',E15.8,' %')
10860 FORMAT (/1X,'Mass Balance')
10870 FORMAT (3X,'Branch ',I0)
10880 FORMAT (5X,A/                                                                                                                &
              5X,'  Spatially integrated mass  [CMBRS] = ',E15.8,1X,A/                                                             &
              5X,'  Temporally integrated mass [CMBRT] = ',E15.8,1X,A/                                                             &
              5X,'  Mass error                         = ',E15.8,1X,A/                                                             &
              5X,'  Percent error                      = ',E15.8,' %')
10890 FORMAT (/1X,A/                                                                                                               &
              3X,'Surface layer [KT] = ',I0/                                                                                       &
              3X,'Elevation   [ELKT] =',F10.3,' m')
10900 FORMAT (/3X,'Current upstream segment [CUS]'/                                                                                &
             (3X,'  Branch ',I0,' = ',I0))     
END PROGRAM CE_QUAL_W2

!***********************************************************************************************************************************
!**                                            F U N C T I O N   D E N S I T Y                                                    **
!***********************************************************************************************************************************

FUNCTION DENSITY (T,TDS,SS)
  USE LOGICC, ONLY: SUSP_SOLIDS, FRESH_WATER, SALT_WATER; USE GLOBAL
                       DENSITY = ((((6.536332E-9*T-1.120083E-6)*T+1.001685E-4)*T-9.09529E-3)*T+6.793952E-2)*T+0.842594
  IF (SUSP_SOLIDS)     DENSITY = DENSITY+6.2E-4*SS
  IF (FRESH_WATER(JW)) DENSITY = DENSITY+TDS*((4.99E-8*T-3.87E-6)*T+8.221E-4)
  IF (SALT_WATER(JW))  DENSITY = DENSITY+TDS*((((5.3875E-9*T-8.2467E-7)*T+7.6438E-5)*T-4.0899E-3)*T+0.824493)                      &
                                 +((-1.6546E-6*T+1.0227E-4)*T-5.72466E-3)*TDS**1.5+4.8314E-4*TDS*TDS
  DENSITY = DENSITY+999.0
END FUNCTION DENSITY

!***********************************************************************************************************************************
!**                                  S U B R O U T I N E   T I M E   V A R Y I N G   D A T A                                      **
!***********************************************************************************************************************************

SUBROUTINE TIME_VARYING_DATA
  USE GLOBAL;  USE SURFHE; USE SCREENC; USE TVDC; USE LOGICC; USE SELWC; USE STRUCTURES; USE NAMESC
  USE KINETIC, ONLY:EXH2O; USE SHADEC 

! Type declaration

  REAL                                   :: NXQWD1, NXQWD2, NXQGT,  NXTVD
  REAL                                   :: NXWSC
  REAL,    ALLOCATABLE, DIMENSION(:)     :: QDTRO,  TDTRO,  ELUHO,  ELDHO,  QWDO,   QTRO,   TTRO,   QINO,   TINO
  REAL,    ALLOCATABLE, DIMENSION(:)     :: TAIRNX, TDEWNX, PHINX,  WINDNX, SRONX,  CLOUDNX,BGTNX
  REAL,    ALLOCATABLE, DIMENSION(:)     :: NXEXT1, NXEXT2, EXTNX,  EXTO   
  REAL,    ALLOCATABLE, DIMENSION(:)     :: TAIRO,  TDEWO,  PHIO,   WINDO,  SROO,   CLOUDO
  REAL,    ALLOCATABLE, DIMENSION(:)     :: QDTRNX, TDTRNX, PRNX,   TPRNX,  ELUHNX, ELDHNX, QWDNX,  QTRNX,  TTRNX,  QINNX,  TINNX
  REAL,    ALLOCATABLE, DIMENSION(:)     :: NXQTR1, NXTTR1, NXCTR1, NXQIN1, NXTIN1, NXCIN1, NXQDT1, NXTDT1, NXCDT1
  REAL,    ALLOCATABLE, DIMENSION(:)     :: NXPR1,  NXTPR1, NXCPR1, NXEUH1, NXTUH1, NXCUH1, NXEDH1, NXTDH1, NXCDH1, NXQOT1, NXMET1
  REAL,    ALLOCATABLE, DIMENSION(:)     :: NXQTR2, NXTTR2, NXCTR2, NXQIN2, NXTIN2, NXCIN2, NXQDT2, NXTDT2, NXCDT2
  REAL,    ALLOCATABLE, DIMENSION(:)     :: NXPR2,  NXTPR2, NXCPR2, NXEUH2, NXTUH2, NXCUH2, NXEDH2, NXTDH2, NXCDH2, NXQOT2, NXMET2
  REAL,    ALLOCATABLE, DIMENSION(:)     :: WSCNX
  REAL,    ALLOCATABLE, DIMENSION(:,:)   :: CTRO,   CINO,   QOUTO,  CDTRO,  TUHO,   TDHO,   QSTRO
  REAL,    ALLOCATABLE, DIMENSION(:,:)   :: CTRNX,  CINNX,  QOUTNX, CDTRNX, CPRNX,  TUHNX,  TDHNX,  QSTRNX
  REAL,    ALLOCATABLE, DIMENSION(:,:,:) :: CUHO,   CDHO,   CUHNX,  CDHNX
  INTEGER                                :: WDQ,    GTQ,    WSH,    SHD    
  INTEGER, ALLOCATABLE, DIMENSION(:)     :: TRQ,    TRT,    TRC,    INQ,    DTQ,    PRE,    UHE,    DHE,    INFT,   DTT
  INTEGER, ALLOCATABLE, DIMENSION(:)     :: PRT,    UHT,    DHT,    INC,    DTC,    PRC,    UHC,    DHC,    OTQ,    MET,    EXT
  LOGICAL, ALLOCATABLE, DIMENSION(:)     :: INFLOW_CONST, TRIB_CONST, DTRIB_CONST, PRECIP_CONST
  SAVE

! Allocation declarations

  ALLOCATE (NXQTR1(NTR), NXTTR1(NTR), NXCTR1(NTR), NXQIN1(NBR), NXTIN1(NBR), NXCIN1(NBR), NXQDT1(NBR), NXTDT1(NBR), NXCDT1(NBR))
  ALLOCATE (NXPR1(NBR),  NXTPR1(NBR), NXCPR1(NBR), NXEUH1(NBR), NXTUH1(NBR), NXCUH1(NBR), NXEDH1(NBR), NXTDH1(NBR), NXCDH1(NBR))
  ALLOCATE (NXQOT1(NBR), NXMET1(NWB), NXQTR2(NTR), NXTTR2(NTR), NXCTR2(NTR), NXQIN2(NBR), NXTIN2(NBR), NXCIN2(NBR), NXQDT2(NBR))
  ALLOCATE (NXTDT2(NBR), NXCDT2(NBR), NXPR2(NBR),  NXTPR2(NBR), NXCPR2(NBR), NXEUH2(NBR), NXTUH2(NBR), NXCUH2(NBR), NXEDH2(NBR))
  ALLOCATE (NXTDH2(NBR), NXCDH2(NBR), NXQOT2(NBR), NXMET2(NWB))
  ALLOCATE (WSCNX(IMX))
  ALLOCATE (QDTRO(NBR),  TDTRO(NBR),  ELUHO(NBR),  ELDHO(NBR),  QWDO(NWD),   QTRO(NTR),   TTRO(NTR),   QINO(NBR),   TINO(NBR))
  ALLOCATE (QDTRNX(NBR), TDTRNX(NBR), PRNX(NBR),   TPRNX(NBR),  ELUHNX(NBR), ELDHNX(NBR), QWDNX(NWD),  QTRNX(NTR),  TTRNX(NTR))
  ALLOCATE (QINNX(NBR),  TINNX(NBR),  SROO(NWB),   TAIRO(NWB),  TDEWO(NWB),  CLOUDO(NWB), PHIO(NWB),   WINDO(NWB),  TAIRNX(NWB))
  ALLOCATE (TDEWNX(NWB), CLOUDNX(NWB),PHINX(NWB),  WINDNX(NWB), SRONX(NWB),  BGTNX(NGT))
  ALLOCATE (TRQ(NTR),    TRT(NTR),    TRC(NTR),    INQ(NBR),    DTQ(NBR),    PRE(NBR),    UHE(NBR),    DHE(NBR),    INFT(NBR))
  ALLOCATE (DTT(NBR),    PRT(NBR),    UHT(NBR),    DHT(NBR),    INC(NBR),    DTC(NBR),    PRC(NBR),    UHC(NBR),    DHC(NBR))
  ALLOCATE (OTQ(NBR),    MET(NWB),    EXT(NWB))
  ALLOCATE (NXEXT1(NWB), NXEXT2(NWB), EXTNX(NWB),  EXTO(NWB))     
  ALLOCATE (CTRO(NCT,NTR),   CINO(NCT,NBR),  QOUTO(KMX,NBR),  CDTRO(NCT,NBR),  TUHO(KMX,NBR),  TDHO(KMX,NBR),  QSTRO(NST,NBR))
  ALLOCATE (CTRNX(NCT,NTR),  CINNX(NCT,NBR), QOUTNX(KMX,NBR), CDTRNX(NCT,NBR), CPRNX(NCT,NBR), TUHNX(KMX,NBR), TDHNX(KMX,NBR))
  ALLOCATE (QSTRNX(NST,NBR))
  ALLOCATE (CUHO(KMX,NCT,NBR), CDHO(KMX,NCT,NBR), CUHNX(KMX,NCT,NBR), CDHNX(KMX,NCT,NBR))
  ALLOCATE (INFLOW_CONST(NBR), TRIB_CONST(NTR),   DTRIB_CONST(NBR),   PRECIP_CONST(NBR))

  NXPR1  = 0.0; NXQTR1 = 0.0; NXTTR1 = 0.0; NXCTR1 = 0.0; NXQIN1 = 0.0; NXTIN1 = 0.0; NXCIN1 = 0.0; NXQDT1 = 0.0; NXTDT1 = 0.0
  NXCDT1 = 0.0; NXTPR1 = 0.0; NXCPR1 = 0.0; NXEUH1 = 0.0; NXTUH1 = 0.0; NXCUH1 = 0.0; NXEDH1 = 0.0; NXTDH1 = 0.0; NXCDH1 = 0.0
  NXQOT1 = 0.0; NXMET1 = 0.0; QSTRNX = 0.0; CDTRNX = 0.0; CTRNX  = 0.0; CINNX  = 0.0; CPRNX  = 0.0; CUHNX  = 0.0; CDHNX  = 0.0
  QINNX  = 0.0; TINNX  = 0.0; CINNX  = 0.0; NXWSC  = 0.0

! Set logical variables

  INFLOW_CONST = CONSTITUENTS .AND. NACIN > 0; TRIB_CONST   = CONSTITUENTS .AND. NACTR > 0
  DTRIB_CONST  = CONSTITUENTS .AND. NACDT > 0; PRECIP_CONST = CONSTITUENTS .AND. NACPR > 0

! Open input files

  NPT = NUNIT
  SHD = NPT; NPT=NPT+1
  OPEN (SHD,FILE=SHDFN,STATUS='OLD')
  READ (SHD,'(///(8X,29F8.0))')       (SHADEI(I),TTLB(I),TTRB(I),CLLB(I),CLRB(I),SRLB1(I),SRLB2(I),SRRB1(I),SRRB2(I),              &
                                      (TOPO(I,J),J=1,IANG),SRFJD1(I),SRFJD2(I),I=1,IMX)
  SHADE = SHADEI 
  WSH = NPT; NPT = NPT+1
  OPEN (WSH,FILE=WSCFN,STATUS='OLD')
  READ (WSH,'(///10F8.0:/(8X,9F8.0))') NXWSC,(WSCNX(I),I=1,IMX)
  WSC = WSCNX
  READ (WSH,'(10F8.0:/(8X,9F8.0))')    NXWSC,(WSCNX(I),I=1,IMX)
  DO JW=1,NWB
    MET(JW) = NPT; NPT = NPT+1
    OPEN (MET(JW),FILE=METFN(JW),STATUS='OLD')
    IF (READ_RADIATION(JW)) THEN
      READ (MET(JW),'(///10F8.0)') NXMET2(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW),SRONX(JW)
      SRONX(JW) = SRONX(JW)*REFL
      SRON(JW)  = SRONX(JW)
      SROO(JW)  = SRON(JW)
    ELSE
      READ (MET(JW),'(///10F8.0)') NXMET2(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW)
    END IF
    TAIR(JW)   = TAIRNX(JW)
    TDEW(JW)   = TDEWNX(JW)
    WIND(JW)   = WINDNX(JW)
    PHI(JW)    = PHINX(JW)
    CLOUD(JW)  = CLOUDNX(JW)
    TAIRO(JW)  = TAIRNX(JW)
    TDEWO(JW)  = TDEWNX(JW)
    WINDO(JW)  = WINDNX(JW)
    PHIO(JW)   = PHINX(JW)
    CLOUDO(JW) = CLOUDNX(JW)
    IF (PHISET > 0) PHI(JW)  = PHISET
    IF (PHISET > 0) PHIO(JW) = PHISET
    IF (READ_RADIATION(JW)) THEN
      READ (MET(JW),'(10F8.0)') NXMET1(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW),SRONX(JW)
      SRONX(JW) = SRONX(JW)*REFL
    ELSE
      READ (MET(JW),'(10F8.0)') NXMET1(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW)
    END IF
    IF (READ_EXTINCTION(JW)) THEN
      EXT(JW) = NPT; NPT = NPT+1
      OPEN (EXT(JW),FILE=EXTFN(JW),STATUS='OLD')     
      READ (EXT(JW),'(///2F8.0)') NXEXT2(JW), EXTNX(JW)
      EXH2O(JW) = EXTNX(JW)
      EXTO(JW)  = EXTNX(JW)
      READ (EXT(JW),'(2F8.0)')    NXEXT1(JW), EXTNX(JW)
    END IF
    DO I=CUS(BS(JW)),DS(BE(JW))
      WIND2(I) = WIND(JW)*WSC(I)*LOG(2.0/0.003)/LOG(WINDH(JW)/0.003)
    END DO
  END DO
  IF (NWD > 0) THEN
    WDQ = NPT; NPT = NPT+1
    OPEN (WDQ,FILE=QWDFN,STATUS='OLD')
    READ (WDQ,'(///10F8.0:/(8X,9F8.0))') NXQWD2,(QWDNX(JW),JW=1,NWD)
    DO JW=1,NWD
      QWD(JW)  = QWDNX(JW)
      QWDO(JW) = QWDNX(JW)
    END DO
    READ (WDQ,'(10F8.0:/(8X,9F8.0))')    NXQWD1,(QWDNX(JW),JW=1,NWD)
  END IF
  IF (TRIBUTARIES) THEN
    DO JT=1,NTR
      TRQ(JT) = NPT; NPT = NPT+1
      TRT(JT) = NPT; NPT = NPT+1
      OPEN (TRQ(JT),FILE=QTRFN(JT),STATUS='OLD')
      OPEN (TRT(JT),FILE=TTRFN(JT),STATUS='OLD')
      READ (TRQ(JT),'(///2F8.0)') NXQTR2(JT),QTRNX(JT)
      READ (TRT(JT),'(///2F8.0)') NXTTR2(JT),TTRNX(JT)
      IF (TRIB_CONST(JT)) THEN
        TRC(JT) = NPT; NPT = NPT+1
        OPEN (TRC(JT),FILE=CTRFN(JT),STATUS='OLD')
        READ (TRC(JT),'(///1000F8.0)') NXCTR2(JT),(CTRNX(TRCN(JAC,JT),JT),JAC=1,NACTR(JT))
      END IF
    END DO
    QTR(1:NTR)    = QTRNX(1:NTR)
    QTRO(1:NTR)   = QTRNX(1:NTR)
    TTR(1:NTR)    = TTRNX(1:NTR)
    TTRO(1:NTR)   = TTRNX(1:NTR)
    CTR(:,1:NTR)  = CTRNX(:,1:NTR)
    CTRO(:,1:NTR) = CTRNX(:,1:NTR)
    DO JT=1,NTR
      READ (TRQ(JT),'(2F8.0)') NXQTR1(JT),QTRNX(JT)
      READ (TRT(JT),'(2F8.0)') NXTTR1(JT),TTRNX(JT)
      IF (TRIB_CONST(JT)) THEN
        READ (TRC(JT),'(1000F8.0)') NXCTR1(JT),(CTRNX(TRCN(JAC,JT),JT),JAC=1,NACTR(JT))
      END IF
    END DO
  END IF
  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      IF (UP_FLOW(JB)) THEN
        IF (.NOT. INTERNAL_FLOW(JB) .AND. .NOT. DAM_INFLOW(JB)) THEN                                                  !TC 08/03/04 RA 1/13/06
          INQ(JB)  = NPT; NPT = NPT+1
          INFT(JB) = NPT; NPT = NPT+1
          OPEN (INQ(JB) ,FILE=QINFN(JB),STATUS='OLD')
          OPEN (INFT(JB),FILE=TINFN(JB),STATUS='OLD')
          READ (INQ(JB), '(///2F8.0)') NXQIN2(JB),QINNX(JB)
          READ (INFT(JB),'(///2F8.0)') NXTIN2(JB),TINNX(JB)
          IF (INFLOW_CONST(JB)) THEN
            INC(JB) = NPT; NPT = NPT+1
            OPEN (INC(JB),FILE=CINFN(JB),STATUS='OLD')
            READ (INC(JB),'(///1000F8.0)') NXCIN2(JB),(CINNX(INCN(JC,JB),JB),JC=1,NACIN(JB))
          END IF
        END IF
        QIN(JB)    = QINNX(JB)
        QIND(JB)   = QINNX(JB)
        QINO(JB)   = QINNX(JB)
        TIN(JB)    = TINNX(JB)
        TIND(JB)   = TINNX(JB)
        TINO(JB)   = TINNX(JB)
        CIN(:,JB)  = CINNX(:,JB)
        CIND(:,JB) = CINNX(:,JB)
        CINO(:,JB) = CINNX(:,JB)
        IF (.NOT. INTERNAL_FLOW(JB) .AND. .NOT. DAM_INFLOW(JB)) THEN                                                  !TC 08/03/04  RA 1/13/06
          READ (INQ(JB), '(2F8.0)') NXQIN1(JB),QINNX(JB)
          READ (INFT(JB),'(2F8.0)') NXTIN1(JB),TINNX(JB)
          IF (INFLOW_CONST(JB)) THEN
            READ (INC(JB),'(1000F8.0)') NXCIN1(JB),(CINNX(INCN(JC,JB),JB),JC=1,NACIN(JB))
          END IF
        END IF
      END IF
      IF (DN_FLOW(JB)) THEN
        IF (NSTR(JB) > 0) THEN
          OTQ(JB) = NPT; NPT = NPT+1
          OPEN (OTQ(JB),FILE=QOTFN(JB),STATUS='OLD')
          READ (OTQ(JB),'(///10F8.0:/(8X,9F8.0))') NXQOT2(JB),(QSTRNX(JS,JB),JS=1,NSTR(JB))
          QSTR(:,JB)  = QSTRNX(:,JB)
          QSTRO(:,JB) = QSTRNX(:,JB)
          READ (OTQ(JB),'(10F8.0:/(8X,9F8.0))')    NXQOT1(JB),(QSTRNX(JS,JB),JS=1,NSTR(JB))
        END IF
      END IF
      IF (PRECIPITATION(JW)) THEN
        PRE(JB) = NPT; NPT = NPT+1
        PRT(JB) = NPT; NPT = NPT+1
        OPEN (PRE(JB),FILE=PREFN(JB),STATUS='OLD')
        OPEN (PRT(JB),FILE=TPRFN(JB),STATUS='OLD')
        READ (PRE(JB),'(///2F8.0)') NXPR2(JB), PRNX(JB)
        READ (PRT(JB),'(///2F8.0)') NXTPR2(JB),TPRNX(JB)
        IF (PRECIP_CONST(JB)) THEN
          PRC(JB) = NPT; NPT = NPT+1
          OPEN (PRC(JB),FILE=CPRFN(JB),STATUS='OLD')
          READ (PRC(JB),'(///1000F8.0)') NXCPR2(JB),(CPRNX(PRCN(JAC,JB),JB),JAC=1,NACPR(JB))
        END IF
        PR(JB)    = PRNX(JB)
        TPR(JB)   = TPRNX(JB)
        CPR(:,JB) = CPRNX(:,JB)
        READ (PRE(JB),'(2F8.0)') NXPR1(JB), PRNX(JB)
        READ (PRT(JB),'(2F8.0)') NXTPR1(JB),TPRNX(JB)
        IF (PRECIP_CONST(JB)) THEN
          READ (PRC(JB),'(1000F8.0)') NXCPR1(JB),(CPRNX(PRCN(JAC,JB),JB),JAC=1,NACPR(JB))
        END IF
      END IF
      IF (DIST_TRIBS(JB)) THEN
        DTQ(JB) = NPT; NPT = NPT+1
        DTT(JB) = NPT; NPT = NPT+1
        OPEN (DTQ(JB),FILE=QDTFN(JB),STATUS='OLD')
        OPEN (DTT(JB),FILE=TDTFN(JB),STATUS='OLD')
        READ (DTQ(JB),'(///2F8.0)') NXQDT2(JB),QDTRNX(JB)
        READ (DTT(JB),'(///2F8.0)') NXTDT2(JB),TDTRNX(JB)
        IF (DTRIB_CONST(JB)) THEN
          DTC(JB) = NPT; NPT = NPT+1
          OPEN (DTC(JB),FILE=CDTFN(JB),STATUS='OLD')
          READ (DTC(JB),'(///1000F8.0)') NXCDT2(JB),(CDTRNX(DTCN(JAC,JB),JB),JAC=1,NACDT(JB))
        END IF
        QDTR(JB)    = QDTRNX(JB)
        QDTRO(JB)   = QDTRNX(JB)
        TDTR(JB)    = TDTRNX(JB)
        TDTRO(JB)   = TDTRNX(JB)
        CDTR(:,JB)  = CDTRNX(:,JB)
        CDTRO(:,JB) = CDTRNX(:,JB)
        READ (DTQ(JB),'(2F8.0)') NXQDT1(JB),QDTRNX(JB)
        READ (DTT(JB),'(2F8.0)') NXTDT1(JB),TDTRNX(JB)
        IF (DTRIB_CONST(JB)) THEN
          READ (DTC(JB),'(1000F8.0)') NXCDT1(JB),(CDTRNX(DTCN(JAC,JB),JB),JAC=1,NACDT(JB))
        END IF
      END IF
      IF (UH_EXTERNAL(JB)) THEN
        UHE(JB) = NPT; NPT = NPT+1
        UHT(JB) = NPT; NPT = NPT+1
        OPEN (UHE(JB),FILE=EUHFN(JB),STATUS='OLD')
        OPEN (UHT(JB),FILE=TUHFN(JB),STATUS='OLD')
        READ (UHE(JB),'(///2F8.0)')              NXEUH2(JB), ELUHNX(JB)      
        READ (UHT(JB),'(///10F8.0:/(8X,9F8.0))') NXTUH2(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
        IF (CONSTITUENTS) THEN
          UHC(JB) = NPT; NPT = NPT+1
          OPEN (UHC(JB),FILE=CUHFN(JB),STATUS='OLD')
          READ (UHC(JB),'(//)')
          DO JAC=1,NAC
            IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (UHC(JB),'(10F8.0:/(8X,9F8.0))') NXCUH2(JB),(CUHNX(K,CN(JAC),JB),     &
                                                              K=2,KB(US(JB)))    
          END DO
        END IF
        ELUH(JB)     = ELUHNX(JB)
        ELUHO(JB)    = ELUHNX(JB)
        TUH(:,JB)    = TUHNX(:,JB)
        TUHO(:,JB)   = TUHNX(:,JB)
        CUH(:,:,JB)  = CUHNX(:,:,JB)
        CUHO(:,:,JB) = CUHNX(:,:,JB)
        READ (UHE(JB),'(2F8.0)')              NXEUH1(JB), ELUHNX(JB)
        READ (UHT(JB),'(10F8.0:/(8X,9F8.0))') NXTUH1(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
        IF (CONSTITUENTS) THEN
          DO JAC=1,NAC
            IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (UHC(JB),'(10F8.0:/(8X,9F8.0))') NXCUH1(JB),(CUHNX(K,CN(JAC),JB),     &
                                                              K=2,KB(US(JB)))
          END DO
        END IF
      END IF
      IF (DH_EXTERNAL(JB)) THEN
        DHE(JB) = NPT; NPT = NPT+1
        DHT(JB) = NPT; NPT = NPT+1
        OPEN (DHE(JB),FILE=EDHFN(JB),STATUS='OLD')
        OPEN (DHT(JB),FILE=TDHFN(JB),STATUS='OLD')
        READ (DHE(JB),'(///10F8.0)')             NXEDH2(JB),ELDHNX(JB)
        READ (DHT(JB),'(///10F8.0:/(8X,9F8.0))') NXTDH2(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
        IF (CONSTITUENTS) THEN
          DHC(JB) = NPT; NPT = NPT+1
          OPEN (DHC(JB),FILE=CDHFN(JB),STATUS='OLD')
          READ (DHC(JB),'(//)')
          DO JAC=1,NAC
            IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (DHC(JB),'(10F8.0:/(8X,9F8.0))') NXCDH2(JB),(CDHNX(K,CN(JAC),JB),     &
                                                              K=2,KB(DS(JB)))
          END DO
        END IF
        ELDH(JB)     = ELDHNX(JB)
        ELDHO(JB)    = ELDHNX(JB)
        TDH(:,JB)    = TDHNX(:,JB)
        TDHO(:,JB)   = TDHNX(:,JB)
        CDH(:,:,JB)  = CDHNX(:,:,JB)
        CDHO(:,:,JB) = CDHNX(:,:,JB)
        READ (DHE(JB),'(10F8.0)')             NXEDH1(JB),ELDHNX(JB)
        READ (DHT(JB),'(10F8.0:/(8X,9F8.0))') NXTDH1(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
        IF (CONSTITUENTS) THEN   
          DO JAC=1,NAC    
            IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (DHC(JB),'(10F8.0:/(8X,9F8.0))') NXCDH1(JB),(CDHNX(K,CN(JAC),JB),     &
                                                              K=2,KB(DS(JB)))    
          END DO
        END IF
      END IF
    END DO
  END DO
  IF (GATES) THEN
    GTQ = NPT; NPT = NPT+1
    OPEN (GTQ,FILE=QGTFN,STATUS='OLD')
    READ (GTQ,'(///1000F8.0)') NXQGT,(BGTNX(JG),JG=1,NGT)
    WHERE (DYNGTC == '     ZGT')
      EGT  = BGTNX
      BGT  = 1.0
      G1GT = 1.0
      G2GT = 1.0
    ELSEWHERE
      BGT = BGTNX
    END WHERE
    READ (GTQ,'(1000F8.0)')    NXQGT,(BGTNX(JG),JG=1,NGT)
  END IF
  NOPEN         = NPT-1
  DYNAMIC_SHADE = SHADEI < 0
RETURN

!***********************************************************************************************************************************
!**                                                  R E A D  I N P U T  D A T A                                                  **
!***********************************************************************************************************************************

ENTRY READ_INPUT_DATA (NXTVD)
  NXTVD = 1.0E10

! Meteorological data

  DO WHILE (JDAY >= NXWSC)
    WSC = WSCNX
    READ (WSH,'(10F8.0:/(8X,9F8.0))') NXWSC,(WSCNX(I),I=1,IMX)
  END DO
  DO JW=1,NWB
    DO WHILE (JDAY >= NXMET1(JW))
      TDEW(JW)   = TDEWNX(JW)
      TDEWO(JW)  = TDEWNX(JW)
      WIND(JW)   = WINDNX(JW)
      WINDO(JW)  = WINDNX(JW)
      PHI(JW)    = PHINX(JW)
      PHIO(JW)   = PHINX(JW)
      IF (PHISET > 0) PHI(JW)  = PHISET
      IF (PHISET > 0) PHIO(JW) = PHISET
      TAIR(JW)   = TAIRNX(JW)
      TAIRO(JW)  = TAIRNX(JW)
      CLOUD(JW)  = CLOUDNX(JW)
      CLOUDO(JW) = CLOUDNX(JW)
      NXMET2(JW) = NXMET1(JW)
      IF (READ_RADIATION(JW)) THEN
        SRON(JW)  = SRONX(JW)
        SROO(JW)  = SRON(JW)
        READ (MET(JW),'(7F8.0)') NXMET1(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW),SRONX(JW)
        SRONX(JW) = SRONX(JW)*REFL
      ELSE
        READ (MET(JW),'(6F8.0)') NXMET1(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW)
      END IF
    END DO
    NXTVD = MIN(NXTVD,NXMET1(JW))
    IF (READ_EXTINCTION(JW)) THEN
      DO WHILE (JDAY >= NXEXT1(JW))
        EXH2O(JW)  = EXTNX(JW)
        EXTO(JW)   = EXTNX(JW)
        NXEXT2(JW) = NXEXT1(JW)
        READ (EXT(JW),'(2F8.0)') NXEXT1(JW),EXTNX(JW)
      END DO     
    END IF     
    DO I=CUS(BS(JW)),DS(BE(JW))
      WIND2(I) = WIND(JW)*WSC(I)*LOG(2.0/0.003)/LOG(WINDH(JW)/0.003)
    END DO     
 END DO

! Withdrawals

  IF (NWD > 0) THEN
    DO WHILE (JDAY >= NXQWD1)
      NXQWD2 = NXQWD1
      DO JWD=1,NWD
        QWD(JWD)  = QWDNX(JWD)
        QWDO(JWD) = QWDNX(JWD)
      END DO
      READ (WDQ,'(10F8.0:/(8X,9F8.0))') NXQWD1,(QWDNX(JWD),JWD=1,NWD)
    END DO
    NXTVD = MIN(NXTVD,NXQWD1)
  END IF

! Tributaries

  IF (TRIBUTARIES) THEN
    DO JT=1,NTR

!**** Inflow

      DO WHILE (JDAY >= NXQTR1(JT))
        QTR(JT)    = QTRNX(JT)
        QTRO(JT)   = QTRNX(JT)
        NXQTR2(JT) = NXQTR1(JT)
        READ (TRQ(JT),'(2F8.0)') NXQTR1(JT),QTRNX(JT)
      END DO
      NXTVD = MIN(NXTVD,NXQTR1(JT))

!**** Inflow temperatures

      IF (JDAY >= NXTTR1(JT)) THEN
        DO WHILE (JDAY >= NXTTR1(JT))
          TTR(JT)    = TTRNX(JT)
          TTRO(JT)   = TTRNX(JT)
          NXTTR2(JT) = NXTTR1(JT)
          READ (TRT(JT),'(2F8.0)') NXTTR1(JT),TTRNX(JT)
        END DO
      END IF
      NXTVD = MIN(NXTVD,NXTTR1(JT))

!**** Inflow constituent concentrations

      IF (TRIB_CONST(JT)) THEN
        DO WHILE (JDAY >= NXCTR1(JT))
          CTR(TRCN(1:NACTR(JT),JT),JT)  = CTRNX(TRCN(1:NACTR(JT),JT),JT)
          CTRO(TRCN(1:NACTR(JT),JT),JT) = CTRNX(TRCN(1:NACTR(JT),JT),JT)
          NXCTR2(JT)                    = NXCTR1(JT)
          READ (TRC(JT),'(1000F8.0)') NXCTR1(JT),(CTRNX(TRCN(JAC,JT),JT),JAC=1,NACTR(JT))
        END DO
        NXTVD = MIN(NXTVD,NXCTR1(JT))
      END IF
    END DO
  END IF

! Branch related inputs

  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)

!**** Inflow

      IF (UP_FLOW(JB)) THEN
        IF (.NOT. INTERNAL_FLOW(JB) .AND. .NOT. DAM_INFLOW(JB)) THEN                                                  !TC 08/03/04 RA 1/13/06
          DO WHILE (JDAY >= NXQIN1(JB))
            QIND(JB)   = QINNX(JB)
            QINO(JB)   = QINNX(JB)
            NXQIN2(JB) = NXQIN1(JB)
            READ (INQ(JB),'(2F8.0)') NXQIN1(JB),QINNX(JB)
          END DO
          NXTVD = MIN(NXTVD,NXQIN1(JB))

!******** Inflow temperature

          DO WHILE (JDAY >= NXTIN1(JB))
            TIND(JB)   = TINNX(JB)
            TINO(JB)   = TINNX(JB)
            NXTIN2(JB) = NXTIN1(JB)
            READ (INFT(JB),'(2F8.0)') NXTIN1(JB),TINNX(JB)
          END DO
          NXTVD = MIN(NXTVD,NXTIN1(JB))

!******** Inflow constituent concentrations

          IF (INFLOW_CONST(JB)) THEN
            DO WHILE (JDAY >= NXCIN1(JB))
              CIND(INCN(1:NACIN(JB),JB),JB) = CINNX(INCN(1:NACIN(JB),JB),JB)
              CINO(INCN(1:NACIN(JB),JB),JB) = CINNX(INCN(1:NACIN(JB),JB),JB)
              NXCIN2(JB)                    = NXCIN1(JB)
              READ (INC(JB),'(1000F8.0)') NXCIN1(JB),(CINNX(INCN(JAC,JB),JB),JAC=1,NACIN(JB))
            END DO
            NXTVD = MIN(NXTVD,NXCIN1(JB))
          END IF
        END IF
      END IF

!**** Outflow

      IF (DN_FLOW(JB) .AND. NSTR(JB) > 0) THEN   
        DO WHILE (JDAY >= NXQOT1(JB))
          QSTR(1:NSTR(JB),JB)  = QSTRNX(1:NSTR(JB),JB)
          QSTRO(1:NSTR(JB),JB) = QSTRNX(1:NSTR(JB),JB)
          NXQOT2(JB)           = NXQOT1(JB)
          READ (OTQ(JB),'(10F8.0:/(8X,9F8.0))') NXQOT1(JB),(QSTRNX(JS,JB),JS=1,NSTR(JB))
        END DO
        NXTVD = MIN(NXTVD,NXQOT1(JB))
      END IF

!**** Distributed tributaries

      IF (DIST_TRIBS(JB)) THEN

!****** Inflow

        DO WHILE (JDAY >= NXQDT1(JB))
          QDTR(JB)   = QDTRNX(JB)
          QDTRO(JB)  = QDTRNX(JB)
          NXQDT2(JB) = NXQDT1(JB)
          READ (DTQ(JB),'(2F8.0)') NXQDT1(JB),QDTRNX(JB)
        END DO
        NXTVD = MIN(NXTVD,NXQDT1(JB))

!****** Temperature

        DO WHILE (JDAY >= NXTDT1(JB))
          TDTR(JB)   = TDTRNX(JB)
          TDTRO(JB)  = TDTRNX(JB)
          NXTDT2(JB) = NXTDT1(JB)
          READ (DTT(JB),'(2F8.0)') NXTDT1(JB),TDTRNX(JB)
        END DO
        NXTVD = MIN(NXTVD,NXTDT1(JB))

!****** Constituent concentrations

        IF (DTRIB_CONST(JB)) THEN
          DO WHILE (JDAY >= NXCDT1(JB))
            CDTR(DTCN(1:NACDT(JB),JB),JB)  = CDTRNX(DTCN(1:NACDT(JB),JB),JB)
            CDTRO(DTCN(1:NACDT(JB),JB),JB) = CDTRNX(DTCN(1:NACDT(JB),JB),JB)
            NXCDT2(JB)                     = NXCDT1(JB)
            READ (DTC(JB),'(1000F8.0)') NXCDT1(JB),(CDTRNX(DTCN(JAC,JB),JB),JAC=1,NACDT(JB))
          END DO
          NXTVD = MIN(NXTVD,NXCDT1(JB))
        END IF
      END IF

!**** Precipitation

      IF (PRECIPITATION(JW)) THEN
        DO WHILE (JDAY >= NXPR1(JB))
          PR(JB)    = PRNX(JB)
          NXPR2(JB) = NXPR1(JB)
          READ (PRE(JB),'(2F8.0)') NXPR1(JB),PRNX(JB)
        END DO
        NXTVD = MIN(NXTVD,NXPR1(JB))

!****** Temperature

        DO WHILE (JDAY >= NXTPR1(JB))
          TPR(JB)    = TPRNX(JB)
          NXTPR2(JB) = NXTPR1(JB)
          READ (PRT(JB),'(2F8.0)') NXTPR1(JB),TPRNX(JB)
        END DO
        NXTVD = MIN(NXTVD,NXTPR1(JB))

!****** Constituent concentrations

        IF (PRECIP_CONST(JB)) THEN
          DO WHILE (JDAY >= NXCPR1(JB))
            CPR(PRCN(1:NACPR(JB),JB),JB) = CPRNX(PRCN(1:NACPR(JB),JB),JB)
            NXCPR2(JB)                   = NXCPR1(JB)
            READ (PRC(JB),'(1000F8.0)') NXCPR1(JB),(CPRNX(PRCN(JAC,JB),JB),JAC=1,NACPR(JB))
          END DO
          NXTVD = MIN(NXTVD,NXCPR1(JB))
        END IF
      END IF

!**** Upstream head conditions

      IF (UH_EXTERNAL(JB)) THEN

!****** Elevations

        DO WHILE (JDAY >= NXEUH1(JB))
          ELUH(JB)   = ELUHNX(JB)
          ELUHO(JB)  = ELUHNX(JB)
          NXEUH2(JB) = NXEUH1(JB)
          READ (UHE(JB),'(2F8.0)') NXEUH1(JB),ELUHNX(JB)
        END DO
        NXTVD = MIN(NXTVD,NXEUH1(JB))

!****** Temperatures

        DO WHILE (JDAY >= NXTUH1(JB))
          DO K=2,KMX-1
            TUH(K,JB)  = TUHNX(K,JB)
            TUHO(K,JB) = TUHNX(K,JB)
          END DO
          NXTUH2(JB) = NXTUH1(JB)
          READ (UHT(JB),'(10F8.0:/(8X,9F8.0))') NXTUH1(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
        END DO
        NXTVD = MIN(NXTVD,NXTUH1(JB))

!****** Constituent concentrations

        IF (CONSTITUENTS) THEN
          DO WHILE (JDAY >= NXCUH1(JB))
            DO K=2,KMX-1
              CUH(K,CN(1:NAC),JB)  = CUHNX(K,CN(1:NAC),JB)
              CUHO(K,CN(1:NAC),JB) = CUHNX(K,CN(1:NAC),JB)
            END DO
            NXCUH2(JB) = NXCUH1(JB)
            DO JAC=1,NAC
              IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (UHC(JB),'(10F8.0:/(8X,9F8.0))') NXCUH1(JB),(CUHNX(K,CN(JAC),JB),   &
                                                                K=2,KB(US(JB)))
            END DO
          END DO
          NXTVD = MIN(NXTVD,NXCUH1(JB))
        END IF
      END IF

!**** Downstream head

      IF (DH_EXTERNAL(JB)) THEN

!****** Elevation

        DO WHILE (JDAY >= NXEDH1(JB))
          ELDH(JB)   = ELDHNX(JB)
          ELDHO(JB)  = ELDHNX(JB)
          NXEDH2(JB) = NXEDH1(JB)
          READ (DHE(JB),'(2F8.0)') NXEDH1(JB),ELDHNX(JB)
        END DO
        NXTVD = MIN(NXTVD,NXEDH1(JB))

!****** Temperature

        DO WHILE (JDAY >= NXTDH1(JB))
          DO K=2,KMX-1
            TDH(K,JB)  = TDHNX(K,JB)
            TDHO(K,JB) = TDHNX(K,JB)
          END DO
          NXTDH2(JB) = NXTDH1(JB)
          READ (DHT(JB),'(10F8.0:/(8X,9F8.0))') NXTDH1(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
        END DO
        NXTVD = MIN(NXTVD,NXTDH1(JB))

!****** Constituents

        IF (CONSTITUENTS) THEN
          DO WHILE (JDAY >= NXCDH1(JB))
            DO K=2,KMX-1
              CDH(K,CN(1:NAC),JB)  = CDHNX(K,CN(1:NAC),JB)
              CDHO(K,CN(1:NAC),JB) = CDHNX(K,CN(1:NAC),JB)
            END DO
            NXCDH2(JB) = NXCDH1(JB)
            DO JAC=1,NAC
              IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (DHC(JB),'(10F8.0:/(8X,9F8.0))') NXCDH1(JB),(CDHNX(K,CN(JAC),JB),   &
                                                                K=2,KB(DS(JB)))
            END DO
          END DO
          NXTVD = MIN(NXTVD,NXCDH1(JB))
        END IF
      END IF
    END DO
  END DO

! Gate height opening

  IF (GATES) THEN
    DO WHILE (JDAY >= NXQGT)
      WHERE (DYNGTC == '     ZGT')
        EGT  = BGTNX  
        BGT  = 1.0
        G1GT = 1.0
        G2GT = 1.0
      ELSEWHERE
        BGT = BGTNX
      ENDWHERE
      READ (GTQ,'(1000F8.0)') NXQGT,(BGTNX(JG),JG=1,NGT)
    END DO
    NXTVD = MIN(NXTVD,NXQGT)
  END IF

! Dead sea case

  DO JW=1,NWB
    IF (NO_INFLOW(JW)) THEN
      QIN(BS(JW):BE(JW))    = 0.0
      QINO(BS(JW):BE(JW))   = 0.0
      QIND(BS(JW):BE(JW))   = 0.0
      QINNX(BS(JW):BE(JW))  = 0.0
      QDTR(BS(JW):BE(JW))   = 0.0
      QDTRO(BS(JW):BE(JW))  = 0.0
      QDTRNX(BS(JW):BE(JW)) = 0.0
      PR(BS(JW):BE(JW))     = 0.0
      PRNX(BS(JW):BE(JW))   = 0.0
    END IF
    IF (NO_OUTFLOW(JW)) THEN
      QSTR(:,BS(JW):BE(JW))   = 0.0
      QSTRO(:,BS(JW):BE(JW))  = 0.0
      QSTRNX(:,BS(JW):BE(JW)) = 0.0
    END IF
  END DO
  WHERE (NO_WIND)
    WIND   = 0.0
    WINDO  = 0.0
    WINDNX = 0.0
  ENDWHERE
  WHERE (READ_RADIATION .AND. NO_HEAT)
    SRON  = 0.0
    SROO  = 0.0
    SRONX = 0.0
  ENDWHERE
  IF (ANY(NO_INFLOW)) THEN
    QTR   = 0.0
    QTRO  = 0.0
    QTRNX = 0.0
    QWD   = 0.0
    QWDO  = 0.0
    QWDNX = 0.0
  END IF
RETURN

!***********************************************************************************************************************************
!**                                              I N T E R P O L A T E  I N P U T S                                               **
!***********************************************************************************************************************************

ENTRY INTERPOLATE_INPUTS

! Meteorological/light extinction data

  DO JW=1,NWB
    IF (INTERP_METEOROLOGY(JW)) THEN
      RATIO     = (NXMET1(JW)-JDAY)/(NXMET1(JW)-NXMET2(JW))
      TDEW(JW)  = (1.0-RATIO)*TDEWNX(JW)+RATIO*TDEWO(JW)
      WIND(JW)  = (1.0-RATIO)*WINDNX(JW)+RATIO*WINDO(JW)
      IF (ABS(PHIO(JW)-PHINX(JW)) > PI) THEN
        PHI(JW) = (1.0-RATIO)*(PHINX(JW)+2.0*PI)+RATIO*PHIO(JW)
      ELSE
        PHI(JW) = (1.0-RATIO)*PHINX(JW)+RATIO*PHIO(JW)
      END IF
      TAIR(JW)  = (1.0-RATIO)*TAIRNX(JW) +RATIO*TAIRO(JW)
      CLOUD(JW) = (1.0-RATIO)*CLOUDNX(JW)+RATIO*CLOUDO(JW)
      IF (READ_RADIATION(JW)) SRON(JW) = (1.0-RATIO)*SRONX(JW)+RATIO*SROO(JW)
    END IF
    IF (INTERP_EXTINCTION(JW)) THEN
      RATIO     = (NXEXT1(JW)-JDAY)/(NXEXT1(JW)-NXEXT2(JW))
      EXH2O(JW) = (1.0-RATIO)*EXTNX(JW)+RATIO*EXTO(JW)
    END IF
  END DO

! Withdrawals

  IF (NWD > 0) THEN
    QRATIO = (NXQWD1-JDAY)/(NXQWD1-NXQWD2)
    DO JWD=1,NWD
      IF (INTERP_WITHDRAWAL(JWD)) QWD(JWD) = (1.0-QRATIO)*QWDNX(JWD)+QRATIO*QWDO(JWD)
    END DO
  END IF

! Tributaries

  IF (NTR > 0) THEN
    DO JT=1,NTR
      IF (INTERP_TRIBS(JT)) THEN
        QRATIO = (NXQTR1(JT)-JDAY)/(NXQTR1(JT)-NXQTR2(JT))
        TRATIO = (NXTTR1(JT)-JDAY)/(NXTTR1(JT)-NXTTR2(JT))
        IF (TRIB_CONST(JT)) CRATIO = (NXCTR1(JT)-JDAY)/(NXCTR1(JT)-NXCTR2(JT))
        QTR(JT)                      = (1.0-QRATIO)*QTRNX(JT)                     +QRATIO*QTRO(JT)
        TTR(JT)                      = (1.0-TRATIO)*TTRNX(JT)                     +TRATIO*TTRO(JT)
        CTR(TRCN(1:NACTR(JT),JT),JT) = (1.0-CRATIO)*CTRNX(TRCN(1:NACTR(JT),JT),JT)+CRATIO*CTRO(TRCN(1:NACTR(JT),JT),JT)
      END IF
    END DO
  END IF

! Branch related inputs

  DO JB=1,NBR

!** Inflow

    IF (UP_FLOW(JB)) THEN
      IF (.NOT. INTERNAL_FLOW(JB) .AND. .NOT. DAM_INFLOW(JB)) THEN                                                    !TC 08/03/04 RA 1/13/06
        IF (INTERP_INFLOW(JB)) THEN
          QRATIO = (NXQIN1(JB)-JDAY)/(NXQIN1(JB)-NXQIN2(JB))
          TRATIO = (NXTIN1(JB)-JDAY)/(NXTIN1(JB)-NXTIN2(JB))
          IF (INFLOW_CONST(JB))  CRATIO = (NXCIN1(JB)-JDAY)/(NXCIN1(JB)-NXCIN2(JB))
          QIND(JB)                      = (1.0-QRATIO)*QINNX(JB)                  +QRATIO*QINO(JB)
          TIND(JB)                      = (1.0-TRATIO)*TINNX(JB)                  +TRATIO*TINO(JB)
          CIND(INCN(1:NACIN(JB),JB),JB) = (1.0-CRATIO)*CINNX(INCN(1:NACIN(JB),JB),JB)+CRATIO*CINO(INCN(1:NACIN(JB),JB),JB)
        END IF
      END IF
    END IF

!** Outflow

    IF (DN_FLOW(JB) .AND. NSTR(JB) > 0) THEN   
      QRATIO = (NXQOT1(JB)-JDAY)/(NXQOT1(JB)-NXQOT2(JB))
      DO JS=1,NSTR(JB)
        IF (INTERP_OUTFLOW(JS,JB)) QSTR(JS,JB) = (1.0-QRATIO)*QSTRNX(JS,JB)+QRATIO*QSTRO(JS,JB)
      END DO
    END IF

!** Distributed tributaries

    IF (DIST_TRIBS(JB)) THEN
      IF (INTERP_DTRIBS(JB)) THEN
        QRATIO = (NXQDT1(JB)-JDAY)/(NXQDT1(JB)-NXQDT2(JB))
        TRATIO = (NXTDT1(JB)-JDAY)/(NXTDT1(JB)-NXTDT2(JB))
        IF (DTRIB_CONST(JB))   CRATIO = (NXCDT1(JB)-JDAY)/(NXCDT1(JB)-NXCDT2(JB))
        QDTR(JB)                      = (1.0-QRATIO)*QDTRNX(JB)                  +QRATIO*QDTRO(JB)
        TDTR(JB)                      = (1.0-TRATIO)*TDTRNX(JB)                  +TRATIO*TDTRO(JB)
        CDTR(DTCN(1:NACDT(JB),JB),JB) = (1.0-CRATIO)*CDTRNX(DTCN(1:NACDT(JB),JB),JB)+CRATIO*CDTRO(DTCN(1:NACDT(JB),JB),JB)
      END IF
    END IF

!** Upstream head

    IF (UH_EXTERNAL(JB)) THEN
      IF (INTERP_HEAD(JB)) THEN
        HRATIO   = (NXEUH1(JB)-JDAY)/(NXEUH1(JB)-NXEUH2(JB))
        TRATIO   = (NXTUH1(JB)-JDAY)/(NXTUH1(JB)-NXTUH2(JB))
        IF (CONSTITUENTS) CRATIO = (NXCUH1(JB)-JDAY)/(NXCUH1(JB)-NXCUH2(JB))
        ELUH(JB) = (1.0-HRATIO)*ELUHNX(JB)+HRATIO*ELUHO(JB)
        DO K=2,KMX-1
          TUH(K,JB)           = (1.0-TRATIO)*TUHNX(K,JB)          +TRATIO*TUHO(K,JB)
          CUH(K,CN(1:NAC),JB) = (1.0-CRATIO)*CUHNX(K,CN(1:NAC),JB)+CRATIO*CUHO(K,CN(1:NAC),JB)
        END DO
      END IF
    END IF

!** Downstream head

    IF (DH_EXTERNAL(JB)) THEN
      IF (INTERP_HEAD(JB)) THEN
        HRATIO = (NXEDH1(JB)-JDAY)/(NXEDH1(JB)-NXEDH2(JB))
        TRATIO = (NXTDH1(JB)-JDAY)/(NXTDH1(JB)-NXTDH2(JB))
        IF (CONSTITUENTS) CRATIO = (NXCDH1(JB)-JDAY)/(NXCDH1(JB)-NXCDH2(JB))
        ELDH(JB) = (1.0-HRATIO)*ELDHNX(JB)+HRATIO*ELDHO(JB)
        DO K=2,KMX-1
          TDH(K,JB)           = (1.0-TRATIO)*TDHNX(K,JB)          +TRATIO*TDHO(K,JB)
          CDH(K,CN(1:NAC),JB) = (1.0-CRATIO)*CDHNX(K,CN(1:NAC),JB)+CRATIO*CDHO(K,CN(1:NAC),JB)
        END DO
      END IF
    END IF
  END DO
RETURN
ENTRY DEALLOCATE_TIME_VARYING_DATA
  DEALLOCATE (NXQTR1, NXTTR1, NXCTR1, NXQIN1, NXTIN1, NXCIN1, NXQDT1, NXTDT1, NXCDT1, NXPR1,  NXTPR1, NXCPR1, NXEUH1, NXTUH1)
  DEALLOCATE (NXCUH1, NXEDH1, NXTDH1, NXCDH1, NXQOT1, NXMET1, NXQTR2, NXTTR2, NXCTR2, NXQIN2, NXTIN2, NXCIN2, NXQDT2, NXTDT2)
  DEALLOCATE (NXCDT2, NXPR2,  NXTPR2, NXCPR2, NXEUH2, NXTUH2, NXCUH2, NXEDH2, NXTDH2, NXCDH2, NXQOT2, NXMET2, WSCNX)
  DEALLOCATE (QDTRO,  TDTRO,  ELUHO,  ELDHO,  QWDO,   QTRO,   TTRO,   QINO,   TINO,   QDTRNX, TDTRNX, PRNX,   TPRNX,  ELUHNX)
  DEALLOCATE (ELDHNX, QWDNX,  QTRNX,  TTRNX,  QINNX,  TINNX,  SROO,   TAIRO,  TDEWO,  CLOUDO, PHIO,   WINDO,  TAIRNX, BGTNX)
  DEALLOCATE (TDEWNX, CLOUDNX,PHINX,  WINDNX, SRONX,  TRQ,    TRT,    TRC,    INQ,    DTQ,    PRE,    UHE,    DHE,    INFT)
  DEALLOCATE (DTT,    PRT,    UHT,    DHT,    INC,    DTC,    PRC,    UHC,    DHC,    OTQ,    MET,    EXT,    EXTNX,  EXTO)
  DEALLOCATE (NXEXT1, NXEXT2, CTRO,   CINO,   QOUTO,  CDTRO,  TUHO,   TDHO,   QSTRO,  CTRNX,  CINNX,  QOUTNX, CDTRNX, CPRNX)
  DEALLOCATE (TUHNX, TDHNX,   QSTRNX, CUHO,   CDHO,   CUHNX,  CDHNX)
  DEALLOCATE (INFLOW_CONST,   TRIB_CONST,     DTRIB_CONST,    PRECIP_CONST)
RETURN
END SUBROUTINE TIME_VARYING_DATA

!***********************************************************************************************************************************
!**                                        S U B R O U T I N E   H E A T  E X C H A N G E                                         **
!***********************************************************************************************************************************

SUBROUTINE HEAT_EXCHANGE
  USE GLOBAL; USE GDAYC; USE SURFHE; USE TVDC; USE SHADEC 

! Type declaration

  REAL :: JDAY, LOCAL, MPS_TO_MPH

! Data declaration

  DATA MPS_TO_MPH          /2.23714/, W_M2_TO_BTU_FT2_DAY /7.60796/, FLUX_BR_TO_FLUX_SI /0.23659/
  DATA BTU_FT2_DAY_TO_W_M2 /0.1314/
  DATA BOWEN_CONSTANT      /0.47/

! Function declaration

  DEG_F(X) =  X*1.8+32.0
  DEG_C(X) = (X-32.0)*5.0/9.0
RETURN

!***********************************************************************************************************************************
!**                                            S H O R T  W A V E  R A D I A T I O N                                              **
!***********************************************************************************************************************************

ENTRY SHORT_WAVE_RADIATION (JDAY)
  LOCAL    =  LONGIT(JW)
  STANDARD =  15.0*INT(LONGIT(JW)/15.0)
  HOUR     = (JDAY-INT(JDAY))*24.0
  IDAY     =  JDAY-((INT(JDAY/365))*365)
  IDAY     =  IDAY+INT(INT(JDAY/365)/4)
  TAUD     = (2*PI*(IDAY-1))/365
  EQTNEW   =  0.170*SIN(4*PI*(IDAY-80)/373)-0.129*SIN(2*PI*(IDAY-8)/355)
  HH(JW)   =  0.261799*(HOUR-(LOCAL-STANDARD)*0.0666667+EQTNEW-12.0)
  DECL(JW) =  0.006918-0.399912*COS(TAUD)+0.070257*SIN(TAUD)-0.006758*COS(2*TAUD)+0.000907*SIN(2*TAUD)-0.002697*COS(3*TAUD)        &
              +0.001480*SIN(3*TAUD)   
  SINAL    =  SIN(LAT(JW)*.0174533)*SIN(DECL(JW))+COS(LAT(JW)*.0174533)*COS(DECL(JW))*COS(HH(JW))
  A00(JW)  =  57.2957795*ASIN(SINAL)
  A0       =  A00(JW)
  IF (A0 > 0.0) THEN
    SRON(JW) = (1.0-0.0065*CLOUD(JW)**2)*24.0*(2.044*A0+0.1296*A0**2-1.941E-3*A0**3+7.591E-6*A0**4)*BTU_FT2_DAY_TO_W_M2
  ELSE
    SRON(JW) = 0.0
  END IF
RETURN

!***********************************************************************************************************************************
!**                                          E Q U I L I B R I U M  T E M P E R A T U R E                                         **
!***********************************************************************************************************************************

ENTRY EQUILIBRIUM_TEMPERATURE

! British units

  TDEW_F   = DEG_F(TDEW(JW))
  TAIR_F   = DEG_F(TAIR(JW))
  SRO_BR   = SRON(JW)*W_M2_TO_BTU_FT2_DAY*SHADE(I)
  WIND_MPH = WIND(JW)*WSC(I)*MPS_TO_MPH
  WIND2M   = WIND_MPH*LOG(2.0/0.003)/LOG(WINDH(JW)/0.003)+NONZERO
  ACONV    = W_M2_TO_BTU_FT2_DAY
  IF (CFW(JW) == 1.0) BCONV = 3.401062
  IF (CFW(JW) == 2.0) BCONV = 1.520411

! Equilibrium temperature and heat exchange coefficient

  ET(I)   =  TDEW_F
  TSTAR   = (ET(I)+TDEW_F)*0.5
  BETA    =  0.255-(8.5E-3*TSTAR)+(2.04E-4*TSTAR*TSTAR)
  FW      =  ACONV*AFW(JW)+BCONV*BFW(JW)*WIND2M**CFW(JW)
  CSHE(I) =  15.7+(0.26+BETA)*FW
  RA      =  3.1872E-08*(TAIR_F+459.67)**4
  ETP     = (SRO_BR+RA-1801.0)/CSHE(I)+(CSHE(I)-15.7)*(0.26*TAIR_F+BETA*TDEW_F)/(CSHE(I)*(0.26+BETA))
  J       =  0
  DO WHILE (ABS(ETP-ET(I)) > 0.05 .AND. J < 10)
    ET(I)   =  ETP
    TSTAR   = (ET(I)+TDEW_F)*0.5
    BETA    =  0.255-(8.5E-3*TSTAR)+(2.04E-4*TSTAR*TSTAR)
    CSHE(I) =  15.7+(0.26+BETA)*FW
    ETP     = (SRO_BR+RA-1801.0)/CSHE(I)+(CSHE(I)-15.7)*(0.26*TAIR_F+BETA*TDEW_F)/(CSHE(I)*(0.26+BETA))
    J       =  J+1
  END DO

! SI units

  ET(I)   = DEG_C(ET(I))
  CSHE(I) = CSHE(I)*FLUX_BR_TO_FLUX_SI/RHOWCP
RETURN

!***********************************************************************************************************************************
!**                                                   S U R F A C E   T E R M S                                                   **
!***********************************************************************************************************************************

ENTRY SURFACE_TERMS (TSUR)

! Partial water vapor pressure of air (mm hg)

  EA = EXP(2.3026*(9.5*TDEW(JW)/(TDEW(JW)+265.5)+0.6609))
  IF (TDEW(JW) > 0.0) EA = EXP(2.3026*(7.5*TDEW(JW)/(TDEW(JW)+237.3)+0.6609))

! Partial water vapor pressure at the water surface

  ES = EXP(2.3026*(9.5*TSUR/(TSUR+265.5)+0.6609))
  IF (TSUR > 0.0) ES = EXP(2.3026*(7.5*TSUR/(TSUR+237.3)+0.6609))

! Wind function

  IF (RH_EVAP(JW)) THEN
    TAIRV = (TAIR(JW)+273.0)/(1.0-0.378*EA/760.0)
    DTV   = (TSUR+273.0)/(1.0-0.378*ES/760.0)-TAIRV
    DTVL  =  0.0084*WIND2(I)**3
    IF (DTV < DTVL) DTV = DTVL
    FW = (3.59*DTV**0.3333+4.26*WIND2(I))
  ELSE
    FW = AFW(JW)+BFW(JW)*WIND2(I)**CFW(JW)
  END IF

! Evaporative flux

  RE(I) = FW*(ES-EA)

! Conductive flux

  RC(I) = FW*BOWEN_CONSTANT*(TSUR-TAIR(JW))

! Back radiation flux

  RB(I) = 5.51E-8*(TSUR+273.15)**4
END SUBROUTINE HEAT_EXCHANGE

!***********************************************************************************************************************************
!**                                                S U B R O U T I N E   S H A D I N G                                            **
!***********************************************************************************************************************************

SUBROUTINE SHADING
  USE SHADEC; USE GLOBAL; USE GDAYC; USE SURFHE; USE GEOMC; USE SCREENC; USE LOGICC
  CHARACTER(1) :: BANK
  REAL         :: LOCAL

! Calculate solar altitude, declination, and local hour angle when short-wave solar radiation is provided as input

  IF (READ_RADIATION(JW)) THEN
    LOCAL    =  LONGIT(JW)
    STANDARD =  15.0*INT(LONGIT(JW)/15.0)
    HOUR     = (JDAY-INT(JDAY))*24.0
    IDAY     =  JDAY-((INT(JDAY/365))*365)
    IDAY     =  IDAY+INT(INT(JDAY/365)/4)
    TAUD     = (2*PI*(IDAY-1))/365
    EQTNEW   =  0.170*SIN(4*PI*(IDAY-80)/373)-0.129*SIN(2*PI*(IDAY-8)/355)
    HH(JW)   =  0.261799*(HOUR-(LOCAL-STANDARD)*0.0666667+EQTNEW-12.0)
    DECL(JW) =  0.006918-0.399912*COS(TAUD)+0.070257*SIN(TAUD)-0.006758*COS(2*TAUD)+0.000907*SIN(2*TAUD)-0.002697*COS(3*TAUD)      &
                +0.001480*SIN(3*TAUD)
    SINAL    =  SIN(LAT(JW)*.0174533)*SIN(DECL(JW))+COS(LAT(JW)*.0174533)*COS(DECL(JW))*COS(HH(JW))
    A00(JW)  =  57.2957795*ASIN(SINAL)
  END IF     

! If the sun is below the horizon, set SHADE(I) to 0.

  IF (A00(JW) < 0.0) THEN
    SHADE(I) = 0.0
  ELSE

!** Calculate solar azimuth angle

    A02 = A00(JW)/57.2957795
    AX  = (SIN(DECL(JW))*COS(LAT(JW)*0.017453)-COS(DECL(JW))*COS(HH(JW))*SIN(LAT(JW)*0.017453))/COS(A02)
    IF (AX >  1.0) AX =  1.0
    IF (AX < -1.0) AX = -1.0
    AZT = ACOS(AX)
    IF (HH(JW) < 0.0) THEN
     AZ00 = AZT
    ELSE
     AZ00 = 2.0*PI-AZT
    END IF
    A0 = A02

!** Interpolate the topographic shade angle

    DO J=1,IANG-1
      IF (AZ00 > ANG(J) .AND. AZ00 <= ANG(J+1)) THEN
        ANG1    =  AZ00-ANG(J)
        ANG2    = (TOPO(I,J+1)-TOPO(I,J))/GAMA                 ! SW 10/17/05
        TOPOANG =  TOPO(I,J)+ANG2*ANG1
      END IF
    END DO
    IF (AZ00 > ANG(IANG) .AND. AZ00 <= 2*PI) THEN
      ANG1    =  AZ00-ANG(IANG)
      ANG2    = (TOPO(I,1)-TOPO(I,IANG))/GAMA                  ! SW 10/17/05
      TOPOANG =  TOPO(I,IANG)+ANG2*ANG1
    END IF

!** Complete topographic shading if solar altitude less than topo angle

    IF (A0 <= TOPOANG) THEN
      SFACT = 0.90
      GO TO 100
    END IF

!** No vegetative shading if azimuth angle is oriented parallel to stream

    IF (AZ00 == PHI0(I) .OR. AZ00 == PHI0(I)+PI .OR. AZ00+PI == PHI0(I)) THEN
      SFACT = 0.0
      GO TO 100
    END IF

!** Bank with the controlling vegetation

    IF (PHI0(I) > 0.0 .AND. PHI0(I) <= PI) THEN
      IF (AZ00 > PHI0(I)     .AND. AZ00 <= PHI0(I)+PI) BANK = 'L'
      IF (AZ00 > 0.0         .AND. AZ00 <= PHI0(I))    BANK = 'R'
      IF (AZ00 > PHI0(I)+PI  .AND. AZ00 <  2.0*PI)     BANK = 'R'
    ELSE IF (PHI0(I) > PI .AND. PHI0(I) <= 2.0*PI) THEN
      IF (AZ00 >= PHI0(I)    .AND. AZ00 < 2.0*PI)      BANK = 'L'
      IF (AZ00 >= 0.0        .AND. AZ00 < PHI0(I)-PI)  BANK = 'L'
      IF (AZ00 >= PHI0(I)-PI .AND. AZ00 < PHI0(I))     BANK = 'R'
    END IF

!** No topographic shading

    WLELEV = EL(KT,I)-Z(I)*COS(ALPHA(JB))
    IF (BANK == 'L') THEN
      IF (TTLB(I) < WLELEV) THEN
        SFACT = 0.0
        GO TO 100
      ELSE
        HT    = TTLB(I)-WLELEV
        CLINE = CLLB(I)
        SRED  = SRLB2(I)
        IF (JDAYG > SRFJD1(I) .AND. JDAYG <= SRFJD2(I)) SRED = SRLB1(I)
      END IF
    ELSE
      IF (TTRB(I) < WLELEV) THEN
        SFACT = 0.0
        GO TO 100
      ELSE
        HT    = TTRB(I)-WLELEV
        CLINE = CLRB(I)
        SRED  = SRRB2(I)
        IF (JDAYG > SRFJD1(I) .AND. JDAYG <= SRFJD2(I)) SRED = SRRB1(I)
      END IF
    END IF
    STLEN = HT/TAN (A0)
    EDGE  = MAX (0.0,CLINE-BI(KT,I)/2.0)

!** Distance from vegetation to water edge on line parallel to azimuth

    EDAZ = EDGE/ABS(SIN(PHI0(I)-AZ00))
    IF (STLEN <= EDAZ) THEN
      SFACT = 0.0
      GO TO 100
    END IF

!** Distance shadow extends over water (perpendicular to segment orientation)

    SN    = MIN (HT*ABS (SIN (ABS (PHI0(I)-AZ00)))/TAN (A0)-EDGE,BI(KT,I))
    SFACT = SRED*SN/BI(KT,I)
100 CONTINUE
    SHADE(I) = MAX (0.0,1-SFACT)
  END IF
  RETURN
END SUBROUTINE SHADING

!***********************************************************************************************************************************
!**                                            S U B R O U T I N E   W I T H D R A W A L                                          **
!***********************************************************************************************************************************

SUBROUTINE WITHDRAWAL
  USE GLOBAL; USE GEOMC; USE TVDC; USE SELWC; USE LOGICC
RETURN

!***********************************************************************************************************************************
!**                                             D O W N S T R E A M   W I T H D R A W A L                                         **
!***********************************************************************************************************************************

ENTRY DOWNSTREAM_WITHDRAWAL (JS)

! Variable initialization

  HSWT = 0.0; HSWB = 0.0; VNORM = 0.0; QNEW = 0.0

! Water surface elevation

  ELR  = SINA(JB)*DLX(ID)*0.5
  WSEL = EL(KT,ID)-Z(ID)*COSA(JB)-ELR

! Structure layer

  DO K=KT,KB(ID)
    IF (EL(K,ID)-ELR < ESTR(JS,JB)) EXIT
  END DO
  KSTR = MAX(K-1,KT)
  KSTR = MIN(KSTR,KB(ID))

! Initial withdrawal limits

  KTOP = MAX(KTSW(JS,JB),KT)
  IF (KSTR < KTOP) KTOP = KSTR
  KBOT = MIN(KBSW(JS,JB),KB(ID))
  IF (KBOT <= KT .AND. KBOT /= KB(ID)) KBOT = KT+1
  IF (KBOT > KB(ID)) KBOT = KB(ID)
  ELSTR = ESTR(JS,JB)
  IF (ESTR(JS,JB) <= EL(KB(ID)+1,ID+1)-ELR) THEN
    KSTR  = KB(ID)
    ELSTR = EL(KB(ID),ID)-ELR
  END IF
  IF (ESTR(JS,JB) > EL(KT,ID)-ELR) ELSTR = WSEL
  IF (KBSW(JS,JB) < KSTR) THEN
    KSTR  = KT
    ELSTR = WSEL
  END IF

! Boundary interference

  COEF = 1.0
  IF ((WSEL-EL(KBOT,ID)-ELR) /= 0.0) THEN
    RATIO = (ELSTR-(EL(KBOT,ID)-ELR))/(WSEL-(EL(KBOT,ID)-ELR))
    IF (RATIO < 0.1 .OR. RATIO > 0.9) COEF = 2.0
  END IF

! Withdrawal zone above structure

  DO K=KSTR-1,KTOP,-1

!** Density frequency

    HT    = (EL(K,ID)-ELR)-ELSTR  
    RHOFT = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/(HT*RHO(KSTR,ID)+NONZERO)*G),NONZERO)

!** Thickness

    IF (POINT_SINK(JS,JB)) THEN
      HSWT = (COEF*QSTR(JS,JB)/RHOFT)**0.333333
    ELSE
      HSWT = SQRT(2.0*COEF*QSTR(JS,JB)/(WSTR(JS,JB)*RHOFT))
    END IF
    IF (HT >= HSWT) THEN
      KTOP = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELSTR+HSWT) < WSEL) THEN
    DLRHOT = ABS(RHO(KSTR,ID)-RHO(KTOP,ID))
  ELSE IF (WSEL == ELSTR) THEN
    DLRHOT = NONZERO
  ELSE
    DLRHOT = ABS(RHO(KSTR,ID)-RHO(KT,ID))*HSWT/(WSEL-ELSTR)
  END IF
  DLRHOT = MAX(DLRHOT,NONZERO)

! Withdrawal zone below structure

  DO K=KSTR+1,KBOT

!** Density frequency

    HB    = ELSTR-(EL(K,ID)-ELR)
    RHOFB = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/(HB*RHO(KSTR,ID)+NONZERO)*G),NONZERO)

!** Thickness

    IF (POINT_SINK(JS,JB)) THEN
      HSWB = (COEF*QSTR(JS,JB)/RHOFB)**0.333333
    ELSE
      HSWB = SQRT(2.0*COEF*QSTR(JS,JB)/(WSTR(JS,JB)*RHOFB))
    END IF
    IF (HB >= HSWB) THEN
      KBOT = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELSTR-HSWB) > EL(KBOT+1,ID)) THEN
    DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))
  ELSE
    DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))*HSWB/(ELSTR-(EL(KBOT+1,ID)-ELR))
  END IF
  DLRHOB = MAX(DLRHOB,NONZERO)
 
! Velocity profile

  VSUM     = 0.0
  DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)
  DO K=KTOP,KBOT
    VNORM(K) = ABS(1.0-((RHO(K,ID)-RHO(KSTR,ID))/DLRHOMAX)**2)*BHR2(K,ID) 
    VSUM     = VSUM+VNORM(K)
  END DO

! Outflows

  DO K=KTOP,KBOT
    QNEW(K)    = (VNORM(K)/VSUM)*QSTR(JS,JB)
    QOUT(K,JB) =  QOUT(K,JB)+QNEW(K)
  END DO

! Inactive layers and total outflow

  IF (JS == NST) THEN
    WHERE (QOUT(:,JB) == 0.0) U(:,ID) = 0.0
  END IF
RETURN

!***********************************************************************************************************************************
!**                                                L A T E R A L   W I T H D R A W A L                                            **
!***********************************************************************************************************************************

ENTRY LATERAL_WITHDRAWAL (JWD)

! Variable initialization

  VNORM = 0.0; QSW(:,JWD) = 0.0; HWDT = 0.0; HWDB = 0.0

! Structure layer

  K = KT
  DO K=KT,KB(I)
    IF (EL(K,I) < EWD(JWD)) EXIT
  END DO
  KWD = MAX(K-1,KT)
  KWD = MIN(KWD,KB(I))

! Initial withdrawal limits

  KTOP = MAX(KTWD(JWD),KT)
  IF (KWD < KTOP) KTOP = KWD
  KBOT = MIN(KBWD(JWD),KB(I))
  IF (KBOT <= KT .AND. KB(I) /= KBOT) KBOT = KT+1
  IF (KBOT > KB(I)) KBOT = KB(I)
  ELWD = EWD(JWD)
  IF (EWD(JWD) <= EL(KB(I)+1,I)) THEN
    KWD  = KB(I)
    ELWD = EL(KB(I),I)
  END IF
  IF (EWD(JWD) > EL(KT,I)) ELWD = EL(KT,I)
  IF (KBWD(JWD) < KWD) THEN
    KWD  = KT
    ELWD = EL(KT,I)
  END IF

! Boundary interference

  COEF = 1.0
  IF (KT /= KBOT) THEN
    RATIO = (ELWD-EL(KBOT,I))/(EL(KT,I)-EL(KBOT,I))
    IF (RATIO < 0.1 .OR. RATIO > 0.9) COEF = 2.0
  END IF

! Withdrawal zone above structure

  DO K=KWD-1,KTOP,-1

!** Density frequency

    HT    = EL(K,I)-ELWD
    RHOFT = MAX(SQRT((ABS(RHO(K,I)-RHO(KWD,I)))/(HT*RHO(KWD,I)+NONZERO)*G),NONZERO)

!** Thickness

    HWDT = (COEF*QWD(JWD)/RHOFT)**0.333333
    IF (HT >= HWDT) THEN
      KTOP = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELWD+HWDT) < EL(KT,I)) THEN
    DLRHOT = ABS(RHO(KWD,I)-RHO(KTOP,I))
  ELSE IF (EL(KT,I) == ELWD) THEN
    DLRHOT = NONZERO
  ELSE
    DLRHOT = ABS(RHO(KWD,I)-RHO(KT,I))*HWDT/(EL(KT,I)-ELWD)
  END IF
  DLRHOT = MAX(DLRHOT,NONZERO)

! Withdrawal zone below structure

  DO K=KWD+1,KBOT

!** Density frequency

    HB    = ELWD-EL(K,I)
    RHOFB = MAX(SQRT((ABS(RHO(K,I)-RHO(KWD,I)))/(HB*RHO(KWD,I)+NONZERO)*G),NONZERO)

!** Thickness

    HWDB = (COEF*QWD(JWD)/RHOFB)**0.333333
    IF (HB >= HWDB) THEN
      KBOT = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELWD-HWDB) > EL(KBOT+1,I)) THEN
    DLRHOB = ABS(RHO(KWD,I)-RHO(KBOT,I))
  ELSE
    DLRHOB = ABS(RHO(KWD,I)-RHO(KBOT,I))*HWDB/(ELWD-EL(KBOT+1,I))
  END IF
  DLRHOB = MAX(DLRHOB,NONZERO)

! Velocity profile

  VSUM     = 0.0
  DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)                                                                             ! SW 1/24/05
  DO K=KTOP,KBOT
    VNORM(K) = ABS(1.0-((RHO(K,I)-RHO(KWD,I))/DLRHOMAX)**2)*BHR2(K,I)
    VSUM     = VSUM+VNORM(K)
  END DO

! Outflows

  DO K=KTOP,KBOT
    QSW(K,JWD) = QSW(K,JWD)+(VNORM(K)/VSUM)*QWD(JWD)
  END DO
  KTW(JWD) = KTOP
  KBW(JWD) = KBOT
END SUBROUTINE WITHDRAWAL

!***********************************************************************************************************************************
!**                                           S U B R O U T I N E   T R A N S P O R T                                             **
!***********************************************************************************************************************************

SUBROUTINE TRANSPORT
  USE GLOBAL; USE GEOMC; USE TVDC; USE TRANS; USE LOGICC; USE STRUCTURES; USE PREC

! Type declarations

  REAL,     SAVE, ALLOCATABLE, DIMENSION(:)     :: RATD,   CURX1,  CURX2,  CURX3
!  REAL,     SAVE, ALLOCATABLE, DIMENSION(:,:)   :: RATZ,   CURZ1,  CURZ2,  CURZ3
  REAL,     SAVE, ALLOCATABLE, DIMENSION(:,:)   :: SF1X,   SF1Z
  REAL,     SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: SF12X,  SF13X
  REAL,     SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: SF2X,   SF3X,   SF4X,   SF5X,   SF6X,   SF7X,   SF8X,   SF9X,   SF10X,  SF11X
  REAL,     SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: SF2Z,   SF3Z,   SF4Z,   SF5Z,   SF6Z,   SF7Z,   SF8Z,   SF9Z,   SF10Z
  REAL(R8), SAVE, ALLOCATABLE, DIMENSION(:,:)   :: DX1,    DX2,    DX3
  REAL(R8), SAVE, ALLOCATABLE, DIMENSION(:,:)   :: AD1X,   AD2X,   AD3X,   AD1Z,   AD2Z,   AD3Z

! Allocation declarations

  ALLOCATE (RATD(IMX),       CURX1(IMX),      CURX2(IMX),       CURX3(IMX))
  ALLOCATE (SF1X(KMX,IMX),   SF1Z(KMX,NWB))
!  ALLOCATE (RATZ(KMX,NWB),   CURZ1(KMX,NWB),  CURZ2(KMX,NWB),   CURZ3(KMX,NWB))
  ALLOCATE (DX1(KMX,IMX),    DX2(KMX,IMX),    DX3(KMX,IMX))
  ALLOCATE (AD1X(KMX,IMX),   AD2X(KMX,IMX),   AD3X(KMX,IMX))
  ALLOCATE (AD1Z(KMX,IMX),   AD2Z(KMX,IMX),   AD3Z(KMX,IMX))
  ALLOCATE (SF2X(KMX,IMX,2), SF3X(KMX,IMX,2), SF4X(KMX,IMX,2),  SF5X(KMX,IMX,2),  SF6X(KMX,IMX,2),  SF7X(KMX,IMX,2))
  ALLOCATE (SF8X(KMX,IMX,2), SF9X(KMX,IMX,2), SF10X(KMX,IMX,2), SF11X(KMX,IMX,2), SF12X(KMX,IMX,2), SF13X(KMX,IMX,2))
  ALLOCATE (SF2Z(KMX,2,NWB), SF3Z(KMX,2,NWB), SF4Z(KMX,2,NWB),  SF5Z(KMX,2,NWB),  SF6Z(KMX,2,NWB),  SF7Z(KMX,2,NWB))
  ALLOCATE (SF8Z(KMX,2,NWB), SF9Z(KMX,2,NWB), SF10Z(KMX,2,NWB))

! Variable initialization

  CT   = 0.0; AT   = 0.0; VT   = 0.0; DT   = 0.0; DX1  = 0.0; DX2  = 0.0; DX3  = 0.0; ADZ  = 0.0; ADX  = 0.0;  AD1X = 0.0
  AD2X = 0.0; AD3X = 0.0; AD1Z = 0.0; AD2Z = 0.0; AD3Z = 0.0
RETURN

!***********************************************************************************************************************************
!**                                        I N T E R P O L A T I O N  M U L T I P L I E R S                                       **
!***********************************************************************************************************************************

ENTRY INTERPOLATION_MULTIPLIERS

! Positive horizontal flows

  DO I=2,IMX-1
    DO K=2,KMX-1
      DLXT = DLX(I-1)
      IF (K > KB(I-1) .OR. INTERNAL_WEIR(K,I)) DLXT = DLX(I)
      DLXMIN       =  MIN(DLX(I+1),DLX(I))
      SF1X(K,I)    = (DLX(I+1)+DLX(I))*0.5
      SF2X(K,I,1)  =  DLX(I)/(DLX(I)+DLX(I+1))
      SF3X(K,I,1)  =  DLX(I)**2
      SF4X(K,I,1)  =  DLX(I+1)/(DLX(I)+DLX(I+1))
      SF5X(K,I,1)  =  0.25*(DLXT+2.0*DLX(I)+DLX(I+1))*(DLXT+DLX(I))
      SF6X(K,I,1)  = -0.25*(         DLX(I)+DLX(I+1))*(DLXT+DLX(I))
      SF7X(K,I,1)  =  0.25*(DLX(I)+DLX(I+1))*(DLXT+2.0*DLX(I)+DLX(I+1))
      SF8X(K,I,1)  =  0.50*(         DLX(I)-DLX(I+1))*DLXMIN
      SF9X(K,I,1)  =  0.50*(DLXT+2.0*DLX(I)-DLX(I+1))*DLXMIN
      SF10X(K,I,1) =  0.50*(DLXT+3.0*DLX(I))         *DLXMIN
      SF11X(K,I,1) =  SF8X(K,I,1) /SF5X(K,I,1)/SF1X(K,I)
      SF12X(K,I,1) =  SF9X(K,I,1) /SF6X(K,I,1)/SF1X(K,I)
      SF13X(K,I,1) =  SF10X(K,I,1)/SF7X(K,I,1)/SF1X(K,I)
    END DO
  END DO

! Negative horizontal flows

  DO I=2,IMX-2
    DO K=2,KMX-1
      DLXT = DLX(I+2)
      IF (K > KB(I+2)) DLXT = DLX(I+1)
      DLXMIN       =  MIN(DLX(I),DLX(I+1))
      SF1X(K,I)    = (DLX(I+1)+DLX(I))*0.5
      SF2X(K,I,2)  =  DLX(I+1)/(DLX(I)+DLX(I+1))
      SF3X(K,I,2)  =  DLX(I+1)**2
      SF4X(K,I,2)  =  DLX(I)/(DLX(I)+DLX(I+1))
      SF5X(K,I,2)  =  0.25*(DLX(I)+2.0*DLX(I+1)+DLXT)*(DLX(I)+DLX(I+1))
      SF6X(K,I,2)  = -0.25*(           DLX(I+1)+DLXT)*(DLX(I)+DLX(I+1))
      SF7X(K,I,2)  =  0.25*(DLX(I)+2.0*DLX(I+1)+DLXT)*(DLX(I+1)+DLXT)
      SF8X(K,I,2)  = -0.50*(       3.0*DLX(I+1)+DLXT)*DLXMIN
      SF9X(K,I,2)  =  0.50*(DLX(I)-2.0*DLX(I+1)-DLXT)*DLXMIN
      SF10X(K,I,2) =  0.50*(DLX(I)-DLX(I+1))*DLXMIN
      SF11X(K,I,2) =  SF8X(K,I,2) /SF5X(K,I,2)/SF1X(K,I)
      SF12X(K,I,2) =  SF9X(K,I,2) /SF6X(K,I,2)/SF1X(K,I)
      SF13X(K,I,2) =  SF10X(K,I,2)/SF7X(K,I,2)/SF1X(K,I)
    END DO
  END DO

! Ultimate multipliers

  IF (ULTIMATE(JW)) THEN
    DO JB=BS(JW),BE(JW)
      DO I=US(JB),DS(JB)
        RATD(I)  =  DLXR(I-1)/DLXR(I)
        CURX1(I) =  2.0*DLX(I)**2/(DLXR(I)+DLXR(I-1))/DLXR(I-1)
        CURX2(I) = -2.0*DLX(I)**2/(DLXR(I)*DLXR(I-1))
        CURX3(I) =  2.0*DLX(I)**2/(DLXR(I)+DLXR(I-1))/DLXR(I)
      END DO
    END DO
  END IF

! Vertical positive flows

  DO K=2,KMX-1
    HT            =  H(K-1,JW)
    HM            =  H(K,JW)
    HB            =  H(K+1,JW)
    HMIN          =  MIN(HB,HM)
    SF1Z(K,JW)    = (HB+HM)*0.5
    SF2Z(K,1,JW)  =  HM**2
    SF3Z(K,1,JW)  =  HM/(HM+HB)
    SF4Z(K,1,JW)  =  HB/(HM+HB)
    SF5Z(K,1,JW)  =  0.25*(HT+2.0*HM+HB)*(HT+HM)
    SF6Z(K,1,JW)  = -0.25*(HM+HB)*(HT+HM)
    SF7Z(K,1,JW)  =  0.25*(HM+HB)*(HT+2.0*HM+HB)
    SF8Z(K,1,JW)  =  0.50*(HM-HB)*HMIN
    SF9Z(K,1,JW)  =  0.50*(HT+2.0*HM-HB)*HMIN
    SF10Z(K,1,JW) =  0.50*(HT+3.0*HM)*HMIN
  END DO

! Vertical negative flows

  DO K=2,KMX-2
    HT            =  H(K,JW)
    HM            =  H(K+1,JW)
    HB            =  H(K+2,JW)
    HMIN          =  MIN(HT,HM)
    SF1Z(K,JW)    = (HM+HT)*0.5
    SF2Z(K,2,JW)  =  HM**2
    SF3Z(K,2,JW)  =  HM/(HT+HM)
    SF4Z(K,2,JW)  =  HT/(HT+HM)
    SF5Z(K,2,JW)  =  0.25*(HT+2.0*HM+HB)*(HT+HM)
    SF6Z(K,2,JW)  = -0.25*(HM+HB)*(HT+HM)
    SF7Z(K,2,JW)  =  0.25*(HT+2.0*HM+HB)*(HM+HB)
    SF8Z(K,2,JW)  = -0.50*(3.0*HM+HB)*HMIN
    SF9Z(K,2,JW)  =  0.50*(HT-2.0*HM-HB)*HMIN
    SF10Z(K,2,JW) =  0.50*(HT-HM)*HMIN
  END DO

! Ultimate multipliers

!  IF (ULTIMATE(JW)) THEN
!    DO K=2,KMX
!      RATZ(K,JW)  =  AVH2(K-1,DS(BE(JW)))/AVH2(K,DS(BE(JW)))                                         ! SW 5/20/05
!      CURZ1(K,JW) =  2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))+AVH2(K,DS(BE(JW))))/AVH2(K-1,DS(BE(JW)))   ! SW 5/20/05
!      CURZ2(K,JW) = -2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))*AVH2(K,DS(BE(JW))))                        ! SW 5/20/05
!      CURZ3(K,JW) =  2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))+AVH2(K,DS(BE(JW))))/AVH2(K,DS(BE(JW)))     ! SW 5/20/05
!    END DO
!  END IF
RETURN

!***********************************************************************************************************************************
!**                                          H O R I Z O N T A L  M U L T I P L I E R S                                           **
!***********************************************************************************************************************************

ENTRY HORIZONTAL_MULTIPLIERS

! Horizontal advection and diffusion multipliers

  IF (UPWIND(JW)) THEN
    DO I=IU,ID-1
      DO K=KT,KB(I)
        IF (U(K,I) >= 0.0) THEN
          C2X      =  COLD(K,I)
          C3X      =  COLD(K,I+1)
          DX2(K,I) = -DX(K,I)/SF1X(K,I)
          DX3(K,I) =  DX(K,I)/SF1X(K,I)
          ADX(K,I) = (DX2(K,I)-U(K,I))*C2X+DX3(K,I)*C3X
        ELSE
          C1X      =  COLD(K,I)
          C2X      =  COLD(K,I+1)
          DX1(K,I) = -DX(K,I)/SF1X(K,I)
          DX2(K,I) =  DX(K,I)/SF1X(K,I)
          ADX(K,I) =  DX1(K,I)*C1X+(DX2(K,I)-U(K,I))*C2X
        END IF
      END DO
    END DO
  ELSE
    DO I=IU,ID-1
      DO K=KT,KB(I)
        COUR = U(K,I)*DLT/DLXR(I)
        IF (U(K,I) >= 0.0) THEN
          C1X = COLD(K,I-1)
          C2X = COLD(K,I)
          C3X = COLD(K,I+1)
          IF (U(K,I-1) <= 0.0 .OR. K > KB(I-1) .OR. INTERNAL_WEIR(K,I-1)) C1X = COLD(K,I)
          IF (INTERNAL_WEIR(K,I)) C3X = COLD(K,I)
          CART      =  C3X
          CALF      =  C1X
          RATS      =  RATD(I)
          CURS1     =  CURX1(I)
          CURS2     =  CURX2(I)
          CURS3     =  CURX3(I)
          DX1(K,I)  =  DX(K,I)*SF11X(K,I,1)
          DX2(K,I)  =  DX(K,I)*SF12X(K,I,1)
          DX3(K,I)  =  DX(K,I)*SF13X(K,I,1)
          ALFA      =  2.0*(DX(K,I)*DLT/(SF1X(K,I)*SF1X(K,I))-(1.0-COUR*COUR)/6.0)*SF3X(K,I,1)
          AD1X(K,I) = (ALFA-COUR*SF8X(K,I,1)*0.5)/SF5X(K,I,1)
          AD2X(K,I) =  SF4X(K,I,1)+(ALFA-COUR*SF9X(K,I,1) *0.5)/SF6X(K,I,1)
          AD3X(K,I) =  SF2X(K,I,1)+(ALFA-COUR*SF10X(K,I,1)*0.5)/SF7X(K,I,1)
        ELSE
          C1X = COLD(K,I)
          C2X = COLD(K,I+1)
          C3X = COLD(K,I+2)
          IF (U(K,I+2) >= 0.0 .OR. K > KB(I+2) .OR. I == ID-1 .OR. INTERNAL_WEIR(K,I+1)) C3X = COLD(K,I+1)
          IF (INTERNAL_WEIR(K,I)) THEN
            C2X = COLD(K,I)
            C3X = COLD(K,I)
          END IF     
          CART      =  C1X
          CALF      =  C3X
          RATS      =  RATD(I+1)
          CURS1     =  CURX1(I+1)
          CURS2     =  CURX2(I+1)
          CURS3     =  CURX3(I+1)
          DX1(K,I)  =  DX(K,I)*SF11X(K,I,2)
          DX2(K,I)  =  DX(K,I)*SF12X(K,I,2)
          DX3(K,I)  =  DX(K,I)*SF13X(K,I,2)
          ALFA      =  2.0*(DX(K,I)*DLT/(SF1X(K,I)*SF1X(K,I))-(1.0-COUR*COUR)/6.0)*SF3X(K,I,2)
          AD1X(K,I) =  SF2X(K,I,2)+(ALFA-COUR*SF8X(K,I,2)*0.5)/SF5X(K,I,2)
          AD2X(K,I) =  SF4X(K,I,2)+(ALFA-COUR*SF9X(K,I,2)*0.5)/SF6X(K,I,2)
          AD3X(K,I) = (ALFA-COUR*SF10X(K,I,2)*0.5)/SF7X(K,I,2)
        END IF
        ADX(K,I) = (DX1(K,I)-U(K,I)*AD1X(K,I))*C1X+(DX2(K,I)-U(K,I)*AD2X(K,I))*C2X+(DX3(K,I)-U(K,I)*AD3X(K,I))*C3X
        IF (ULTIMATE(JW)) THEN
          RATDI = 1.0/RATS
          DELC  = RATS*C3X+(RATDI-RATS)*C2X-RATDI*C1X
          DELC  = SIGN(1.0,U(K,I))*DELC
          ADELC = ABS(DELC)
          ACURZ = ABS(CURS3*C3X+CURS2*C2X+CURS1*C1X)
          IF (ACURZ <= 0.6*ADELC) THEN
            FLUX = AD1X(K,I)*C1X+AD2X(K,I)*C2X+AD3X(K,I)*C3X
          ELSE IF (ACURZ >= ADELC) THEN
            FLUX = C2X
          ELSE IF (ABS(COUR) > 0.0) THEN
            FTEMP = AD1X(K,I)*C1X+AD2X(K,I)*C2X+AD3X(K,I)*C3X
            CREF  = CALF+(C2X-CALF)/ABS(COUR)
            IF (DELC > 0.0) THEN
              CMAX1 = MIN(CREF,CART)
              IF (CREF < C2X) CMAX1 = CART
              FLUX = 0.5*(C2X+CMAX1)
              IF (FTEMP <= CMAX1 .AND. FTEMP >= C2X) FLUX = FTEMP
            ELSE
              CMIN1 = MAX(CREF,CART)
              IF (CREF > C2X) CMIN1 = CART
              IF (FTEMP >= CMIN1 .AND. FTEMP <= C2X) THEN
                FLUX = FTEMP
              ELSE IF (FTEMP > 0.0) THEN
                FLUX = 0.5*(C2X+CMIN1)
              ELSE
                FLUX = 0.0
              END IF
            END IF
          ELSE
            FLUX = 0.0
          END IF
          ADX(K,I) = (DX1(K,I)*C1X+DX2(K,I)*C2X+DX3(K,I)*C3X)-U(K,I)*FLUX
        END IF
      END DO
    END DO
  END IF
RETURN

!***********************************************************************************************************************************
!**                                            V E R T I C A L  M U L T I P L I E R S                                             **
!***********************************************************************************************************************************

ENTRY VERTICAL_MULTIPLIERS

! Vertical advection multipliers

  IF (UPWIND(JW)) THEN
    DO I=IU,ID
      DO K=KT,KB(I)-1
        C2Z = COLD(K+1,I)
        IF (W(K,I) >= 0.0) C2Z = COLD(K,I)
        ADZ(K,I) = -W(K,I)*C2Z
      END DO
    END DO
  ELSE
    DO I=IU,ID
      DO K=KT,KB(I)-1
        IF (W(K,I) >= 0.0) THEN
          C1Z   = COLD(K-1,I)
          C2Z   = COLD(K,I)
          C3Z   = COLD(K+1,I)
          CART  = C3Z
          CALF  = C1Z
          RATS  = RATZ(K,JW)
          CURS1 = CURZ1(K,JW)
          CURS2 = CURZ2(K,JW)
          CURS3 = CURZ3(K,JW)
          IF (K <= KT+1) THEN
            C1Z   =  COLD(KT,I)
            HT    =  H1(KT,I)
            HM    =  H1(K,I)
            HB    =  H1(K+1,I)
            CALF  =  C1Z
            RATS  =  AVH1(KT,I)/AVH1(K,I)
            CURS1 =  2.0*HM*HM/(AVH1(KT,I)+AVH1(K,I))/AVH1(KT,I)
            CURS2 = -2.0*HM*HM/(AVH1(KT,I)*AVH1(K,I))
            CURS3 =  2.0*HM*HM/(AVH1(KT,I)+AVH1(K,I))/AVH1(K,I)
            IF (K == KT) THEN
              HM    =  H1(KT,I)
              RATS  =  1.0
              CURS3 =  1.0
              CURS2 = -2.0
              CURS1 =  1.0
            END IF
            HMIN          =  MIN(HB,HM)
            SF1Z(K,JW)    = (HB+HM)*0.5
            SF2Z(K,1,JW)  =  HM**2
            SF3Z(K,1,JW)  =  HM/(HM+HB)
            SF4Z(K,1,JW)  =  HB/(HM+HB)
            SF5Z(K,1,JW)  =  0.25*(HT+2.0*HM+HB)*(HT+HM)
            SF6Z(K,1,JW)  = -0.25*(HM+HB)*(HT+HM)
            SF7Z(K,1,JW)  =  0.25*(HM+HB)*(HT+2.0*HM+HB)
            SF8Z(K,1,JW)  =  0.5*(HM-HB)*HMIN
            SF9Z(K,1,JW)  =  0.5*(HT+2.0*HM-HB)*HMIN
            SF10Z(K,1,JW) =  0.5*(HT+3.0*HM)*HMIN
          END IF
          COUR      =  W(K,I)*DLT/SF1Z(K,JW)
          ALFA      =  2.0*(DZQ(K,I)*DLT/(SF1Z(K,JW)*SF1Z(K,JW))-(1.0-COUR*COUR)/6.0)*SF2Z(K,1,JW)
          AD1Z(K,I) = (ALFA-COUR*SF8Z(K,1,JW)*0.5)/SF5Z(K,1,JW)
          AD2Z(K,I) =  SF4Z(K,1,JW)+(ALFA-COUR*SF9Z(K,1,JW) *0.5)/SF6Z(K,1,JW)
          AD3Z(K,I) =  SF3Z(K,1,JW)+(ALFA-COUR*SF10Z(K,1,JW)*0.5)/SF7Z(K,1,JW)
        ELSE
          C1Z = COLD(K,I)
          C2Z = COLD(K+1,I)
          C3Z = COLD(K+2,I)
          IF (K == KB(I)-1) C3Z = COLD(K+1,I)
          CART  = C1Z
          CALF  = C3Z
          CURS3 = CURZ3(K+1,JW)
          CURS2 = CURZ2(K+1,JW)
          CURS1 = CURZ1(K+1,JW)
          RATS  = AVH1(K,I)/AVH1(K+1,I)
          IF (K == KT) THEN
            HT            =  H1(KT,I)
            HM            =  H1(KT+1,I)
            HB            =  H1(KT+2,I)
            HMIN          =  MIN(HT,HM)
            RATS          =  AVH1(KT,I)/AVH1(K,I)
            CURS1         =  2.0*HM*HM/(AVH1(KT,I)+AVH1(K,I))/AVH1(KT,I)
            CURS2         = -2.0*HM*HM/(AVH1(KT,I)*AVH1(K,I))
            CURS3         =  2.0*HM*HM/(AVH1(KT,I)+AVH1(K,I))/AVH1(K,I)
            SF1Z(K,JW)    = (HM+HT)*0.5
            SF2Z(K,2,JW)  =  HM**2
            SF3Z(K,2,JW)  =  HM/(HT+HM)
            SF4Z(K,2,JW)  =  HT/(HT+HM)
            SF5Z(K,2,JW)  =  0.25*(HT+2.0*HM+HB)*(HT+HM)
            SF6Z(K,2,JW)  = -0.25*(HM+HB)*(HT+HM)
            SF7Z(K,2,JW)  =  0.25*(HT+2.0*HM+HB)*(HM+HB)
            SF8Z(K,2,JW)  = -0.5*(3.0*HM+HB)*HMIN
            SF9Z(K,2,JW)  =  0.5*(HT-2.0*HM-HB)*HMIN
            SF10Z(K,2,JW) =  0.5*(HT-HM)*HMIN
          END IF
          COUR      =  W(K,I)*DLT/SF1Z(K,JW)
          ALFA      =  2.0*(DZQ(K,I)*DLT/(SF1Z(K,JW)*SF1Z(K,JW))-(1.0-COUR*COUR)/6.0)*SF2Z(K,2,JW)
          AD1Z(K,I) =  SF3Z(K,2,JW)+(ALFA-COUR*SF8Z(K,2,JW)*0.5)/SF5Z(K,2,JW)
          AD2Z(K,I) =  SF4Z(K,2,JW)+(ALFA-COUR*SF9Z(K,2,JW)*0.5)/SF6Z(K,2,JW)
          AD3Z(K,I) = (ALFA-COUR*SF10Z(K,2,JW)*0.5)/SF7Z(K,2,JW)
        END IF
        ADZ(K,I) = -W(K,I)*(AD1Z(K,I)*C1Z+AD2Z(K,I)*C2Z+AD3Z(K,I)*C3Z)
        IF (ULTIMATE(JW)) THEN
          RATZI = 1.0/RATS
          DELC  = RATS*C3Z+(RATZI-RATS)*C2Z-RATZI*C1Z
          DELC  = SIGN(1.0,W(K,I))*DELC
          ADELC = ABS(DELC)
          ACURZ = ABS(CURS3*C3Z+CURS2*C2Z+CURS1*C1Z)
          IF (ACURZ <= 0.6*ADELC) THEN
            FLUX = AD1Z(K,I)*C1Z+AD2Z(K,I)*C2Z+AD3Z(K,I)*C3Z
          ELSE IF (ACURZ >= ADELC) THEN
            FLUX = C2Z
          ELSE IF (ABS(COUR) > 0.0) THEN
            FTEMP = AD1Z(K,I)*C1Z+AD2Z(K,I)*C2Z+AD3Z(K,I)*C3Z
            CREF  = CALF+(C2Z-CALF)/ABS(COUR)
            IF (DELC > 0.0) THEN
              CMAX1 = CART
              IF (CREF >= C2Z) CMAX1 = MIN(CREF,CART)
              FLUX = 0.5*(C2Z+CMAX1)
              IF (FTEMP <= CMAX1 .AND. FTEMP >= C2Z) FLUX = FTEMP
            ELSE
              CMIN1 = MAX(CREF,CART)
              IF (CREF > C2Z) CMIN1 = CART
              IF (FTEMP >= CMIN1 .AND. FTEMP <= C2Z) THEN
                FLUX = FTEMP
              ELSE IF (FTEMP > 0.0) THEN
                FLUX = 0.5*(C2Z+CMIN1)
              ELSE
                FLUX = 0.0
              END IF
            END IF
          ELSE
            FLUX = 0.0
          END IF
          ADZ(K,I) = -W(K,I)*FLUX
        END IF
      END DO
    END DO
  END IF
RETURN

!***********************************************************************************************************************************
!**                                           H O R I Z O N T A L  T R A N S P O R T                                              **
!***********************************************************************************************************************************

ENTRY HORIZONTAL_TRANSPORT
  IF (CONSTITUENTS) THEN
    DO I=IU,ID
      DO K=KT,KB(I)
        CNEW(K,I) = (COLD(K,I)*BH2(K,I)/DLT+(ADX(K,I)*BHR1(K,I)-ADX(K,I-1)*BHR1(K,I-1))/DLX(I)+(1.0-THETA(JW))                     &
                    *(ADZ(K,I)*BB(K,I)-ADZ(K-1,I)*BB(K-1,I))+SSB(K,I)/DLX(I))*DLT/BH1(K,I)+SSK(K,I)*DLT
      END DO
    END DO
  ELSE
    DO I=IU,ID
      DO K=KT,KB(I)
        CNEW(K,I) = (COLD(K,I)*BH2(K,I)/DLT+(ADX(K,I)*BHR1(K,I)-ADX(K,I-1)*BHR1(K,I-1))/DLX(I)+(1.0-THETA(JW))                     &
                    *(ADZ(K,I)*BB(K,I)-ADZ(K-1,I)*BB(K-1,I))+SSB(K,I)/DLX(I))*DLT/BH1(K,I)
      END DO
    END DO
  END IF
RETURN
ENTRY DEALLOCATE_TRANSPORT
  DEALLOCATE (RATD, CURX1, CURX2, CURX3, SF1X,  SF1Z, RATZ, CURZ1, CURZ2, CURZ3, CT,   AT,   VT,   DT)     
  DEALLOCATE (DX1,  DX2,   DX3,   AD1X,  AD2X,  AD3X,  AD1Z, AD2Z, AD3Z,  SF2X,  SF3X, SF4X, SF5X, SF6X, SF7X)
  DEALLOCATE (SF8X, SF9X,  SF10X, SF11X, SF12X, SF13X, SF2Z, SF3Z, SF4Z,  SF5Z,  SF6Z, SF7Z, SF8Z, SF9Z, SF10Z)
RETURN
END SUBROUTINE TRANSPORT

!***********************************************************************************************************************************
!*                                              S U B R O U T I N E    T R I D I A G                                              **
!***********************************************************************************************************************************

SUBROUTINE TRIDIAG(A,V,C,D,S,E,N,U)
  USE PREC
  INTEGER,                             INTENT(IN)  :: S, E, N
  REAL(R8),              DIMENSION(:), INTENT(IN)  :: A(E),V(E),C(E),D(E)
  REAL,                  DIMENSION(:), INTENT(OUT) :: U(N)
  REAL(R8), ALLOCATABLE, DIMENSION(:)              :: BTA, GMA
  ALLOCATE (BTA(N),GMA(N))

  BTA(S) = V(S)
  GMA(S) = D(S)
  DO I=S+1,E
    BTA(I) = V(I)-A(I)/BTA(I-1)*C(I-1)
    GMA(I) = D(I)-A(I)/BTA(I-1)*GMA(I-1)
  END DO
  U(E) = GMA(E)/BTA(E)
  DO I=E-1,S,-1
    U(I) = (GMA(I)-C(I)*U(I+1))/BTA(I)
  END DO
  Deallocate (BTA, GMA)                                                                                             ! SW 10/17/05
END SUBROUTINE TRIDIAG

!***********************************************************************************************************************************
!**                                        S U B R O U T I N E   C A L C U L A T E  A Z                                           **
!***********************************************************************************************************************************

RECURSIVE SUBROUTINE CALCULATE_AZ
  USE GEOMC; USE GLOBAL; USE TRANS; USE PREC; USE EDDY; USE KINETIC; use macrophytec; USE LOGICC; SAVE     !v3.5

  REAL, DIMENSION (2) :: SIG
  DATA SIG  /1.0, 1.3/, TKEMIN1 /1.25E-7/, TKEMIN2 /1.0E-9/

  DO K=KT,KBMIN(I)-1
    VSH(K,I) = ((U(K+1,I)-U(K,I))/((AVH2(K,I)+AVH2(K,I+1))*0.5))**2
  END DO
  IF (AZC(JW) == '     TKE') THEN 
    CALL CALCULATE_TKE
  ELSE
    DO K=KT,KBMIN(I)-1
      CALL CALCULATE_AZ0
      BUOY    = (RHO(K+1,I)-RHO(K,I)+RHO(K+1,I+1)-RHO(K,I+1))/(2.0*(AVH2(K,I)+AVH2(K,I+1))*0.5)            !AVH2(K,I)       SW 10/20/07
      RIAZ0   =  LOG(AZ0/AZMAX(JW))/1.5
      RI      =  G*BUOY/(RHOW*VSH(K,I)+NONZERO)
      RIAZ1   =  MAX(RI,RIAZ0)
      RIAZ1   =  MIN(RIAZ1,10.0)
      EXPAZ   =  EXP(-1.5*RIAZ1)
      AZ(K,I) =  MAX(AZMIN,AZ0*EXPAZ+AZMIN*(1.0-EXPAZ))    ! AZ computed at lower edge of cell
      DZT(K,I) =  MAX(DZMIN,FRAZDZ*(AZ0*EXPAZ+DZMIN*(1.0-EXPAZ)))   ! DZ computed at lower edge of cell - later averaged to cell center lower edge
    END DO
  END IF
RETURN

ENTRY CALCULATE_AZ0
  IF (K == KT) THEN
    DEPTHL = (EL(KT,I)  -Z(I)  *COSA(JB)-EL(KB(I),I)    +H2(KB(I),I)  *COSA(JB))/COSA(JB)
    DEPTHR = (EL(KT,I+1)-Z(I+1)*COSA(JB)-EL(KB(I+1),I+1)+H2(KB(I+1),I)*COSA(JB))/COSA(JB)
    ZDLR   =  DEPTHL-H2(KT,I)+DEPTHR-H2(KT,I+1)
    IF (AZC(JW) == '     RNG' .OR. AZC(JW) == '   PARAB') THEN
      DEPTH = (DEPTHR+DEPTHL)*0.5
      USTAR =  SQRT(G*DEPTH*SLOPE(JB))
      IF (SLOPE(JB) == 0.0) THEN
        USTAR = 0.0
        DO KL=KT,KBMIN(I)
          USTAR = USTAR+SQRT(AZ(KL-1,I)*SQRT(VSH(KL-1,I))/RHO(KL,I))
        END DO
        USTAR = USTAR/(KBMIN(I)-KT+1)
      END IF
    END IF
  ELSE
    ZDLR = (EL(K-1,I)-EL(KB(I),I)+H2(KB(I),I)*COSA(JB))/COSA(JB)+(EL(K-1,I+1)-EL(KB(I),I+1)+H2(KB(I+1),I)*COSA(JB))/COSA(JB)
  END IF
  ZD = ZDLR/(DEPTHL+DEPTHR)
  IF (AZC(JW) == '    NICK') THEN
    SLM = (DEPTH*(0.14-0.08*(1.0-ZD)**2-0.06*(1.0-ZD)**4))**2
    AZ0 =  MAX(AZMIN,SLM*SQRT(VSH(K,I)))
  ELSE IF (AZC(JW) == '     RNG') THEN
    VISCK = EXP((T2(K,I)+495.691)/(-37.3877))
    IF (T2(K,I) > 30.0)  VISCK = EXP((T2(K,I)+782.190)/(-57.7600))
    VISCF = MAX(0.0,0.08477*((ZDLR*0.5*USTAR/VISCK)**3)*((1.0-ZDLR*0.5/DEPTH)**3)-100.0)
    VISCF = (1.0+VISCF)**0.33333333333333
    AZ0   = MAX(AZMIN,VISCK*VISCF)
  ELSE IF (AZC(JW) == '   PARAB') THEN
    AZ0 = MAX(AZMIN,0.41*USTAR*ZDLR*0.5*(1.0-ZD))
  ELSE
    SLM = HMAX2
    IF (AZC(JW) == '     W2N') SLM = ((DEPTHR+DEPTHL)*0.5*(0.14-0.08*(1.0-ZD)**2-0.06*(1.0-ZD)**4))**2
    AZ0 = 0.4*SLM*SQRT(VSH(K,I)+((FRICBR(K,I)+WSHY(I)*DECAY(K,I))/(AZ(K,I)+NONZERO))**2)+AZMIN
  END IF
RETURN

ENTRY CALCULATE_TKE
  USTAR = SQRT(1.25*CZ(I)*WIND10(I)**2/RHO(KT,I))
  IF (MANNINGS_N(JW)) THEN
   HRAD = BH1(KT,I)/(B(KTI(I),I)-B(KT+1,I)+2.*AVH1(KT,I))      ! HRAD = BHR1(KT,I)/(BR(KTI(I),I)-BR(KT+1,I)+2.*AVH1(KT,I))  SW 10/5/07  These calculations are at the segment centers and vertical center of a layer
    if(macrophyte_on.and.mannings_n(jw))then
      call macrophyte_friction(hrad,fric(i),effric,kt,i)
      gc2=g*effric*effric/hrad**0.33333333
    else if(.not.macrophyte_on.and.mannings_n(jw))then
      gc2=g*fric(i)*fric(i)/hrad**0.33333333
    end if
  ELSE
    GC2 = 0.0
    IF (FRIC(I) /= 0.0) GC2 = G/(FRIC(I)*FRIC(I))
  END IF
  USTARBKT     = SQRT(GC2)*ABS(0.5*(U(KT,I)+U(KT,I-1)))                                                 ! SG 10/4/07 
  TKE(KT,I,1) = (3.33*(USTAR*USTAR+USTARBKT*USTARBKT))*(BH2(KT,I)/BH1(KT,I))                            ! SG 10/4/07
  TKE(KT,I,2) = (USTAR*USTAR*USTAR+USTARBKT*USTARBKT*USTARBKT)*5.0/H1(KT,I)*(BH2(KT,I)/BH1(KT,I))       ! SG 10/4/07
  DO K=KT+1,KB(I)-1                                                                                     
    BOUK = MAX(AZ(K,I)*G*(RHO(K+1,I)-RHO(K,I))/(H(K,JW)*RHOW),0.0)
    PRDK = AZ(K,I)*(0.5*(U(K,I)+U(K,I-1)-U(K+1,I)-U(K+1,I-1))/(H(K,JW)*0.5 + H(K+1,JW)*0.5))**2.0       ! SG 10/4/07
    PRHE = 10.0*GC2**1.25*ABS(0.5*(U(K,I)+U(K,I-1)))**4.0/(0.5*B(K,I))**2.0
    IF (MANNINGS_N(JW)) THEN
      HRAD = BH(K,I)/(B(K,I)-B(K+1,I)+2.0*H(K,JW))          ! HRAD = BHR(K,I)/(BR(K,I)-BR(K+1,I)+2.0*H(K,JW))  SW 10/5/07
      if(macrophyte_on.and.mannings_n(jw))then
        call macrophyte_friction(hrad,fric(i),effric,k,i)
        gc2=g*effric*effric/hrad**0.33333333
      else if (.not.macrophyte_on.and.mannings_n(jw))then
        gc2=g*fric(i)*fric(i)/hrad**0.33333333
      end if
    END IF
    PRHK       = GC2/(0.5*B(K,I))*ABS(0.5*(U(K,I)+U(K,I-1)))**3.0
    UNST       = PRDK-TKE(K,I,2)
    UNSE       = 1.44*TKE(K,I,2)/TKE(K,I,1)*PRDK-1.92*(TKE(K,I,2)/TKE(K,I,1)*TKE(K,I,2))
    TKE(K,I,1) = TKE(K,I,1)+DLT*(UNST+PRHK-BOUK)
    TKE(K,I,2) = TKE(K,I,2)+DLT*(UNSE+PRHE)
  END DO
  USTARB         = SQRT(GC2)*ABS(0.5*(U(KB(I),I)+U(KB(I),I-1)))                                 
  TKE(KB(I),I,1) = 0.5*(3.33*USTARB*USTARB+TKE(KB(I),I,1))                                      
  TKE(KB(I),I,2) = 0.5*(USTARB*USTARB*USTARB*5.0/H(KB(I),JW)+TKE(KB(I),I,2))                   
    
  DO J=1,2                                                    ! SG 10/4/07 Series of bug fixes for TKE
	  K=KT			
	  AT(K,I) = 0.0
	  CT(K,I) = 0.0
	  VT(K,I) = 1.0
	  DT(K,I) = TKE(K,I,J)	
	  DO K = KT+1,KB(I)-1                                                                         
		   AT(K,I) = -DLT/BH1(K,I)*BB(K-1,I)/SIG(J)*AZ(K-1,I)/AVH1(K-1,I)
		   CT(K,I) = -DLT/BH1(K,I)*BB(K,I)/SIG(J)*AZ(K,I)/AVH1(K,I)
		   VT(K,I) =  1.0-AT(K,I)-CT(K,I)
		   DT(K,I) =  TKE(K,I,J)
	  END DO
	  K = KB(I)                                                                                   
	  AT(K,I) = 0.0
	  CT(K,I) = 0.0
	  VT(K,I) = 1.0
	  DT(K,I) = TKE(K,I,J)
	  CALL TRIDIAG(AT(:,I),VT(:,I),CT(:,I),DT(:,I),KT,KB(I),KMX,TKE(:,I,J))
	END DO

  DO K=KT,KB(I)                                                                                  
    TKE(K,I,1) = MAX(TKE(K,I,1),TKEMIN1)
    TKE(K,I,2) = MAX(TKE(K,I,2),TKEMIN2)
    AZT(K,I)     = 0.09*TKE(K,I,1)*TKE(K,I,1)/TKE(K,I,2) 
  END DO

  DO K=KT,KB(I)-1                                                                             
	AZ(K,I) = 0.5*(AZT(K,I)+AZT(K+1,I))                                          ! SW 10/18/07
	AZ(K,I) = MAX(AZMIN,AZ(K,I))
	AZ(K,I) = MIN(AZMAX(JW),AZ(K,I))
	DZ(K,I) = MAX(DZMIN,FRAZDZ*AZ(K,I))     ! No further averaging required for DZ, AZ is at lower edge of cell at center, in main program averged to edge
  END DO
  AZ(KB(I),I) = AZMIN                                                                            

END SUBROUTINE CALCULATE_AZ

!***********************************************************************************************************************************
!**                                            S U B R O U T I N E   K I N E T I C S                                              **
!***********************************************************************************************************************************

SUBROUTINE KINETICS
  USE SCREENC; USE GLOBAL; USE KINETIC; USE GEOMC; USE TVDC; USE LOGICC; USE SURFHE
  use macrophytec; use zooplanktonc  !v3.5

! Type declarations

  REAL                                :: LAM1,   LAM2,   NH4PR,  NO3PR,  LIMIT,  LIGHT,  L, L0, L1
  REAL                                :: KW,     INCR,   OH,     K1,     K2
  real                                :: ltcoefm, lavg,  macext  !v3.5
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: OMTRM,  SODTRM, NH4TRM, NO3TRM
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: DOM,    POM,    PO4BOD, NH4BOD, TICBOD
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: lam2m  ! v3.5
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ATRM,   ATRMR,  ATRMF
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ETRM,   ETRMR,  ETRMF
  SAVE   

! Allocation declarations

  ALLOCATE (OMTRM(KMX,IMX),    SODTRM(KMX,IMX),    NH4TRM(KMX,IMX),    NO3TRM(KMX,IMX), DOM(KMX,IMX), POM(KMX,IMX))
  ALLOCATE (PO4BOD(KMX,IMX),   NH4BOD(KMX,IMX),    TICBOD(KMX,IMX))
  ALLOCATE (ATRM(KMX,IMX,NAL), ATRMR(KMX,IMX,NAL), ATRMF(KMX,IMX,NAL))
  ALLOCATE (ETRM(KMX,IMX,NEP), ETRMR(KMX,IMX,NEP), ETRMF(KMX,IMX,NEP))
  ALLOCATE (lam2m(KMX,kmx))  !v3.5
RETURN

!***********************************************************************************************************************************
!**                                      T E M P E R A T U R E  R A T E  M U L T I P L I E R S                                    **
!***********************************************************************************************************************************

ENTRY TEMPERATURE_RATES
  DO I=IU,ID
    DO K=KT,KB(I)
      LAM1        = FR(T1(K,I),NH4T1(JW),NH4T2(JW),NH4K1(JW),NH4K2(JW))
      NH4TRM(K,I) = LAM1/(1.0+LAM1-NH4K1(JW))
      LAM1        = FR(T1(K,I),NO3T1(JW),NO3T2(JW),NO3K1(JW),NO3K2(JW))
      NO3TRM(K,I) = LAM1/(1.0+LAM1-NO3K1(JW))
      LAM1        = FR(T1(K,I),OMT1(JW),OMT2(JW),OMK1(JW),OMK2(JW))
      OMTRM(K,I)  = LAM1/(1.0+LAM1-OMK1(JW))
      LAM1        = FR(T1(K,I),SODT1(JW),SODT2(JW),SODK1(JW),SODK2(JW))
      SODTRM(K,I) = LAM1/(1.0+LAM1-SODK1(JW))
      DO JA=1,NAL
        LAM1          = FR(T1(K,I),AT1(JA),AT2(JA),AK1(JA),AK2(JA))
        LAM2          = FF(T1(K,I),AT3(JA),AT4(JA),AK3(JA),AK4(JA))
        ATRMR(K,I,JA) = LAM1/(1.0+LAM1-AK1(JA))
        ATRMF(K,I,JA) = LAM2/(1.0+LAM2-AK4(JA))
        ATRM(K,I,JA)  = ATRMR(K,I,JA)*ATRMF(K,I,JA)
      END DO
      DO JE=1,NEP
        LAM1          = FR(T1(K,I),ET1(JE),ET2(JE),EK1(JE),EK2(JE))
        LAM2          = FF(T1(K,I),ET3(JE),ET4(JE),EK3(JE),EK4(JE))
        ETRMR(K,I,JE) = LAM1/(1.0+LAM1-EK1(JE))
        ETRMF(K,I,JE) = LAM2/(1.0+LAM2-EK4(JE))
        ETRM(K,I,JE)  = ETRMR(K,I,JE)*ETRMF(K,I,JE)
      END DO
! v3.5 start
      do m=1,nmc
        LAM1    = FR(T1(K,I),mT1(m),mT2(m),mK1(m),mk2(m))
        LAM2    = FF(T1(K,I),mT3(m),mT4(m),mK3(m),mk4(m))
        MACTRMR(K,I,m) = LAM1/(1.0+LAM1-mK1(m))
        MACTRMF(K,I,m) = LAM2/(1.0+LAM2-mK4(m))
        MACTRM(K,I,m)  = macTRMR(K,I,m)*macTRMF(K,I,m)
      end do
      if(zooplankton_calc)then
	    DO JZ = 1, NZP
          LAM1       = FR(T1(K,I),zt1(jz),zt2(jz),zk1(jz),zk2(jz))
          LAM2       = FF(T1(K,I),zt3(jz),zt4(jz),zk3(jz),zk4(jz))
          zoormr(k,i,jz)= lam1/(1.+lam1-zk1(jz))
          zoormf(k,i,jz)= lam2/(1.+lam2-zk4(jz))
          zoorm(k,i,jz) = zoormr(k,i,jz)*zoormf(k,i,jz)
        END DO
	  end if
! v3.5 end
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                 K I N E T I C   R A T E S                                                     **
!***********************************************************************************************************************************

ENTRY KINETIC_RATES

! Decay rates

  DO I=IU,ID
    DO K=KT,KB(I)
      do1(k,i)          = o2(k,i)/(o2(k,i)+kdo)                  !v3.5      
      do2(k,i)          = 1.0 - o2(k,i)/(o2(k,i)+kdo)           !v3.5
      DO3(K,I)          = (1.0+SIGN(1.0,O2(K,I)-1.E-10)) *0.5
! v3.5 start      
      SEDD(K,I)         =   SODTRM(K,I) *SDKv(k,i)   *SED(K,I) *DO3(K,I)   !cb 10/22/06      
      SEDDp(K,I)         =  SODTRM(K,I) *SDKv(k,i)   *SEDp(K,I) *DO3(K,I)
      SEDDn(K,I)         =  SODTRM(K,I) *SDKv(k,i)   *SEDn(K,I) *DO3(K,I)
      SEDDc(K,I)         =  SODTRM(K,I) *SDKv(k,i)   *SEDc(K,I) *DO3(K,I)      
      sedbr(k,i)         =  sedb(jw)    *sed(k,i)                           !cb 11/30/06
      sedbrp(k,i)        =  sedb(jw)    *sedp(k,i)                          !cb 11/30/06
      sedbrn(k,i)        =  sedb(jw)    *sedn(k,i)                          !cb 11/30/06
      sedbrc(k,i)        =  sedb(jw)    *sedc(k,i)                          !cb 11/30/06      
! v3.5 end
      NH4D(K,I)         =  NH4TRM(K,I) *NH4DK(JW) *NH4(K,I) *DO1(K,I)
      NO3D(K,I)         =  NO3TRM(K,I) *NO3DK(JW) *NO3(K,I) *DO2(K,I)
      LDOMD(K,I)        =  OMTRM(K,I)  *LDOMDK(JW)*LDOM(K,I)*DO3(K,I)
      RDOMD(K,I)        =  OMTRM(K,I)  *RDOMDK(JW)*RDOM(K,I)*DO3(K,I)
      LPOMD(K,I)        =  OMTRM(K,I)  *LPOMDK(JW)*LPOM(K,I)*DO3(K,I)
      RPOMD(K,I)        =  OMTRM(K,I)  *RPOMDK(JW)*RPOM(K,I)*DO3(K,I)
      LRDOMD(K,I)       =  OMTRM(K,I)  *LRDDK(JW) *LDOM(K,I)*DO3(K,I)
      LRPOMD(K,I)       =  OMTRM(K,I)  *LRPDK(JW) *LPOM(K,I)*DO3(K,I)
      CBODD(K,I,1:NBOD) =  KBOD(1:NBOD)*TBOD(1:NBOD)**(T1(K,I)-20.0)*DO3(K,I)
      if(k == kb(i))then     ! SW 4/18/07
	  SODD(K,I)         =  SOD(I)/BH2(K,I)*SODTRM(K,I)*BI(K,I)
	  else
      SODD(K,I)         =  SOD(I)/BH2(K,I)*SODTRM(K,I)*(BI(K,I)-BI(K+1,I))
	  endif
      ALGEX = 0.0; SSEXT = 0.0
      IF (.NOT. READ_EXTINCTION(JW)) THEN
        DO JA=1,NAL
          ALGEX = ALGEX+EXA(JA)*ALG(K,I,JA)
        END DO
        DO JS=1,NSS
          SSEXT = SSEXT+EXSS(JW)*SS(K,I,JS)
        END DO
! v3.5 start
        totmac=0.0                    
        do m=1,nmc
          if(macrophyte_calc(jw,m))then
            jt=kti(i)
            je=kb(i)
            do jj=jt,je
              totmac = exm(m)*macrm(jj,K,I,m)+totmac
            end do
          end if
        end do       
        macext=totmac/(bh2(k,i)*dlx(i))
        zooext = 0.0
	    do jz = 1,nzp
	      zooext = zooext + zoo(k,i,jz)*exz(jz) 
	    end do
        if(kticol(i))then
          jt=kti(i)
        else
          jt=kti(i)+1
        end if        
        je=kb(i)        
        do jj=jt,je
          totmac=0.0
          do m=1,nmc
            if(macrophyte_calc(jw,m))then            
              totmac = exm(m)*macrm(jj,K,I,m)+totmac
            end if
          end do
          if(cw(jj,i).gt.0.0)then
            macext=totmac/(cw(jj,i)*dlx(i)*h2(k,i))
          else
            macext=0.0
          end if          
          GAMMAj(jj,K,I) = EXH2O(JW)+SSEXT+EXOM(JW)*(LPOM(K,I)+RPOM(K,I))+ALGEX+macext+zooext
        end do          
      END IF
      GAMMA(K,I) = EXH2O(JW)+SSEXT+EXOM(JW)*(LPOM(K,I)+RPOM(K,I))+ALGEX+macext+zooext      
! v3.5 end
    END DO
  END DO

! Inorganic suspended solids settling rates

  DO I=IU,ID
    DO K=KT,KB(I)
      FPSS(K,I) = PARTP(JW)         /(PARTP(JW)*TISS(K,I)+PARTP(JW)*FE(K,I)*DO1(K,I)+1.0)
      FPFE(K,I) = PARTP(JW)*FE(K,I) /(PARTP(JW)*TISS(K,I)+PARTP(JW)*FE(K,I)*DO1(K,I)+1.0)
      SSSI(K,I) = SSSO(K-1,I)
      TOTSS0    = 0.0
      DO JS=1,NSS
        TOTSS0 = TOTSS0+SSS(JS)*FPSS(K,I)*SS(K,I,JS)
      END DO
      SSSO(K,I) = (TOTSS0+FES(JW)*FPFE(K,I))*BI(K,I)/BH2(K,I)*DO1(K,I)                           ! SW 11/08/07
      FPSS(K,I) =  FPSS(K,I)*TISS(K,I)
    END DO
  END DO

! Algal rates

  DO JA=1,NAL
    DO I=IU,ID

!**** Limiting factor

      LIGHT = (1.0-BETA(JW))*SRON(JW)*SHADE(I)/ASAT(JA)
      LAM1  =  LIGHT
      LAM2  =  LIGHT
      DO K=KT,KB(I)

!****** Limiting factor

        LAM1           = LAM2
        LAM2           = LAM1*EXP(-GAMMA(K,I)*H2(K,I))
        FDPO4          = 1.0-FPSS(K,I)-FPFE(K,I)
        ALLIM(K,I,JA)  = 2.718282*(EXP(-LAM2)-EXP(-LAM1))/(GAMMA(K,I)*H2(K,I))
        IF (AHSP(JA)  /= 0.0) APLIM(K,I,JA) =  FDPO4*PO4(K,I)/(FDPO4*PO4(K,I)+AHSP(JA)+NONZERO)
        IF (AHSN(JA)  /= 0.0) ANLIM(K,I,JA) = (NH4(K,I)+NO3(K,I))/(NH4(K,I)+NO3(K,I)+AHSN(JA)+NONZERO)
        IF (AHSSI(JA) /= 0.0) ASLIM(K,I,JA) =  DSI(K,I)/(DSI(K,I)+AHSSI(JA)+NONZERO)
        LIMIT          = MIN(APLIM(K,I,JA),ANLIM(K,I,JA),ASLIM(K,I,JA),ALLIM(K,I,JA))

!****** Algal rates

        AGR(K,I,JA) =  ATRM(K,I,JA)*AG(JA)*LIMIT
        ARR(K,I,JA) =  ATRM(K,I,JA)*AR(JA)*DO3(K,I)
        AMR(K,I,JA) = (ATRMR(K,I,JA)+1.0-ATRMF(K,I,JA))*AM(JA)
        AER(K,I,JA) =  MIN((1.0-ALLIM(K,I,JA))*AE(JA)*ATRM(K,I,JA),AGR(K,I,JA))
        IF (AS(JA) >= 0.0) THEN
          ASR(K,I,JA) =  AS(JA)*(ALG(K-1,I,JA)-ALG(K,I,JA))*BI(K,I)/BH2(K,I)
        ELSE
          ASR(K,I,JA) = -AS(JA)*(ALG(K+1,I,JA)*BI(K+1,I)/BH2(K,I)-ALG(K,I,JA)*BI(K,I)/BH2(K,I))             !SP 8/27/07
        END IF
      END DO
    END DO
    IF (AS(JA) < 0) ASR(KT,:,JA) = -AS(JA)*ALG(KT+1,:,JA)*BI(KT+1,:)*DLX(:)/VOL(KT,:)
  END DO

! v3.5 start
!****** OM stoichiometry  
    DO I=IU,ID
      DO K=KT,KB(I)
        orgpld(k,i)=0.0
        orgprd(k,i)=0.0
        orgplp(k,i)=0.0
        orgprp(k,i)=0.0
        orgnld(k,i)=0.0
        orgnrd(k,i)=0.0
        orgnlp(k,i)=0.0
        orgnrp(k,i)=0.0
        if(CAC(NLDOMP) == '      ON')then
          if(ldom(k,i).gt.0.0)orgpld(k,i)=ldomp(k,i)/ldom(k,i)
        else
          orgpld(k,i)=orgp(jw)
        end if
        if(CAC(NRDOMP) == '      ON')then
          if(rdom(k,i).gt.0.0)orgprd(k,i)=rdomp(k,i)/rdom(k,i)
        else
          orgprd(k,i)=orgp(jw)
        end if
        if(CAC(NLPOMP) == '      ON')then
          if(lpom(k,i).gt.0.0)orgplp(k,i)=lpomp(k,i)/lpom(k,i)
        else
          orgplp(k,i)=orgp(jw)
        end if
        if(CAC(NRPOMP) == '      ON')then
          if(rpom(k,i).gt.0.0)orgprp(k,i)=rpomp(k,i)/rpom(k,i)
        else
          orgprp(k,i)=orgp(jw)
        end if
        if(CAC(NLDOMN) == '      ON')then
          if(ldom(k,i).gt.0.0)orgNld(k,i)=ldomN(k,i)/ldom(k,i)
        else
          orgNld(k,i)=orgN(jw)
        end if
        if(CAC(NRDOMN) == '      ON')then
          if(rdom(k,i).gt.0.0)orgNrd(k,i)=rdomN(k,i)/rdom(k,i)
        else
          orgNrd(k,i)=orgN(jw)
        end if
        if(CAC(NLPOMN) == '      ON')then
          if(lpom(k,i).gt.0.0)orgnlp(k,i)=lpomn(k,i)/lpom(k,i)
        else
          orgnlp(k,i)=orgn(jw)
        end if
        if(CAC(NRPOMP) == '      ON')then
          if(rpom(k,i).gt.0.0)orgnrp(k,i)=rpomn(k,i)/rpom(k,i)
        else
          orgnrp(k,i)=orgn(jw)
        end if
      end do    
    end do

  do m=1,nmc
    DO I=IU,ID        
      LTCOEFm = (1.0-BETA(jw))*SRON(jw)*SHADE(I)            
      if(kticol(i))then
        jt=kti(i)
      else
        jt=kti(i)+1
      end if        
      je=kb(i)        
      do jj=jt,je                  
        lam1=ltcoefm
        lam2m(jj,kt)=lam1*exp(-gammaj(jj,kt,i)*h2(kt,i))
        lavg=(lam1-lam2m(jj,kt))/(GAMMAj(jj,kt,i)*H2(kt,i))              
        mLLIM(jj,kt,I,m) = lavg/(lavg+msat(m))            
        IF (mHSP(m)  /= 0.0.and.psed(m) < 1.0)then
          mPLIM(kt,I,m) =  FDPO4*PO4(kt,I)/(FDPO4*PO4(kt,I)+mHSP(m)+nonzero)
        else
          mPLIM(kt,I,m)=1.0
        end if
        IF (mHSN(m)  /= 0.0.and.nsed(m) < 1.0)then
          mNLIM(kt,I,m) = NH4(kt,I)/(NH4(kt,I)+mHSN(m)+nonzero)     
        else
          mNLIM(kt,I,m)=1.0
        end if
        IF (mHSc(m) /= 0.0)then
          mcLIM(kt,i,m) = co2(kt,I)/(co2(kt,I)+mHSc(m)+NONZERO)
        end if
        LIMIT          = MIN(mPLIM(kt,I,m),mNLIM(kt,I,m),mcLIM(kt,I,m),mLLIM(jj,kt,I,m))

!************* sources/sinks

        mGR(jj,Kt,I,m) = macTRM(Kt,I,m)*mG(m)*LIMIT
 
      end do

      mRR(Kt,I,m) = macTRM(Kt,I,m)*mR(m)*DO3(Kt,I)
      mMR(Kt,I,m) = (macTRMR(Kt,I,m)+1.0-mAcTRMF(Kt,I,m))*mM(m)      

      DO K=KT+1,KB(I)                
        jt=k
        je=kb(i)        
        do jj=jt,je                              
          lam1=lam2m(jj,k-1)
          lam2m(jj,k)=lam1*exp(-gammaj(jj,k,i)*h2(k,i))
          lavg=(lam1-lam2m(jj,k))/(GAMMAj(jj,k,i)*H2(k,i))              
          mLLIM(jj,K,I,m) = lavg/(lavg+msat(m))            
          IF (mHSP(m)  /= 0.0.and.psed(m) < 1.0)then
            mPLIM(K,I,m) =  FDPO4*PO4(K,I)/(FDPO4*PO4(K,I)+mHSP(m)+nonzero)
          else
            mPLIM(K,I,m)=1.0
          end if
          IF (mHSN(m)  /= 0.0.and.nsed(m) < 1.0)then
            mNLIM(K,I,m) = NH4(K,I)/(NH4(K,I)+mHSN(m)+nonzero)     
          else
             mNLIM(K,I,m)=1.0
          end if
          IF (mHSc(m) /= 0.0)then
            mcLIM(k,i,m) = co2(K,I)/(co2(K,I)+mHSc(m)+NONZERO)
          end if
          LIMIT          = MIN(mPLIM(K,I,m),mNLIM(K,I,m),mcLIM(K,I,m),mLLIM(jj,K,I,m))

!************* sources/sinks

          mGR(jj,K,I,m) = macTRM(K,I,m)*mG(m)*LIMIT
 
        end do

        mRR(K,I,m) = macTRM(K,I,m)*mR(m)*DO3(K,I)
        mMR(K,I,m) = (macTRMR(K,I,m)+1.0-mAcTRMF(K,I,m))*mM(m)
      end do     
    END DO
  END DO

  DO Jz=1,nzp
    DO I=IU,ID
      do k=kt,kb(i)
!****** zooplankton rates
        tgraze(k,i,jz)=PREFP(jz)*lpom(k,i)
        do jjz = 1, nzp
          tgraze(k,i,jz) = tgraze(k,i,jz) + prefz(jjz,jz)*zoo(k,i,jjz)          !cb 5/17/2007
	    end do
        do ja=1,nal
          tgraze(k,i,jz)=prefa(ja,jz)*alg(k,i,ja)+tgraze(k,i,jz)
        end do        
        ZMINFAC  = (1.0+SIGN(1.0,ZOO(K,I,jz)-ZOOMIN(jz)))*0.5
        ZRt(K,I,jz) =  ZooRMR(K,I,jz)*zr(jz)*ZMINFAC*DO3(K,I)
        IF (tgraze(K,I,jz) <= 0.0 .OR. O2(K,I) < 2.0) THEN
          Zmu(K,I,jz)       = 0.0                        
          Agz(K,I,1:nal,jz) = 0.0 
		  ZGZ(k,i,jz,:) = 0.0                       
          IF (O2(K,I) < 2.0) ZMINFAC = 2*ZMINFAC      
        ELSE                                          
          Zmu(K,I,jz) = MAX(ZooRM(K,I,jz)*zg(jz)*(tgraze(K,I,jz)-ZooMIN(jz))/(tgraze(K,I,jz)+ZS2p(jz)), 0.0)
          DO JA=1,NAL                                                           
!            Agz(K,I,JA,jz) = Zmu(K,I,jz)*ZOO(K,I,jz)*prefA(JA,jz)/tgraze(K,I,jz)     
          Agz(K,I,JA,jz) = Zmu(K,I,jz)*ZOO(K,I,jz)*(alg(k,i,ja)*prefA(JA,jz)/tgraze(K,I,jz))                      !  kv 5/26/2007
          END DO 
		  
          do jjz = 1,nzp ! omnivorous zooplankton 
!		    ZGZ(k,i,jjz,jz)       = Zmu(K,I,jz)*ZOO(K,I,jz)*prefZ(jjz,jz)/tgraze(K,I,jz)      !kv 5/9/2007
          ZGZ(k,i,jjz,jz)  = Zmu(K,I,jz)*zoo(k,i,jz)*(ZOO(K,I,jjz)*prefZ(jjz,jz)/tgraze(K,I,jz))         !kv 5/26/2007
          end do
         END IF                                                                  
         Zmt(K,I,jz) = MAX(1.0-ZooRMF(K,I,jz),0.02)*zm(jz)*ZMINFAC        
      end do
    END DO
  end do
! v3.5 end



RETURN

!***********************************************************************************************************************************
!**                                             G E N E R I C   C O N S T I T U E N T                                             **
!***********************************************************************************************************************************

ENTRY GENERIC_CONST (JG)
  IF (CGQ10(JG) /= 0.0) THEN
    DO I=IU,ID
      DO K=KT,KB(I)
        CGSS(K,I,JG) = -CG0DK(JG)-CG1DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)*CG(K,I,JG)+CGS(JG)*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I)
      END DO
    END DO
  ELSE
    DO I=IU,ID
      DO K=KT,KB(I)
        CGSS(K,I,JG) = -CG0DK(JG)-CG1DK(JG)*CG(K,I,JG)+CGS(JG)*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I)
      END DO
    END DO
  END IF
RETURN

!***********************************************************************************************************************************
!**                                               S U S P E N D E D   S O L I D S                                                 **
!***********************************************************************************************************************************

ENTRY SUSPENDED_SOLIDS (J)
  DO I=IU,ID
    SSR = 0.0
    IF (SEDIMENT_RESUSPENSION(J)) THEN
      FETCH = FETCHD(I,JB)
      IF (COS(PHI(JW)-PHI0(I)) < 0.0) FETCH = FETCHU(I,JB)
      FETCH = MAX(FETCH,BI(KT,I),DLX(I))
      U2    = WIND(JW)*WSC(I)*WIND(JW)*WSC(I)+NONZERO
      COEF1 = 0.53  *(G*DEPTHB(KT,I)/U2)**0.75
      COEF2 = 0.0125*(G*FETCH/U2)**0.42
      COEF3 = 0.833* (G*DEPTHB(KT,I)/U2)**0.375
      COEF4 = 0.077* (G*FETCH/U2)**0.25
      HS    = 0.283 *U2/G*0.283*TANH(COEF1)*TANH(COEF2/TANH(COEF1))
      TS    = 2.0*PI*U2/G*1.2*  TANH(COEF3)*TANH(COEF4/TANH(COEF3))
      L0    = G*TS*TS/(2.0*PI)
    END IF
    SSSS(KT,I,J) = -SSS(J)*SS(KT,I,J)*BI(KT,I)/BH2(KT,I)+SSR
    DO K=KT,KB(I)
      IF (SEDIMENT_RESUSPENSION(J)) THEN
        L1 = L0
        L  = L0*TANH(2.0*PI*DEPTHB(K,I)/L1)
        DO WHILE (ABS(L-L1) > 0.001)
          L1 = L
          L  = L0*TANH(2.0*PI*DEPTHB(K,I)/L1)
        END DO
        COEF = MIN(710.0,2.0*PI*DEPTHB(K,I)/L)
        UORB = PI*HS/TS*100.0/SINH(COEF)
        TAU  = 0.003*UORB*UORB
        IF (TAU-TAUCR(J) > 0.0) EPSILON = MAX(0.0,0.008/49.0*(TAU-TAUCR(J))**3*10000.0/DLT)
		if(k == kb(i))then   ! SW 4/18/07
		SSR = EPSILON*DLX(I)*BI(K,I)/VOL(K,I)
		else
        SSR = EPSILON*DLX(I)*(BI(K,I)-BI(K+1,I))/VOL(K,I)
		endif
      END IF
      SSSS(K,I,J) = SSS(J)*(SS(K-1,I,J)-SS(K,I,J))*BI(K,I)/BH2(K,I)+SSR
    END DO
    IF (SEDIMENT_RESUSPENSION(J)) SSR = EPSILON*DLX(I)*BI(KB(I),I)/VOL(KB(I),I)
    SSSS(KB(I),I,J) = SSS(J)*(SS(KB(I)-1,I,J)-SS(KB(I),I,J))/H(KB(I),JW)+SSR
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      P H O S P H O R U S                                                      **
!***********************************************************************************************************************************

ENTRY PHOSPHORUS
  PO4AR(:,IU:ID) = 0.0; PO4AG(:,IU:ID) = 0.0; PO4ER(:,IU:ID) = 0.0; PO4EG(:,IU:ID) = 0.0; PO4BOD(:,IU:ID) = 0.0
  po4mr(:,iu:id) = 0.0; po4mg(:,iu:id) = 0.0; po4zr(:,iu:id)=0.0   !v3.5

  DO I=IU,ID
    DO K=KT,KB(I)
      DO JCB=1,NBOD
        PO4BOD(K,I) = PO4BOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODP(JCB)
      END DO
      DO JA=1,NAL
        PO4AG(K,I) = PO4AG(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*AP(JA)
        PO4AR(K,I) = PO4AR(K,I)+ARR(K,I,JA)*ALG(K,I,JA)*AP(JA)
      END DO
      DO JE=1,NEP
        PO4EG(K,I) = PO4EG(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*EP(JE)
        PO4ER(K,I) = PO4ER(K,I)+ERR(K,I,JE)*EPC(K,I,JE)*EP(JE)
      END DO
      PO4EP(K,I)  = PO4ER(K,I)-PO4EG(K,I)
      PO4AP(K,I)  = PO4AR(K,I)-PO4AG(K,I)
! v3.5 start
      PO4POM(K,I) = ORGPLP(k,i)*LPOMD(K,I)+orgprp(k,i)*RPOMD(K,I)
      PO4DOM(K,I) = ORGPLD(k,i)*LDOMD(K,I)+orgprd(k,i)*RDOMD(K,I)
! v3.5 end
      PO4OM(K,I)  = PO4POM(K,I)+PO4DOM(K,I)
! v3.5 start!           
      PO4SD(K,I)  = SEDDp(K,I)      
! v3.5 end
      PO4SR(K,I)  = PO4R(JW)*SODD(K,I)*DO2(K,I)
      PO4NS(K,I)  = SSSI(K,I)*PO4(K-1,I)-SSSO(K,I)*PO4(K,I)
! v3.5 start      
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          if(k.eq.kt)then
            jt=kti(i)
          else
            jt=k
          end if
          je=kb(i)
          do jj=jt,je
            po4mg(k,i)= po4mg(k,i)+mgr(jj,k,i,m)*macrm(jj,k,i,m)*mp(m)*(1.0-psed(m))     
            po4mr(k,i)= po4mr(k,i)+mrr(k,i,m)*macrm(jj,k,i,m)*mp(m)
          end do
        end if
      end do      
      po4mr(k,i)=po4mr(k,i)/(dlx(i)*bh(k,i))
      po4mg(k,i)=po4mg(k,i)/(dlx(i)*bh(k,i))      
      do jz = 1,nzp
        PO4zR(K,I) = po4zr(k,i) + zrt(K,I,jz)*zoo(K,I,jz)*zp(jz)
	  end do
! v3.5 end     

      PO4SS(K,I)  = PO4AP(K,I)+PO4EP(K,I)+PO4OM(K,I)+PO4SD(K,I)+PO4SR(K,I)+PO4NS(K,I)+PO4BOD(K,I)  &
                    +po4mr(k,i)-po4mg(k,i) +po4zr(k,i)     !v3.5

    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                        A M M O N I U M                                                        **
!***********************************************************************************************************************************

ENTRY AMMONIUM
  NH4AG(:,IU:ID) = 0.0; NH4AR(:,IU:ID) = 0.0; NH4ER(:,IU:ID) = 0.0; NH4EG(:,IU:ID) = 0.0; NH4BOD(:,IU:ID) = 0.0
  nh4mg(:,iu:id) = 0.0; nh4mr(:,iu:id) = 0.0; nh4zr(:,iu:id)=0.0   !v3.5
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JCB=1,NBOD
        NH4BOD(K,I) =  NH4BOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODN(JCB)
      END DO
      DO JA=1,NAL
        NH4PR = NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO)
        IF (ANEQN(JA).EQ.2) NH4PR      = NH4(K,I)*NO3(K,I)/((ANPR(JA)+NH4(K,I))*(ANPR(JA)+NO3(K,I)))+NH4(K,I)*ANPR(JA)/((NO3(K,I)  &
                                        +NH4(K,I)+NONZERO)*(ANPR(JA)+NO3(K,I)))
        IF (AHSN(JA) > 0.0) NH4AG(K,I) = NH4AG(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*AN(JA)*NH4PR
        NH4AR(K,I) = NH4AR(K,I)+ARR(K,I,JA)*ALG(K,I,JA)*AN(JA)
      END DO
      DO JE=1,NEP
        NH4PR = NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO)
        IF (ENEQN(JE) == 2) NH4PR      = NH4(K,I)*NO3(K,I)/((ENPR(JE)+NH4(K,I))*(ENPR(JE)+NO3(K,I)))+NH4(K,I)*ENPR(JE)/((NO3(K,I)  &
                                        +NH4(K,I)+NONZERO)*(ENPR(JE)+NO3(K,I)))
        IF (EHSN(JE) > 0.0) NH4EG(K,I) = NH4EG(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*EN(JE)*NH4PR
        NH4ER(K,I) = NH4ER(K,I)+ERR(K,I,JE)*EPC(K,I,JE)*EN(JE)
      END DO
      NH4EP(K,I)  =  NH4ER(K,I) -NH4EG(K,I)
      NH4AP(K,I)  =  NH4AR(K,I) -NH4AG(K,I)
! v3.5 start
      NH4DOM(K,I) = LDOMD(K,I)*orgnld(k,i) +RDOMD(K,I)*ORGNrd(k,i)
      NH4POM(K,I) = LPOMD(K,I)*orgnlp(k,i) +RPOMD(K,I)*ORGNrp(k,i)
! v3.5 end
      NH4OM(K,I)  =  NH4DOM(K,I)+NH4POM(K,I)
! v3.5 start      
      NH4SD(K,I)  =  SEDDn(K,I)      
! v3.5 end
      NH4SR(K,I)  =  NH4R(JW) *SODD(K,I)*DO2(K,I)
! v3.5 start
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          if(k.eq.kt)then
            jt=kti(i)
          else
            jt=k
          end if
          je=kb(i)
          do jj=jt,je
            nh4mr(k,i)= nh4mr(k,i)+mrr(k,i,m)*macrm(jj,k,i,m)*mn(m)     
            nh4mg(k,i)= nh4mg(k,i)+mgr(jj,k,i,m)*macrm(jj,k,i,m)*mn(m)*(1.0-nsed(m))        
          end do
        end if
      end do            
      nh4mr(k,i)=nh4mr(k,i)/(dlx(i)*bh(k,i))
      nh4mg(k,i)=nh4mg(k,i)/(dlx(i)*bh(k,i))      
	  do jz = 1,nzp
	    nh4zr(k,i) = nh4zr(k,i) + zrt(k,i,jz)*zoo(k,i,jz)*zn(jz) ! Baulk
	  end do
      NH4SS(K,I)  =  NH4AP(K,I)+NH4EP(K,I)+NH4OM(K,I)+NH4SD(K,I)+NH4SR(K,I)+NH4BOD(K,I)-NH4D(K,I)  &
         +nh4mr(k,i)-nh4mg(k,i) +nh4zr(k,i)	     !v3.5
! v3.5 end
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                          N I T R A T E                                                        **
!***********************************************************************************************************************************

ENTRY NITRATE
  NO3AG(:,IU:ID) = 0.0; NO3EG(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        NO3PR = 1.0-NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO)
        IF (ANEQN(JA).EQ.2)  NO3PR      = 1.0-(NH4(K,I)*NO3(K,I)/((ANPR(JA)+NH4(K,I))*(ANPR(JA)+NO3(K,I)))+NH4(K,I)*ANPR(JA)       &
                                          /((NO3(K,I)+NH4(K,I)+NONZERO)*(ANPR(JA)+NO3(K,I))))
        IF (AHSN(JA).GT.0.0) NO3AG(K,I) = NO3AG(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*NO3PR*AN(JA)
      END DO
      DO JE=1,NEP
        NO3PR = 1.0-NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO)
        IF (ENEQN(JE).EQ.2)  NO3PR      = 1.0-(NH4(K,I)*NO3(K,I)/((ENPR(JE)+NH4(K,I))*(ENPR(JE)+NO3(K,I)))+NH4(K,I)*ENPR(JE)       &
                                          /((NO3(K,I)+NH4(K,I)+NONZERO)*(ENPR(JE)+NO3(K,I))))
        IF (EHSN(JE).GT.0.0) NO3EG(K,I) = NO3EG(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*NO3PR*EN(JE)
      END DO
      if(k == kb(i)) then      ! SW 4/18/07
      NO3SED(K,I) = NO3(K,I)*NO3S(JW)*NO3TRM(K,I)*(BI(K,I))/BH2(K,I)
	  else
      NO3SED(K,I) = NO3(K,I)*NO3S(JW)*NO3TRM(K,I)*(BI(K,I)-BI(K+1,I))/BH2(K,I)
	  endif
      NO3SS(K,I)  = NH4D(K,I)-NO3D(K,I)-NO3AG(K,I)-NO3EG(K,I)-NO3SED(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  D I S S O L V E D   S I L I C A                                              **
!***********************************************************************************************************************************

ENTRY DISSOLVED_SILICA
  DSIAG(:,IU:ID) = 0.0; DSIEG(:,IU:ID) = 0.0; DSIBOD = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        DSIAG(K,I) = DSIAG(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*ASI(JA)
      END DO
      DO JE=1,NEP
        DSIEG(K,I) = DSIEG(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*ESI(JE)
      END DO
      DSID(K,I)  =  PSIDK(JW)*PSI(K,I)
      DSISD(K,I) =  SEDD(K,I)*ORGSI(JW)
      DSISR(K,I) =  DSIR(JW)*SODD(K,I)*DO2(K,I)
      DSIS(K,I)  = (SSSI(K,I)*DSI(K-1,I)-SSSO(K,I)*DSI(K,I))*PARTSI(JW)
      DSISS(K,I) =  DSID(K,I)+DSISD(K,I)+DSISR(K,I)+DSIS(K,I)+DSIBOD-DSIAG(K,I)-DSIEG(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                P A R T I C U L A T E   S I L I C A                                            **
!***********************************************************************************************************************************

ENTRY PARTICULATE_SILICA
  PSIAM(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        PSIAM(K,I) = PSIAM(K,I)+AMR(K,I,JA)*PSI(K,I)*ASI(JA)
      END DO
      PSID(K,I)  = PSIDK(JW)*PSI(K,I)
      PSINS(K,I) = PSIS(JW)*(PSI(K-1,I)*DO1(K-1,I)-PSI(K,I)*DO1(K,I))*BI(K,I)/BH2(K,I)
      PSISS(K,I) = PSIAM(K,I)-PSID(K,I)+PSINS(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                            I R O N                                                            **
!***********************************************************************************************************************************

ENTRY IRON
  DO I=IU,ID
    DO K=KT,KB(I)
      FENS(K,I) = FES(JW)*(FE(K-1,I)*DO1(K-1,I)-FE(K,I)*DO1(K,I))*BI(K,I)/BH2(K,I)
      FESR(K,I) = FER(JW)*SODD(K,I)*DO2(K,I)
      FESS(K,I) = FESR(K,I)+FENS(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                       L A B I L E   D O M                                                     **
!***********************************************************************************************************************************

ENTRY LABILE_DOM
  LDOMAP(:,IU:ID) = 0.0; LDOMEP(:,IU:ID) = 0.0; ldommac(:,iu:id)= 0.0  !v3.5
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        LDOMAP(K,I) = LDOMAP(K,I)+(AER(K,I,JA)+(1.0-APOM(JA))*AMR(K,I,JA))*ALG(K,I,JA)
      END DO
      DO JE=1,NEP
        LDOMEP(K,I) = LDOMEP(K,I)+(EER(K,I,JE)+(1.0-EPOM(JE))*EMR(K,I,JE))*EPC(K,I,JE)
      END DO
! v3.5 start      
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          if(k.eq.kt)then
            jt=kti(i)
          else
            jt=k
          end if
          je=kb(i)
          do jj=jt,je
            ldommac(k,i)=ldommac(k,i)+(1.0-mpom(m))*mmr(k,i,m)*macrm(jj,k,i,m)     
          end do
        end if
      end do          
      ldommac(k,i)=ldommac(k,i)/(dlx(i)*bh(k,i))          
      LDOMSS(K,I) = LDOMAP(K,I)+LDOMEP(K,I)-LDOMD(K,I)-LRDOMD(K,I)+ldommac(k,i)
! v3.5 end
      
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   D O M                                                  **
!***********************************************************************************************************************************

ENTRY REFRACTORY_DOM
  DO I=IU,ID
    DO K=KT,KB(I)
      RDOMSS(K,I) = LRDOMD(K,I)-RDOMD(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      L A B I L E   P O M                                                      **
!***********************************************************************************************************************************

ENTRY LABILE_POM
!  LPOMAP(:,IU:ID) = 0.0
  LPOMAP(:,IU:ID) = 0.0; lpomep(:,iu:id) = 0.0;   lpommac(:,iu:id) = 0.0; lpzooin(:,iu:id)=0.0;lpzooout(:,iu:id)=0.0  !v3.5  ! cb 5/19/06
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        LPOMAP(K,I) = LPOMAP(K,I)+APOM(JA)*(AMR(K,I,JA)*ALG(K,I,JA))
      END DO
      DO JE=1,NEP                                                          ! cb 5/19/06
        LPOMEP(K,I) = LPOMEP(K,I)+EPOM(JE)*(EMR(K,I,JE)*EPC(K,I,JE))       ! cb 5/19/06            
      END DO                                                               ! cb 5/19/06
      LPOMNS(K,I) = POMS(JW)*(LPOM(K-1,I)-LPOM(K,I))*BI(K,I)/BH2(K,I)
! v3.5 start      
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          jt=k
          je=kb(i)
          do jj=jt,je
            lpommac(k,i)=lpommac(k,i)+mpom(m)*lrpmac(m)*mmr(k,i,m)*macrm(jj,k,i,m)     
          end do
        end if
      end do            
      lpommac(k,i)=lpommac(k,i)/(dlx(i)*bh(k,i))
      do jz = 1,nzp
        if(tgraze(k,i,jz) > 0.0)then      
          lpzooout(k,i)=lpzooout(k,i)+zoo(k,i,jz)*(zmt(k,i,jz)+(zmu(k,i,jz)-(zmu(k,i,jz)*zeff(jz))))
          lpzooin(k,i)=lpzooin(k,i)+zoo(k,i,jz)*prefp(jz)*zmu(k,i,jz)*LPOM(k,i)/tgraze(k,i,jz)
        else
          lpzooout(k,i)=lpzooout(k,i)+zoo(k,i,jz)*(zmt(k,i,jz)+(zmu(k,i,jz)-(zmu(k,i,jz)*zeff(jz))))
          lpzooin(k,i)=0.0
        end if
      end do
      LPOMSS(K,I) = LPOMAP(K,I)+lpomep(k,i)-LPOMD(K,I)+LPOMNS(K,I)-LRPOMD(K,I)+lpommac(k,i)     &
                  + lpzooout(k,i)-lpzooin(k,i)
! v3.5 end      
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   P O M                                                  **
!***********************************************************************************************************************************

ENTRY REFRACTORY_POM
  rpommac(:,iu:id) = 0.0 !v3.5
  DO I=IU,ID
    DO K=KT,KB(I)
      RPOMNS(K,I) = POMS(JW)*(RPOM(K-1,I)-RPOM(K,I))*BI(K,I)/BH2(K,I)
! v3.5 start            
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          jt=k
          je=kb(i)
          tmac=0.0
          do jj=jt,je
            rpommac(k,i)=rpommac(k,i)+mpom(m)*(1.0-lrpmac(m))*mmr(k,i,m)*macrm(jj,k,i,m)     
          end do
        end if
      end do
      rpommac(k,i)=rpommac(k,i)/(dlx(i)*bh(k,i))
      RPOMSS(K,I) = LRPOMD(K,I)+RPOMNS(K,I)-RPOMD(K,I)+rpommac(k,i)
! v3.5 end
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                         A L G A E                                                             **
!***********************************************************************************************************************************

ENTRY ALGAE (J)      
  agzt(:,iu:id,j) = 0.0 !v3.5      
  DO I=IU,ID
    DO K=KT,KB(I)
! v3.5 start      
      do jz = 1,nzp
!	   agzt(k,i,j) = agzt(k,i,j) + agz(k,i,j,jz)*ALG(K,I,J)                                                                         !MM 8/13/06
	  agzt(k,i,j) = agzt(k,i,j) + agz(k,i,j,jz)                       ! cb 5/26/07
	  end do
      ASS(K,I,J) = ASR(K,I,J)+(AGR(K,I,J)-AER(K,I,J)-AMR(K,I,J)-ARR(K,I,J))*ALG(K,I,J)-agzt(k,i,j)	  
! v3.5 end
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                            B I O C H E M I C A L   O 2   D E M A N D                                          **
!***********************************************************************************************************************************

ENTRY BIOCHEMICAL_O2_DEMAND(JBOD)
  if(jbod == 1)cbodns(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
! v3.5 start
      cbodset = cbodS(jbod)*(cbod(K-1,I,jbod)-cbod(K,I,jbod))*BI(K,I)/BH2(K,I)
      cbodns(k,i)=cbodns(k,i)+cbodset
      CBODSS(K,I,JBOD) = -CBODD(K,I,JBOD)*CBOD(K,I,JBOD)+cbodset
! v3.5 end
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                D I S S O L V E D   O X Y G E N                                                **
!***********************************************************************************************************************************

ENTRY DISSOLVED_OXYGEN
  DOAP(:,IU:ID) = 0.0; DOAR(:,IU:ID) = 0.0; DOEP(:,IU:ID) = 0.0; DOER(:,IU:ID) = 0.0; DOBOD(:,IU:ID) = 0.0
  domp(:,iu:id) = 0.0; domr(:,iu:id) = 0.0; dozr(:,iu:id)=0.0   !v3.5

  DO I=IU,ID
    DOSS(KT,I) = 0.0
    DO K=KT,KB(I)
      DO JCB=1,NBOD
        DOBOD(K,I) = DOBOD(K,I)+RBOD(JCB)*CBODD(K,I,JCB)*CBOD(K,I,JCB)
      END DO
      DO JA=1,NAL
        DOAP(K,I) = DOAP(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*O2AG(JA)
        DOAR(K,I) = DOAR(K,I)+ARR(K,I,JA)*ALG(K,I,JA)*O2AR(JA)
      END DO
      DO JE=1,NEP    
        DOEP(K,I) = DOEP(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*O2EG(JE)
        DOER(K,I) = DOER(K,I)+ERR(K,I,JE)*EPC(K,I,JE)*O2ER(JE)
      END DO
! v3.5 start               
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          if(k.eq.kt)then
            jt=kti(i)
          else
            jt=k
          end if
          je=kb(i)
          do jj=jt,je
            domp(k,i)=domp(k,i)+mgr(jj,k,i,m)*macrm(jj,k,i,m)*o2mg(m)
            domr(k,i)=domr(k,i)+mrr(k,i,m)*macrm(jj,k,i,m)*o2mr(m)
          end do
        end if
      end do              
      domp(k,i)=domp(k,i)/(dlx(i)*bh(k,i))
      domr(k,i)=domr(k,i)/(dlx(i)*bh(k,i))        
! v3.5 end
      DOPOM(K,I) = (LPOMD(K,I)+RPOMD(K,I))*O2OM(JW)
      DODOM(K,I) = (LDOMD(K,I)+RDOMD(K,I))*O2OM(JW)
      DOOM(K,I)  =  DOPOM(K,I)+DODOM(K,I)+DOBOD(K,I)
      DONIT(K,I) =  NH4D(K,I)*O2NH4(JW)
      DOSED(K,I) =  SEDD(K,I)*O2OM(JW)
      DOSOD(K,I) =  SODD(K,I)*DO3(K,I)
      ! v3.5 start
     do jz = 1, nzp
      dozr(k,i)  = dozr(k,i)+zrt(k,i,jz)*zoo(k,i,jz)*o2zr(jz)
	 end do
    DOSS(K,I)  =  DOAP(K,I)+DOEP(K,I)-DOAR(K,I)-DOER(K,I)-DOOM(K,I)-DONIT(K,I)-DOSOD(K,I)-DOSED(K,I)  &
                    +domp(k,i)-domr(k,i)-dozr(k,i)     
! v3.5 end
    END DO
    DOSAT = SATO(T1(KT,I),TDS(KT,I),PALT(I),SALT_WATER(JW))
    IF (.NOT. ICE(I)) THEN
      CALL GAS_TRANSFER
      O2EX       =  REAER(I)
      DOAE(KT,I) = (DOSAT-O2(KT,I))*O2EX*BI(KT,I)/BH2(KT,I)
      DOSS(KT,I) =  DOSS(KT,I)+DOAE(KT,I)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                              I N O R G A N I C   C A R B O N                                                  **
!***********************************************************************************************************************************

ENTRY INORGANIC_CARBON
  TICAP(:,IU:ID) = 0.0; TICEP(:,IU:ID) = 0.0; TICBOD(:,IU:ID) = 0.0
  ticmc(:,iu:id) = 0.0; ticzr(:,iu:id)=0.0  !v3.5
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JCB=1,NBOD
        TICBOD(K,I) = TICBOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODC(JCB)
      END DO
      DO JA=1,NAL
        TICAP(K,I) = TICAP(K,I)+AC(JA)*(ARR(K,I,JA)-AGR(K,I,JA))*ALG(K,I,JA)
      END DO
      DO JE=1,NEP
        TICEP(K,I) = TICEP(K,I)+EC(JE)*(ERR(K,I,JE)-EGR(K,I,JE))*EPC(K,I,JE)
      END DO
! v3.5 start
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          if(k.eq.kt)then
            jt=kti(i)
          else
            jt=k
          end if
          je=kb(i)
          do jj=jt,je
            ticmc(k,i)=ticmc(k,i)+(mrr(k,i,m)-mgr(jj,k,i,m))*macrm(jj,k,i,m)*mc(m)     
          end do
        end if
      end do      
      ticmc(k,i)=ticmc(k,i)/(dlx(i)*bh(k,i))          
      do jz = 1,nzp
        ticzr(k,i)=ticzr(k,i)+zrt(k,i,jz)*zoo(k,i,jz)*zc(jz) !mlm
	  end do      
      TICSS(K,I) = TICAP(K,I)+TICEP(K,I)+seddc(k,i)+ORGC(JW)*(LPOMD(K,I)+RPOMD(K,I)+LDOMD(K,I)+RDOMD(K,I)) &
                   +CO2R(JW)*SODD(K,I)*DO3(K,I)+TICBOD(K,I)+ticmc(k,i)+ticzr(k,i)      	  
! v3.5 end      
    END DO
    IF (.NOT. ICE(I)) THEN
      IF (REAER(I) == 0.0) CALL GAS_TRANSFER
      CO2EX       = REAER(I)*0.923
      TICSS(KT,I) = TICSS(KT,I)+CO2EX*(0.286*EXP(-0.0314*(T2(KT,I))*PALT(I))-CO2(KT,I))*BI(KT,I)/BH2(KT,I)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T                                                          **
!***********************************************************************************************************************************

ENTRY SEDIMENT
  SEDAS(:,IU:ID) = 0.0; LPOMEP(:,IU:ID) = 0.0; SEDCB(:,IU:ID) = 0.0 !v3.5
! v3.5 start
  DO I=IU,ID
    sedsi=0.0
    DO K=KT,KB(I)
      DO JA=1,NAL
	    if(k == kb(i))then    ! SW 4/18/07
        SEDAS(K,I) = SEDAS(K,I)+MAX(AS(JA),0.0)*ALG(K,I,JA)*BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
		else
        SEDAS(K,I) = SEDAS(K,I)+MAX(AS(JA),0.0)*ALG(K,I,JA)*BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
		endif
      END DO
      sedem = 0.0   ! cb 5/19/06
      DO JE=1,NEP
!        LPOMEP(K,I) = LPOMEP(K,I)+EPOM(JE)*(EMR(K,I,JE)*EPC(K,I,JE))
        sedem = sedem+ebr(k,i,je)/h1(k,i)*EPC(K,I,JE)    ! cb 5/19/06
      END DO
! v3.5 start
      do jd=1,nbod
        SEDcb(K,I) = SEDcb(K,I)+MAX(cbods(jd),0.0)*cbod(K,I,Jd)*BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      end do
! v3.5 end
      SEDOMS(K,I) = pomS(JW)*(LPOM(K,I)+RPOM(K,I))*BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))  !cb 10/22/06
      SEDSO       = sedS(JW)*SED(K,I)*BI(K+1,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      sedbr(K,I)  = sedb(jw)*sed(k,i)                                                                                ! SW 10/15/07
      SEDNS(K,I)  = SEDSI-SEDSO
      SEDSI       = SEDSO
!      SED(K,I)    = MAX(SED(K,I)+(LPOMEP(K,I)+SEDAS(K,I)+SEDOMS(K,I)+SEDNS(K,I)-SEDD(K,I))*DLT,0.0)
      SED(K,I)    = MAX(SED(K,I)+(sedem+SEDAS(K,I)+sedcb(k,i)+SEDOMS(K,I)+SEDNS(K,I)-SEDD(K,I)-sedbr(k,i))*DLT,0.0)   ! cb 11/30/06
    END DO
  END DO
RETURN

! v3.5 start

!***********************************************************************************************************************************
!**                                                      S E D I M E N T   P H O S P H O R U S                                    **
!***********************************************************************************************************************************

ENTRY SEDIMENTP
  SEDASp(:,IU:ID) = 0.0; LPOMEPp(:,IU:ID) = 0.0; sedcbp(:,IU:ID) = 0.0
  DO I=IU,ID
    sedsip=0.0
    DO K=KT,KB(I)
      DO JA=1,NAL
        SEDASp(K,I) = SEDASp(K,I)+MAX(AS(JA),0.0)*ap(ja)*ALG(K,I,JA)*BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      END DO
      DO JE=1,NEP
        LPOMEPp(K,I) = LPOMEPp(K,I)+EPOM(JE)*ep(je)*(EMR(K,I,JE)*EPC(K,I,JE))
      END DO
      do jd=1,nbod
        sedcbp(k,i)=sedcbp(k,i)+MAX(cbods(jd),0.0)*bodp(jd)*cbod(K,I,jd)*BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      end do
      SEDOMSp(K,I) = poms(JW)*(LPOMp(K,I)+RPOMp(K,I))*BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))  !cb 10/22/06
      SEDSOp       = seds(JW)*SEDp(K,I)*BI(K+1,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      SEDNSp(K,I)  = SEDSIp-SEDSOp
      SEDSIp       = SEDSOp
      SEDp(K,I)    = MAX(SEDp(K,I)+(LPOMEPp(K,I)+SEDASp(K,I)+SEDOMSp(K,I)+sedcbp(k,i)+SEDNSp(K,I)-SEDDp(K,I)   &
                     -sedbrp(k,i))*DLT,0.0)                                                                 !cb 11/30/06
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T   N I T R O G E N                                        **
!***********************************************************************************************************************************

ENTRY SEDIMENTn
  SEDASn(:,IU:ID) = 0.0; LPOMEPn(:,IU:ID) = 0.0; sedcbn(:,IU:ID) = 0.0
  DO I=IU,ID
    sedsin=0.0
    DO K=KT,KB(I)
      DO JA=1,NAL
        SEDASn(K,I) = SEDASn(K,I)+MAX(AS(JA),0.0)*an(ja)*ALG(K,I,JA)*BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      END DO
      DO JE=1,NEP
        LPOMEPn(K,I) = LPOMEPn(K,I)+EPOM(JE)*en(je)*(EMR(K,I,JE)*EPC(K,I,JE))
      END DO
      do jd=1,nbod
        sedcbn(k,i)=sedcbn(k,i)+MAX(cbods(jd),0.0)*bodn(jd)*cbod(K,I,jd)*BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      end do
      SEDOMSn(K,I) = poms(JW)*(LPOMn(K,I)+RPOMn(K,I))*BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))  !cb 10/22/06
      SEDSOn       = seds(JW)*SEDn(K,I)*BI(K+1,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      if(k == kb(i)) then      ! SW 12/16/07
      sedNO3(K,I)  = NO3(K,I)*NO3S(JW)*NO3TRM(K,I)*(BI(K,I))/BH2(K,I)
	  else
      sedNO3(K,I)  = NO3(K,I)*NO3S(JW)*NO3TRM(K,I)*(BI(K,I)-BI(K+1,I))/BH2(K,I)
	  endif
      SEDNSn(K,I)  = SEDSIn-SEDSOn
      SEDSIn       = SEDSOn
      SEDn(K,I)    = MAX(SEDn(K,I)+(LPOMEPn(K,I)+SEDASn(K,I)+SEDOMSn(K,I)+sedcbn(k,i)+SEDNSn(K,I)+sedno3(k,i)-   &
                     SEDDn(K,I)-sedbrn(k,i))*DLT,0.0)  !cb 11/30/06
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T   C A R B O N                                            **
!***********************************************************************************************************************************

ENTRY SEDIMENTC
  SEDASc(:,IU:ID) = 0.0; LPOMEPc(:,IU:ID) = 0.0; sedcbc(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      sedsip=0.0
      DO JA=1,NAL
        SEDASc(K,I) = SEDASc(K,I)+MAX(AS(JA),0.0)*ac(ja)*ALG(K,I,JA)*BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      END DO
      DO JE=1,NEP
        LPOMEPc(K,I) = LPOMEPc(K,I)+EPOM(JE)*ec(je)*(EMR(K,I,JE)*EPC(K,I,JE))
      END DO
      do jd=1,nbod
        sedcbc(k,i)=sedcbc(k,i)+MAX(cbods(jd),0.0)*bodc(jd)*cbod(K,I,jd)*BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      end do
      SEDOMSc(K,I) = poms(JW)*orgc(jw)*(LPOM(K,I)+RPOM(K,I))*BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))   !cb 10/22/06
      SEDSOc       = seds(JW)*SEDc(K,I)*BI(K+1,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      SEDNSc(K,I)  = SEDSIc-SEDSOc
      SEDSIc       = SEDSOc
      SEDc(K,I)    = MAX(SEDc(K,I)+(LPOMEPc(K,I)+SEDASc(K,I)+SEDOMSc(K,I)+sedcbc(k,i)+SEDNSc(K,I)-SEDDc(K,I)    &
                     -sedbrc(k,i))*DLT,0.0)                                                                   !cb 11/30/06
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T   D E C A Y    R A T E                                   **
!***********************************************************************************************************************************

ENTRY SEDIMENT_decay_rate  
  DO I=IU,ID
    sedsidk=0.0
    DO K=KT,KB(I)
      sedsum=0.0
      DO JA=1,NAL
        SEDsum = sedsum + MAX(AS(JA),0.0)*ALG(K,I,JA)*BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I)) * lpomdk(jw)
      END DO
      DO JE=1,NEP
        SEDsum = sedsum + EPOM(JE)*(EMR(K,I,JE)*EPC(K,I,JE)) * lpomdk(jw)
      END DO
      do jd=1,nbod
        sedsum = sedsum+MAX(cbods(jd),0.0)*cbod(K,I,Jd)*BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      end do
      sedsum = sedsum - sedb(jw) * sed(k,i) * sdkv(k,i)              !cb 11/30/06
      sedsum = sedsum+poms(JW)*(LPOM(K,I)*lpomdk(jw)+RPOM(K,I)*rpomdk(jw))*BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))  ! cb 10/22/06
      SEDsodk       = seds(JW)*SED(K,I)*BI(K+1,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))*sdkv(k,i)      
      sednsdk       = sedsidk-sedsodk
      sedsidk       = sedsodk
      sedsum        = sedsum + sednsdk
      SEDdktot(K,I)    = MAX(SEDdktot(K,I)+sedsum*DLT,0.0)
      if(sed(k,i).ge.0.001)sdkv(k,i)    = max(seddktot(k,i)/(sed(k,i)+nonzero),0.0)  ! cb 10/23/06
    END DO
  END DO
RETURN

! v3.5 end

!***********************************************************************************************************************************
!*                                                         E P I P H Y T O N                                                      **
!***********************************************************************************************************************************

ENTRY EPIPHYTON (J)
  DO I=IU,ID

!** Limiting factor

    LIGHT = (1.0-BETA(JW))*SRON(JW)*SHADE(I)/ESAT(J)
    LAM2  =  LIGHT
    LAM1  =  LIGHT
    DO K=KT,KB(I)

!**** Limiting factor

      LAM1          = LAM2
      LAM2          = LAM1*EXP(-GAMMA(K,I)*H1(K,I))
      FDPO4         = 1.0-FPSS(K,I)-FPFE(K,I)
      ELLIM(K,I,J)  = 2.718282*(EXP(-LAM2)-EXP(-LAM1))/(GAMMA(K,I)*H1(K,I))
      IF (EHSP(J)  /= 0.0) EPLIM(K,I,J) =  FDPO4*PO4(K,I)/(FDPO4*PO4(K,I)+EHSP(J)+NONZERO)
      IF (EHSN(J)  /= 0.0) ENLIM(K,I,J) = (NH4(K,I)+NO3(K,I))/(NH4(K,I)+NO3(K,I)+EHSN(J)+NONZERO)
      IF (EHSSI(J) /= 0.0) ESLIM(K,I,J) =  DSI(K,I)/(DSI(K,I)+EHSSI(J)+NONZERO)
      LIMIT         =  MIN(EPLIM(K,I,J),ENLIM(K,I,J),ESLIM(K,I,J),ELLIM(K,I,J))
      BLIM          =  1.0-EPD(K,I,J)/(EPD(K,I,J)+EHS(J))

!**** Sources/sinks

      EGR(K,I,J) =  MIN(ETRM(K,I,J)*EG(J)*LIMIT*BLIM,PO4(K,I)/(EP(J)*DLT*EPD(K,I,J)/H1(KT,I)+NONZERO),(NH4(K,I)+NO3(K,I))/(EN(J)   &
                    *DLT*EPD(K,I,J)/H1(K,I)+NONZERO))
      ERR(K,I,J) =  ETRM(K,I,J)*ER(J)*DO3(K,I)
      EMR(K,I,J) = (ETRMR(K,I,J)+1.0-ETRMF(K,I,J))*EM(J)
      EER(K,I,J) =  MIN((1.0-ELLIM(K,I,J))*EE(J)*ETRM(K,I,J),EGR(K,I,J))
!      EPD(K,I,J) =  MAX(EPD(K,I,J)+EPD(K,I,J)*(EGR(K,I,J)-ERR(K,I,J)-EMR(K,I,J)-EER(K,I,J)-EBR(K,I,J)/(H1(K,I)*0.0025))*DLT,0.0)
      EPD(K,I,J) =  MAX(EPD(K,I,J)+EPD(K,I,J)*(EGR(K,I,J)-ERR(K,I,J)-EMR(K,I,J)-EER(K,I,J)-EBR(K,I,J)/H1(K,I))*DLT,0.00)   ! cb 5/18/06
      if(k == kb(i)) then      ! SW 12/16/07
      EPM(K,I,J) =  EPD(K,I,J)*(BI(K,I)+2.0*H1(K,I))*DLX(I)
	  else
      EPM(K,I,J) =  EPD(K,I,J)*(BI(K,I)-BI(K+1,I)+2.0*H1(K,I))*DLX(I)
	  endif
      EPC(K,I,J) =  EPM(K,I,J)/VOL(K,I)
    END DO
  END DO
RETURN

! v3.5 start
!***********************************************************************************************************************************
!**                                                       L A B I L E   D O M   P H O S P H O R U S                               **
!***********************************************************************************************************************************

ENTRY LABILE_DOM_P
  LDOMpAP(:,IU:ID) = 0.0; LDOMpEP(:,IU:ID) = 0.0; ldompmp(:,iu:id)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        LDOMpAP(K,I) = LDOMpAP(K,I)+(AER(K,I,JA)+(1.0-APOM(JA))*AMR(K,I,JA))*ALG(K,I,JA)*ap(ja)
      END DO
      DO JE=1,NEP
        LDOMpEP(K,I) = LDOMpEP(K,I)+(EER(K,I,JE)+(1.0-EPOM(JE))*EMR(K,I,JE))*EPC(K,I,JE)*ep(je)
      END DO
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          if(k.eq.kt)then
            jt=kti(i)
          else
            jt=k
          end if
          je=kb(i)
          do jj=jt,je
            ldompmp(k,i)=ldompmp(k,i)+(1.0-mpom(m))*mmr(k,i,m)*macrm(jj,k,i,m)*mp(m)
          end do
        end if
      end do          
      ldompmp(k,i)=ldompmp(k,i)/(dlx(i)*bh(k,i))
      LDOMpSS(K,I) = LDOMpAP(K,I)+LDOMpEP(K,I)+ldompmp(k,i)-(LDOMD(K,I)+LRDOMD(K,I))*orgpld(k,i)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   D O M   P H O S P H O R U S                            **
!***********************************************************************************************************************************

ENTRY REFRACTORY_DOM_P
  DO I=IU,ID
    DO K=KT,KB(I)
      RDOMpSS(K,I) = LRDOMD(K,I)*orgpld(k,i)-RDOMD(K,I)*orgprd(k,i)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      L A B I L E   P O M   P H O S P H O R U S                                **
!***********************************************************************************************************************************

ENTRY LABILE_POM_p
  LPOMpAP(:,IU:ID) = 0.0;lpompmp(:,iu:id)=0.0;lpzooinp(:,iu:id)=0.0; lpzoooutp(:,iu:id)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        LPOMpAP(K,I) = LPOMpAP(K,I)+APOM(JA)*(AMR(K,I,JA)*ALG(K,I,JA))*ap(ja)
      END DO
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          jt=k
          je=kb(i)
          do jj=jt,je
            lpompmp(k,i)=lpompmp(k,i)+mpom(m)*lrpmac(m)*mmr(k,i,m)*macrm(jj,k,i,m)*mp(m)
          end do
        end if
      end do            
      lpompmp(k,i)=lpompmp(k,i)/(dlx(i)*bh(k,i))
	do jz = 1,nzp
      if(tgraze(k,i,jz) > 0.0)then      
        lpzoooutp(k,i)=lpzoooutp(k,i) + zoo(k,i,jz)*(zmt(k,i,jz)+(zmu(k,i,jz)-(zmu(k,i,jz)*zeff(jz))))*zp(jz)
        lpzooinp(k,i)=lpzooinp(k,i) + zoo(k,i,jz)*zmu(k,i,jz)*prefp(jz)*LPOM(k,i)/tgraze(k,i,jz)*zp(jz)
      else
        lpzoooutp(k,i)=lpzoooutp(k,i)+zoo(k,i,jz)*(zmt(k,i,jz)+(zmu(k,i,jz)-(zmu(k,i,jz)*zeff(jz))))*zp(jz)
        lpzooinp(k,i)=0.0
      end if
    end do
      LPOMpNS(K,I) = POMS(JW)*(LPOM(K-1,I)*orgplp(k-1,i)-LPOM(K,I)*orgplp(k,i))*BI(K,I)/BH2(K,I)
      LPOMpSS(K,I) = LPOMpAP(K,I)+lpompmp(k,i)-LPOMD(K,I)*orgplp(k,i)+LPOMpNS(K,I)-LRPOMD(K,I)*orgplp(k,i)
	  do jz = 1,nzp
	   LPOMpSS(K,I) =LPOMpSS(K,I) + lpzoooutp(k,i)-lpzooinp(k,i)
	  end do

	END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   P O M   P H O S P H O R U S                            **
!***********************************************************************************************************************************

ENTRY REFRACTORY_POM_P
  rpompmp(:,iu:id)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)    
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          jt=k
          je=kb(i)
          tmac=0.0
          do jj=jt,je
            rpompmp(k,i)=rpompmp(k,i)+mpom(m)*(1.0-lrpmac(m))*mmr(k,i,m)*macrm(jj,k,i,m)*mp(m)
          end do
        end if
      end do
      rpompmp(k,i)=rpompmp(k,i)/(dlx(i)*bh(k,i))
      RPOMpNS(K,I) = POMS(JW)*(RPOM(K-1,I)*orgprp(k-1,i)-RPOM(K,I)*orgprp(k,i))*BI(K,I)/BH2(K,I)
      RPOMpSS(K,I) = LRPOMD(K,I)*orgplp(k,i)+RPOMpNS(K,I)-RPOMD(K,I)*orgprp(k,i)+rpompmp(k,i)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                       L A B I L E   D O M   N I T R O G E N                                   **
!***********************************************************************************************************************************

ENTRY LABILE_DOM_n
  LDOMnAP(:,IU:ID) = 0.0; LDOMnEP(:,IU:ID) = 0.0; ldomnmp(:,iu:id)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        LDOMnAP(K,I) = LDOMnAP(K,I)+(AER(K,I,JA)+(1.0-APOM(JA))*AMR(K,I,JA))*ALG(K,I,JA)*an(ja)
      END DO
      DO JE=1,NEP
        LDOMnEP(K,I) = LDOMnEP(K,I)+(EER(K,I,JE)+(1.0-EPOM(JE))*EMR(K,I,JE))*EPC(K,I,JE)*en(je)
      END DO
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          if(k.eq.kt)then
            jt=kti(i)
          else
            jt=k
          end if
          je=kb(i)
          do jj=jt,je
            ldomnmp(k,i)=ldomnmp(k,i)+(1.0-mpom(m))*mmr(k,i,m)*macrm(jj,k,i,m)*mn(m)
          end do
        end if
      end do
      ldomnmp(k,i)=ldomnmp(k,i)/(dlx(i)*bh(k,i))
      LDOMnSS(K,I) = LDOMnAP(K,I)+LDOMnEP(K,I)+ldomnmp(k,i)-(LDOMD(K,I)+LRDOMD(K,I))*orgnld(k,i)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   D O M   N I T R O G E N                                **
!***********************************************************************************************************************************

ENTRY REFRACTORY_DOM_n
  DO I=IU,ID
    DO K=KT,KB(I)
      RDOMnSS(K,I) = LRDOMD(K,I)*orgnld(k,i)-RDOMD(K,I)*orgnrd(k,i)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      L A B I L E   P O M   N I T R O G E N                                    **
!***********************************************************************************************************************************

ENTRY LABILE_POM_n
  LPOMnAP(:,IU:ID) = 0.0;lpomnmp(:,iu:id)=0.0;lpzooinn(:,iu:id)=0.0; lpzoooutn(:,iu:id)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        LPOMnAP(K,I) = LPOMnAP(K,I)+APOM(JA)*(AMR(K,I,JA)*ALG(K,I,JA))*an(ja)
      END DO
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          jt=k
          je=kb(i)
          do jj=jt,je
            lpomnmp(k,i)=lpomnmp(k,i)+mpom(m)*lrpmac(m)*mmr(k,i,m)*macrm(jj,k,i,m)*mn(m)
          end do
        end if
      end do            
      lpomnmp(k,i)=lpomnmp(k,i)/(dlx(i)*bh(k,i))
	do jz = 1,nzp
      if(tgraze(k,i,jz) > 0.0)then      
        lpzoooutn(k,i)=lpzoooutn(k,i)+zoo(k,i,jz)*(zmt(k,i,jz)+(zmu(k,i,jz)-(zmu(k,i,jz)*zeff(jz))))*zn(jz)
        lpzooinn(k,i)=lpzooinn(k,i)+zoo(k,i,jz)*prefp(jz)*zmu(k,i,jz)*LPOM(k,i)/tgraze(k,i,jz)*zn(jz)
      else
        lpzoooutn(k,i)=lpzoooutn(k,i)+zoo(k,i,jz)*(zmt(k,i,jz)+(zmu(k,i,jz)-(zmu(k,i,jz)*zeff(jz))))*zn(jz)
        lpzooinn(k,i)=0.0
      end if
	end do
      LPOMnNS(K,I) = POMS(JW)*(LPOM(K-1,I)*orgnlp(k-1,i)-LPOM(K,I)*orgnlp(k,i))*BI(K,I)/BH2(K,I)
      LPOMnSS(K,I) = LPOMnAP(K,I)+lpomnmp(k,i)-LPOMD(K,I)*orgnlp(k,i)+LPOMnNS(K,I)-LRPOMD(K,I)*orgnlp(k,i) &
            + lpzoooutn(k,i)-lpzooinn(k,i)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   P O M   N I T R O G E N                                **
!***********************************************************************************************************************************

ENTRY REFRACTORY_POM_n
  rpomnmp(:,iu:id)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          jt=k
          je=kb(i)
          tmac=0.0
          do jj=jt,je
            rpomnmp(k,i)=rpomnmp(k,i)+mpom(m)*(1.0-lrpmac(m))*mmr(k,i,m)*macrm(jj,k,i,m)*mn(m)
          end do
        end if
      end do
      rpomnmp(k,i)=rpomnmp(k,i)/(dlx(i)*bh(k,i))
      RPOMnNS(K,I) = POMS(JW)*(RPOM(K-1,I)*orgnrp(k-1,i)-RPOM(K,I)*orgnrp(k,i))*BI(K,I)/BH2(K,I)
      RPOMnSS(K,I) = LRPOMD(K,I)*orgnlp(k,i)+RPOMnNS(K,I)-RPOMD(K,I)*orgnrp(k,i)+rpomnmp(k,i)
    END DO
  END DO
RETURN


!************************************************************************
!**                          M A C R O P H Y T E                       **
!************************************************************************

ENTRY macrophyte(LLM)
  m=LLM
  DO I=IU,ID         
    wlelev=(el(kt,i)-z(i)*cos(alpha(jb)))    
    if(kticol(i))then
      jt=kti(i)
    else
      jt=kti(i)+1
    end if
    je=kb(i)          
    do jj=jt,je
      if(jj.lt.kt)then                
        colb=el(jj+1,i)
      else
        colb=el(kt+1,i)
      end if
      coldep=wlelev-colb
      if(macrc(jj,Kt,I,m).gt.mmax(m))then
        mgr(jj,kt,i,m)=0.0
      end if                        
      macss(jj,Kt,I,m) = (mGR(jj,Kt,I,m)-mMR(Kt,I,m)-mRR(Kt,I,m))*macrc(jj,Kt,I,m)     
      macrm(jj,kt,i,m)   = macrm(jj,kt,i,m)+macss(jj,kt,i,m)*dlt*coldep*cw(jj,i)*DLX(I)                                        
    end do      

    DO K=KT+1,KB(I)       
      jt=k            
      je=kb(i)        
      do jj=jt,je
        if(macrc(jj,K,I,m).gt.mmax(m))then
          mgr(jj,k,i,m)=0.0
        end if
        macss(jj,K,I,m) = (mGR(jj,K,I,m)-mMR(K,I,m)-mRR(K,I,m))*macrc(jj,K,I,m)                 
        if(mact(jj,k,i).gt.mbmp(m).and.mact(jj,k-1,i).lt.mbmp(m).and.macss(jj,k,i,m).gt.0.0)then     
          if(k-1.eq.kt)then
            bmass=macss(jj,k,i,m)*dlt*h2(k,i)*cw(jj,i)*DLX(I)
            macrm(jj,k-1,i,m)=macrm(jj,k-1,i,m)+bmass            
            colb=el(kt+1,i)      
            coldep=wlelev-colb              
            macss(jj,k-1,i,m)=bmass/dlt/(coldep*cw(jj,i)*DLX(I)) + macss(jj,k-1,i,m)                   
          else
            bmass=macss(jj,k,i,m)*dlt*h2(k,i)*cw(jj,i)*DLX(I)
            macrm(jj,k-1,i,m)=macrm(jj,k-1,i,m)+bmass                 
            macss(jj,k-1,i,m)=bmass/dlt/(h2(k-1,i)*cw(jj,i)*DLX(I))+ macss(jj,k-1,i,m)                       
          end if
          macss(jj,k,i,m)=0.0                
        else
          bmasstest=macrm(jj,k,i,m)+macss(jj,k,i,m)*dlt*H2(k,i)*cw(jj,i)*DLX(I)
          if(bmasstest.ge.0.0)then
            macrm(jj,k,i,m)   = bmasstest
          else
            macss(jj,k,i,m)=-macrm(jj,k,i,m)/dlt/(H2(k,i)*cw(jj,i)*DLX(I))
            macrm(jj,k,i,m)=0.0
          end if                
        end if                      
      end do                                
    END DO
  END DO
  DO I=IU,ID
    tmac=0.0
    cvol=0.0          
    if(kticol(i))then
      jt=kti(i)
    else
      jt=kti(i)+1
    end if
    je=kb(i)
    wlelev=(el(kt,i)-z(i)*cos(alpha(jb)))
         
    do jj=jt,je
      if(jj.lt.kt)then                
        colb=el(jj+1,i)
      else
        colb=el(kt+1,i)
      end if
      coldep=wlelev-colb
      if(cw(jj,i).gt.0.0)then
        macrc(jj,kt,i,m)=macrm(jj,kt,i,m)/(cw(jj,i)*coldep*dlx(i))
      else
        macrc(jj,kt,i,m)=0.0
      end if
      tmac=tmac+macrm(jj,kt,i,m)
      cvol=cvol+cw(jj,i)*coldep*dlx(i)
    end do
        
    mac(kt,i,m)=tmac/cvol
                
    DO K=KT+1,KB(I)            
      jt=k              
      je=kb(i)
      tmac=0.0
      cvol=0.0
      do jj=jt,je
        if(cw(jj,i).gt.0.0)then
          macrc(jj,k,i,m)=macrm(jj,k,i,m)/(cw(jj,i)*h2(k,i)*dlx(i))
        else
          macrc(jj,k,i,m)=0.0
        end if
        tmac=tmac+macrm(jj,k,i,m)
        cvol=cvol+cw(jj,i)*h2(k,i)*dlx(i)
      end do
      mac(k,i,m)=tmac/cvol
    end do
  end do

  DO I=IU,ID
    tmac=0.0
    cvol=0.0
    do k=kt,kb(i)
      if(k.eq.kt)then
        jt=kti(i)
      else
        jt=k
      end if
      je=kb(i)
      do jj=jt,je
        mact(jj,k,i)=0.0
        do mi=1,nmc
          if(macrophyte_calc(jw,mi))then
            mact(jj,k,i)=macrc(jj,k,i,mi)+mact(jj,k,i)
          end if
        end do
      end do
    end do
  end do
      

  RETURN

! v3.5 end

!***********************************************************************************************************************************
!*                                                  K I N E T I C   F L U X E S                                                   **
!***********************************************************************************************************************************

ENTRY KINETIC_FLUXES
  DO JAF=1,NAF(JW)
    DO I=CUS(BS(JW)),DS(BE(JW))
      DO K=KT,KB(I)
        KFS(K,I,KFCN(JAF,JW)) = KFS(K,I,KFCN(JAF,JW))+KF(K,I,KFCN(JAF,JW))*VOL(K,I)*DLT
      END DO
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                       p H   C O 2                                                             **
!***********************************************************************************************************************************

ENTRY PH_CO2

! pH and carbonate species

  DO I=IU,ID
    DO K=KT,KB(I)
      CART = TIC(K,I)/12000.0
      ALKT = ALK(K,I)/5.0E+04
      T1K  = T1(K,I)+273.15

!**** Ionic strength

      IF (FRESH_WATER(JW)) S2 = 2.5E-05*TDS(K,I)
      IF (SALT_WATER(JW))  S2 = 1.47E-3+1.9885E-2*TDS(K,I)+3.8E-5*TDS(K,I)*TDS(K,I)

!**** Debye-Huckel terms and activity coefficients

      SQRS2  =  SQRT(S2)
      DH1    = -0.5085*SQRS2/(1.0+1.3124*SQRS2)+4.745694E-03+4.160762E-02*S2-9.284843E-03*S2*S2
      DH2    = -2.0340*SQRS2/(1.0+1.4765*SQRS2)+1.205665E-02+9.715745E-02*S2-2.067746E-02*S2*S2
      H2CO3T =  10.0**(0.0755*S2)
      HCO3T  =  10.0**DH1
      CO3T   =  10.0**DH2
      OH     =  HCO3T

!**** Temperature adjustment

      KW = 10.0**(-283.971-0.05069842*T1K+13323.0/T1K+102.24447*LOG10(T1K)-1119669.0/(T1K*T1K))/OH
      K1 = 10.0**(-3404.71/T1K+14.8435-0.032786*T1K)*H2CO3T/HCO3T
      K2 = 10.0**(-2902.39/T1K+ 6.4980-0.023790*T1K)*HCO3T/CO3T

!**** pH evaluation

      PHT = -PH(K,I)-2.1
      IF (PH(K,I) <= 0.0) PHT = -14.0
      INCR = 10.0
      DO N=1,3
        F    = 1.0
        INCR = INCR/10.0
        ITER = 0
        DO WHILE (F > 0.0 .AND. ITER < 12)
          PHT    = PHT+INCR
          HION   = 10.0**PHT
          BICART = CART*K1*HION/(K1*HION+K1*K2+HION*HION)
          F      = BICART*(HION+2.0*K2)/HION+KW/HION-ALKT-HION/OH
          ITER   = ITER+1
        END DO
        PHT = PHT-INCR
      END DO

!**** pH, carbon dioxide, bicarbonate, and carbonate concentrations

      HION      =  10.0**PHT
      PH(K,I)   = -PHT
      CO2(K,I)  =  TIC(K,I)/(1.0+K1/HION+K1*K2/(HION*HION))
      HCO3(K,I) =  TIC(K,I)/(1.0+HION/K1+K2/HION)
      CO3(K,I)  =  TIC(K,I)/((HION*HION)/(K1*K2)+HION/K2+1.0)
    END DO
  END DO
RETURN

! v3.5 start
!**********************************************************
!**           SUBROUTINE ZOOPLANKTON                     **
!**********************************************************

ENTRY zooplankton
  DO I=IU,ID
    DO K=KT,KB(I) 
	  do jz = 1, nzp
            zgztot=0.0                                                                                                   ! kv 5/9/2007
	        do jjz = 1,nzp
!             zgztot=zgztot+zgz(k,i,jz,jjz)*zoo(k,i,jz)                                                                   ! kv 5/9/2007          
            zgztot=zgztot+zgz(k,i,jz,jjz)                                                                             ! cb 5/26/07
            end do
        zooss(k,i,jz)= (zmu(k,i,jz)*zeff(jz)-zrt(k,i,jz)-zmt(k,i,jz))*zoo(k,i,jz) - zgztot   ! omnivorous zooplankton    ! kv 5/9/2007
	  end do
    END DO
  END DO
RETURN
! v3.5 end

!***********************************************************************************************************************************
!**                                              D E R I V E D   C O N S T I T U E N T S                                          **
!***********************************************************************************************************************************

ENTRY DERIVED_CONSTITUENTS
  APR = 0.0; ATOT = 0.0; TOTSS = 0.0; CHLA = 0.0
  DO JW=1,NWB
    KT = KTWB(JW)
    DO JB=BS(JW),BE(JW)
      DO I=CUS(JB),DS(JB)
        DO K=KT,KB(I)
          DO JA=1,NAL
            APR(K,I) = APR(K,I)+(AGR(K,I,JA)-ARR(K,I,JA))*ALG(K,I,JA)*H2(K,I)*DAY
          END DO
        END DO
        DO K=KT,KB(I)
          CBODC = 0.0; CBODN = 0.0; CBODP = 0.0; BODTOT = 0.0; ALGP = 0.0; ALGN = 0.0
          DO JA=1,NAL
            ATOT(K,I) = ATOT(K,I)+ALG(K,I,JA)
          END DO
          DO IBOD=1,NBOD   
            CBODC  = CBODC+CBOD(K,I,IBOD)*BODC(IBOD)
            CBODN  = CBODN+CBOD(K,I,IBOD)*BODN(IBOD)
            CBODP  = CBODP+CBOD(K,I,IBOD)*BODP(IBOD)
            BODTOT = BODTOT+CBOD(K,I,IBOD)
          END DO
          DOM(K,I) = LDOM(K,I)+RDOM(K,I)
          POM(K,I) = LPOM(K,I)+RPOM(K,I)
          DOC(K,I) = DOM(K,I)*ORGC(JW)+CBODC
          POC(K,I) = POM(K,I)*ORGC(JW)
          DO JA=1,NAL
            POC(K,I) = POC(K,I)+ALG(K,I,JA)*AC(JA)
            ALGP     = ALGP+ALG(K,I,JA)*AP(JA)
            ALGN     = ALGN+ALG(K,I,JA)*AN(JA)
          END DO
! v3.5 start
          do jz=1,nzp
            poc(k,i)=poc(k,i)+zc(jz)*zoo(k,i,jz) !mlm baulk
            zoop=zoo(k,i,jz)*zp(jz) !mlm baulk
            zoon=zoo(k,i,jz)*zn(jz) !mlm baulk
	      end do
          TOC(K,I)   = DOC(K,I)+POC(K,I)
          DOP(K,I)   = LDOM(K,I)*ORGPLD(k,i)+RDOM(k,i)*orgprd(k,i)+CBODP
          DON(K,I)   = LDOM(K,I)*ORGNLD(k,i)+RDOM(k,i)*orgnrd(k,i)+CBODN
          POP(K,I)   = LPOM(K,I)*ORGPLP(k,i)+RPOM(k,i)*orgprp(k,i)+ALGP+zoop
          PON(K,I)   = LPOM(K,I)*ORGNLP(k,i)+RPOM(k,i)*orgnrp(k,i)+ALGN+zoop
! v3.5 end           
          TOP(K,I)   = DOP(K,I)+POP(K,I)
          TON(K,I)   = DON(K,I)+PON(K,I)
          TKN(K,I)   = TON(K,I)+NH4(K,I)
! v3.5 start
          CBODU(K,I) = O2OM(JW)*(DOM(K,I)+POM(K,I)+ATOT(K,I))+BODTOT 
		  do jz = 1,nzp
		    cbodu(k,i) = cbodu(k,i) + o2om(jw)*zoo(k,i,jz)
		  end do
! v3.5 end
          TPSS       = 0.0
          DO JS=1,NSS    
            TPSS = TPSS+SS(K,I,JS)*PARTP(JW)
          END DO
          TP(K,I)   =  TOP(K,I)+PO4(K,I)+TPSS
          TN(K,I)   =  TON(K,I)+NH4(K,I)+NO3(K,I)
          O2DG(K,I) = (O2(K,I)/SATO(T1(K,I),TDS(K,I),PALT(I),SALT_WATER(JW)))*100.0
          DO JA=1,NAL
            CHLA(K,I)  = CHLA(K,I) +ALG(K,I,JA)/ACHLA(JA)
            TOTSS(K,I) = TOTSS(K,I)+ALG(K,I,JA)
          END DO
          TOTSS(K,I) = TOTSS(K,I)+TISS(K,I)+POM(K,I)
        END DO
      END DO
    END DO
  END DO
RETURN
ENTRY DEALLOCATE_KINETICS
  DEALLOCATE (OMTRM,  SODTRM, NH4TRM, NO3TRM, DOM, POM, PO4BOD, NH4BOD, TICBOD, ATRM,   ATRMR,  ATRMF, ETRM,   ETRMR,  ETRMF)
! v3.5 start
  deALLOCATE (lam2m)  
! v3.5 end
RETURN
END SUBROUTINE KINETICS

!***********************************************************************************************************************************
!**                                          S U B R O U T I N E   G A S   T R A N S F E R                                        **
!***********************************************************************************************************************************

SUBROUTINE GAS_TRANSFER
  USE GLOBAL; USE GEOMC; USE KINETIC
  REAL, PARAMETER :: THETA_REAERATION = 1.024, M_TO_FT = 3.2808

  IF (REAERC(JW) == '   RIVER') THEN

!** Average depth in ft

    AREA = 0.0
    DO K=KT,KBMIN(I)
      AREA = AREA+BHR1(K,I)
    END DO
    ADEPTH = AREA/BR(KTI(I),I)*M_TO_FT

!** Average velocity in feet/second

    UAVG = ABS(QC(I))/AREA*M_TO_FT

!** Reaeration factor

    IF (NEQN(JW) == 0) THEN
      IF (ADEPTH <= 2.0) THEN
        REAER(I) = 21.64*UAVG**0.67/ADEPTH**1.85
      ELSE IF (UAVG <= 1.8) THEN
        REAER(I) = 12.96*SQRT(UAVG)/ADEPTH**1.5
      ELSE
        HDEPTH = -11.875*UAVG+23.375
        IF (HDEPTH >= ADEPTH) THEN
          REAER(I) = 12.96*SQRT(UAVG)/ADEPTH**1.5
        ELSE
          REAER(I) = 11.57*UAVG**0.969/ADEPTH**1.673
        END IF
      END IF
    ELSE IF (NEQN(JW) == 1) THEN                                                                              !O'connor-Dobbins
      REAER(I) = 12.96*SQRT(UAVG)/ADEPTH**1.5
    ELSE IF (NEQN(JW) == 2) THEN                                                                              !Churchill
      REAER(I) = 11.57*UAVG**0.969/ADEPTH**1.673
    ELSE IF (NEQN(JW) == 3) THEN                                                                              !Tsivoglou
      S = SLOPE(JB)*5280.0
      IF (ABS(QC(I))*35.5 >= 10.0) THEN
        REAER(I) = 0.88*S*UAVG
      ELSE
        REAER(I) = 1.8*S*UAVG
      END IF
    ELSE IF (NEQN(JW) == 4) THEN                                                                              !Owens
      REAER(I) = 21.64*UAVG**0.67/ADEPTH**1.85
    ELSE IF (NEQN(JW) == 5) THEN                                                                              !Thackston and Krenkel
      USTAR    = SQRT(ADEPTH*SLOPE(JB)*32.2)                                                                  ! SR 5/10/05
      REAER(I) = 24.88*(1.0+SQRT(0.176*UAVG/SQRT(ADEPTH)))*USTAR/ADEPTH                                       ! SR 5/10/05
    ELSE IF (NEQN(JW) == 6) THEN                                                                              !Langbien and Durum
      REAER(I) = 7.60*UAVG/ADEPTH**1.33
    ELSE IF (NEQN(JW) == 7) THEN                                                                              !Melching and Flores 
      UAVG = UAVG/M_TO_FT
      IF (QC(I) == 0.0) THEN   
        REAER(I) = 0.0      
      ELSE IF (ABS(QC(I)) < 0.556) THEN   
        REAER(I) = 517.0*((UAVG*SLOPE(JB))**0.524)*ABS(QC(I))**(-0.242)
      ELSE
        REAER(I) = 596.0*((UAVG*SLOPE(JB))**0.528)*ABS(QC(I))**(-0.136)
      END IF
    ELSE IF (NEQN(JW) == 8) THEN                                                                              !Melching and Flores
      UAVG   = UAVG/M_TO_FT
      ADEPTH = ADEPTH/M_TO_FT
      IF (ABS(QC(I)) < 0.556) THEN
        REAER(I) = 88.0*((UAVG*SLOPE(JB))**0.313)*ADEPTH**(-0.353)
      ELSE
        REAER(I) = 142.0*((UAVG*SLOPE(JB))**0.333)*ADEPTH**(-0.66)*BI(KT,I)**(-0.243)
      END IF
    ELSE IF (NEQN(JW) == 9) THEN                                                                              !User defined SI units
      UAVG     = UAVG/M_TO_FT
      ADEPTH   = ADEPTH/M_TO_FT
      REAER(I) = RCOEF1(JW)*(UAVG**RCOEF2(JW))*(ADEPTH**RCOEF3(JW))*(SLOPE(JB)**RCOEF4(JW))
    ELSE IF (NEQN(JW) == 10) THEN                                                                             ! Thackston and Krenkel - updated
      USTAR    = SQRT(ADEPTH*SLOPE(JB)*32.2)                                                                  ! SR 5/10/05     
      REAER(I) = 4.99*(1.0+9.0*(0.176*UAVG/SQRT(ADEPTH))**0.25)*USTAR/ADEPTH                                  ! SR 5/10/05  
    END IF
    REAER(I) = REAER(I)*ADEPTH/M_TO_FT
  ELSE IF (REAERC(JW) == '    LAKE') THEN
    IF (NEQN(JW) == 1) THEN                                                                                   !Broecker
      REAER(I) = 0.864*WIND10(I)
    ELSE IF (NEQN(JW) == 2) THEN
      IF (WIND10(I) <= 3.5) THEN                                                                              !Gelda
        A     = 0.2
        BCOEF = 1.0
      ELSE
        A     = 0.057
        BCOEF = 2.0
      END IF
      REAER(I) = A*WIND10(I)**BCOEF
    ELSE IF (NEQN(JW) == 3) THEN                                                                              !Banks & Herrera
      REAER(I) = (0.728*SQRT(WIND10(I))-0.317*WIND10(I)+0.0372*WIND10(I)**2)
    ELSE IF (NEQN(JW) == 4) THEN                                                                              !Wanninkhof
      REAER(I) = 0.0986*WIND10(I)**1.64
    ELSE IF (NEQN(JW) == 5) THEN                                                                              !Chen & Kanwisher
      DMO2     = 2.04E-9
      REAER(I) = DAY*DMO2/((200.0-60.0*SQRT(MIN(WIND10(I),11.0)))*1.E-6)
    ELSE IF (NEQN(JW) == 6) THEN                                                                              !Cole & Buchak
      REAER(I) = (0.5+0.05*WIND10(I)*WIND10(I))
    ELSE IF (NEQN(JW) == 7) THEN                                                                              !Banks
      IF (WIND10(I) <= 5.5) THEN
        REAER(I) = 0.362*SQRT(WIND10(I))
      ELSE
        REAER(I) = 0.0277*WIND10(I)**2
      END IF
    ELSE IF (NEQN(JW) == 8) THEN                                                                              !Smith
      REAER(I) = 0.64+0.128*WIND10(I)**2
    ELSE IF (NEQN(JW) == 9) THEN                                                                              !Liss
      IF (WIND10(I) <= 4.1) THEN
        REAER(I) = 0.156*WIND10(I)**0.63
      ELSE
        REAER(I) = 0.0269*WIND10(I)**1.9
      END IF
    ELSE IF (NEQN(JW) == 10) THEN                                                                             !Downing and Truesdale
      REAER(I) = 0.0276*WIND10(I)**2
    ELSE IF (NEQN(JW) == 11) THEN                                                                             !Kanwisher
      REAER(I) = 0.0432*WIND10(I)**2
    ELSE IF (NEQN(JW) == 12) THEN                                                                             !Yu, et al
      REAER(I) = 0.319*WIND10(I)
    ELSE IF (NEQN(JW) == 13) THEN                                                                             !Weiler
      IF (WIND10(I) <= 1.6) THEN
        REAER(I) = 0.398
      ELSE
        REAER(I) = 0.155*WIND10(I)**2
      END IF
    ELSE IF (NEQN(JW) == 14) THEN                                                                             !User defined
      REAER(I) = RCOEF1(JW)+RCOEF2(JW)*WIND10(I)**RCOEF3(JW)
    END IF
  ELSE IF (REAERC(JW) == ' ESTUARY') THEN
    AREA = 0.0
    DO K=KT,KBMIN(I)
      AREA = AREA+BHR1(K,I)
    END DO
    ADEPTH = AREA/BR(KTI(I),I)*M_TO_FT
    UAVG   = ABS(QC(I))/AREA*M_TO_FT

!** Reaeration factor

    IF (NEQN(JW) == 0) THEN
      IF (ADEPTH <= 2.0) THEN
        REAER(I) = 21.64*UAVG**0.67/ADEPTH**1.85
      ELSE IF (UAVG <= 1.8) THEN
        REAER(I) = 12.96*SQRT(UAVG)/ADEPTH**1.5
      ELSE
        HDEPTH = -11.875*UAVG+23.375
        IF (HDEPTH >= ADEPTH) THEN
          REAER(I) = 12.96*SQRT(UAVG)/ADEPTH**1.5
        ELSE
          REAER(I) = 11.57*UAVG**0.969/ADEPTH**1.673
        END IF
      END IF
    ELSE IF (NEQN(JW) == 1) THEN                                                                           !Thomann and Fitzpatrick
      REAER(I) = (0.728*SQRT(WIND10(I))-0.317*WIND10(I)+0.0372*WIND10(I)**2)+3.93*SQRT(UAVG/M_TO_FT)/(ADEPTH/M_TO_FT)**0.5
    END IF
  END IF
  IF (REAER(I) <= 0.6) REAER(I) = 0.6
  REAER(I) = REAER(I)*THETA_REAERATION**(T1(KT,I)-20.0)
  REAER(I) = REAER(I)/DAY
END SUBROUTINE GAS_TRANSFER

!***********************************************************************************************************************************
!**                                       S U B R O U T I N E   T O T A L  D I S S O L V E D  G A S                               **
!***********************************************************************************************************************************

SUBROUTINE TOTAL_DISSOLVED_GAS (P,NSG,N,T,C)
  USE TDGAS; USE STRUCTURES; USE GLOBAL

  SAT = EXP(7.7117-1.31403*(LOG(T+45.93)))*P
  IF (NSG == 0) THEN
    IF (EQSP(N) == 1) THEN
      TDG = AGASSP(N)*.035313*QSP(N)+BGASSP(N)
      IF (TDG > 145.0) TDG = 145.0
      C = SAT    
      IF (TDG >= 100.0) C = TDG*SAT/100.0
    ELSE IF (EQSP(N) == 2) THEN
      TDG = AGASSP(N)+BGASSP(N)*EXP(0.035313*QSP(N)*CGASSP(N))
      IF (TDG > 145.0)  TDG = 145.0
      C = SAT    
      IF (TDG >= 100.0) C = TDG*SAT/100.0
    ELSE
      DB = SAT-C
      DA = DB*(1.0+0.38*AGASSP(N)*BGASSP(N)*CGASSP(N)*(1.0-0.11*CGASSP(N))*(1.0+0.046*T))
      C  = SAT-DA
    END IF
  ELSE IF (EQGT(N) == 1) THEN
    TDG = AGASGT(N)*0.035313*QGT(N)+BGASGT(N)
    IF (TDG > 145.0) TDG = 145.0
    C = SAT    
    IF (TDG >= 100.0) C = TDG*SAT/100.0
  ELSE IF (EQGT(N) == 2) THEN
    TDG = AGASGT(N)+BGASGT(N)*EXP(.035313*QGT(N)*CGASGT(N))
    IF (TDG > 145.0) TDG = 145.0
    C = SAT    
    IF (TDG >= 100.0) C = TDG*SAT/100.0
  ELSE
    DB = SAT-C
    DA = DB*(1.0+0.38*AGASGT(N)*BGASGT(N)*CGASGT(N)*(1.0-0.11*CGASGT(N))*(1.0+0.046*T))
    C  = SAT-DA
  END IF
END SUBROUTINE TOTAL_DISSOLVED_GAS

!***********************************************************************************************************************************
!**                                           S U B R O U T I N E   O U T P U T                                                   **
!***********************************************************************************************************************************

SUBROUTINE OUTPUT (JDAY,IUPR,IDPR,KBR,ISNP,BL,NBL)
  USE GLOBAL; USE GDAYC;  USE GEOMC;  USE KINETIC; USE TVDC; USE NAMESC; USE LOGICC
  use macrophytec !v3.5

! Type declaration

  REAL                        :: JDAY, LIMIT
  INTEGER, DIMENSION(IMX)     :: BL
  INTEGER, DIMENSION(IMX,NWB) :: ISNP
  LOGICAL                     :: NEW_PAGE
  CHARACTER(8)                :: LFAC
  CHARACTER(10)               :: BLANK

! Data declaration

  DATA BLANK /'          '/

! Variable initialization

  NEW_PAGE = .TRUE.

! Blank inactive cells

  NBL  = 1
  JB   = 1
  IUPR = 1
  DO I=1,IDPR-1
    IF (CUS(JB) > ISNP(I,JW)) THEN
      BL(NBL) = I
      NBL     = NBL+1
      IF (JB == 1) IUPR = I+1
    END IF
    IF (ISNP(I+1,JW) > DS(JB)) JB = JB+1
  END DO
  NBL = NBL-1

! Water surface elevation, water surface deviation, ice cover, and sediment oxygen demand

  CONV(1,:) = BLANK
  DO I=IUPR,IDPR
    DO JJB=1,NBR
      IF (ISNP(I,JW) >= US(JJB)-1 .AND. ISNP(I,JW) <= DS(JJB)+1) EXIT
    END DO
    WRITE (CONV(1,I),'(F10.3)') EL(KTWB(JW),ISNP(I,JW))-Z(ISNP(I,JW))*COSA(JJB)
  END DO
  WRITE (SNP(JW),'(/A//2X,1000I10)') '          Water Surface, m',(ISNP(I,JW),I=IUPR,IDPR)
  WRITE (SNP(JW),'(2X,1000A10/)') (CONV(1,I),I=IUPR,IDPR)
  DO I=IUPR,IDPR
    WRITE (CONV(1,I),'(F10.4)') SNGL(Z(ISNP(I,JW)))
  END DO
  WRITE (SNP(JW),'(/A//2X,1000I10)') '          Water Surface Deviation (positive downwards), m',(ISNP(I,JW),I=IUPR,IDPR)
  WRITE (SNP(JW),'(2X,1000A10/)') (CONV(1,I),I=IUPR,IDPR)
  IF (ICE_CALC(JW)) THEN
    DO I=IUPR,IDPR
      WRITE (CONV(1,I),'(F10.3)') ICETH(ISNP(I,JW))
    END DO
    WRITE (SNP(JW),'(/A//3X,1000A10)') '          Ice Thickness, m',(CONV(1,I),I=IUPR,IDPR)
  END IF
  IF (CONSTITUENTS) THEN
    DO I=IUPR,IDPR
      WRITE (CONV(1,I),'(F10.3)') SOD(ISNP(I,JW))*DAY
    END DO
    IF (OXYGEN_DEMAND) THEN
      WRITE (SNP(JW),'(/A//3X,1000A10/)') '          Sediment Oxygen Demand, g/m^2/day',(CONV(1,I),I=IUPR,IDPR)
    END IF
  END IF

! Hydrodynamic variables and temperature

  CONV = BLANK
  DO JH=1,NHY
    L = LEN_TRIM(FMTH(JH))
    IF (PRINT_HYDRO(JH,JW)) THEN
      DO I=IUPR,IDPR
        IF (JH == 1) THEN
          DO K=KTWB(JW),KB(ISNP(I,JW))
            WRITE (CONV(K,I),FMTH(JH)) INT(HYD(K,ISNP(I,JW),JH))
          END DO
        ELSE IF (JH > 6) THEN
          DO K=KTWB(JW),KB(ISNP(I,JW))
            WRITE (CONV(K,I),FMTH(JH)(1:L)) HYD(K,ISNP(I,JW),JH)*DLT
          END DO
        ELSE
          DO K=KTWB(JW),KB(ISNP(I,JW))
            WRITE (CONV(K,I),FMTH(JH)(1:L)) HYD(K,ISNP(I,JW),JH)*HMULT(JH)
          END DO
        END IF
      END DO
      IF (NEW_PAGE) THEN
        WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
        NLINES = KMX-KTWB(JW)+14
      END IF
      NLINES   = NLINES+KMX-KTWB(JW)+11
      NEW_PAGE = NLINES > 72
      WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian day = ',INT(JDAY),' days ',                     &
                                                  (JDAY-INT(JDAY))*24.0,' hours   '//HNAME(JH)
      WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
      DO K=KTWB(JW),KBR
        WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)
      END DO
    END IF
  END DO

! Constituent concentrations

  IF (CONSTITUENTS) THEN
    DO JAC=1,NAC
      JC = CN(JAC)
      L  = LEN_TRIM(FMTC(JC))
      IF (PRINT_CONST(JC,JW)) THEN
        DO I=IUPR,IDPR
          DO K=KTWB(JW),KB(ISNP(I,JW))
            WRITE (CONV(K,I),FMTC(JC)(1:L)) C2(K,ISNP(I,JW),JC)*CMULT(JC)
          END DO
        END DO
        IF (NEW_PAGE) THEN
          WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
          NLINES = KMX-KTWB(JW)+14
        END IF
        NLINES   = NLINES+KMX-KTWB(JW)+11
        NEW_PAGE = NLINES > 72
        WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))    &
                                                     *24.0,' hours   '//CNAME(JC)
        WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
        DO K=KTWB(JW),KBR
          WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)
        END DO
      END IF
    END DO

!** Derived constituent concentrations

    DO JD=1,NDC
      L = LEN_TRIM(FMTCD(JD))
      IF (PRINT_DERIVED(JD,JW)) THEN
        DO I=IUPR,IDPR
          DO K=KTWB(JW),KB(ISNP(I,JW))
            WRITE (CONV(K,I),FMTCD(JD)(1:L)) CD(K,ISNP(I,JW),JD)*CDMULT(JD)
          END DO
        END DO
        IF (NEW_PAGE) THEN
          WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
          NLINES = KMX-KTWB(JW)+14
        END IF
        NLINES   = NLINES+KMX-KTWB(JW)+11
        NEW_PAGE = NLINES > 72
        WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))    &
                                                     *24.0,' hours    '//CDNAME(JD)
        WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
        DO K=KTWB(JW),KBR
          WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)
        END DO
      END IF
    END DO

!** Sediment

    IF (PRINT_SEDIMENT(JW)) THEN
      DO I=IUPR,IDPR
        DO K=KTWB(JW),KB(ISNP(I,JW))
          WRITE (CONV(K,I),'(F10.2)') SED(K,ISNP(I,JW))
        END DO
      END DO
      IF (NEW_PAGE) THEN
        WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
        NLINES = KMX-KTWB(JW)+14
      END IF
      NLINES   = NLINES+KMX-KTWB(JW)+11
      NEW_PAGE = NLINES > 72
      WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))*24.0,&
                                                 ' hours     Organic sediments, g/m^3'
      WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
      DO K=KTWB(JW),KBR
        WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)
      END DO
    END IF

! v3.5 start
    IF (PRINT_SEDIMENT(JW)) THEN
      DO I=IUPR,IDPR
        DO K=KTWB(JW),KB(ISNP(I,JW))
          WRITE (CONV(K,I),'(F10.2)') SEDp(K,ISNP(I,JW))
        END DO
      END DO
      IF (NEW_PAGE) THEN
        WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
        NLINES = KMX-KTWB(JW)+14
      END IF
      NLINES   = NLINES+KMX-KTWB(JW)+11
      NEW_PAGE = NLINES > 72
      WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))*24.0,&
                                                 ' hours     Organic phosphorus sediments, g/m^3'
      WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
      DO K=KTWB(JW),KBR
        WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)
      END DO
    END IF

    IF (PRINT_SEDIMENT(JW)) THEN
      DO I=IUPR,IDPR
        DO K=KTWB(JW),KB(ISNP(I,JW))
          WRITE (CONV(K,I),'(F10.2)') SEDn(K,ISNP(I,JW))
        END DO
      END DO
      IF (NEW_PAGE) THEN
        WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
        NLINES = KMX-KTWB(JW)+14
      END IF
      NLINES   = NLINES+KMX-KTWB(JW)+11
      NEW_PAGE = NLINES > 72
      WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))*24.0,&
                                                 ' hours     Organic nitrogen sediments, g/m^3'
      WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
      DO K=KTWB(JW),KBR
        WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)
      END DO
    END IF

    IF (PRINT_SEDIMENT(JW)) THEN
      DO I=IUPR,IDPR
        DO K=KTWB(JW),KB(ISNP(I,JW))
          WRITE (CONV(K,I),'(F10.2)') SEDc(K,ISNP(I,JW))
        END DO
      END DO
      IF (NEW_PAGE) THEN
        WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
        NLINES = KMX-KTWB(JW)+14
      END IF
      NLINES   = NLINES+KMX-KTWB(JW)+11
      NEW_PAGE = NLINES > 72
      WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))*24.0,&
                                                 ' hours     Organic carbon sediments, g/m^3'
      WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
      DO K=KTWB(JW),KBR
        WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)
      END DO
    END IF

! v3.5 end


!** Epiphyton

    DO JE=1,NEP
      IF (PRINT_EPIPHYTON(JW,JE)) THEN
        DO I=IUPR,IDPR
          DO K=KTWB(JW),KB(ISNP(I,JW))
            WRITE (CONV(K,I),'(F10.2)') EPD(K,ISNP(I,JW),JE)
          END DO
        END DO
        IF (NEW_PAGE) THEN
          WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
          NLINES = KMX-KTWB(JW)+14
        END IF
        NLINES   = NLINES+KMX-KTWB(JW)+11
        NEW_PAGE = NLINES > 72
        WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',                    &
                                                    (JDAY-INT(JDAY))*24.0,' hours     Epiphyton, g/m^2'
        WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
        DO K=KTWB(JW),KBR
          WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(ADJUSTR(CONV(K,I)),I=IUPR,IDPR)
        END DO
      END IF
    END DO

! v3.5 start

!********* macrophytes
    do L=1,nmc
      CONV = BLANK
      IF (PRINT_macrophyte(jw,L)) THEN
        DO I=IUPR,IDPR
          DO K=ktwb(jw),KB(ISNP(I,jw))
            WRITE (CONV(K,I),'(F10.2)') mac(K,ISNP(I,jw),L)
          END DO
        END DO        
        IF (NEW_PAGE) THEN
          WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
          NLINES = KMX-KTWB(JW)+14
        END IF
        NLINES   = NLINES+KMX-KTWB(JW)+11
        NEW_PAGE = NLINES > 72
        WRITE (SNP(jw),'(/1X,3(A,1X,I0),A,F0.2,A/)')  MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',  &
                                                      (JDAY-INT(JDAY))*24.0,' hours   ','  Macrophyte, g/m^3'          
        WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
        DO K=ktwb(jw),KBR
          WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)     
        END DO
        do i=iupr,idpr
          conv2=blank
          DO K=ktwb(jw),KB(ISNP(I,jw))
            if(k.eq.ktwb(jw))then
              jt=kti(isnp(i,jw))
            else
              jt=k
            end if
            je=kb(isnp(i,jw))
            DO jj=jt,je              
              WRITE (CONV2(K,jj),'(F10.2)') macrc(jj,K,ISNP(I,jw),L)
            END DO
          END DO
          IF (NEW_PAGE) THEN
            WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
            NLINES = KMX-KTWB(JW)+14
          END IF
          NLINES   = NLINES+KMX-KTWB(JW)+11
          NEW_PAGE = NLINES > 72
          WRITE (SNP(jw),'(/1X,3(A,1X,I0),A,F0.2,A/)')  MONTH,GDAY,',',YEAR,'    Julian Date ', &
                                                        INT(JDAY),' days ',(JDAY-INT(JDAY))     &
                                                     *24.0,' hours   ','  Macrophyte Columns g/m^3'                                
          write (SNP(jw),3052)isnp(i,jw)
3052      format(7x,"Segment   ",i8)
          jt=kti(isnp(i,jw))
          je=kb(isnp(i,jw))
          WRITE (SNP(jw),'(A,200I10)') ' Layer  Depth',(jj,jj=jt,je)
          DO K=ktwb(jw),KB(isnp(i,jw))                                    
            WRITE (SNP(jw),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV2(K,jj),jj=jt,je)     
          END DO
          write (SNP(jw),3053)isnp(i,jw)
3053      format(7x,"Macrophyte Limitation Segment   ",i8)            
          WRITE (SNP(jw),'(A,200I10)')' Layer  Depth', (jj,jj=jt,je)
          mlfpr=blank
          DO K=ktwb(jw),KB(isnp(i,jw))
            if(k.eq.ktwb(jw))then
              jt=kti(isnp(i,jw))
            else
              jt=k
            end if
            do jj=jt,je
              LIMIT = MIN(mPLIM(K,isnp(i,jw),L),mNLIM(K,isnp(i,jw),L),mcLIM(K,isnp(i,jw),L),mLLIM(jj,K,isnp(i,jw),L))
              IF (LIMIT == mPLIM(K,isnp(i,jw),L)) THEN
                WRITE (LFAC,'(F8.4)') mPLIM(K,isnp(i,jw),L)
                mLFPR(jj,K,isnp(i,jw),L) = ' P'//LFAC
              ELSE IF (LIMIT == mNLIM(K,I,L)) THEN
                WRITE (LFAC,'(F8.4)') mNLIM(K,isnp(i,jw),L)
                mLFPR(jj,K,isnp(i,jw),L) = ' N'//LFAC
              ELSE IF (LIMIT == mcLIM(K,I,L)) THEN
                WRITE (LFAC,'(F8.4)') mcLIM(K,isnp(i,jw),L)
                mLFPR(jj,K,isnp(i,jw),L) = ' C'//LFAC
              ELSE IF (LIMIT == mLLIM(jj,K,isnp(i,jw),L)) THEN
                WRITE (LFAC,'(F8.4)') mLLIM(jj,K,isnp(i,jw),L)
                mLFPR(jj,K,isnp(i,jw),L) = ' L'//LFAC
              END IF     
            end do                                           
          END DO
          jt=kti(isnp(i,jw))
          je=kb(isnp(i,jw))
          DO K=ktwb(jw),KB(isnp(i,jw))
            WRITE (SNP(jw),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(mlfpr(jj,k,isnp(i,jw),L),jj=jt,je)     
          end do
        end do
      end if
    end do

! v3.5 end

!** Algal nutrient limitations

    DO JA=1,NAL
      IF (LIMITING_FACTOR(JA)) THEN
        IF (NEW_PAGE) THEN
          WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
          NLINES = KMX-KTWB(JW)+14
        END IF
        NLINES   = NLINES+KMX-KTWB(JW)+11
        NEW_PAGE = NLINES > 72
        DO I=IUPR,IDPR     !mlm 6/30/06
          DO K=KTWB(JW),KB(ISNP(I,JW))   !mlm  6/30/06
            LIMIT = MIN(APLIM(K,ISNP(I,JW),JA),ANLIM(K,ISNP(I,JW),JA),ASLIM(K,ISNP(I,JW),JA),ALLIM(K,ISNP(I,JW),JA))
            IF (LIMIT == APLIM(K,ISNP(I,JW),JA)) THEN
              WRITE (LFAC,'(F8.4)') APLIM(K,ISNP(I,JW),JA)
              LFPR(K,ISNP(I,JW)) = ' P'//LFAC 
            ELSE IF (LIMIT == ANLIM(K,ISNP(I,JW),JA)) THEN
              WRITE (LFAC,'(F8.4)') ANLIM(K,ISNP(I,JW),JA)
              LFPR(K,ISNP(I,JW)) = ' N'//LFAC
            ELSE IF (LIMIT == ASLIM(K,ISNP(I,JW),JA)) THEN
              WRITE (LFAC,'(F8.4)') ASLIM(K,ISNP(I,JW),JA)
              LFPR(K,ISNP(I,JW)) = ' S'//LFAC
            ELSE IF (LIMIT == ALLIM(K,ISNP(I,JW),JA)) THEN
              WRITE (LFAC,'(F8.4)') ALLIM(K,ISNP(I,JW),JA)
              LFPR(K,ISNP(I,JW)) = ' L'//LFAC
            END IF
          END DO
        END DO
        WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A,I0,A/)') MONTH,GDAY,',',YEAR,'    Julian Date',INT(JDAY),' days ',                &
                                                         (JDAY-INT(JDAY))*24.0,' hours    Algal group ',JA,' limiting factor'
        WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
        DO K=KTWB(JW),KBR
          WRITE (SNP(JW),'(1X,I4,F8.2,1000A)') K,DEPTHM(K,DS(BS(JW))),(LFPR(K,ISNP(I,JW)),I=IUPR,IDPR)
        END DO
      END IF
    END DO

!** Epiphyton nutrient limitations

    DO JE=1,NEP
      IF (PRINT_EPIPHYTON(JW,JE)) THEN
        IF (NEW_PAGE) THEN
          WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
          NLINES = KMX-KTWB(JW)+14
        END IF
        NLINES   = NLINES+KMX-KTWB(JW)+11
        NEW_PAGE = NLINES > 72
        DO I=IUPR,IDPR                                            !mlm   6/30/2006
          DO K=KTWB(JW),KB(ISNP(I,JW))                            !mlm  6/30/2006
            LIMIT = MIN(EPLIM(K,ISNP(I,JW),JE),ENLIM(K,ISNP(I,JW),JE),ESLIM(K,ISNP(I,JW),JE),ELLIM(K,ISNP(I,JW),JE))
            IF (LIMIT == EPLIM(K,ISNP(I,JW),JE)) THEN
              WRITE (LFAC,'(F8.4)') EPLIM(K,ISNP(I,JW),JE)
              LFPR(K,ISNP(I,JW)) = ' P'//LFAC 
            ELSE IF (LIMIT == ENLIM(K,ISNP(I,JW),JE)) THEN
              WRITE (LFAC,'(F8.4)') ENLIM(K,ISNP(I,JW),JE)
              LFPR(K,ISNP(I,JW)) = ' N'//LFAC
            ELSE IF (LIMIT == ESLIM(K,ISNP(I,JW),JE)) THEN
              WRITE (LFAC,'(F8.4)') ESLIM(K,ISNP(I,JW),JE)
              LFPR(K,ISNP(I,JW)) = ' S'//LFAC
            ELSE IF (LIMIT == ELLIM(K,ISNP(I,JW),JE)) THEN
              WRITE (LFAC,'(F8.4)') ELLIM(K,ISNP(I,JW),JE)
              LFPR(K,ISNP(I,JW)) = ' L'//LFAC
            END IF
          END DO
        END DO
        WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A,I0,A/)') MONTH,GDAY,',',YEAR,'    Julian Date',INT(JDAY),' days ',                &
                                                         (JDAY-INT(JDAY))*24.0,' hours    Epiphyton group ',JE,' limiting factor'
        WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
        DO K=KTWB(JW),KBR
          WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(LFPR(K,ISNP(I,JW)),I=IUPR,IDPR) 
        END DO
      END IF
    END DO
  END IF
END SUBROUTINE OUTPUT

!***********************************************************************************************************************************
!*                                       S U B R O U T I N E    R E S T A R T   O U T P U T                                       **
!***********************************************************************************************************************************

SUBROUTINE RESTART_OUTPUT (RSOFN)
!  USE GLOBAL; USE SCREENC; USE RSTART; USE GDAYC; USE GEOMC; USE KINETIC, ONLY:EPM, EPD; USE TVDC, ONLY:QSUM; USE KINETIC, ONLY:SED
  USE GLOBAL; USE SCREENC; USE RSTART; USE GDAYC; USE GEOMC; USE KINETIC, ONLY:EPM,EPD,SEDC,SEDN,SEDP,PH, sdkv; USE TVDC, ONLY:QSUM !MLM 6/10/07;10/06
  USE KINETIC, ONLY:SED ; USE ZOOPLANKTONC, ONLY: ZOO; USE EDDY, ONLY:TKE                      !MLM 10/06  ! SW 10/5/07
  CHARACTER(*) :: RSOFN
  OPEN  (RSO,FILE=RSOFN,FORM='UNFORMATTED',STATUS='UNKNOWN')
  WRITE (RSO) NIT,    NV,     KMIN,   IMIN,   NSPRF,  CMBRT,  ZMIN,   IZMIN,  START,  CURRENT
  WRITE (RSO) DLTDP,  SNPDP,  TSRDP,  VPLDP,  PRFDP,  CPLDP,  SPRDP,  RSODP,  SCRDP,  FLXDP,  WDODP
  WRITE (RSO) JDAY,   YEAR,   ELTM,   ELTMF,  DLT,    DLTAV,  DLTS,   MINDLT, JDMIN,  CURMAX
  WRITE (RSO) NXTMSN, NXTMTS, NXTMPR, NXTMCP, NXTMVP, NXTMRS, NXTMSC, NXTMSP, NXTMFL, NXTMWD
  WRITE (RSO) VOLIN,  VOLOUT, VOLUH,  VOLDH,  VOLPR,  VOLTRB, VOLDT,  VOLWD,  VOLEV,  VOLSBR, VOLTR, VOLSR
  WRITE (RSO) TSSEV,  TSSPR,  TSSTR,  TSSDT,  TSSWD,  TSSIN,  TSSOUT, TSSS,   TSSB,   TSSICE
  WRITE (RSO) TSSUH,  TSSDH,  TSSUH2, TSSDH2, CSSUH2, CSSDH2, VOLUH2, VOLDH2, QUH1
  WRITE (RSO) ESBR,   ETBR,   EBRI
  WRITE (RSO) Z,      SZ,     ELWS,   SAVH2,  SAVHR,  H2
  WRITE (RSO) KTWB,   KTI,    SKTI,   SBKT   
  WRITE (RSO) ICE,    ICETH,  CUF,    QSUM
  WRITE (RSO) U,      W,      SU,     SW,     AZ,     SAZ,    DLTLIM
  WRITE (RSO) T1,     T2,     C1,     C2,     C1S,    EPD,    SED,    KFS,    CSSK
  WRITE (RSO) SEDC, SEDN, SEDP, ZOO, CD  ! mlm 10/06
  write(rso) sdkv                        ! mlm 6/10/07
  WRITE(RSO) TKE                         ! SW 10/5/07
  CLOSE (RSO)
END SUBROUTINE RESTART_OUTPUT

!***********************************************************************************************************************************
!**                                     S U B R O U T I N E    G R E G O R I A N   D A T E                                        **
!***********************************************************************************************************************************

SUBROUTINE GREGORIAN_DATE
  USE GDAYC

! Determine if new year (regular or leap) and increment year

  DO WHILE (JDAYG >= 366)
    IF (.NOT. LEAP_YEAR .AND. JDAYG >= 366) THEN
      JDAYG     = JDAYG-365
      YEAR      = YEAR+1
      LEAP_YEAR = MOD(YEAR,4) == 0
    ELSE IF (JDAYG >= 367) THEN
      JDAYG     = JDAYG-366
      YEAR      = YEAR+1
      LEAP_YEAR = MOD(YEAR,4) == 0
    ELSE
      EXIT
    END IF
  END DO
  INCR = 0
  IF (LEAP_YEAR) INCR = 1

! Determine month and day of year

  IF (JDAYG >= 1 .AND. JDAYG < 32) THEN
    GDAY  = JDAYG
    DAYM  = 31.0
    MONTH = '  January'
    M     = 1
  ELSE IF (JDAYG >= 32 .AND. JDAYG < 60+INCR) THEN
    GDAY  = JDAYG-31
    DAYM  = 29.0
    MONTH = ' February'
    M     = 2
  ELSE IF (JDAYG >= 60 .AND. JDAYG < 91+INCR) THEN
    GDAY  = JDAYG-59-INCR
    DAYM  = 31.0
    MONTH = '    March'
    M     = 3
  ELSE IF (JDAYG >= 91 .AND. JDAYG < 121+INCR) THEN
    GDAY  = JDAYG-90-INCR
    DAYM  = 30.0
    MONTH = '    April'
    M     = 4
  ELSE IF (JDAYG >= 121 .AND. JDAYG < 152+INCR) THEN
    GDAY  = JDAYG-120-INCR
    DAYM  = 31.0
    MONTH = '      May'
    M     = 5
  ELSE IF (JDAYG >= 152 .AND. JDAYG < 182+INCR) THEN
    GDAY  = JDAYG-151-INCR
    DAYM  = 30.0
    MONTH = '     June'
    M     = 6
  ELSE IF (JDAYG >= 182 .AND. JDAYG < 213+INCR) THEN
    GDAY  = JDAYG-181-INCR
    DAYM  = 31.0
    MONTH = '     July'
    M     = 7
  ELSE IF (JDAYG >= 213 .AND. JDAYG < 244+INCR) THEN
    GDAY  = JDAYG-212-INCR
    DAYM  = 31.0
    MONTH = '   August'
    M     = 8
  ELSE IF (JDAYG >= 244 .AND. JDAYG < 274+INCR) THEN
    GDAY  = JDAYG-243-INCR
    DAYM  = 30.0
    MONTH = 'September'
    M     = 9
  ELSE IF (JDAYG >= 274 .AND. JDAYG < 305+INCR) THEN
    GDAY  = JDAYG-273-INCR
    DAYM  = 31.0
    MONTH = '  October'
    M     = 10
  ELSE IF (JDAYG >= 305 .AND. JDAYG < 335+INCR) THEN
    GDAY  = JDAYG-304-INCR
    DAYM  = 30.0
    MONTH = ' November'
    M     = 11
  ELSE IF (JDAYG >= 335 .AND. JDAYG < 366+INCR) THEN
    GDAY  = JDAYG-334-INCR
    DAYM  = 31.0
    MONTH = ' December'
    M     = 12
  END IF
END SUBROUTINE GREGORIAN_DATE

!***********************************************************************************************************************************
!**                                            S U B R O U T I N E    W A T E R B O D Y                                           **
!***********************************************************************************************************************************

SUBROUTINE WATERBODY
  USE GLOBAL; USE GEOMC; USE TVDC; USE LOGICC

! Type declarations

  REAL, SAVE, ALLOCATABLE, DIMENSION(:)   :: ELL,    ELR,    CL
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: QU,     QD
  REAL                                    :: C(KMX,IMX), SS(KMX,IMX)

! Allocation declarations

  ALLOCATE (ELL(KMX), ELR(KMX), CL(NCT), QU(KMX,IMX), QD(KMX,IMX))

! Variable initialization

  ELL = 0.0; ELR = 0.0; CL = 0.0; QU = 0.0; QD = 0.0

! Debug variable
!  ncount=0
! End debug

RETURN

!***********************************************************************************************************************************
!**                                               U P S T R E A M   V E L O C I T Y                                               **
!***********************************************************************************************************************************

ENTRY UPSTREAM_VELOCITY
  DO JJB=1,NBR
    IF (UHS(JB) >= US(JJB) .AND. UHS(JB) <= DS(JJB)) EXIT
  END DO
  DO JJW=1,NWB
    IF (JJB >=BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  DO K=KTWB(JJW),KB(UHS(JB))+1
    ELL(K) = EL(K,UHS(JB))-SINA(JJB)*DLX(UHS(JB))*0.5
  END DO
  DO K=KT,KB(IU)+1
    ELR(K) = EL(K,IU)+SINA(JB)*DLX(IU)*0.5
  END DO
  ELWS = EL(KTWB(JJW),UHS(JB))-Z(UHS(JB))*COSA(JJB)
  EL1  = ELWS-SINA(JJB)*DLX(UHS(JB))*0.5
  ELWS = EL1
  KL   = KTWB(JJW)+1
  DO K=KT+1,KB(IU)+1
    IF (ELR(K) >= ELL(KL)) THEN
      IF (KL == KTWB(JJW)+1) THEN
        Q1 = U(KL-1,UHS(JB))*BHR1(KTWB(JJW),UHS(JB))
        IF (KL == KB(UHS(JB))+1 .AND. ELL(KL) < ELR(KB(IU)+1)) THEN
          HT = ELWS-ELR(KB(IU)+1)
        ELSE
          HT = H1(KTWB(JJW),UHS(JB))
        END IF
      ELSE
        Q1 = U(KL-1,UHS(JB))*BHR1(KL-1,UHS(JB))
        HT = H1(KL-1,UHS(JB))
      END IF
      IF (K == KT+1) THEN
        EL1         = EL(KTWB(JJW),UHS(JB))-Z(UHS(JB))*COSA(JJB)-SINA(JJB)*DLX(UHS(JB))*0.5
        U(K-1,IU-1) = Q1*((EL1-ELR(K))/HT)/BHR1(KT,IU-1)
      ELSE
        U(K-1,IU-1) = Q1*((EL1-ELR(K))/HT)/BHR1(K-1,IU-1)
      END IF
      EL1 = ELR(K)
      IF (ELR(K) == ELL(KL)) KL = KL+1
    ELSE
      Q1 = 0.0
      DO WHILE (ELR(K) <= ELL(KL))
        IF (KL == KTWB(JJW)+1 .AND. K == KT+1) THEN
          Q1 = Q1+U(KL-1,UHS(JB))*BHR1(KTWB(JJW),UHS(JB))
        ELSE IF (KL == KTWB(JJW)+1) THEN   
          Q1 = Q1+U(KL-1,UHS(JB))*BHR1(KTWB(JJW),UHS(JB))*(EL1-ELL(KL))/H1(KTWB(JJW),UHS(JB))
        ELSE   
          FRAC = (EL1-ELL(KL))/H1(KL-1,UHS(JB))
          Q1   =  Q1+U(KL-1,UHS(JB))*BHR1(KL-1,UHS(JB))*FRAC
        END IF
        EL1 = ELL(KL)
        KL  = KL+1
        IF (KL > KB(UHS(JB))) EXIT
      END DO
      IF (K == KT+1) THEN
        BRTOT = BHR1(KT,IU-1)
      ELSE
        BRTOT = BHR1(K-1,IU-1)
      END IF
      FRAC = 0.0
!      IF (KL < KMX) THEN                             ! SW 6/29/06
        HT   =  H1(KL-1,UHS(JB))                       
        FRAC = (EL1-ELR(K))/HT
        IF (KB(UHS(JB)) >= KB(IU-1) .AND. K > KB(IU-1)) FRAC = (EL1-ELL(KL))/HT     ! SW 6/29/06
        IF (KL >= KMX .and. ELR(K) < ELL(KL)) FRAC=1.0          ! SW 6/29/06
        Q1   =  Q1+U(KL-1,UHS(JB))*BHR1(KL-1,UHS(JB))*FRAC
!      ELSE                                           ! SW 6/29/06
!        Q1 = Q1+U(KL-1,UHS(JB))*BHR1(KL-1,UHS(JB))   ! SW 6/29/06
!      END IF                                         ! SW 6/29/06
      U(K-1,IU-1) = Q1/BRTOT
      IF (KL > KB(UHS(JB))) THEN
        IF (FRAC < 1.0 .AND. FRAC /= 0.0) THEN
          IF (K == KB(IU)+1) THEN
            U(K-1,IU-1) = (Q1+U(KL-1,UHS(JB))*BHR1(KL-1,UHS(JB))*(1.0-FRAC))/BRTOT
          ELSE
            U(K,IU-1) = U(KL-1,UHS(JB))*BHR1(KL-1,UHS(JB))*(1.0-FRAC)/BHR(K,IU-1)
          END IF
        END IF
        GO TO 100
      END IF
      EL1 = ELR(K)
    END IF
  END DO
100 CONTINUE

! Debug 
!    QL=0.0; QR=0.0
!  do k=ktwb(jjw),kb(uhs(jb))
!    QL=QL+u(k,uhs(jb))*BHR1(k,uhs(jb))
!  enddo
!  do k=kt,kb(iu-1)
!    QR=QR+u(k,iu-1)*BHR1(k,iu-1)
!  enddo
!  if((QR-QL)/QL > 0.02)then
!     if(ncount.eq.0)open(299,file='debug_out.txt',status='unknown')
!     ncount=ncount+1
!     write(299,*)'***QL=',QL,' QR=',QR
!     write(299,*)'UHS(JB)=',uhs(jb),'  iu-1=',iu-1,' ELWS=', elws
!     write(299,*)'LEFT: K      U     BHR1    ELEV      H1     BHR    ELL'
!     do k=ktwb(jjw),kb(uhs(jb))
!       write(299,'(i8,f8.3,f8.3,5f8.3)')k,u(k,uhs(jb)),BHR1(k,uhs(jb)),el(k,uhs(jb)),h1(k,uhs(jb)),bhr(k,uhs(jb)),ell(k)
!     enddo
!     write(299,*)'RIGHT: K     U     BHR1    ELEV      H1     BHR    ELR'
!     do k=kt,kb(iu-1)
!       write(299,'(i8,f8.3,f8.3,5f8.3)')k,u(k,iu-1),BHR1(k,iu-1),el(k,iu-1),h1(k,iu-1),bhr(k,iu-1),elr(k)
!     enddo
!  end if
! End debug

RETURN

!***********************************************************************************************************************************
!**                                              U P S T R E A M   W A T E R B O D Y                                              **
!***********************************************************************************************************************************

ENTRY UPSTREAM_WATERBODY
  DO JJB=1,NBR
    IF (UHS(JB) >= US(JJB) .AND. UHS(JB) <= DS(JJB)) EXIT
  END DO
  DO JJW=1,NWB
    IF (JJB >=BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  DO K=KTWB(JJW),KB(UHS(JB))+1
    ELL(K) = EL(K,UHS(JB))-SINA(JJB)*DLX(UHS(JB))*0.5
  END DO
  DO K=KT,KB(IU)+1
    ELR(K) = EL(K,IU)+SINA(JB)*DLX(IU)*0.5
  END DO
  ELWS = EL(KTWB(JJW),UHS(JB))-Z(UHS(JB))*COSA(JJB)
  EL1  = ELWS-SINA(JJB)*DLX(UHS(JB))*0.5
  KL   = KTWB(JJW)+1
  DO K=KT+1,KB(IU)+1
    IF (ELR(K) >= ELL(KL)) THEN
      T1(K-1,IU-1)            = T1(KL-1,UHS(JB))
      T2(K-1,IU-1)            = T2(KL-1,UHS(JB))
      C1S(K-1,IU-1,CN(1:NAC)) = C1S(KL-1,UHS(JB),CN(1:NAC))
      C1(K-1,IU-1,CN(1:NAC))  = C1S(KL-1,UHS(JB),CN(1:NAC))
      C2(K-1,IU-1,CN(1:NAC))  = C1S(KL-1,UHS(JB),CN(1:NAC))
      EL1                     = ELR(K)
      IF (ELR(K) == ELL(KL)) KL = KL+1
    ELSE
      BRTOT = 0.0
      CL    = 0.0
      T1L   = 0.0
      T2L   = 0.0
      DO WHILE (ELR(K) <= ELL(KL))
        IF (KL == KTWB(JJW)+1 .AND. K == KT+1) THEN
          B1 = BH2(KTWB(JJW),UHS(JB))
        ELSE
          B1 = B(KL-1,UHS(JB))*(EL1-ELL(KL))
        END IF
        BRTOT         = BRTOT+B1
        T1L           = T1L+B1*T1(KL-1,UHS(JB))
        T2L           = T2L+B1*T2(KL-1,UHS(JB))
        CL(CN(1:NAC)) = CL(CN(1:NAC))+B1*C1S(KL-1,UHS(JB),CN(1:NAC))
        EL1 = ELL(KL)
        KL  = KL+1
        IF (KL > KB(UHS(JB))+1) EXIT
      END DO
      IF (KL <= KB(UHS(JB))+1) THEN
        B1    = B(KL-1,UHS(JB))*(EL1-ELR(K))
        BRTOT = BRTOT+B1
        IF (BRTOT > 0.0) THEN   
          T1(K-1,IU-1)            = (T1L+B1*T1(KL-1,UHS(JB)))/BRTOT
          T2(K-1,IU-1)            = (T2L+B1*T2(KL-1,UHS(JB)))/BRTOT
          C1S(K-1,IU-1,CN(1:NAC)) = (CL(CN(1:NAC))+B1*C1S(KL-1,UHS(JB),CN(1:NAC)))/BRTOT
          C1(K-1,IU-1,CN(1:NAC))  =  C1S(K-1,IU-1,CN(1:NAC))
          C2(K-1,IU-1,CN(1:NAC))  =  C1S(K-1,IU-1,CN(1:NAC))
        ELSE   
          T1(K-1,IU-1)            = T1(KL-1,UHS(JB))
          T2(K-1,IU-1)            = T2(KL-1,UHS(JB))
          C1S(K-1,IU-1,CN(1:NAC)) = C1S(KL-1,UHS(JB),CN(1:NAC))
          C1(K-1,IU-1,CN(1:NAC))  = C1S(K-1,IU-1,CN(1:NAC))
          C2(K-1,IU-1,CN(1:NAC))  = C1S(K-1,IU-1,CN(1:NAC))
        END IF     
      ELSE
        IF (BRTOT > 0.0) THEN   
          T1(K-1,IU-1)            = T1L          /BRTOT
          T2(K-1,IU-1)            = T2L          /BRTOT
          C1S(K-1,IU-1,CN(1:NAC)) = CL(CN(1:NAC))/BRTOT
          C1(K-1,IU-1,CN(1:NAC))  = C1S(K-1,IU-1,CN(1:NAC))
          C2(K-1,IU-1,CN(1:NAC))  = C1S(K-1,IU-1,CN(1:NAC))
        ELSE   
          T1(K-1,IU-1)            = T1(KL-1,UHS(JB))
          T2(K-1,IU-1)            = T2(KL-1,UHS(JB))
          C1S(K-1,IU-1,CN(1:NAC)) = C1S(KL-1,UHS(JB),CN(1:NAC))
          C1(K-1,IU-1,CN(1:NAC))  = C1S(K-1,IU-1,CN(1:NAC))
          C2(K-1,IU-1,CN(1:NAC))  = C1S(K-1,IU-1,CN(1:NAC))
        END IF
        EXIT
      END IF
      EL1 = ELR(K)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                            D O W N S T R E A M   W A T E R B O D Y                                            **
!***********************************************************************************************************************************

ENTRY DOWNSTREAM_WATERBODY
  DO JJB=1,NBR
    IF (CDHS(JB) >= CUS(JJB) .AND. CDHS(JB) <= DS(JJB)) EXIT
  END DO
  DO JJW=1,NWB
    IF (JJB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  KR = KTWB(JJW)+1
  DO K=KT,KB(ID)+1
    ELL(K) = EL(K,ID)-SINA(JB)*DLX(ID)*0.5
  END DO
  DO K=KTWB(JJW),KB(CDHS(JB))+1
    ELR(K) = EL(K,CDHS(JB))+SINA(JJB)*DLX(CDHS(JB))*0.5
  END DO
  ELWS = EL(KTWB(JW),ID)-Z(ID)*COSA(JB)
  EL1  = ELWS-SINA(JB)*DLX(ID)*0.5
  IDT  = ID+1
  DO K=KT+1,KB(ID)+1
    IF (ELL(K) >= ELR(KR)) THEN
      T1(K-1,IDT)            = T1(KR-1,CDHS(JB))
      T2(K-1,IDT)            = T2(KR-1,CDHS(JB))
      C1S(K-1,IDT,CN(1:NAC)) = C1S(KR-1,CDHS(JB),CN(1:NAC))
      C1(K-1,IDT,CN(1:NAC))  = C1S(KR-1,CDHS(JB),CN(1:NAC))
      C2(K-1,IDT,CN(1:NAC))  = C1S(KR-1,CDHS(JB),CN(1:NAC))
      EL1                    = ELL(K)
      IF (ELL(K) == ELR(KR)) KR = KR+1
      IF (KR > KB(CDHS(JB))+1) EXIT
    ELSE
      BRTOT = 0.0
      CL    = 0.0
      T1L   = 0.0
      T2L   = 0.0
      DO WHILE (ELL(K) <= ELR(KR))
        IF (KR == KTWB(JJW)+1 .AND. K == KT+1) THEN
          B1    = BH2(KTWB(JJW),CDHS(JB))
          BRTOT = BRTOT+B1
        ELSE
          B1    = B(KR-1,CDHS(JB))*(EL1-ELR(KR))
          BRTOT = BRTOT+B1
        END IF
        T1L           = T1L+B1*T1(KR-1,CDHS(JB))
        T2L           = T2L+B1*T2(KR-1,CDHS(JB))
        CL(CN(1:NAC)) = CL(CN(1:NAC))+B1*C1S(KR-1,CDHS(JB),CN(1:NAC))
        EL1           = ELR(KR)
        KR            = KR+1
        IF (KR > KB(CDHS(JB))+1) EXIT
      END DO
      IF (KR <= KB(CDHS(JB)+1)) THEN
        B1 = B(KR-1,CDHS(JB))*(EL1-ELL(K))
      ELSE   
        B1 = 0.0
      END IF
      BRTOT = BRTOT+B1
      IF (BRTOT == 0.0) EXIT
      T1(K-1,IDT)            = (T1L+B1*T1(KR-1,CDHS(JB)))/BRTOT
      T2(K-1,IDT)            = (T2L+B1*T2(KR-1,CDHS(JB)))/BRTOT
      C1S(K-1,IDT,CN(1:NAC)) = (CL(CN(1:NAC))+B1*C1S(KR-1,CDHS(JB),CN(1:NAC)))/BRTOT
      C1(K-1,IDT,CN(1:NAC))  =  C1S(K-1,IDT,CN(1:NAC))
      C2(K-1,IDT,CN(1:NAC))  =  C1S(K-1,IDT,CN(1:NAC))
      EL1                    =  ELL(K)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                 U P S T R E A M   B R A N C H                                                 **
!***********************************************************************************************************************************

ENTRY UPSTREAM_BRANCH
  DO JJW=1,NWB
    IF (JJB >=BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  DO K=KT,KB(I)+1
    ELL(K) = EL(K,I)
  END DO
  DO K=KTWB(JWUH(JB)),KB(CUS(JJB))+1
    ELR(K) = EL(K,CUS(JJB))+SINA(JJB)*DLX(CUS(JJB))*0.5
  END DO
  ELWS = EL(KTWB(JWUH(JB)),CUS(JJB))-Z(CUS(JJB))*COSA(JJB)
  EL1  = ELWS+SINA(JJB)*DLX(CUS(JJB))*0.5
  KR   = KTWB(JWUH(JB))+1
  DO K=KT+1,KB(I)+1
    IF (ELL(K) >= ELR(KR)) THEN
      Q1 = VOLUH2(KR-1,JJB)/DLT
      B1 = BHR(KR-1,CUS(JJB)-1)
      IF (KR == KTWB(JWUH(JB))+1) B1 = BHR2(KT,CUS(JJB)-1)
      U1 = U(KR-1,CUS(JJB)-1)*B1
      HT = H2(KR-1,JWUH(JB))
      IF (KR == KTWB(JWUH(JB))+1) HT = H2(KT,CUS(JJB))
      Q1 = Q1*((EL1-ELL(K))/HT)
      B2 = BHR(K-1,I)
      IF (K == KTWB(JW)+1) B2 = BHR2(KT,I)
      UXBR(K-1,I) = UXBR(K-1,I)+(ABS((U1/B2)*COS(BETABR)*Q1)/DLX(I))
      UYBR(K-1,I) = UYBR(K-1,I)+ABS(Q1*SIN(BETABR))
      EL1         = ELL(K)
      IF (ELL(K) == ELR(KR)) KR = KR+1
    ELSE
      U1    = 0.0
      Q1    = 0.0
      BRTOT = 0.0
      DO WHILE (ELL(K) <= ELR(KR))
        IF (KR /= KTWB(JWUH(JB))+1) THEN
          FRAC = (EL1-ELR(KR))/(H(KR-1,JWUH(JB)))
          B1   =  BHR(KR-1,CUS(JJB)-1)
        ELSE
          FRAC = (EL1-ELR(KR))/H2(KT,CUS(JJB))
          B1   =  BHR2(KT,CUS(JJB)-1)
        END IF
        U1  = U1+U(KR-1,CUS(JJB)-1)*B1*FRAC
        Q1  = Q1+VOLUH2(KR-1,JJB)/DLT*FRAC
        EL1 = ELR(KR)
        KR  = KR+1
        IF (KR > KB(CUS(JJB)+1)) EXIT
      END DO
      IF (K == KTWB(JW)+1) THEN
        B2 = BHR2(KT,I)
      ELSE   
        B2 = BHR(K-1,I)
      END IF
      IF (H(KR-1,JWUH(JB)) /= 0.0) THEN
        IF (KR-1 == KTWB(JWUH(JB))) THEN
          HT = H2(KT,CUS(JJB))
        ELSE
          HT = H2(KR-1,JWUH(JB))
        END IF     
        FRAC        = (EL1-ELL(KR-1))/HT
        Q1          =  Q1+FRAC*VOLUH2(KR-1,JJB)/DLT
        UXBR(K-1,I) =  UXBR(K-1,I)+(ABS((U1/B2)*COS(BETABR)*Q1)/DLX(I))
        UYBR(K-1,I) =  UYBR(K-1,I)+ABS(Q1*SIN(BETABR))
      END IF
      IF (KR > KB(CUS(JJB)+1)) EXIT
      EL1 = ELL(K)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                               D O W N S T R E A M   B R A N C H                                               **
!***********************************************************************************************************************************

ENTRY DOWNSTREAM_BRANCH
  DO JJW=1,NWB
    IF (JJB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  DO K=KT,KB(I)+1
    ELR(K) = EL(K,I)
  END DO
  DO K=KTWB(JJW),KB(DS(JJB))+1
    ELL(K) = EL(K,DS(JJB))+SINA(JJB)*DLX(DS(JJB))*0.5
  END DO
  ELWS = EL(KTWB(JJW),DS(JJB))-Z(DS(JJB))*COSA(JJB)
  EL1  = ELWS-SINA(JJB)*DLX(DS(JJB))*0.5
  KL   = KTWB(JJW)+1
  DO K=KT+1,KB(I)+1
    IF (ELR(K) >= ELL(KL)) THEN
      Q1 = VOLDH2(KL-1,JJB)/DLT
      B1 = BHR(KL-1,DS(JJB)-1)
      IF (KL == KTWB(JJW)+1) B1 = BHR2(KT,DS(JJB))
      U1 = U(KL-1,DS(JJB)-1)*B1
      HT = H2(KL-1,JJW)
      IF (KL == KTWB(JJW)+1) HT = H2(KT,DS(JJB))
      Q1 = Q1*((EL1-ELR(K))/HT)
      B2 = BHR(K-1,I)
      IF (K == KTWB(JW)+1) B2 = BHR2(KT,I)
      UXBR(K-1,I) = UXBR(K-1,I)+(ABS((U1/B2)*COS(BETABR)*Q1)/DLX(I))
      UYBR(K-1,I) = UYBR(K-1,I)+ABS(Q1*SIN(BETABR))
      EL1 = ELR(K)
      IF (ELR(K) == ELL(KL)) KL = KL+1
    ELSE
      U1    = 0.0
      Q1    = 0.0
      BRTOT = 0.0
      DO WHILE (ELR(K) <= ELL(KL))
        IF (KL /= KTWB(JJW)+1) THEN
          FRAC = (EL1-ELL(KL))/H(KL-1,JJW)
          B1   = BHR(KL-1,DS(JJB)-1)
        ELSE
          FRAC = (EL1-ELL(KL))/H2(KT,DS(JJB))
          B1   = BHR2(KT,DS(JJB))
        END IF
        U1  = U1+U(KL-1,DS(JJB))*B1*FRAC
        Q1  = Q1+(VOLDH2(KL-1,JJB)/DLT)*FRAC
        EL1 = ELL(KL)
        KL  = KL+1
        IF (KL > KB(DS(JJB)+1)) EXIT
      END DO
      IF (K == KTWB(JW)+1) THEN
        B2 = BHR2(KT,I)
      ELSE
        B2 = BHR(K-1,I)
      END IF
      IF (H(KL-1,JJW) /= 0.0) THEN
        IF (KL-1 == KTWB(JJW)) THEN
          HT = H2(KT,DS(JJB))
        ELSE
          HT = H2(KL-1,JJW)
        END IF     
        FRAC        = (EL1-ELL(KL-1))/HT
        Q1          =  Q1+FRAC*(VOLDH2(KL-1,JJB)/DLT)
        UXBR(K-1,I) =  UXBR(K-1,I)+(ABS((U1/B2)*COS(BETABR)*Q1)/DLX(I))
        UYBR(K-1,I) =  UYBR(K-1,I)+ABS(Q1*SIN(BETABR))
      END IF
      IF (KL > KB(DS(JJB)+1)) EXIT
      EL1 = ELR(K)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                   U P S T R E A M   F L O W                                                   **
!***********************************************************************************************************************************

ENTRY UPSTREAM_FLOW
  DO JJB=1,NBR
    IF (UHS(JB) >= US(JJB) .AND. UHS(JB) <= DS(JJB)) EXIT
  END DO
  DO JJW=1,NWB
    IF (JJB >=BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  DO K=KTWB(JWUH(JB)), KB(UHS(JB))+1
    ELL(K) = EL(K,UHS(JB))
  END DO
  DO K=KT,KB(IU)+1
    ELR(K) = EL(K,IU)+SINA(JB)*DLX(IU)*0.5
  END DO
  ELWS = EL(KTWB(JW),IU)-Z(IU)*COSA(JB)
  EL1  = ELWS+SINA(JB)*DLX(IU)*0.5
  KR   = KT+1
  DO K=KTWB(JWUH(JB))+1,KB(UHS(JB))+1
    IF (ELL(K) >= ELR(KR)) THEN
      Q1 = VOLUH2(KR-1,JB)/DLT
      HT = H2(KR-1,JW)
      IF (KR == KTWB(JW)+1) THEN
        QU(K-1,UHS(JB))  = Q1
        QSS(K-1,UHS(JB)) = QSS(K-1,UHS(JB))-Q1
      ELSE
        QU(K-1,UHS(JB))  = Q1*((EL1-ELL(K))/HT)
        QSS(K-1,UHS(JB)) = QSS(K-1,UHS(JB))-QU(K-1,UHS(JB))
      END IF
      EL1 = ELL(K)
      IF (ELL(K) == ELR(KR)) KR = KR+1
    ELSE
      Q1 = 0.0
      DO WHILE (ELL(K) <= ELR(KR))
        IF (KR /= KTWB(JW)+1) THEN
          FRAC = (EL1-ELR(KR))/(H(KR-1,JW))
        ELSE 
          FRAC = (EL1-ELR(KR))/(H2(KT,IU))
        END IF
        Q1  = Q1+VOLUH2(KR-1,JB)/DLT*FRAC
        EL1 = ELR(KR)
        KR  = KR+1
        IF (KR > KB(IU)+1) EXIT
      END DO
      IF (H(KR-1,JW) /= 0.0) THEN
        FRAC = (EL1-ELL(K))/H(KR-1,JW)
        QU(K-1,UHS(JB)) = Q1+FRAC*VOLUH2(KR-1,JB)/DLT
      ELSE
        QU(K-1,UHS(JB)) = Q1
      END IF
      QSS(K-1,UHS(JB)) = QSS(K-1,UHS(JB))-QU(K-1,UHS(JB))
      IF (KR > KB(IU)+1) EXIT
      EL1 = ELL(K)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                 D O W N S T R E A M   F L O W                                                 **
!***********************************************************************************************************************************

ENTRY DOWNSTREAM_FLOW
  DO JJB=1,NBR
    IF (CDHS(JB) >= US(JJB) .AND. CDHS(JB) <= DS(JJB)) EXIT
  END DO
  DO JJW=1,NWB
    IF (JJB >=BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  DO K=KTWB(JJW),KB(CDHS(JB))+1
    ELR(K) = EL(K,CDHS(JB))
  END DO
  DO K=KT,KB(ID)+1
    ELL(K) = EL(K,ID)-SINA(JB)*DLX(ID)*0.5
  END DO
  ELWS = EL(KTWB(JW),ID)-Z(ID)*COSA(JB)
  EL1  = ELWS-SINA(JB)*DLX(ID)*0.5
  KL   = KT+1
  DO K=KTWB(JJW)+1,KB(CDHS(JB))+1
    IF (ELR(K) >= ELL(KL)) THEN
      Q1 = VOLDH2(KL-1,JB)/DLT
      HT = H2(KL-1,JW)
      IF (KL == KTWB(JW)+1) HT = H2(KT,ID)
      QD(K-1,CDHS(JB))  = Q1*((EL1-ELR(K))/HT)
      QSS(K-1,CDHS(JB)) = QSS(K-1,CDHS(JB))+QD(K-1,DHS(JB))
      EL1               = ELR(K)
      IF (ELR(K) == ELL(KL)) KL = KL+1
    ELSE
      Q1 = 0.0
      DO WHILE (ELR(K) <= ELL(KL))
        IF (KL /= KTWB(JW)+1) THEN
          FRAC = (EL1-ELL(KL))/(H(KL-1,JW))
        ELSE   
          FRAC = (EL1-ELL(KL))/(H2(KT,ID))
        END IF     
        Q1  = Q1+(VOLDH2(KL-1,JB)/DLT)*FRAC
        EL1 = ELL(KL)
        KL  = KL+1
        IF (KL > KB(ID)+1) EXIT
      END DO
      IF (H(KL-1,JW) /= 0.0) THEN
        FRAC             = (EL1-ELR(K))/H(KL-1,JW)
        QD(K-1,CDHS(JB)) =  Q1+FRAC*(VOLDH2(KL-1,JB)/DLT)
      ELSE
        QD(K-1,CDHS(JB)) = Q1
      END IF
      QSS(K-1,CDHS(JB)) = QSS(K-1,CDHS(JB))+QD(K-1,CDHS(JB))
      IF (KL > KB(ID)+1) EXIT
      EL1 = ELR(K)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                            U P S T R E A M   C O N S T I T U E N T                                            **
!***********************************************************************************************************************************

ENTRY UPSTREAM_CONSTITUENT(C,SS)
  DO JJB=1,NBR
    IF (UHS(JB) >= US(JJB) .AND. UHS(JB) <= DS(JJB)) EXIT
  END DO
  DO JJW=1,NWB
    IF (JJB >=BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  DO K=KTWB(JWUH(JB)),KB(UHS(JB))+1
    ELL(K) = EL(K,UHS(JB))
  END DO
  DO K=KT,KB(IU)+1
    ELR(K) = EL(K,IU)+SINA(JB)*DLX(IU)*0.5
  END DO
  ELWS = EL(KTWB(JW),IU)-Z(IU)*COSA(JB)
  EL1  = ELWS+SINA(JB)*DLX(IU)*0.5
  KR   = KT+1
  DO K=KTWB(JWUH(JB))+1,KB(UHS(JB))+1
    IUT = IU
    IF (QU(K-1,UHS(JB)) >= 0.0) IUT = IU-1
    IF (ELL(K) >= ELR(KR)) THEN
      T1L             = C(KR-1,IUT)
      SS(K-1,UHS(JB)) = SS(K-1,UHS(JB))-T1L*QU(K-1,UHS(JB))
      EL1             = ELL(K)
      IF (ELL(K) == ELR(KR)) KR = KR+1
    ELSE
      T1L   = 0.0
      BRTOT = 0.0
      DO WHILE (ELL(K) <= ELR(KR))
        IF (KR == KT+1 .AND. K == KTWB(JWUH(JB))+1) THEN
          B1    = BH2(KT,IU)
          BRTOT = BRTOT+B1
        ELSE
          B1    = B(KR-1,IU)*(EL1-ELR(KR))
          BRTOT = BRTOT+B1
        END IF
        IUT = IU
        IF (QU(K-1,UHS(JB)) >= 0.0) IUT = IU-1
        T1L = T1L+B1*C(KR-1,IUT)
        EL1 = ELR(KR)
        KR  = KR+1
        IF (KR > KB(IU)+1) EXIT
      END DO
      IUT = IU     
      IF (QU(K-1,UHS(JB)) >= 0.0) IUT = IU-1
      B1              =  B(KR-1,IU)*(EL1-ELL(K))
      BRTOT           =  BRTOT+B1
      T1L             = (T1L+B1*C(KR-1,IUT))/BRTOT
      SS(K-1,UHS(JB)) = TSS(K-1,UHS(JB))-T1L*QU(K-1,UHS(JB))
      IF (KR > KB(IU)+1) EXIT
      EL1 = ELL(K)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                          D O W N S T R E A M   C O N S T I T U E N T                                          **
!***********************************************************************************************************************************

ENTRY DOWNSTREAM_CONSTITUENT (C,SS)
  DO JJB=1,NBR
    IF (CDHS(JB) >= US(JJB) .AND. CDHS(JB) <= DS(JJB)) EXIT
  END DO
  DO JJW=1,NWB
    IF (JJB >=BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  DO K=KTWB(JJW),KB(CDHS(JB))+1
    ELR(K) = EL(K,CDHS(JB))
  END DO
  DO K=KT,KB(ID)+1
    ELL(K) = EL(K,ID)+SINA(JB)*DLX(ID)*0.5
  END DO
  ELWS = EL(KTWB(JW),ID)-Z(ID)*COSA(JB)
  EL1  = ELWS+SINA(JB)*DLX(ID)*0.5
  KL   = KT+1
  DO K=KTWB(JJW)+1,KB(CDHS(JB))+1
    IDT = ID+1      
    IF (QD(K-1,CDHS(JB)) >= 0.0) IDT = ID
    IF (ELR(K) >= ELL(KL)) THEN
      T1L              = C(KL-1,IDT)
      SS(K-1,CDHS(JB)) = SS(K-1,CDHS(JB))+T1L*QD(K-1,CDHS(JB))
      EL1              = ELR(K)
      IF (ELR(K) == ELL(KL)) KL = KL+1
    ELSE
      T1L   = 0.0
      BRTOT = 0.0
      DO WHILE (ELR(K) <= ELL(KL))
        IF (KL == KTWB(JW)+1 .AND. K == KTWB(JJW)+1) THEN
          B1    = BH2(KT,ID)
          BRTOT = BRTOT+B1
        ELSE
          B1    = B(KL-1,ID)*(EL1-ELL(KL))
          BRTOT = BRTOT+B1
        END IF
        T1L = T1L+B1*C(KL-1,IDT)
        EL1 = ELL(KL)
        KL  = KL+1
        IF (KL > KB(ID)) THEN
          B1               =  B(KL-1,ID)*(EL1-ELR(K))
          BRTOT            =  BRTOT+B1
          T1L              = (T1L+B1*C(KL-1,IDT))/BRTOT
          SS(K-1,CDHS(JB)) =  SS(K-1,CDHS(JB))+T1L*QD(K-1,CDHS(JB))
          GO TO 200
        END IF
      END DO
      B1               =  B(KL-1,ID)*(EL1-ELR(K))
      BRTOT            =  BRTOT+B1
      T1L              = (T1L+B1*C(KL-1,IDT))/BRTOT
      SS(K-1,CDHS(JB)) =  SS(K-1,CDHS(JB))+T1L*QD(K-1,CDHS(JB))
      EL1              =  ELL(KL)
    END IF
  END DO
200 CONTINUE
RETURN
ENTRY DEALLOCATE_WATERBODY
  DEALLOCATE (ELL, ELR, CL, QU, QD)
RETURN
END SUBROUTINE WATERBODY

!***********************************************************************************************************************************
!**                                          S U B R O U T I N E   G A T E  F L O W                                               **
!***********************************************************************************************************************************

SUBROUTINE GATE_FLOW
  USE STRUCTURES; USE GLOBAL; USE GEOMC

  DO JG=1,NGT
   if(dyngtc(jg) == '    FLOW')then    ! V3.5   
    QGT(JG) = BGT(JG)
   else  

    ELIU =  EL(KTWB(JWUGT(JG)),IUGT(JG))-Z(IUGT(JG))*COSA(BS(JWUGT(JG)))
    IF (LATERAL_GATE(JG)) THEN   
      ELIU =  EL(KTWB(JWUGT(JG)),IUGT(JG))-Z(IUGT(JG))*COSA(BS(JWUGT(JG)))
    ELSE   
      ELIU=EL(KTWB(JWUGT(JG)),IUGT(JG))-Z(IUGT(JG))*COSA(BS(JWUGT(JG)))-SINA(JBUGT(JG))*DLX(IUGT(JG))*0.5
    END IF  
    IF (IDGT(JG) /= 0)THEN
      IF (US(JBDGT(JG)) /= IDGT(JG)) THEN
        ELID = EL(KTWB(JWDGT(JG)),IDGT(JG))-Z(IDGT(JG))*COSA(BS(JWDGT(JG)))
      ELSE   
        ELID = EL(KTWB(JWDGT(JG)),IDGT(JG))-Z(IDGT(JG))*COSA(BS(JWDGT(JG)))+SINA(JBDGT(JG))*DLX(IDGT(JG))*0.5
      END IF
    ELSE   
      ELID = -100.0
    END IF
    IF (BGT(JG) /= 0.0) THEN
      IF (ELID > EGT(JG) .OR. ELIU > EGT(JG)) THEN
        ISUB = 0
        IF (A2GT(JG) /= 0.0 .AND. IDGT(JG) /= 0.0) THEN
          HTAIL =  ELID-EGT(JG)                                                          ! SW 5/10/05
          IF (HTAIL > 0) THEN
            HENERGY = (U(KTWB(JWUGT(JG)),IUGT(JG))**2)/(2.0*G)+ELIU-EGT(JG)              ! SW 5/10/05
            IF (HTAIL/HENERGY > 0.67) ISUB = 1
          END IF
        END IF
        IGT = 0
        IF (BGT(JG) >= 0.8*(ELIU-EGT(JG)) .AND. GTA1(JG) /= 0.0) IGT = 1
        IF (IGT == 0) THEN
          IF (ISUB == 0) THEN
            DLEL = ELIU-EGT(JG)
            IF (A2GT(JG) == 0.0 .AND. G2GT(JG) /= 0.0) DLEL = ELIU-G2GT(JG)
            IF (DLEL < 0.0) THEN
              DLEL    = -DLEL
              QGT(JG) = -A1GT(JG)*(DLEL**B1GT(JG))*BGT(JG)**G1GT(JG)
            ELSE
              QGT(JG) =  A1GT(JG)*(DLEL**B1GT(JG))*BGT(JG)**G1GT(JG)
            END IF
          ELSE IF (ELID > ELIU) THEN
            DLEL    =  ELID-ELIU
            QGT(JG) = -A2GT(JG)*DLEL**B2GT(JG)*BGT(JG)**G2GT(JG)
          ELSE
            DLEL    =  ELIU-ELID
            QGT(JG) =  A2GT(JG)*DLEL**B2GT(JG)*BGT(JG)**G2GT(JG)
          END IF
        ELSE IF (ISUB == 0) THEN
          DLEL = ELIU-EGT(JG)
          IF (ELID > EGT(JG)) DLEL = ELIU-ELID
          IF (DLEL < 0.0) THEN
            DLEL    = -DLEL
            QGT(JG) = -GTA1(JG)*DLEL**GTB1(JG)
          ELSE
            QGT(JG) =  GTA1(JG)*DLEL**GTB1(JG)
          END IF
        ELSE IF (ELID > ELIU) THEN
          DLEL    =  ELID-ELIU
          QGT(JG) = -GTA2(JG)*DLEL**GTB2(JG)
        ELSE
          DLEL    =  ELIU-ELID
          QGT(JG) =  GTA2(JG)*DLEL**GTB2(JG)
        END IF
      ELSE
        QGT(JG) = 0.0
      END IF
    ELSE
      QGT(JG) = 0.0
    END IF
   endif
  END DO
END SUBROUTINE GATE_FLOW

!***********************************************************************************************************************************
!**                                         S U B R O U T I N E   S P I L L W A Y  F L O W                                        **
!***********************************************************************************************************************************

SUBROUTINE SPILLWAY_FLOW
  USE STRUCTURES; USE GLOBAL; USE GEOMC

  DO JS=1,NSP
    IF (LATERAL_SPILLWAY(JS)) THEN
       ELIU =  EL(KTWB(JWUSP(JS)),IUSP(JS))-Z(IUSP(JS))*COSA(BS(JWUSP(JS)))
    ELSE
       ELIU =  EL(KTWB(JWUSP(JS)),IUSP(JS))-Z(IUSP(JS))*COSA(BS(JWUSP(JS)))-SINA(JBUSP(JS))*DLX(IUSP(JS))*0.5
    END IF  
    IF (IDSP(JS) /= 0) THEN
      IF (US(JBDSP(JS)) /= IDSP(JS)) THEN
         ELID = EL(KTWB(JWDSP(JS)),IDSP(JS))-Z(IDSP(JS))*COSA(BS(JWDSP(JS)))
      ELSE   
         ELID = EL(KTWB(JWDSP(JS)),IDSP(JS))-Z(IDSP(JS))*COSA(BS(JWDSP(JS)))+SINA(JBDSP(JS))*DLX(IDSP(JS))*0.5
      END IF
    ELSE   
      ELID = -1.0
    END IF
    IF (ELID >= ESP(JS) .OR. ELIU >= ESP(JS)) THEN
      ISUB = 0
      IF (A2SP(JS) /= 0.0 .AND. IDSP(JS) /= 0) THEN
        HTAIL   =  ELID-ESP(JS)                                                 ! SW 5/10/05
        IF (HTAIL > 0) THEN
          HENERGY = (U(KTWB(JWUSP(JS)),IUSP(JS))**2)/(2.0*G)+ELIU-ESP(JS)       ! SW 5/10/05
          IF (HTAIL/HENERGY > 0.67) ISUB = 1
        END IF     
      END IF
      IF (ISUB == 0) THEN
        DLEL = ELIU-ESP(JS)
        IF (DLEL < 0.0) THEN
          DLEL    = -DLEL
          QSP(JS) = -A1SP(JS)*DLEL**B1SP(JS)
        ELSE
          QSP(JS) =  A1SP(JS)*DLEL**B1SP(JS)
        END IF
      ELSE IF (ELID > ELIU) THEN
        DLEL    =  ELID-ELIU
        QSP(JS) = -A2SP(JS)*DLEL**B2SP(JS)
      ELSE
        DLEL    =  ELIU-ELID
        QSP(JS) =  A2SP(JS)*DLEL**B2SP(JS)
      END IF
    ELSE
      QSP(JS) = 0.0
    END IF
  END DO
END SUBROUTINE SPILLWAY_FLOW

!***********************************************************************************************************************************
!**                                            S U B R O U T I N E   P I P E  F L O W                                             **
!***********************************************************************************************************************************

SUBROUTINE PIPE_FLOW_INITIALIZE
  USE GLOBAL; USE GEOMC; USE STRUCTURES
  REAL :: JDAY
  SAVE

  ALLOCATE (BEGIN(NPI), WLFLAG(NPI), VMAX(NPI))
  QOLD   =  0.01;  VMAX   =  0.01
  BEGIN  = .TRUE.; WLFLAG = .TRUE.
RETURN

ENTRY PIPE_FLOW (NIT,JDAY)
  DTQ = DLT/10.0
  DO JP=1,NPI
    DIA   = WPI(JP)
    CLEN  = DLXPI(JP)
    FMAN  = FPI(JP)
    CLOSS = FMINPI(JP)
    UPIE  = EUPI(JP)
    DNIE  = EDPI(JP)
    DLTX  = CLEN/(REAL(NC-1)*0.5)
    IF (LATERAL_PIPE(JP)) THEN   
      EL1   = EL(KTWB(JWUPI(JP)),IUPI(JP))-Z(IUPI(JP))*COSA(JBUPI(JP))
    ELSE
      EL1   = EL(KTWB(JWUPI(JP)),IUPI(JP))-Z(IUPI(JP))*COSA(JBUPI(JP))-SINA(JBDPI(JP))*DLX(IUPI(JP))*0.5
    END IF
    IF (IDPI(JP) /= 0) THEN
      IF (US(JBDPI(JP)) /= IDPI(JP)) THEN
        EL2   = EL(KTWB(JWDPI(JP)),IDPI(JP))-Z(IDPI(JP))*COSA(JBDPI(JP))
      ELSE
        EL2   = EL(KTWB(JWDPI(JP)),IDPI(JP))-Z(IDPI(JP))*COSA(JBDPI(JP))+SINA(JBDPI(JP))*DLX(IDPI(JP))*0.5
      END IF
    ELSE
      EL2 = -1.0
    END IF
    HIE = MAX(UPIE,DNIE)
    IF (DIA == 0.0) THEN
      QPI(JP)    =  0.0
      WLFLAG(JP) = .TRUE.
      GO TO 140
    END IF
    EPS = 0.001
    IF ((HIE+EPS) >= EL1 .AND. (HIE+EPS) >= EL2) THEN
      QPI(JP)    =  0.0
      WLFLAG(JP) = .TRUE.
      GO TO 140
    END IF
    IF (EL1 > EL2) THEN
      DCHECK = EL1-UPIE
    ELSE
      DCHECK = EL2-DNIE
    END IF
    IF (DCHECK < 0.02) THEN
      QPI(JP)    =  0.0
      WLFLAG(JP) = .TRUE.
      GO TO 140
    END IF
    IF (ABS(QOLD(JP)) < 0.001) QOLD(JP) = 0.001
    IF (EL1 >= (UPIE+DIA) .AND. EL2 >= (DNIE+DIA)) THEN
      D1 = EL1
      D2 = EL2
      GO TO 120
    END IF
    IF (EL1 > EL2) THEN
      DTEST = EL2-DNIE
    ELSE
      DTEST = EL1-UPIE
    END IF
    DCRIT = DEPTHCRIT(ABS(QOLD(JP)))
    IF (DTEST <= DCRIT) THEN
      IF (EL1 <= EL2) THEN
        D1 = UPIE+DCRIT
        D2 = EL2
      ELSE
        D1 = EL1
        D2 = DNIE+DCRIT
      END IF
      VTOT = 0.0
      TOTT = 0.0
110   CONTINUE
      IF (NIT /= 0) THEN
        DTQ = OMEGA*DLTX/VMAX(JP)
        IF (DTQ > (DLT-TOTT)) THEN
          DTQ = DLT-TOTT
        ELSE IF ((2.0*DTQ) > (DLT-TOTT)) THEN
          DTQ = (DLT-TOTT)*0.5
        END IF
      END IF
      CALL OPEN_CHANNEL (D1,D2,QPI(JP),JP,DTQ,JDAY)
      DCRIT = DEPTHCRIT(ABS(QPI(JP)))
      IF (EL1 <= EL2) THEN
        D1 = UPIE+DCRIT
      ELSE
        D2 = DNIE+DCRIT
      END IF
      VTOT = VTOT+DTQ*QPI(JP)
      TOTT = DTQ+TOTT
      IF (TOTT < (DLT-EPS2)) GO TO 110
      QPI(JP) = VTOT/DLT
      GO TO 140
    END IF
    D1 = EL1
    D2 = EL2
120 CONTINUE
    TOTT = 0.0
    VTOT = 0.0
130 CONTINUE
    IF (NIT /= 0) THEN
      DTQ = OMEGA*DLTX/VMAX(JP)
      IF (DTQ > (DLT-TOTT)) THEN
        DTQ = DLT-TOTT
      ELSE IF ((2.0*DTQ) > (DLT-TOTT)) THEN
        DTQ = (DLT-TOTT)*0.5
      END IF
    END IF
    CALL OPEN_CHANNEL (D1,D2,QPI(JP),JP,DTQ,JDAY)
    VTOT = VTOT+DTQ*QPI(JP)
    TOTT = DTQ+TOTT
    IF (TOTT < (DLT-EPS2)) GO TO 130
    QPI(JP) = VTOT/DLT
140 CONTINUE
    QOLD(JP) = QPI(JP)
    IF (QPI(JP) == 0.0) WLFLAG(JP) = .TRUE.
  END DO
RETURN
ENTRY DEALLOCATE_PIPE_FLOW
  DEALLOCATE (BEGIN, WLFLAG, VMAX)
RETURN
END SUBROUTINE PIPE_FLOW_INITIALIZE

!***********************************************************************************************************************************
!**                                        S U B R O U T I N E   O P E N  C H A N N E L                                           **
!***********************************************************************************************************************************

SUBROUTINE OPEN_CHANNEL_INITIALIZE
  USE GLOBAL; USE STRUCTURES
  REAL, PARAMETER :: THETA=0.55

! Type declarations

  REAL                                 :: JDAY
  REAL,    ALLOCATABLE, DIMENSION(:)   :: Y,   D,  B,   V,   CAREA, TOPW,  BELEV, Q, VOLD, YOLD
  REAL,    ALLOCATABLE, DIMENSION(:)   :: YT,  VT, VPR, YPR, TAREA, TOPWT, RT
  REAL,    ALLOCATABLE, DIMENSION(:,:) :: DAA, AL     
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: INDX
  LOGICAL                              :: SMOOTH_WATER_LEVELS, OPENWRN
  SAVE

! Allocation declarations

  ALLOCATE (Y(NN),    V(NN),     CAREA(NN),  TOPW(NN),   BELEV(NN),  Q(NN),     VOLD(NN), YOLD(NN), D(NN), B(NN))
  ALLOCATE (YT(NN),   VT(NN),    VPR(NN),    YPR(NN),    TAREA(NN),  TOPWT(NN), RT(NN),   INDX(NN))     
  ALLOCATE (AL(NN,2), DAA(NN,NN))     
RETURN

ENTRY OPEN_CHANNEL (EL1,EL2,QOUT,IC,DT,JDAY)

! Variable initializtion

  B     = 0.0; Y     = 0.0; V = 0.0; VT = 0.0; YT = 0.0; RT = 0.0; DAA = 0.0; YPR = 0.0; VPR = 0.0; TOPW = 0.0; TOPWT = 0.0 
  CAREA = 0.0; TAREA = 0.0
  BELEV(1)  = UPIE
  BELEV(NC) = DNIE
  PHI       = ASIN((UPIE-DNIE)/CLEN)
  DLTX      = CLEN/(REAL(NC-1)*0.5)
  DO J=2,NC-1
    DLTX2    =  DLTX*0.5
    SLOPE    = (UPIE-DNIE)/CLEN
    DIST     = (REAL(J-1)*DLTX2)
    BELEV(J) =  UPIE-SLOPE*DIST
  END DO
  BEPR1 =  UPIE+SLOPE*DLTX2
  BEPR2 =  DNIE-SLOPE*DLTX2
  BC1   = (EL1-BEPR1)*COS(PHI)
  IF (BC1 <= 0.0) BC1 = EL1-UPIE
  BC2 = (EL2-BEPR2)*COS(PHI)
  IF (BC2 <= 0.0) BC2 = EL2-DNIE
  IF (.NOT. BEGIN(IC)) THEN
    IF (WLFLAG(IC)) THEN
      DO J=2,NC-1,2
        WLSLOPE = ((BC1-BC2)/(CLEN+DLTX))*COS(PHI)
        DIST    = (REAL(J-1)*0.5*DLTX)+DLTX2
        Y(J)    =  BC1-WLSLOPE*DIST
        YT(J)   =  Y(J)
        DTP(IC) =  DT
      END DO
    ELSE
      DO I=2,NC-1,2
        Y(I)  = YS(I,IC)
        YT(I) = YST(I,IC)
      END DO
    END IF
  END IF
  DO I=1,NC,2
    V(I)  = VS(I,IC)
    VT(I) = VST(I,IC)
  END DO
  IF (BEGIN(IC)) THEN
    BEGIN(IC) = .FALSE.
    DO J=2,NC-1,2
      WLSLOPE = ((BC1-BC2)/(CLEN+DLTX))*COS(PHI)
      DIST    = (REAL(J-1)*0.5*DLTX)+DLTX2
      Y(J)    =  BC1-WLSLOPE*DIST
      YT(J)   =  Y(J)
      DTP(IC) =  DT
    END DO
    DO J=1,NC,2
      V(J)  = 0.0
      VT(J) = V(J)
    END DO
    OPENWRN = .TRUE.
  END IF
  SMOOTH_WATER_LEVELS = .FALSE.
  DO N=1,NC,2
    IF (N == NC) THEN
      BAR1 = BAREA(BC2,DIA)
      RAD1 = BAR1/WETPER(BC2,DIA)
    ELSE
      BAR1 = BAREA(Y(N+1),DIA)
      RAD1 = BAR1/WETPER(Y(N+1),DIA)
    END IF
    IF (N == 1) THEN
      BAR2 = BAREA(BC1,DIA)
      RAD2 = BAR2/WETPER(BC1,DIA)
    ELSE
      BAR2 = BAREA(Y(N-1),DIA)
      RAD2 = BAR2/WETPER(Y(N-1),DIA)
    END IF
    RT(N) = (RAD1+RAD2)*0.5
  END DO
  DO N=2,NC-1,2
    TAREA(N) = BAREA(Y(N),DIA)
    TOPWT(N) = TWIDTH(Y(N),DIA)
    CAREA(N) = BAREA(Y(N),DIA)
  END DO

! Projected water levels and velocities

  DO J=1,NC,2
    VPR(J) = V(J)+DT*(V(J)-VT(J))/DTP(IC)
  END DO
  DO J=2,NC-1,2
    YPR(J) = Y(J)+DT*(Y(J)-YT(J))/DTP(IC)
  END DO

! Matrix setup

  VTOT = 0.0
  DO J=1,NC,2
    VTOT = VTOT+V(J)
  END DO
  VAVG = VTOT/(REAL(NC-1)*0.5)

! Continuity

  DO N=2,NC-1,2
    VPR(N) = (VPR(N-1)+VPR(N+1))*0.5
    V(N)   = (V(N-1)+V(N+1))*0.5
    IF (N /= 2) DAA(N,N-2) = -THETA*(DT/DLTX)*(VPR(N)*0.5)
    DAA(N,N-1) = -THETA*(DT/DLTX)*(TAREA(N)/TOPWT(N))
    DAA(N,N)   =  1.0
    DAA(N,N+1) =  THETA*(DT/DLTX)*(TAREA(N)/TOPWT(N))
    IF (N /= NC-1) DAA(N,N+2) = THETA*(DT/DLTX)*(VPR(N)*0.5)
    IF (N == 2) THEN
      B(N) = Y(N)-(1.0-THETA)*(DT/DLTX)*(TAREA(N)/TOPWT(N))*(V(N+1)-V(N-1))-(1.0-THETA)*(DT/DLTX)*(V(N)*0.5)*(Y(N+2)-BC1)          &
             +THETA*(DT/DLTX)*(VPR(N)*0.5)*BC1
    ELSE IF (N == NC-1) THEN
      B(N) = Y(N)-(1.0-THETA)*(DT/DLTX)*(TAREA(N)/TOPWT(N))*(V(N+1)-V(N-1))-(1.0-THETA)*(DT/DLTX)*(V(N)*0.5)*(BC2-Y(N-2))          &
             -THETA*(DT/DLTX)*(VPR(N)*0.5)*BC2
    ELSE
      B(N) = Y(N)-(1.0-THETA)*(DT/DLTX)*(TAREA(N)/TOPWT(N))*(V(N+1)-V(N-1))-(1.0-THETA)*(DT/DLTX)*(V(N)*0.5)*(Y(N+2)-Y(N-2))
    END IF
  END DO
  IF (VAVG > 0.0 .OR. (VAVG == 0.0 .AND. EL1 > EL2)) THEN

!** Momentum 

    DO N=1,NC,2
      IF (N /= 1) THEN
        DAA(N,N-2) = -THETA*(DT/DLTX)*VPR(N)
        DAA(N,N-1) = -THETA*(DT/DLTX)*G*COS(PHI)
      END IF
      DAA(N,N) = 1.0+THETA*DT*G*(FMAN**2)*ABS(VPR(N))/(RT(N)**(4.0/3.0))+THETA*(DT/DLTX)*VPR(N)+THETA*(CLOSS*0.5)*(DT/CLEN)        &
                 *ABS(VPR(N))
      IF (N /= NC) DAA(N,N+1) = THETA*(DT/DLTX)*G*COS(PHI)
      IF (N == 1) THEN
        B(N) = V(N)-(1.0-THETA)*(DT/DLTX)*G*(Y(N+1)-BC1)*COS(PHI)-(1.0-THETA)*V(N)*(DT/DLTX)*V(N)-(1.0-THETA)*DT*G*(FMAN**2)       &
               /(RT(N)**(4.0/3.0))*V(N)*ABS(V(N))+DT*G*SIN(PHI)-(1.0-THETA)*(DT/CLEN)*(CLOSS*0.5)*V(N)*ABS(V(N))+THETA*(DT/DLTX)   &
               *G*COS(PHI)*BC1
      ELSE IF (N == NC) THEN
        B(N) = V(N)-(1.0-THETA)*(DT/DLTX)*G*(BC2-Y(N-1))*COS(PHI)-(1.0-THETA)*V(N)*(DT/DLTX)*(V(N)-V(N-2))-(1.0-THETA)             &
               *DT*G*(FMAN**2)/(RT(N)**(4.0/3.0))*V(N)*ABS(V(N))+DT*G*SIN(PHI)-(1.0-THETA)*(DT/CLEN)*(CLOSS*0.5)*V(N)*ABS(V(N))    &
               -THETA*(DT/DLTX)*G*COS(PHI)*BC2
      ELSE
        B(N) = V(N)-(1.0-THETA)*(DT/DLTX)*G*(Y(N+1)-Y(N-1))*COS(PHI)-(1.0-THETA)*V(N)*(DT/DLTX)*(V(N)-V(N-2))-(1.0-THETA)          &
               *DT*G*(FMAN**2)/(RT(N)**(4.0/3.0))*V(N)*ABS(V(N))+DT*G*SIN(PHI)-(1.0-THETA)*(DT/CLEN)*(CLOSS*0.5)*V(N)*ABS(V(N))
      END IF
    END DO
  ELSE
    DO N=1,NC,2
      IF (N /= NC) THEN
        DAA(N,N+2) = THETA*(DT/DLTX)*VPR(N)
        DAA(N,N+1) = THETA*(DT/DLTX)*G*COS(PHI)
      END IF
      DAA(N,N) = 1.0+THETA*DT*G*(FMAN**2)*ABS(VPR(N))/(RT(N)**(4.0/3.0))-THETA*(DT/DLTX)*VPR(N)+THETA*(CLOSS*0.5)*(DT/CLEN)        &
                 *ABS(VPR(N))
      IF (N /= 1) DAA(N,N-1) = -THETA*(DT/DLTX)*G*COS(PHI)
      IF (N == NC) THEN
        B(N) = V(N)-(1.0-THETA)*(DT/DLTX)*G*(BC2-Y(N-1))*COS(PHI)-(1.0-THETA)*V(N)*(DT/DLTX)*(-V(N))-(1.0-THETA)*DT*G*(FMAN**2)    &
               /(RT(N)**(4.0/3.0))*V(N)*ABS(V(N))+DT*G*SIN(PHI)-(1.0-THETA)*(DT/CLEN)*(CLOSS*0.5)*V(N)*ABS(V(N))-THETA*(DT/DLTX)   &
               *G*COS(PHI)*BC2
      ELSE IF (N == 1) THEN
        B(N) = V(N)-(1.0-THETA)*(DT/DLTX)*G*(Y(N+1)-BC1)*COS(PHI)-(1.0-THETA)*V(N)*(DT/DLTX)*(V(N+2)-V(N))-(1.0-THETA)             &
               *DT*G*(FMAN**2)/(RT(N)**(4.0/3.0))*V(N)*ABS(V(N))+DT*G*SIN(PHI)-(1.0-THETA)*(DT/CLEN)*(CLOSS*0.5)*V(N)*ABS(V(N))    &
               +THETA*(DT/DLTX)*G*COS(PHI)*BC1
      ELSE
        B(N) = V(N)-(1.0-THETA)*(DT/DLTX)*G*(Y(N+1)-Y(N-1))*COS(PHI)-(1.0-THETA)*V(N)*(DT/DLTX)*(V(N+2)-V(N))-(1.0-THETA)          &
               *DT*G*(FMAN**2)/(RT(N)**(4.0/3.0))*V(N)*ABS(V(N))+DT*G*SIN(PHI)-(1.0-THETA)*(DT/CLEN)*(CLOSS*0.5)*V(N)*ABS(V(N))
      END IF
    END DO
  END IF
  NP = NN
  CALL LUDCMP (DAA,NC,NP,INDX,D)
  CALL LUBKSB (DAA,NC,NP,INDX,B)
  DO I=2,NC-1,2
    YOLD(I)   = Y(I)
    YST(I,IC) = Y(I)
  END DO
  DO I=2,NC-1,2
    Y(I) = B(I)
  END DO

! Smooth water levels

  DO I=2,NC-1,2
    IF (Y(I) <= 0.0) THEN
      IF (OPENWRN) THEN
        OPEN (391,FILE='culvert.wrn',STATUS='unknown')
        OPENWRN = .FALSE.
      END IF
      SMOOTH_WATER_LEVELS = .TRUE.
    END IF
  END DO
  IF (SMOOTH_WATER_LEVELS) THEN
    DO J=2,NC-1,2
      WLSLOPE = ((BC1-BC2)/(CLEN+DLTX))*COS(PHI)
      DIST    = (REAL(J-1)*0.5*DLTX)+DLTX2
      Y(J)    =  BC1-WLSLOPE*DIST
    END DO
    WRITE (391,10010) IC, JDAY
    SMOOTH_WATER_LEVELS = .FALSE.
  END IF

! Flows

  NQCNT = 0
  QSUM  = 0.0
  DO I=1,NC,2
    VOLD(I)   = V(I)
    VST(I,IC) = V(I)
    V(I)      = B(I)
    IF (I == NC) THEN
      BAR1 = BAREA(BC2,DIA)
    ELSE
      BAR1 = BAREA(Y(I+1),DIA)
    END IF
    IF (I == 1) THEN
      BAR2 = BAREA(BC1,DIA)
    ELSE
      BAR2 = BAREA(Y(I-1),DIA)
    END IF
    CAREA(I) = (BAR1+BAR2)*0.5
    Q(I)     =  V(I)*CAREA(I)
    NQCNT    =  NQCNT+1
    QSUM     =  QSUM+Q(I)
  END DO
  QAVG = QSUM/REAL(NQCNT)
  DO I=2,NC-1,2
    YS(I,IC) = Y(I)
  END DO
  VMAX(IC) = 0.0
  DO I=1,NC,2
    VS(I,IC) = V(I)
    VMAX(IC) = MAX(ABS(V(I)),VMAX(IC))
  END DO
  DTP(IC)    =  DT
  QOUT       =  QAVG
  QOLD(IC)   =  QOUT
  WLFLAG(IC) = .FALSE.
10010 FORMAT ('water levels for culvert ',I3,' on Julian Day ',F10.3,' are <= 0 - predictions have been smoothed')
RETURN
ENTRY DEALLOCATE_OPEN_CHANNEL
  DEALLOCATE (Y, V, CAREA, TOPW, BELEV, Q, VOLD, YOLD, D, B, YT, VT, VPR, YPR, TAREA, TOPWT, RT, INDX, AL, DAA)
RETURN
END SUBROUTINE OPEN_CHANNEL_INITIALIZE

!***********************************************************************************************************************************
!**                                             S U B R O U T I N E   G R I D  A R E A 1                                          **
!***********************************************************************************************************************************

SUBROUTINE GRID_AREA1 (EL1,EL2,DIFF,BTOP)
  USE GLOBAL; USE GEOMC

! Difference in areas for trapezoidal geometry

  DO K=2,KB(I)
    IF (EL(K,I) <= EL1) THEN
      K1 = K
      EXIT
    END IF
  END DO
  DO K=2,KB(I)
    IF (EL(K,I) <= EL2) THEN
      K2 = K     
      EXIT
    END IF
  END DO
  BAREA1 = 0.0
  BAREA2 = 0.0
  DO K=KB(I),K1,-1   
    BAREA1 = BAREA1+BH(K,I)   
  END DO
  DIST = EL1-EL(K1,I)    
  IF (H(K1-1,JW)/2.0 < DIST) THEN     
    DIST1  =  H(K1-1,JW)*0.5
    SLPE   = (B(K1-1,I)-BB(K1-1,I))/(0.5*H(K1-1,JW))    
    BAREA1 =  BAREA1+BB(K1-1,I)*DIST1+0.5*SLPE*DIST1*DIST1     
    DIST2  =  DIST-H(K1-1,JW)*0.5    
    SLOPE  = (BB(K1-2,I)-B(K1-1,I))/(0.5*H(K1-1,JW))    
    BAREA1 =  BAREA1+B(K1-1,I)*DIST2+0.5*SLPE*DIST2*DIST2     
    BTOP   =  B(K1-1,I)+DIST2*SLPE
  ELSE
    SLPE   = (B(K1-1,I)-BB(K1-1,I))/(0.5*H(K1-1,JW))    
    BAREA1 =  BAREA1+BB(K1-1,I)*DIST+0.5*SLPE*DIST*DIST  
    BTOP   =  BB(K1-1,I)+DIST*SLPE
  END IF
  DO K=KB(I),K2,-1   
    BAREA2 = BAREA2+BH(K,I)   
  END DO
  DIST = EL2-EL(K2,I)   
  IF (H(K2-1,JW)/2. < DIST) THEN     
    DIST1  =  H(K2-1,JW)*0.5
    SLPE   = (B(K2-1,I)-BB(K2-1,I))/(0.5*H(K2-1,JW))    
    BAREA2 =  BAREA2+BB(K2-1,I)*DIST1+0.5*SLPE*DIST1*DIST1     
    DIST2  =  DIST-H(K2-1,JW)*0.5   
    SLPE   = (BB(K2-2,I)-B(K2-1,I))/(0.5*H(K2-1,JW))    
    BAREA2 =  BAREA2+B(K2-1,I)*DIST2+0.5*SLPE*DIST2*DIST2    
  ELSE
    SLPE   = (B(K2-1,I)-BB(K2-1,I))/(0.5*H(K2-1,JW))    
    BAREA2 =  BAREA2+BB(K2-1,I)*DIST+0.5*SLPE*DIST*DIST     
  END IF
  DIFF = BAREA1-BAREA2
  RETURN
END SUBROUTINE

!***********************************************************************************************************************************
!**                                             S U B R O U T I N E   G R I D  A R E A 2                                          **
!***********************************************************************************************************************************

SUBROUTINE GRID_AREA2
  USE GLOBAL; USE GEOMC; USE RSTART

  AREA   = (EL(KT,I)-SZ(I)-(EL(KT,I)-Z(I)))*BI(KT,I)  
  SL     = (B(KT,I)-BB(KT,I))/(0.5*H(KT,JW))
  A_COEF = -1.0
  B_COEF =  SZ(I)*2.+BI(KT,I)/(0.5*SL)
  C_COEF = -AREA/(0.5*SL)-SZ(I)**2-BI(KT,I)*2.*SZ(I)/SL
  Z(I)   = (-B_COEF+SQRT(B_COEF**2-4.*A_COEF*C_COEF))/(2.0*A_COEF)
  KTI(I) = 2   
  DO K=2,KB(I)
    IF (EL(K,I) <= EL(KT,I)-Z(I)) THEN
      KTI(I) = K-1
      EXIT
    END IF
  END DO
  RETURN
END SUBROUTINE

!***********************************************************************************************************************************
!**                                              S U B R O U T I N E   L U D C M P                                                **
!***********************************************************************************************************************************

SUBROUTINE LUDCMP (A,N,NP,INDX,D)
  REAL            :: A(NP,NP)
  REAL            :: VV(500)
  REAL, PARAMETER :: TINY=1.0E-20
  INTEGER         :: INDX(NP)

  D = 1.0
  DO I=1,N
    AAMAX = 0.0
    DO J=1,N
      IF (ABS(A(I,J)) > AAMAX) AAMAX = ABS(A(I,J))
    END DO
    VV(I) = 1.0/AAMAX
  END DO
  DO J=1,N
    DO I=1,J-1
      SUM = A(I,J)
      DO K=1,I-1
        SUM = SUM-A(I,K)*A(K,J)
      END DO
      A(I,J) = SUM
    END DO
    AAMAX = 0.0
    DO I=J,N
      SUM = A(I,J)
      DO K=1,J-1
        SUM = SUM-A(I,K)*A(K,J)
      END DO
      A(I,J) = SUM
      DUM = VV(I)*ABS(SUM)
      IF (DUM >= AAMAX) THEN
        IMAX  = I
        AAMAX = DUM
      END IF
    END DO
    IF (J /= IMAX) THEN
      DO K=1,N
        DUM       = A(IMAX,K)
        A(IMAX,K) = A(J,K)
        A(J,K)    = DUM
      END DO
      D        = -D
      VV(IMAX) =  VV(J)
    END IF
    INDX(J) = IMAX
    IF (A(J,J) == 0.0) A(J,J) = TINY
    IF (J /= N) THEN
      DUM = 1.0/A(J,J)
      DO I=J+1,N
        A(I,J) = A(I,J)*DUM
      END DO
    END IF
  END DO
END SUBROUTINE LUDCMP

!***********************************************************************************************************************************
!**                                              S U B R O U T I N E   L U B K S B                                                **
!***********************************************************************************************************************************

SUBROUTINE LUBKSB (A,N,NP,INDX,B)
  REAL    :: A(NP,NP), B(N)
  INTEGER :: N, NP, INDX(NP)
  INTEGER :: I, II, J, LL

  II = 0
  DO I=1,N
    LL    = INDX(I)
    SUM   = B(LL)
    B(LL) = B(I)
    IF (II /= 0) THEN
      DO J=II,I-1
        SUM = SUM-A(I,J)*B(J)
      END DO
    ELSE IF (SUM /= 0.0) THEN
      II = I
    END IF
    B(I) = SUM
  END DO
  DO I=N,1,-1
    SUM = B(I)
    DO J=I+1,N
      SUM = SUM-A(I,J)*B(J)
    END DO
    B(I) = SUM/A(I,I)
  END DO
END SUBROUTINE LUBKSB

!***********************************************************************************************************************************
!**                                                  F U N C T I O N   B A R E A                                                  **
!***********************************************************************************************************************************

FUNCTION BAREA (DEPTH,DIA)
  PARAMETER (PI=3.14159265359)
  IF (DEPTH < DIA) THEN
    BAREA = (DEPTH-DIA*0.5)*SQRT(DEPTH*DIA-DEPTH**2)+(DIA**2*0.25)*ASIN((2.0/DIA)*(DEPTH-DIA*0.5))+(PI*DIA**2)/8.0
  ELSE
    BAREA = (PI*DIA**2)*0.25
  END IF
END FUNCTION BAREA

!***********************************************************************************************************************************
!**                                                  F U N C T I O N   T W I D T H                                                **
!***********************************************************************************************************************************

FUNCTION TWIDTH (DEPTH,DIA)
  IF (DEPTH < DIA) THEN
    TWIDTH = 2.0*SQRT((DIA*DEPTH)-DEPTH**2)
  ELSE
    TWIDTH = 0.005*DIA
  END IF
END FUNCTION TWIDTH

!***********************************************************************************************************************************
!**                                                  F U N C T I O N   W E T P E R                                                **
!***********************************************************************************************************************************

FUNCTION WETPER (DEPTH,DIA)
  REAL, PARAMETER :: PI=3.14159265359
  IF (DEPTH < DIA) THEN
    WETPER = DIA*(ASIN((2.0/DIA)*(DEPTH-DIA*0.5))+PI*0.5)
  ELSE
    WETPER = PI*DIA
  END IF
END FUNCTION WETPER

!***********************************************************************************************************************************
!**                                                F U N C T I O N   D E P T H C R I T                                            **
!***********************************************************************************************************************************

FUNCTION DEPTHCRIT (FLOW)
  USE STRUCTURES
  EXTERNAL CDFUNC
  X1        = DIA/1.0E7
  X2        = DIA
  TOL       = 0.001
  DEPTHCRIT = ZBRENT1(CDFUNC,X1,X2,TOL,FLOW)
END FUNCTION DEPTHCRIT

!***********************************************************************************************************************************
!**                                                  F U N C T I O N   C D F U N C                                                **
!***********************************************************************************************************************************

FUNCTION CDFUNC (DEPTH,FLOW)
  USE STRUCTURES
  CDFUNC = (FLOW**2*TWIDTH(DEPTH,DIA))/(BAREA(DEPTH,DIA)**3*9.81)-1.0
END FUNCTION CDFUNC

!***********************************************************************************************************************************
!**                                                  F U N C T I O N   Z B R E N T                                                **
!***********************************************************************************************************************************

FUNCTION ZBRENT1 (FUNC,X1,X2,TOL,BARG)
  EXTERNAL   FUNC
  PARAMETER (FACTOR=0.1,NTRY=50,ITMAX=100,EPS=3.E-8)
  F1 = FUNC(X1,BARG)
  F2 = FUNC(X2,BARG)
  IF (F1 <= 0.0) THEN
    DO I=1,40
      X1 = X1/10.0
      F1 = FUNC(X1,BARG)
      IF (F1 > 0.0) EXIT
    END DO
  END IF
  DO J=1,NTRY
    IF (F1*F2 < 0.0) EXIT
    IF (ABS(F1) < ABS(F2)) THEN
      X1 = X1+FACTOR*(X1-X2)
      F1 = FUNC(X1,BARG)
    ELSE
      X2 = X2+FACTOR*(X2-X1)
      F2 = FUNC(X2,BARG)
    END IF
  END DO
  BA = X1
  B  = X2
  FA = FUNC(BA,BARG)
  FB = FUNC(B,BARG)
  FC = FB
  DO ITER=1,ITMAX
    IF (FB*FC > 0.0) THEN
      C  = BA
      FC = FA
      D  = B-BA
      E  = D
    END IF
    IF (ABS(FC) < ABS(FB)) THEN
      BA = B
      B  = C
      C  = BA
      FA = FB
      FB = FC
      FC = FA
    END IF
    TOL1 = 2.0*EPS*ABS(B)+0.5*TOL
    XM   = 0.5*(C-B)
    IF (ABS(XM) <= TOL1 .OR. FB == 0.0) THEN
      ZBRENT1 = B; EXIT
    END IF
    IF (ABS(E) >= TOL1 .AND. ABS(FA) > ABS(FB)) THEN
      S = FB/FA
      IF (BA == C) THEN
        P = 2.0*XM*S
        Q = 1.0-S
      ELSE
        Q =  FA/FC
        R =  FB/FC
        P =  S*(2.*XM*Q*(Q-R)-(B-BA)*(R-1.0))
        Q = (Q-1.0)*(R-1.0)*(S-1.0)
      END IF
      IF (P > 0.0) Q = -Q
      P = ABS(P)
      IF (2.0*P < MIN(3.0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
        E = D
        D = P/Q
      ELSE
        D = XM
        E = D
      END IF
    ELSE
      D = XM
      E = D
    END IF
    BA = B
    FA = FB
    IF (ABS(D) > TOL1) THEN
      B = B+D
    ELSE
      B = B+SIGN(TOL1,XM)
    END IF
    FB = FUNC(B,BARG)
  END DO
  ZBRENT1 = B
END FUNCTION ZBRENT1

!v3.5 start
!************************************************************************
!**               S U B R O U T I N E    porosity                      **
!************************************************************************

      SUBROUTINE porosity

      use geomc;use global;use macrophytec;use porosityc;use LOGICC


      use screenc

    if(nit.eq.0)then

 1040 FORMAT((8X,i8,3F8.0))

      DO Jw=1,Nwb
        KT = KTwb(Jw)
        DO JB=BS(Jw),BE(Jw)
          IU   = CUS(JB)
          ID   = DS(JB)
          DO I=IU,ID            
            DO K=2,KB(I)
              VOLi(K,I) = BH(K,I)*DLX(I)
            END DO
            VOLi(KT,I)    = BH2(KT,I)*DLX(I)
          end do
        end do
      end do

    end if

    do jb=1,nbr
      cosa(jb)=cos(alpha(jb))
    end do

!c  calculating # of macrophyte stems in each cell

    DO Jw=1,Nwb
      KT = KTwb(Jw)
      DO JB=BS(Jw),BE(Jw)
        IU   = CUS(JB)
        ID   = DS(JB)
        DO I=IU,ID
          HkTi  = H(KT,jw)-Z(I)
          if(kt.eq.kti(i))then
            volkti(i)=hkti*bic(kt,i)*dlx(i)
          else
            volkti(i) = bic(KTI(I),I)*(EL(KT,i)-EL(KTI(I)+1,i)-Z(I)*cosa(jb))/cosa(jb)*dlx(i)
          end if
          DO K=KTI(I)+1,KT
            volkti(I) = volkti(I)+voli(k,i)
          END DO

          do m=1,nmc
            vstemkt(i,m)=(mac(kt,i,m)*volkti(i))/dwv(m)    !cb 6/29/06
          end do
      
          do k=kt+1,kb(i)
            do m=1,nmc
              vstem(k,i,m)=(mac(k,i,m)*voli(k,i))/dwv(m)   !cb 6/29/06
            end do
          end do
        end do
      end do
    end do

    por=1.0
    DO Jw=1,Nwb
      KT = KTwb(Jw)
      DO JB=BS(Jw),BE(Jw)
        IU = cUS(JB)
        ID = DS(JB)

        do i=iu,id
          do k=kt,kb(i)
            if(k.eq.kt)then
              vstot=0.0
              do m=1,nmc
                  vstot=vstot+vstemkt(i,m)
              end do
              por(kt,i)=(VOLkti(I)-vstot)/volkti(i)
            else
              vstot=0.0
              do m=1,nmc
                vstot=vstot+vstem(k,i,m)
              end do
              por(k,i)=(VOLi(K,I)-vstot)/voli(k,i)
            end if
          end do
        end do

        do i=iu,id
          do k=kti(i),kb(i)
            if(k.le.kt)then
              b(k,i)=por(kt,i)*bic(k,i)
            else
              b(k,i)=por(k,i)*bic(k,i)
            end if

          end do
        end do

      end do
    end do



! Boundary widths

  DO JW=1,NWB
    KT = KTWB(JW)
    DO JB=BS(JW),BE(JW)
      IU = US(JB)
      ID = DS(JB)
      DO I=IU-1,ID+1
        B(1,I) = B(2,I)
        DO K=KB(I)+1,KMX
          B(K,I) = B(KB(I),I)
        END DO
      END DO
    END DO
  END DO
  DO JW=1,NWB
    KT = KTWB(JW)
    DO JB=BS(JW),BE(JW)
      IU    = US(JB)
      ID    = DS(JB)
      IEXIT = 0
      DO K=1,KMX-1
        B(K,IU-1) = B(K,IU)
        IF (UH_INTERNAL(JB) .OR. HEAD_FLOW(JB)) THEN   
          IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
            B(K,IU-1) = B(K,UHS(JB))
          ELSE
            ELR = EL(K,IU)+SINA(JB)*DLX(IU)*0.5
            ELL = EL(2,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
            IF (ELR >= ELL) THEN
              B(K,IU-1) = B(2,UHS(JB))
            ELSE
              DO KUP=2,KMX-1
                ELL1 = EL(KUP,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                ELL2 = EL(KUP+1,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                IF (ELL1 > ELR .AND. ELL2 <= ELR) THEN
                  IF (KUP > KB(UHS(JB))) THEN
                    KB(IU-1)    = K-1
                    KBMIN(IU-1) = MIN(KB(IU),KB(IU-1))
                    IEXIT       = 1
                    EXIT
                  END IF
                  ELR2 = EL(K+1,IU)+SINA(JB)*DLX(IU)*0.5
                  IF (ELR2 >= ELL2) THEN
                    B(K,IU-1) = B(KUP,UHS(JB)); EXIT
                  ELSE
                    K1 = KUP+1
                    IF (K1 > KMX) EXIT
                    B11 = 0.0
                    EL1 = ELR
                    EL2 = EL(K1,UHS(JB))-SINA(JBUH(JB))*DLX(IU)*0.5
                    DO WHILE (ELR2 <= EL2)
                      B11 = B11+(EL1-EL2)*B(K1-1,UHS(JB))
                      EL1 = EL2
                      K1  = K1+1
                      IF (K1 >= KMX+1 .OR. EL2 == ELR2) EXIT
                      EL2 = EL(K1,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                      IF (EL2 <= ELR2) EL2 = ELR2   
                    END DO
                    B(K,IU-1) = B11/H(K,JW); EXIT
                  END IF
                END IF
              END DO
              IF (EL(KMX,UHS(JB)) > EL(K,IU)) B(K,IU-1) = B(K-1,IU-1)
              IF (B(K,IU-1) == 0.0) B(K,IU-1) = B(K-1,IU-1)
              IF (IEXIT == 1) EXIT
            END IF
          END IF
        END IF
      END DO
      IEXIT = 0 
      DO K=1,KMX-1
        B(K,ID+1) = B(K,ID)
        IF (DH_INTERNAL(JB)) THEN
          IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
            B(K,ID+1) = B(K,DHS(JB))
          ELSE
            ELL = EL(K,ID)-SINA(JB)*DLX(ID)*0.5
            ELR = EL(2,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
            IF (ELL >= ELR) THEN
              B(K,ID+1) = B(2,DHS(JB))
            ELSE
              DO KDN=2,KMX-1
                ERR1 = EL(KDN,DHS(JB))  +SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                ERR2 = EL(KDN+1,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                IF (ERR1 >= ELL .AND. ERR2 < ELL) THEN
                  IF (KDN > KB(DHS(JB))) THEN
                    KB(ID+1)  = K-1
                    KBMIN(ID) = MIN(KB(ID),KB(ID+1))
                    IEXIT     = 1
                    EXIT
                  END IF
                  ELL2 = EL(K+1,ID)-SINA(JB)*DLX(ID)*0.5
                  IF (ELL2 >= ERR2) THEN
                    B(K,ID+1) = B(KDN,DHS(JB)); EXIT
                  ELSE
                    K1  = KDN+1
                    IF (K1 > KMX) EXIT
                    B11 = 0.0
                    EL2 = ELL
                    EL1 = EL(K1,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                    DO WHILE (ELL2 <= EL1)      
                      B11 = B11+(EL2-EL1)*B(K1-1,DHS(JB))
                      EL2 = EL1
                      K1  = K1+1
                      IF (K1 >= KMX+1 .OR. EL1 == ELL2) EXIT
                      EL1 = EL(K1,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                      IF (EL1 <= ELL2) EL1 = ELL2   
                    END DO
                    B(K,ID+1) = B11/H(K,JW); EXIT
                  END IF
                END IF
              END DO
              IF (EL(KMX,DHS(JB)) > EL(K,ID)) B(K,ID+1) = B(K-1,ID+1)
              IF (B(K,ID+1) == 0.0) B(K,ID+1) = B(K-1,ID+1)
              IF (IEXIT == 1) EXIT
            END IF
          END IF
        END IF
      END DO

!**** Areas and bottom widths

      IF (.NOT. TRAPEZOIDAL(JW)) THEN                                                                                  !SW 07/16/04
        DO I=IU-1,ID+1
          DO K=1,KMX-1
            BH2(K,I) = B(K,I)*H(K,JW)
            BH(K,I)  = B(K,I)*H(K,JW)
            BB(K,I)  = B(K,I)-(B(K,I)-B(K+1,I))/(0.5*(H(K,JW)+H(K+1,JW)))*H(K,JW)*0.5                                  !SW 08/02/04
          END DO
          BH(KMX,I) = BH(KMX-1,I)
        END DO
!****** Derived geometry

        DO I=IU-1,ID+1
          BH2(KT,I) = B(KTI(I),I)*(EL(KT,I)-EL(KTI(I)+1,I)-Z(I)*COSA(JB))/COSA(JB)
          IF (KT == KTI(I)) BH2(KT,I) = H2(KT,I)*B(KT,I)
          DO K=KTI(I)+1,KT
            BH2(KT,I) = BH2(KT,I)+BH(K,I)
          END DO
          BKT(I)   = BH2(KT,I)/H2(KT,I)
          BI(KT,I) = B(KTI(I),I)
        END DO
      ELSE                                                                                                             !SW 07/16/04
        DO I=IU-1,ID+1
          DO K=1,KMX-1
            BB(K,I)  = B(K,I)-(B(K,I)-B(K+1,I))/(0.5*(H(K,JW)+H(K+1,JW)))*H(K,JW)*0.5 
          END DO
          BB(KB(I),I) = B(KB(I),I)*0.5
          BH2(1,I)    = B(1,I)*H(1,JW)  
          BH(1,I)     = BH2(1,I)   
          DO K=2,KMX-1
            BH2(K,I) = 0.25*H(K,JW)*(BB(K-1,I)+2.*B(K,I)+BB(K,I))
            BH(K,I)  = BH2(K,I)
          END DO
          BH(KMX,I) = BH(KMX-1,I)
        END DO
        DO I=IU-1,ID+1
          ELT = EL(KT,I)-Z(I)
          CALL GRID_AREA1(ELT,EL(KT+1,I),BH2(KT,I),BI(KT,I))
          BKT(I) = BH2(KT,I)/H2(KT,I)
        END DO
      END IF
      DO I=IU-1,ID
        DO K=1,KMX-1
          AVH2(K,I) = (H2(K,I) +H2(K+1,I)) *0.5
          AVHR(K,I) =  H2(K,I)+(H2(K,I+1)-H2(K,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                  !SW 07/29/04
        END DO
        AVH2(KMX,I) = H2(KMX,I)
        DO K=1,KMX
          BR(K,I)   = B(K,I)  +(B(K,I+1)  -B(K,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                !SW 07/29/04
          BHR(K,I)  = BH(K,I) +(BH(K,I+1) -BH(K,I)) /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                !SW 07/29/04
          BHR2(K,I) = BH2(K,I)+(BH2(K,I+1)-BH2(K,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                !SW 07/29/04
        END DO
      END DO
      DO K=1,KMX-1  
        AVH2(K,ID+1) = (H2(K,ID+1)+H2(K+1,ID+1))*0.5
        BR(K,ID+1)   =   B(K,ID+1)
        BHR(K,ID+1)  =   BH(K,ID+1)
      END DO
      AVH2(KMX,ID+1) = H2(KMX,ID+1)
      AVHR(KT,ID+1)  = H2(KT,ID+1)
      BHR2(KT,ID+1)  = BH2(KT,ID+1)
      IUT = IU
      IF (UP_HEAD(JB)) IUT = IU-1
      DO I=IUT,ID
        DO K=1,KMX-1
          VOL(K,I) = B(K,I)*H2(K,I)*DLX(I)
        END DO
        VOL(KT,I)    = BH2(KT,I)*DLX(I)
        DEPTHB(KT,I) = H2(KT,I)
        DEPTHM(KT,I) = H2(KT,I)*0.5
        DO K=KT+1,KMX
          DEPTHB(K,I) = DEPTHB(K-1,I)+ H2(K,I)
          DEPTHM(K,I) = DEPTHM(K-1,I)+(H2(K-1,I)+H2(K,I))*0.5
        END DO
      END DO
    END DO
  END DO

10 continue

    return
    end

!************************************************************************
!**               S U B R O U T I N E    macrophyte_friction           **
!************************************************************************

      SUBROUTINE macrophyte_friction(hrad,bedfr,effric,k,ii)

    use geomc;use global;use macrophytec;use porosityc  
              

  do m=1,nmc
    savolrat=dwv(m)/dwsa(m)     !cb 6/29/2006
    if(k.eq.kt)then
!      sarea(m)=vstemkt(ii,m)*savolrat/pi
      sarea(m)=vstemkt(ii,m)*savolrat*anorm(m)     !cb 6/29/2006
    else
!      sarea(m)=vstem(k,ii,m)*savolrat/pi
      sarea(m)=vstem(k,ii,m)*savolrat*anorm(m)     !cb 6/29/2006
    end if
  end do
  xsarea=bh2(k,ii)

  tsarea=0.0
  artot=0.0
  sctot=0.0
  do m=1,nmc
    artot=artot+sarea(m)
    sctot=sctot+cddrag(m)*sarea(m)
    tsarea=tsarea+sarea(m)
  end do

  if(artot.gt.0.0)then
    cdavg=sctot/artot
    frin=cdavg*tsarea*hrad**(4./3.)/(2.0*g*xsarea*dlx(ii)*bedfr**2)     
    effric=bedfr*sqrt(1.0+frin)
  else
    effric=bedfr
  end if

  return
  end

! v3.5
