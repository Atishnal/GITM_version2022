! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine calc_chemistry(iBlock)

  use ModSizeGitm
  use ModGITM
  use ModPlanet
  use ModRates
  use ModEUV
  use ModSources
  use ModInputs, only: iDebugLevel, UseIonChemistry, UseNeutralChemistry, f107, DoCheckForNans, UseReactionRatePerturbations
  use ModConstants
  use ModReactionRatePerturb !! Atishnal to use module ModReactionRatePerturb
  
  implicit none

  integer, intent(in) :: iBlock

  real :: IonSources(nIons), NeutralSources(nSpeciesTotal)
  real :: IonLosses(nIons), NeutralLosses(nSpeciesTotal)
  real :: DtSub, DtOld, DtTotal, DtMin, DtAve, Source, Reaction, tr, tr3, rr
  real :: te3, ti, tn, tn1, tn06, dtsubtmp, losstmp, dentmp, l, t,m1,m2,y1,y2,k1,k2
  real :: Ions(nIons), Neutrals(nSpeciesTotal)
  real :: tli(nIons), tsi(nIons), tln(nSpeciesTotal), tsn(nSpeciesTotal)
  real :: szap

  integer :: iLon, iLat, iAlt, iIon, nIters, iNeutral
  
  real :: lon, lat, alt
  real :: LonDeg, LatDeg, AltKm, LST
  real,dimension(7)    :: AP  

  real :: ChemicalHeatingSub,percent,o2ptotal
  real :: ChemicalHeatingSubI, ChemicalHeatingSubE
  real :: Emission(nEmissions), EmissionTotal(nEmissions)
  real, dimension(nLons,nLats,nAlts) :: &
       tr3d, tr33d,te12d, tr3m0443d, tr3m083d, tr3m043d, &
       tn3d, ti3d, ti33d, ti103d,ti93d, ti153d, ti3m0443d, ti3m0243d, &
       ti3m0233d, ti3m1163d, ti10m0673d, ti3m0453d, ti10m2123d, ti3m0873d, &
       ti3m0523d, ti9m0923d, ti3m0553d, ti15m0233d, &
       te3m0393d, te3m0853d, te33d, te3m053d,te3m073d,&
       te12m0563d,te227d, te3m0913d, te3m0813d, te073d

  real, dimension(nLons,nLats,nAlts) :: &
       teffective_n2, teffective_o2, teffective_no, u2, mb, mbb, &
       k1_n2, k2_o2, k3_no

  real :: k1_n2_point, k2_o2_point, k3_no_point

  real :: te3m05,te3m07,te12m056, tr3m044, tr3m04, tr3m08, te07
  real :: ti3m044, ti3m024, ti3m023, ti3m116, ti10m067, ti3m045
  real :: ti10m212, ti3m087, ti3m052, ti9m092, ti3m055, ti15m023
  real :: te3m039, te3m085, rr_opn2, te22m05, te3m091, te3m081
  real :: ionso, ionlo, neuso, neulo

  logical :: UseNeutralConstituent(nSpeciesTotal)
  logical :: UseIonConstituent(nIons)
  
  !---------------------------------------------------------------------------
  !! Atishnal passing perturbed reaction rates to GITM chemical scheme
  if (UseReactionRatePerturbations) then
         !call get_reaction_rate('o+o+m=>o2+m+5p12ev', rr_o_p_o_p_m__o2_p_m_5p12ev)  ! R1
         call get_reaction_rate('o+2d+n2=>n2++o+1p35ev', rr_op2d_p_n2__n2p_p_1p35ev)  ! R2
         !call get_reaction_rate('o+2p+n2=>n2++o+3p05ev', rr_op2p_p_n2__n2p_p_o_p_3p05ev)  ! R3
         call get_reaction_rate('n2++o2=>o2++n2+3p53ev_lt', rr_n2p_p_o2__o2p_p_n2_p_3p53ev_lt)  ! R4
         call get_reaction_rate('n2++o2=>o2++n2+3p53ev_gt', rr_n2p_p_o2__o2p_p_n2_p_3p53ev_gt)  ! R5
         call get_reaction_rate('n2++o=>no++n2d+op70ev=>no++n4s+3p08ev_lt', rr_n2p_p_o__nop_p_n2d_p_0p70ev__nop_p_n4s_p_3p08ev_lt)  ! R6
         call get_reaction_rate('n2++o=>no++n2d+op70ev=>no++n4s+3p08ev_gt', rr_n2p_p_o__nop_p_n2d_p_0p70ev__nop_p_n4s_p_3p08ev_gt)  ! R7
         call get_reaction_rate('n2++e=>2n2d+1p04ev=>2n4s+5p77ev', rr_n2p_p_e__2n2d_p_1p04ev__2n4s_p_5p77ev)  ! R8
         !call get_reaction_rate('n2++n4s=>n2+n+2p48ev', rr_n2p_p_n4s__n2_p_np_p_2p48ev)  ! R9
         !call get_reaction_rate('n2++o=>op4s+n2+1p96ev', rr_n2p_p_o__op4s_p_n2_p_1p96ev)  ! R10
         call get_reaction_rate('n2++no=>no++n2+6p33', rr_n2p_p_no__nop_p_p_n2_p_6p33ev)  ! R11
         !call get_reaction_rate('o+4s+o2=>o2++o+1p55ev_lt', rr_op4s_p_o2__o2p_p_o_p_1p55ev_lt)  ! R12
         !call get_reaction_rate('o+4s+o2=>o2++o+1p55ev_gt', rr_op4s_p_o2__o2p_p_o_p_1p55ev_gt)  ! R13
         !call get_reaction_rate('o+2d+o2=>o2++o+4p865ev', rr_op2d_p_o2__o2p_p_o_p_4p865ev)  ! R14
         !call get_reaction_rate('o+2p+o2=>o2++o+6p54ev', rr_op2p_p_o2__o2p_p_o_p_6p54ev)  ! R15
         !call get_reaction_rate('o+2p+o2=>o+4s+o2+5p016ev', rr_op2p_p_o2__op4s_p_o2_p_5p016ev)  ! R16
         !call get_reaction_rate('n++o2=>o2++n4s+2p5ev_lt', rr_np_p_o2__o2p_p_n4s_p_2p5ev_lt)  ! R17
         !call get_reaction_rate('n++o2=>o2++n4s+2p5ev_gt', rr_np_p_o2__o2p_p_n4s_p_2p5ev_gt)  ! R18
         !call get_reaction_rate('n++o2=>o2++n2d+0p1ev_lt', rr_np_p_o2__o2p_p_n2d_p_0p1ev_lt)  ! R19
         !call get_reaction_rate('n++o2=>o2++n2d+0p1ev_gt', rr_np_p_o2__o2p_p_n2d_p_0p1ev_gt)  ! R20
         !call get_reaction_rate('o2++e=>o1d+o1d+3p06ev_lt', rr_o2p_p_e__o1d_p_o1d_p_3p06ev_lt)  ! R21
         !call get_reaction_rate('o2++e=>o1d+o1d+3p06ev_gt', rr_o2p_p_e__o1d_p_o1d_p_3p06ev_gt)  ! R22
         !call get_reaction_rate('o2++n4s=>no++o+4p21ev', rr_o2p_p_n4s__nop_p_o_4p21ev)  ! R23
         !call get_reaction_rate('o2++n2d=>no++o+6p519ev', rr_o2p_p_n2d__nop_p_o_p_6p519ev)  ! R24
         !call get_reaction_rate('o2++n2p=>o2++n4s+3p565ev', rr_o2p_p_n2p__o2p_p_n4s_p_3p565ev)  ! R25
         call get_reaction_rate('o2++no=>no++o2+2p813ev', rr_o2p_p_no_nop_p_o2_p_2p813ev)  ! R26
         !call get_reaction_rate('o+2d+o=>o+4s+o3p+3p31ev=>o+4s+o1d+1p35ev', rr_op2d_p_o_op4s_p_o3p_p_3p31ev__op4s_p_o1d_p_1p35ev)  ! R27
         !call get_reaction_rate('o+2d+e=>o+4s+e+3p31ev', rr_op2d_p_e__op4s_p_e_p_3p31ev)  ! R28
         !call get_reaction_rate('o+2p+o=>o+4s+o+5p0ev', rr_op2p_p_o__op4s_p_o_p_5p0ev)  ! R29
         !call get_reaction_rate('o+2p+e=>o+4s+e+5p0ev', rr_op2p_p_e__op4s_p_e_p_5p0ev)  ! R30
         !call get_reaction_rate('o+2p=>o+4s+2470A', rr_op2p__op4s_p_2470a)  ! R31
         call get_reaction_rate('n++o2=>o+4s+no+2p31ev_lt', rr_np_p_o2__op4s_p_no_p_2p31ev_lt)  ! R32
         call get_reaction_rate('n++o2=>o+4s+no+2p31ev_gt', rr_np_p_o2__op4s_p_no_p_2p31ev_gt)  ! R33
         call get_reaction_rate('o+4s+n2=>no++n4s+1p10ev_lt', rr_op4s_p_n2__nop_p_n4S_p_1p10ev_lt)  ! R34
         call get_reaction_rate('o+4s+n2=>no++n4s+1p10ev_gt', rr_op4s_p_n2__nop_p_n4S_p_1p10ev_gt)  ! R35
         call get_reaction_rate('o+4s+no=>no++o+4p36ev', rr_op4s_p_no__nop_p_o_p_4p36ev)  ! R36
         call get_reaction_rate('o+4s+n2d=>n++o+1p45ev', rr_op4s_p_n2d__np_p_o_p_1p45ev)  ! R37
         !call get_reaction_rate('o+2p+e=>o+2d+e+1p69ev', rr_op2p_p_e__op2d_p_e_p_1p69ev)  ! R38
         !call get_reaction_rate('o+2d+n2=>no++n+4p41ev', rr_op2d_p_n2__nop_p_n_p_4p41ev)  ! R39
         call get_reaction_rate('o+2d+no=>no++o+4p37ev', rr_op2d_p_no__nop_p_o_p_4p37ev)  ! R40
         !call get_reaction_rate('o+2p=>o+2d+7320a', rr_op2p__op2d_p_7320a)  ! R41
         !call get_reaction_rate('o+2d=>o+4s+3726a', rr_op2d__op4s_p_3726a)  ! R42
         !call get_reaction_rate('he++n2=>n++n+he+0p28ev', rr_hep_p_n2__np_p_n_p_he_p_0p28ev)  ! R43
         !call get_reaction_rate('he++n2=>n2++he', rr_hep_p_n2__n2p_p_he)  ! R44
         !call get_reaction_rate('he++o2=>o+o+he', rr_hep_p_o2__o_p_o_p_he)  ! R45
         !call get_reaction_rate('he++e-=>he', rr_hep_p_em__he)  ! R46
         !call get_reaction_rate('o+2p+n=>n++o+2p7ev', rr_op2p_p_n__np_p_o_p_2p7ev)  ! R47
         call get_reaction_rate('n++no=>n2++o+2p2ev', rr_np_p_no__n2p_p_o_p_2p2ev)  ! R48
         call get_reaction_rate('n++no=>no++n4s+3p4ev', rr_np_p_no__nop_p_n4s_p_3p4ev)  ! R49
         !call get_reaction_rate('n++o2=>no++o3p+6p67ev_lt', rr_np_p_o2__nop_p_o3p_p_6p67ev_lt)  ! R50
         !call get_reaction_rate('n++o2=>no++o3p+6p67ev_gt', rr_np_p_o2__nop_p_o3p_p_6p67ev_gt)  ! R51
         !call get_reaction_rate('o+2d+n=>n++o+1p0ev', rr_op2d_p_n__np_p_o_p_1p0ev)  ! R52
         call get_reaction_rate('n++o2=>no++o1d+4p71ev_lt', rr_np_p_o2__nop_p_o1d_p_4p71ev_lt)  ! R53
         call get_reaction_rate('n++o2=>no++o1d+4p71ev_gt', rr_np_p_o2__nop_p_o1d_p_4p71ev_gt)  ! R54
         !call get_reaction_rate('n++o=>o++n+0p93ev', rr_np_p_o__op_p_n_p_0p93ev)  ! R55
         call get_reaction_rate('no++e=>o+n2d+0p38ev', rr_nop_p_e__o_p_n2d_p_0p38ev)  ! R56
         call get_reaction_rate('no++e=>o+n4s+2p77ev', rr_nop_p_e__o_p_n4s_p_2p77ev)  ! R57
         call get_reaction_rate('n2d+e=>n4s+e+2p38ev', rr_n2d_p_e__n4s_p_e_p_2p38ev)  ! R58
         call get_reaction_rate('n2d+o=>n4s+o3p+2p38ev=>n4s+o1d+0p42ev', rr_n2d_p_o__n4s_p_o3p_p_2p38ev__n4s_p_o1d_p_0p42ev)  ! R59
         !call get_reaction_rate('n2d=>n4s+5200a', rr_n2d__n4s_p_5200a)  ! R60
         call get_reaction_rate('no=>n4s+o', rr_no__n4s_p_o)  ! R61
         call get_reaction_rate('n4s+o2=>no+o+1p385ev', rr_n4s_p_o2__no_p_o_p_1p385ev)  ! R62
         call get_reaction_rate('n4s+no=>n2+o+3p25ev', rr_n4s_p_no__n2_p_o_p_3p25ev)  ! R63
         !call get_reaction_rate('n2p=>n2d+10400a', rr_n2p__n2d_p_10400a)  ! R64
         call get_reaction_rate('n2d+o2=>no+o3p+3p76ev=>no+o1d+1p80ev', rr_n2d_p_o2__no_p_o3p_p_3p76ev__no_p_o1d_p_1p80ev)  ! R65
         call get_reaction_rate('n2d+no=>n2+o+5p63ev', rr_n2d_p_no__n2_p_o_p_5p63ev)  ! R66
         !call get_reaction_rate('o1d=>o3p+6300a', rr_o1d__o3p_p_6300a)  ! R67
         !call get_reaction_rate('o1d=>o3p+6364a', rr_o1d__o3p_p_6364a)  ! R68
         !call get_reaction_rate('o1d+e=>o3p+e+1p96ev', rr_o1d_p_e__o3p_p_e_p_1p96ev)  ! R69
         !call get_reaction_rate('o1d+n2=>o3p+n2+1p96ev', rr_o1d_p_n2__o3p_p_n2_p_1p96ev)  ! R70
         !call get_reaction_rate('o1d+o2=>o3p+o2+1p96ev', rr_o1d_p_o2__o3p_p_o2_p_1p96ev)  ! R71
         !call get_reaction_rate('o1d+o3p=>o3p+o3p+1p96ev', rr_o1d_p_o3p__o3p_p_o3p_p_1p96ev)  ! R72
         call get_reaction_rate('no=>no++e', rr_no__nop_p_e)  ! R73
  endif
  
  !! End Atishnal

  if (iDebugLevel > 3) then
     do iIon = 1, nIons
        write(*,*) "====> start calc_chemistry: Max Ion Density: ", iIon, &
             maxval(IDensityS(1:nLons,1:nLats,(nAlts*4)/5,iIon,iBlock))
     enddo
  endif

  ChemicalHeatingRate = 0.0
  ChemicalHeatingRateIon = 0.0
  ChemicalHeatingRateEle = 0.0
  ChemicalHeatingSpecies = 0.0
 
  UseNeutralConstituent = .true.
  UseIonConstituent     = .true.

!  UseNeutralConstituent(iO_1D_) = .false.
!  UseIonConstituent(iO_2PP_) = .false.
!  UseIonConstituent(iO_2DP_) = .false.
!  
!  UseNeutralConstituent(iN_4S_) = .false.
!  UseNeutralConstituent(iN_2D_) = .false.
!  UseNeutralConstituent(iO2_) = .false.
!
!  UseNeutralConstituent = .false.
!  UseIonConstituent = .false.
!  UseIonConstituent(1) = .true.
!  UseIonConstituent(2) = .true.
!  UseIonConstituent(3) = .true.
!
!  UseNeutralConstituent(1) = .true.
!  UseNeutralConstituent(2) = .true.
!  UseNeutralConstituent(3) = .true.

!  open(unit=95,file='data.dat')
  DtMin = Dt

  Emissions(:,:,:,:,iBlock) = 0.0
  
  if (.not.UseIonChemistry) return

  call report("Chemistry",2)
  call start_timing("calc_chemistry")

  DtAve = 0.0

  nIters=0

!  AuroralIonRateS = 0.0

  if (DoCheckForNans) then
     call check_for_nans_ions('before chemistry')
     call check_for_nans_neutrals('before chemistry')
     call check_for_nans_temps('before chemistry')
  endif

  u2 = IVelocityPerp(1:nLons,1:nLats,1:nAlts,iEast_,iBlock)**2 + &
       IVelocityPerp(1:nLons,1:nLats,1:nAlts,iNorth_,iBlock)**2 + &
       IVelocityPerp(1:nLons,1:nLats,1:nAlts,iUp_,iBlock)**2

  mb  = 0.0
  mbb = 0.0

  ! This is from Schunk and Nagy 2nd Ed, formula 12.13 (pg 416)
  do iNeutral = 1, nSpeciesTotal
     ! Collisions should be better defined
     mbb = mbb + &
          (Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)) / &
          (mass(iNeutral) + MassI(iO_4SP_))
     mb  = mb + &
          (mass(iNeutral) * Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)) / &
          (mass(iNeutral) + MassI(iO_4SP_))
  enddo

  mb = mb/mbb

  teffective_n2 = iTemperature(1:nLons,1:nLats,1:nAlts,iBlock) + &
       MassI(iO_4SP_)/(MassI(iO_4SP_) + Mass(iN2_)) * &
       (Mass(iN2_) - mb)/(3*Boltzmanns_Constant) * u2
       
  teffective_o2 = iTemperature(1:nLons,1:nLats,1:nAlts,iBlock) + &
       MassI(iO_4SP_)/(MassI(iO_4SP_) + Mass(iO2_)) * &
       (Mass(iO2_) - mb)/(3*Boltzmanns_Constant) * u2
       
  teffective_no = iTemperature(1:nLons,1:nLats,1:nAlts,iBlock) + &
       MassI(iO_4SP_)/(MassI(iO_4SP_) + Mass(iNO_)) * &
       (Mass(iNO_) - mb)/(3*Boltzmanns_Constant) * u2

  where (teffective_n2 < 350.0)
     teffective_n2 = 350.0
  endwhere

  where (teffective_n2 > 6000.0)
     teffective_n2 = 6000.0
  endwhere

  where (teffective_o2 < 350.0)
     teffective_o2 = 350.0
  endwhere

  where (teffective_o2 > 6000.0)
     teffective_o2 = 6000.0
  endwhere

  where (teffective_no < 320.0)
     teffective_no = 320.0
  endwhere

  where (teffective_no > 6000.0)
     teffective_no = 6000.0
  endwhere

  where (teffective_n2 <= 1700.0)
     k1_n2 =   1.533e-18 &
          - 5.920e-19*(teffective_n2/300.0) &
          + 8.600e-20*(teffective_n2/300.0)**2
  endwhere
  
  where (teffective_n2 > 1700.0)
     k1_n2 =   2.730e-18 &
          - 1.155e-18*(teffective_n2/300.0) &
          + 1.483e-19*(teffective_n2/300.0)**2
  endwhere
  
  k2_o2 =   2.820e-17 &
       - 7.740e-18*(teffective_o2/300.0) &
       + 1.073e-18*(teffective_o2/300.0)**2 &
       - 5.170e-20*(teffective_o2/300.0)**3 &
       + 9.650e-22*(teffective_o2/300.0)**4
  
  where (teffective_no <= 1500.0)
     k3_no =   8.360e-19 &
          - 2.020e-19*(teffective_no/300.0) &
          + 6.950e-20*(teffective_no/300.0)**2
  endwhere
  
  where (teffective_no > 1500.0)
     k3_no =   5.330e-19 &
          - 1.640e-20*(teffective_no/300.0) &
          + 4.720e-20*(teffective_no/300.0)**2 &
          - 7.050e-22*(teffective_no/300.0)**3
  endwhere
  
   tr3d = (iTemperature(1:nLons,1:nLats,1:nAlts,iBlock) &
          + Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
          TempUnit(1:nLons,1:nLats,1:nAlts)) / 2.0

  tn3d = Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
                TempUnit(1:nLons,1:nLats,1:nAlts)
				
  ti3d = 	iTemperature(1:nLons,1:nLats,1:nAlts,iBlock)			
		 
		 
  ti33d = ti3d/300.0
  ti103d = ti3d/1000.0
  ti93d = ti3d/900.0
  ti153d = ti3d/1500.0
  te33d = eTemperature(1:nLons,1:nLats,1:nAlts,iBlock)/300.0
  te12d = eTemperature(1:nLons,1:nLats,1:nAlts,iBlock)/1200.0
  te227d = -22740.0/eTemperature(1:nLons,1:nLats,1:nAlts,iBlock)
  te073d = (250.0/eTemperature(1:nLons,1:nLats,1:nAlts,iBlock))**0.7

  te3m073d    = te33d**(-0.7)
  te12m0563d  = te12d**(-0.56)
  te3m053d    = te33d**(-0.5)
  te3m0393d   = te33d**(-0.39)
  te3m0853d   = te33d**(-0.85)
  te3m0913d   = te33d**(0.91)
  te3m0813d   = te33d**(0.81)
  ti3m0443d   = ti33d**(-0.44)
  ti3m0243d   = ti33d**(-0.24)
  ti3m0233d   = ti33d**(-0.23)
  ti3m1163d   = ti33d**(-1.16)
  ti10m0673d  = ti103d**(0.67)
  ti3m0453d   = ti33d**(-0.45)
  ti10m2123d  = ti103d**(2.12)
  ti3m0873d   = ti33d**(0.87)
  ti3m0523d   = ti33d**(-0.52)
  ti9m0923d   = ti93d**(0.92)
  ti3m0553d   = ti33d**(0.55)
  ti15m0233d  = ti153d**(0.2)

  m1 = ALOG(1.0/100000.0)/(115.0-100.0)
  k1 = 100000.0*exp(-m1*100.0)
  m2 = ALOG(1.0/1000.0)/(180.0-100.0)
  k2 = 1000.0*exp(-m2*100.0)


  do iLon = 1, nLons
     do iLat = 1, nLats

        szap = cos(sza(iLon, iLat,iBlock))
        if (szap < 0.0) szap = 0.0

        ChemicalHeating2d(iLon, iLat) = 0.0
        
        do iAlt = 1, nAlts

           y1 = max(1.0,k1*exp(m1*altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0))
           y2 = max(1.0,k2*exp(m2*altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0))
           NeutralSourcesTotal = 0.0
           NeutralLossesTotal = 0.0
           tr       = tr3d(iLon,iLat,iAlt)
           tr3      = tr33d(iLon,iLat,iAlt)
           te3      = te33d(iLon,iLat,iAlt)
           te3m05   = te3m053d(iLon,iLat,iAlt)
           te3m07   = te3m073d(iLon,iLat,iAlt)
           te07     = te073d(iLon,iLat,iAlt)
           te3m085  = te3m0853d(iLon,iLat,iAlt)
           te3m091  = te3m0913d(iLon,iLat,iAlt)
           te3m081  = te3m0813d(iLon,iLat,iAlt)
           te3m039  = te3m0393d(iLon,iLat,iAlt)
           te12m056 = te12m0563d(iLon,iLat,iAlt)
           tr3m044  = tr3m0443d(iLon,iLat,iAlt)
           tr3m04   = tr3m043d(iLon,iLat,iAlt)
           ti3m044  = ti3m0443d (iLon,iLat,iAlt)
           ti3m024  = ti3m0243d (iLon,iLat,iAlt)
           ti3m023  = ti3m0233d (iLon,iLat,iAlt)
           ti3m116  = ti3m1163d (iLon,iLat,iAlt)
           ti10m067 = ti10m0673d (iLon,iLat,iAlt)
           ti3m045  = ti3m0453d (iLon,iLat,iAlt)
           ti10m212 = ti10m2123d (iLon,iLat,iAlt)
           ti3m087  = ti3m0873d (iLon,iLat,iAlt)
           ti3m052  = ti3m0523d (iLon,iLat,iAlt)
           ti9m092  = ti9m0923d (iLon,iLat,iAlt)
           ti3m055  = ti3m0553d (iLon,iLat,iAlt)
           ti15m023 = ti15m0233d (iLon,iLat,iAlt)
           te22m05  = eTemperature(iLon,iLat,iAlt,iBlock)**(.5) * &
                exp(te227d(iLon,iLat,iAlt))
           tr3m08  = tr3m083d(iLon,iLat,iAlt)
           ti = iTemperature(iLon,iLat,iAlt,iBlock)
           tn = Temperature(iLon,iLat,iAlt,iBlock)*&
                TempUnit(iLon,iLat,iAlt)

           tn1 = exp(107.8/tn)
           tn06 = exp(67.5/tn)

           rr_opn2 = min(5.0e-19,4.5e-20*tr3**2)

           k1_n2_point = k1_n2(iLon,iLat,iAlt)
           k2_o2_point = k2_o2(iLon,iLat,iAlt)
           k3_no_point = k3_no(iLon,iLat,iAlt)

           DtTotal = 0.0
           EmissionTotal = 0.0

           Ions = IDensityS(iLon,iLat,iAlt,:,iBlock)

           Neutrals = NDensityS(iLon,iLat,iAlt,:,iBlock)

           niters = 0
           o2ptotal = 0

           do while (DtTotal < Dt)
              
              ChemicalHeatingSub = 0.0
              ChemicalHeatingSubI = 0.0
              ChemicalHeatingSubE = 0.0
              ChemicalHeatingS = 0
              Emission = 0.0

              DtSub = Dt - DtTotal

              IonSources = 0.0
              NeutralSources = 0.0
              IonLosses  = 0.0
              NeutralLosses = 0.0

              ! ----------------------------------------------------------
              ! O2 -> 2O
              ! ----------------------------------------------------------
              Reaction = EuvDissRateS(iLon,iLat,iAlt,iO2_,iBlock)

              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 2*Reaction

              ! -----------------------------------------------------------
              ! O+O+M -> O2+M + 5.12 eV
              ! -----------------------------------------------------------
              rr=rr_o_p_o_p_m__o2_p_m_5p12ev *((300./tn)**(2)) ! 4.7e-33 R1
              rr= rr*1.e-12  !cm6s-1-->m6s-1

              Reaction = rr * Neutrals(iO_3P_)**2 *&
                   (Neutrals(iO2_)+ &
                    Neutrals(iO_3P_)+ &
                    Neutrals(iN2_))

              NeutralLosses(iO_3P_) = NeutralLosses(iO_3P_) + 2*Reaction
              NeutralSources(iO2_) = NeutralSources(iO2_) + Reaction
			  
              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 5.12 

              ! ----------------------------------------------------------
              ! N2 -> 2N
              ! ----------------------------------------------------------

              Reaction = EuvDissRateS(iLon,iLat,iAlt,iN2_,iBlock)

              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction
              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + .25*Reaction
              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + .60*Reaction

              ! Solar EUV

              ! ----------------------------------------------------------
              ! N2+
              ! ----------------------------------------------------------

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iN2P_,iBlock)

              IonSources(iN2P_)   = IonSources(iN2P_)   + Reaction
              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction

              ! Aurora

              Reaction = AuroralIonRateS(iLon,iLat,iAlt,iN2_, iBlock) + &
                   IonPrecipIonRateS(iLon,iLat,iAlt,iN2_, iBlock)

              IonSources(iN2P_)   = IonSources(iN2P_) + Reaction
              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction

              ! O+(2D) + N2 -> N2+ + O + 1.35 eV
			  
			  rr = rr_op2d_p_n2__n2p_p_1p35ev * ti3m055  ! 1.5e-16 R2
              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Neutrals(iN2_)

              IonSources(iN2P_)    = IonSources(iN2P_)   + Reaction
              NeutralSources(iO_3P_)  = NeutralSources(iO_3P_) + Reaction
              IonLosses(iO_2DP_)   = IonLosses(iO_2DP_)  + Reaction
              NeutralLosses(iN2_)  = NeutralLosses(iN2_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + Reaction * 0.859
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 0.491
              ChemicalHeatingS(iop2d_n2) =  &
                   ChemicalHeatingS(iop2d_n2) + &
                   Reaction * 1.33

              ! O+(2P) + N2 -> N2+ + O + 3.05 eV

              rr = rr_op2p_p_n2__n2p_p_o_p_3p05ev * ti3m055  ! 2.0e-16 R3
              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Neutrals(iN2_)

              IonSources(iN2P_)    = IonSources(iN2P_)   + Reaction
              NeutralSources(iO_3P_)  = NeutralSources(iO_3P_) + Reaction
              IonLosses(iO_2PP_)   = IonLosses(iO_2PP_)  + Reaction
              NeutralLosses(iN2_)  = NeutralLosses(iN2_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + Reaction * 1.941
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.109

               ChemicalHeatingS(iop2p_n2) =  &
                   ChemicalHeatingS(iop2p_n2) + &
                   Reaction * 3.02

              ! N2+ + O2 -> O2+ + N2 + 3.53 eV

               if (ti<=1000.0) then
                  rr = rr_n2p_p_o2__o2p_p_n2_p_3p53ev_lt * ti3m116  ! 5.1e-17 R4
               else 
                  rr = rr_n2p_p_o2__o2p_p_n2_p_3p53ev_gt * ti10m067 ! 1.26e-17 R5
               endif

               Reaction = &
                    rr * &
                    Ions(iN2P_) * &
                    Neutrals(iO2_)

              IonSources(iO2P_)    = IonSources(iO2P_)   + Reaction
              NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
              IonLosses(iN2P_)     = IonLosses(iN2P_)  + Reaction
              NeutralLosses(iO2_)  = NeutralLosses(iO2_) + Reaction

              o2ptotal = o2ptotal + reaction
              ChemicalHeatingSub = &
                   ChemicalHeatingSub + Reaction * 1.822
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.708
              
              ChemicalHeatingS(in2p_o2) =  &
                   ChemicalHeatingS(in2p_o2) + &
                   Reaction * 3.53

              ! N2+ + O -> NO+ + N(2D) + 0.70 eV
              !!!!!             -> NO+ + N(4S) + 3.08 eV

              if (ti<=1500.0) then
                 rr = rr_n2p_p_o__nop_p_n2d_p_0p70ev__nop_p_n4s_p_3p08ev_lt * ti3m044  ! 1.33e-16 R6
              else 
                 rr = rr_n2p_p_o__nop_p_n2d_p_0p70ev__nop_p_n4s_p_3p08ev_gt * ti15m023 ! 6.55e-17 R7
              endif

              Reaction = &
                   rr * &
                   Ions(iN2P_) * &
                   Neutrals(iO_3P_)

              IonSources(iNOP_)      = IonSources(iNOP_)      + Reaction
              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction
              IonLosses(iN2P_)       = IonLosses(iN2P_)       + Reaction
              NeutralLosses(iO_3P_)     = NeutralLosses(iO_3P_)     + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.477 

              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 0.223
                   
              ChemicalHeatingS(in2p_o) =  &
                   ChemicalHeatingS(in2p_o) + &
                   Reaction * 0.70
                   
              ! N2+ + e -> 2 N(2D) + 1.04 eV
              ! N2+ + e -> 2 N(4S) + 5.77 eV

              rr = rr_n2p_p_e__2n2d_p_1p04ev__2n4s_p_5p77ev * te3m039 ! 2.2e-13 R8

              Reaction = &
                   rr * &
                   Ions(iN2P_) * &
                   Ions(ie_)

              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + 2*Reaction*0.56
              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + 2*Reaction*0.44
              IonLosses(iN2P_)       = IonLosses(iN2P_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.04 *0.56 + &
                   Reaction * 5.77 * 0.44

              ChemicalHeatingS(in2p_e) =  &
                   ChemicalHeatingS(in2p_e) + &
                   Reaction * 3.44

              ! N2+ + N(4S) -> N2 + N+ + 2.48 eV

              rr = rr_n2p_p_n4s__n2_p_np_p_2p48ev ! 1.0e-17 R9

              Reaction = &
                   rr * &
                   Ions(iN2P_) * &
                   Neutrals(iN_4S_)

              NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
              IonSources(iNP_) = IonSources(iNP_) + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              IonLosses(iN2P_)       = IonLosses(iN2P_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.827
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.653
                  

              ! N2+ + O -> O+(4S) + N2 + 1.96 eV

              rr = rr_n2p_p_o__op4s_p_n2_p_1p96ev * ti3m023  ! 7.0e-18 R10

              Reaction = &
                   rr * &
                   Ions(iN2P_) * &
                   Neutrals(iO_3P_)

              NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
              IonSources(iO_4SP_)  = IonSources(iO_4SP_)  + Reaction
              NeutralLosses(iO_3P_)   = NeutralLosses(iO_3P_) + Reaction
              IonLosses(iN2P_)     = IonLosses(iN2P_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.712
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.248

              ChemicalHeatingS(in2p_o) =  &
                   ChemicalHeatingS(in2p_o) + &
                   Reaction * 1.96

              ! N2+ + NO -> NO+ + N2 + 6.33 eV

              rr = rr_n2p_p_no__nop_p_p_n2_p_6p33ev  ! Richards 3.6e-16  ! R11

              Reaction = &
                   rr * &
                   Ions(iN2P_) * &
                   Neutrals(iNO_)

              NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
              IonSources(iNOP_)    = IonSources(iNOP_)     + Reaction
              NeutralLosses(iNO_)  = NeutralLosses(iNO_)  + Reaction
              IonLosses(iN2P_)     = IonLosses(iN2P_)    + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.274
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 3.056

              ! ----------------------------------------------------------
              ! O2+
              ! ----------------------------------------------------------

              ! -----------
              ! Solar EUV
              ! -----------

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO2P_,iBlock)
           
              IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction

              o2ptotal = o2ptotal + reaction

              ! -----------
              ! Aurora
              ! -----------

              Reaction = AuroralIonRateS(iLon,iLat,iAlt,iO2_, iBlock) + &
                   IonPrecipIonRateS(iLon,iLat,iAlt,iO2_, iBlock)

              IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction

              o2ptotal = o2ptotal + reaction

              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction

              ! -----------
              ! O+(4S) + O2 -> O2+ + O + 1.55 eV
              ! -----------

              if (ti<=900.0) then
                 rr = rr_op4s_p_o2__o2p_p_o_p_1p55ev_lt * ti3m052 ! 1.6e-17 R12
              else
                 rr = rr_op4s_p_o2__o2p_p_o_p_1p55ev_gt * ti9m092  ! 9e-18 R13
              endif

              Reaction = &
                   rr * &
                   Ions(iO_4SP_) * &
                   Neutrals(iO2_)

              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction
              o2ptotal = o2ptotal + reaction

              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
              IonLosses(iO_4SP_)  = IonLosses(iO_4SP_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.033
				  
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 0.517

              ChemicalHeatingS(iop_o2) =  &
                   ChemicalHeatingS(iop_o2) + &
                   Reaction * 1.55

              ! -----------
              ! O+(2D) + O2 -> O2+ + O + 4.865 eV
              ! -----------

              rr = rr_op2d_p_o2__o2p_p_o_p_4p865ev  ! 7.0e-16 R14

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Neutrals(iO2_)

              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction
              o2ptotal = o2ptotal + reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction/1.0
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.243
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.622
              
              ! -----------
              ! O+(2P) + O2 -> O2+ + O + 6.54 eV
              ! -----------

              rr = rr_op2p_p_o2__o2p_p_o_p_6p54ev  ! 1.3e-16 R15
              
              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Neutrals(iO2_)

              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction
              o2ptotal = o2ptotal + reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 4.36
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 2.18

              ! -----------
              ! O+(2P) + O2 -> O+(4S) + O2 + 5.016 eV
              ! -----------

              rr = rr_op2p_p_o2__op4s_p_o2_p_5p016ev ! 1.3e-16 R16

              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Neutrals(iO2_)

              NeutralSources(iO2_) = NeutralSources(iO2_) + Reaction
              IonSources(iO_4SP_)   = IonSources(iO_4SP_)   + Reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.672
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 3.344

              ! -----------
              ! N+ + O2 -> O2+ + N(4S) + 2.5 eV
              ! -----------

              if (ti<=1000.0) then
                 rr = rr_np_p_o2__o2p_p_n4s_p_2p5ev_lt * ti3m045 ! 1.925e-16 R17
              else 
                 rr = rr_np_p_o2__o2p_p_n4s_p_2p5ev_gt ! 3.325e-16 R18
              endif

             Reaction = &
                  rr * &
                  Ions(iNP_) * &
                  Neutrals(iO2_)

             NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
             IonSources(iO2P_)      = IonSources(iO2P_)      + Reaction
              o2ptotal = o2ptotal + reaction
             NeutralLosses(iO2_)    = NeutralLosses(iO2_)    + Reaction
             IonLosses(iNP_)        = IonLosses(iNP_)       + Reaction

             ChemicalHeatingSub = &
                  ChemicalHeatingSub + &
                  Reaction * 1.739
				  
             ChemicalHeatingSubI = &
                  ChemicalHeatingSubI + Reaction * 0.761
             
             ChemicalHeatingS(inp_o2) =  &
                  ChemicalHeatingS(inp_o2) + &
                  Reaction * 2.486
             
             ! -----------
             ! N+ + O2 -> O2+ + N(2D) + 0.1 eV
             ! -----------

             if (ti<=1000.0) then
                rr = rr_np_p_o2__o2p_p_n2d_p_0p1ev_lt  * ti3m045 ! 0.825e-16 R19
             else
                rr = rr_np_p_o2__o2p_p_n2d_p_0p1ev_gt  ! 1.425e-16 R20
             endif

             Reaction = &
                  rr * &
                  Ions(iNP_) * &
                  Neutrals(iO2_)

             NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction

             IonSources(iO2P_)      = IonSources(iO2P_)      + Reaction
             o2ptotal               = o2ptotal + reaction
             NeutralLosses(iO2_)    = NeutralLosses(iO2_)    + Reaction
             IonLosses(iNP_)        = IonLosses(iNP_)       + Reaction

             ChemicalHeatingSub = &
                  ChemicalHeatingSub + &
                  Reaction * 0.07
				  
             ChemicalHeatingSubI = &
                  ChemicalHeatingSubI + Reaction * 0.03

             ChemicalHeatingS(inp_o2) =  &
                  ChemicalHeatingS(inp_o2) + &
                  Reaction * 0.1

             ! -----------
             ! O2+ + e -> O(1D) + O(1D) + 3.06 eV
             ! O2+ + e -> O(3P) + O(1D) + 5.02 eV
             ! O2+ + e -> O(3P) + O(3P) + 6.99 eV
             ! -----------

             if (ti<=1200.0) then
                rr = rr_o2p_p_e__o1d_p_o1d_p_3p06ev_lt * te3m07 ! 1.95e-13 R21
             else
                rr = rr_o2p_p_e__o1d_p_o1d_p_3p06ev_gt * te12m056 ! 7.39e-14 R22
             endif

             Reaction = &
                  rr * &
                  Ions(iO2P_) * &
                  Ions(ie_)
           
             if (UseNeutralConstituent(iO_1D_)) then
                NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 0.22*Reaction * 2.0
                NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 0.42*Reaction
                NeutralSources(iO_1D_) = NeutralSources(iO_1D_) + 0.42*Reaction
                NeutralSources(iO_1D_) = NeutralSources(iO_1D_) + 0.31*Reaction * 2.0
                ! This really should be 0.05 to O(1D) and O(1S)
                NeutralSources(iO_1D_) = NeutralSources(iO_1D_) + 0.05*Reaction
                NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 0.05*Reaction

                ChemicalHeatingSub = &
                     ChemicalHeatingSub + &
                     Reaction * 6.99 * 0.22 + &
                     Reaction * 5.02 * 0.42 + &
                     Reaction * 3.06 * 0.31
                 
                ChemicalHeatingS(io2p_e) = &
                     ChemicalHeatingS(io2p_e) + &
                     Reaction * 6.99 * 0.22 + &
                     Reaction * 5.02 * 0.42 + &
                     Reaction * 3.06 * 0.31

             else
                 
                NeutralSources(iO_3P_)    = &
                     NeutralSources(iO_3P_) + Reaction * 2.0

                ChemicalHeatingSub = &
                     ChemicalHeatingSub + &
                     Reaction * 5.0
                 
                ChemicalHeatingS(io2p_e) = &
                     ChemicalHeatingS(io2p_e) + &
                     Reaction * 5.0

             endif
              
             IonLosses(iO2P_)      = IonLosses(iO2P_)   + Reaction

             ! -----------
             ! O2+ + N(4S) -> NO+ + O + 4.21 eV
             ! -----------

             rr = rr_o2p_p_n4s__nop_p_o_4p21ev ! 1.0e-16 Richards R23

             Reaction = &
                  rr * &
                  Ions(iO2P_) * &
                  Neutrals(iN_4S_)

             NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
             IonSources(iNOP_)     = IonSources(iNOP_)     + Reaction
             NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
             IonLosses(iO2P_)      = IonLosses(iO2P_)      + Reaction
             ChemicalHeatingSub = &
                  ChemicalHeatingSub + &
                  Reaction * 2.746
				   
             ChemicalHeatingSubI = &
                  ChemicalHeatingSubI + Reaction * 1.464

             ChemicalHeatingS(io2p_n) =  &
                  ChemicalHeatingS(io2p_n) + &
                  Reaction * 4.25

             ! -----------
             ! O2+ + N(2D) -> NO+ + O + 6.519 eV
             ! -----------

             rr = rr_o2p_p_n4s__nop_p_o_4p21ev !1.8e-16 ! Richards R24

             Reaction = &
                  rr * &
                  Ions(iO2P_) * &
                  Neutrals(iN_2D_)

             NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
             IonSources(iNOP_)     = IonSources(iNOP_)     + Reaction
             NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
             IonLosses(iO2P_)      = IonLosses(iO2P_)      + Reaction
             
             ChemicalHeatingSub = &
                  ChemicalHeatingSub + &
                  Reaction * 4.252
				   
             ChemicalHeatingSubI = &
                  ChemicalHeatingSubI + Reaction * 2.267

             ! -----------
             ! O2+ + N(2P) -> O2+ + N(4S) + 3.565 eV
             ! -----------

             rr = rr_o2p_p_n2p__o2p_p_n4s_p_3p565ev !2.2e-17 ! Richards R25

             Reaction = &
                  rr * &
                  Ions(iO2P_) * &
                  Neutrals(iN_2P_)

             NeutralSources(iN_4S_)   = NeutralSources(iN_4S_)   + Reaction
             IonSources(iO2P_)     = IonSources(iO2P_)     + Reaction
             NeutralLosses(iN_2P_) = NeutralLosses(iN_2P_) + Reaction
             IonLosses(iO2P_)      = IonLosses(iO2P_)      + Reaction
              
             ChemicalHeatingSub = &
                  ChemicalHeatingSub + &
                  Reaction * 2.479

             ChemicalHeatingSubI = &
                  ChemicalHeatingSubI + Reaction * 1.086

             ! -----------
             ! O2+ + NO -> NO+ + O2 + 2.813 eV
             ! -----------

              rr = rr_o2p_p_no_nop_p_o2_p_2p813ev ! 4.5e-16 ! schunk and nagy !R26

              Reaction = &
                   rr * &
                   Ions(iO2P_) * &
                   Neutrals(iNO_)

              NeutralSources(iO2_) = NeutralSources(iO2_) + Reaction
              IonSources(iNOP_)    = IonSources(iNOP_)    + Reaction
              NeutralLosses(iNO_)  = NeutralLosses(iNO_)  + Reaction
              IonLosses(iO2P_)     = IonLosses(iO2P_)     + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.361
				   
				   ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.452

              ChemicalHeatingS(io2p_no) =  &
                   ChemicalHeatingS(io2p_no) + &
                    Reaction * 2.813

!!! Temp change to stop crash
!
!              ! -----------
!              ! O2+ + N2 -> NO+ + NO + 0.9333 eV
!              ! -----------
!               
!              rr = 5.0e-22
!               
!              Reaction = &
!                   rr * &
!                   Ions(iO2P_) * &
!                   Neutrals(iN2_)
!
!              NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction
!
!              IonSources(iNOP_)    = IonSources(iNOP_)    + Reaction
!              NeutralLosses(iN2_)  = NeutralLosses(iN2_)  + Reaction
!              IonLosses(iO2P_)     = IonLosses(iO2P_)     + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 0.9333
!
!              ChemicalHeatingS(io2p_n2) =  &
!                   ChemicalHeatingS(io2p_n2) + &
!                   Reaction * 0.9333

              ! ----------------------------------------------------------
              ! O(4S)+
              ! ----------------------------------------------------------

              ! Solar EUV

!              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO_4SP_,iBlock) * &
!                   Neutrals(iO_3P_)

              rr=EuvIonRateS(iLon,iLat,iAlt,iO_4SP_,iBlock)
              Reaction = rr

              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! Aurora

              ! Aurora goes 0.4, 0.4, 0.2 into O(4S), O(2D) and O(2P) respectively
              Reaction = &
                   0.4 * AuroralIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock) + &
                   0.4 * IonPrecipIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock)

              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! Aurora goes 0.4, 0.4, 0.2 into O(3P), O(2D) and O(2P) respectively
              Reaction = &
                   0.4 * AuroralIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock) + &
                   0.4 * IonPrecipIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock)

              IonSources(iO_2DP_) = IonSources(iO_2DP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! Aurora goes 0.4, 0.4, 0.2 into O(3P), O(2D) and O(2P) respectively
              Reaction = &
                   0.2 * AuroralIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock) + &
                   0.2 * IonPrecipIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock)

              IonSources(iO_2PP_) = IonSources(iO_2PP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! -----------
              ! O+(2D) + O -> O+(4S) + O(3P) + 3.31 eV
              ! O+(2D) + O -> O+(4S) + O(1D) + 1.35 eV
              ! -----------

              rr = rr_op2d_p_o_op4s_p_o3p_p_3p31ev__op4s_p_o1d_p_1p35ev ! 1.0e-17 ! R27

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Neutrals(iO_3P_)

              ! We create and loose the same amount of O 
              ! (when O(1D) is not used...
             
              if (UseNeutralConstituent(iO_1D_)) then
                 
                 NeutralSources(iO_3P_)  = &
                      NeutralSources(iO_3P_) + 0.5 * Reaction
                 NeutralSources(iO_1D_)  = &
                      NeutralSources(iO_1D_) + 0.5 * Reaction
                 
                 NeutralLosses(iO_3P_)   = NeutralSources(iO_3P_)  + Reaction

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 1.655 * 0.5 + &
                      Reaction * 0.675 * 0.5
					  
                 ChemicalHeatingSubI = &
                      ChemicalHeatingSubI + &
                      Reaction * 1.655 * 0.5 + &
                      Reaction * 0.675 * 0.5 
              else

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 1.655
					  
                 ChemicalHeatingSubI = &
                      ChemicalHeatingSubI + Reaction * 1.655

              endif
              
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction
              
              ! -----------
              ! O+(2D) + e -> O+(4S) + e + 3.31 eV
              ! -----------

              rr = rr_op2d_p_e__op4s_p_e_p_3p31ev * te3m05 ! 6.03e-14 R28

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Ions(ie_)

              ! We create and loose the same amount of e
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.31
				   
              ChemicalHeatingS(iop2d_e) =  &
                   ChemicalHeatingS(iop2d_e) + &
                   Reaction * 0.0

!!! Temp change to stop crash
!               ! -----------
!               ! O+(2D) + N2 -> O+(4S) + N2 + 3.31 eV
!               ! -----------
!
!               rr = 8.0e-16
!
!               Reaction = &
!                    rr * &
!                    Ions(iO_2DP_) * &
!                    Neutrals(iN2_)
!
!               ! We create and loose the same amount of N2
!               IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
!               IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction
!
!               ChemicalHeatingSubI = &
!                    ChemicalHeatingSubI + &
!                    Reaction * 3.31
!
!               ChemicalHeatingS(iop2d_n2) =  &
!                    ChemicalHeatingS(iop2d_n2) + &
!                    Reaction * 3.31

              ! -----------
              ! O+(2P) + O -> O+(4S) + O + 5.0 eV
              ! -----------

              rr = rr_op2p_p_o__op4s_p_o_p_5p0ev ! 4.0e-16  ! R29

              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Neutrals(iO_3P_)

              ! We create and loose the same amount of O
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.5
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 2.5

              ChemicalHeatingS(iop2p_o) =  &
                   ChemicalHeatingS(iop2p_o) + &
                   Reaction * 5.0

              ! -----------
              ! O+(2P) + e -> O+(4S) + e + 5.0 eV
              ! -----------

              rr = rr_op2p_p_e__op4s_p_e_p_5p0ev * te3m05 ! 3.03e-14 R30

              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Ions(ie_)

              ! We create and loose the same amount of e
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 5.0

              ! maybe should be 0 because the energy should go the electrons
              ChemicalHeatingS(iop2p_e) =  &
                   ChemicalHeatingS(iop2p_e) + &
                   Reaction*5 

              ! -----------
              ! O+(2P) -> O+(4S) + 2470A
              ! -----------

              rr = rr_op2p__op4s_p_2470a ! 0.047  ! R31

              Reaction = &
                   rr * &
                   Ions(iO_2PP_)

              ! We create and loose the same amount of e
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              Emission(iE2470_) = Emission(iE2470_) + Reaction

              ! -----------
              ! N+ + O2 -> O+(4S) + NO + 2.31 eV
              ! -----------
              
              if (ti<=1000.0) then
                 rr = rr_np_p_o2__op4s_p_no_p_2p31ev_lt * ti3m045 ! 0.275e-16 ! R32
              else 
                 rr = rr_np_p_o2__op4s_p_no_p_2p31ev_gt ! 0.475e-16  ! R33
              endif
			 
              Reaction = &
                   rr * &
                   Ions(iNP_) * &
                   Neutrals(iO2_)

              NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction
              IonSources(iO_4SP_)  = IonSources(iO_4SP_)  + Reaction
              NeutralLosses(iO2_)  = NeutralLosses(iO2_)  + Reaction
              IonLosses(iNP_)      = IonLosses(iNP_)      + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.803
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.507

              ChemicalHeatingS(inp_o2) =  &
                   ChemicalHeatingS(inp_o2) + &
                   Reaction * 2.31

              ! -----------
              ! O+(4S) + N2 -> NO+ + N(4S) + 1.10 eV
              ! -----------

              if (ti<=1000.0) then
                 rr = rr_op4s_p_n2__nop_p_n4S_p_1p10ev_lt * ti3m045 ! 1.2e-18 ! R34
              else
                 rr = rr_op4s_p_n2__nop_p_n4S_p_1p10ev_gt * ti10m212 ! 7.0e-19 ! R35
              endif

              Reaction = &
                   rr * &
                   Ions(iO_4SP_) * &
                   Neutrals(iN2_)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              IonSources(iNOP_)      = IonSources(iNOP_)      + Reaction
              NeutralLosses(iN2_)    = NeutralLosses(iN2_)    + Reaction
              IonLosses(iO_4SP_)     = IonLosses(iO_4SP_)     + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.75

              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 0.35

              ChemicalHeatingS(iop_n2) =  &
                   ChemicalHeatingS(iop_n2) + &
                   Reaction * 1.10

              ! -----------
              ! O+(4S) + NO -> NO+ + O + 4.36 eV
              ! -----------

              rr = rr_op4s_p_no__nop_p_o_p_4p36ev * ti3m087  ! 7.0e-19 ! R36

              Reaction = &
                   rr * &
                   Ions(iO_4SP_) * &
                   Neutrals(iNO_)

              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              IonSources(iNOP_)   = IonSources(iNOP_)   + Reaction
              NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction
              IonLosses(iO_4SP_)  = IonLosses(iO_4SP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.844
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.516

              ! -----------
              ! O+(4S) + N(2D) -> N+ + O + 1.45 eV
              ! -----------

              rr = rr_op4s_p_n2d__np_p_o_p_1p45ev ! 1.3e-16  ! R37

              Reaction = &
                   rr * &
                   Ions(iO_4SP_) * &
                   Neutrals(iN_2D_)

              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              IonSources(iNP_)      = IonSources(iNP_)      + Reaction
              NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              IonLosses(iO_4SP_)    = IonLosses(iO_4SP_)    + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.677
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 0.773
              
              ! ----------------------------------------------------------
              ! O(2D)+
              ! ----------------------------------------------------------

              ! Solar EUV

!              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO_2DP_,iBlock) * &
!                   Neutrals(iO_3P_)

              rr=EuvIonRateS(iLon,iLat,iAlt,iO_2DP_,iBlock)
              Reaction = rr 

              IonSources(iO_2DP_) = IonSources(iO_2DP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! -----------
              ! O+(2P) + e -> O+(2D) + e + 1.69 eV
              ! -----------
              rr = rr_op2p_p_e__op2d_p_e_p_1p69ev * te3m05  ! 1.84e-13 ! R38

              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Ions(ie_)

              ! We create and loose the same amount of e
              IonSources(iO_2DP_) = IonSources(iO_2DP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              ChemicalHeatingSubI = &
                   ChemicalHeatingSub + &
                   Reaction * 1.69
				   
              ChemicalHeatingS(iop2p_e) =  &
                   ChemicalHeatingS(iop2p_e) + &
                   Reaction * 1.69
                   
               ! -----------
              ! O+(2D) + N2 -> NO+ + N + 4.41 eV
              ! -----------

              rr = rr_op2d_p_n2__nop_p_n_p_4p41ev ! 2.5e-17  ! R39

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Neutrals(iN2_)

              NeutralLosses(iN2_)  = NeutralLosses(iN2_)  + Reaction
              IonSources(iNOP_) = IonSources(iNOP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction
              NeutralSources(iN_4S_)  = NeutralSources(iN_4S_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.007
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.403
                   
              ! -----------
              ! O+(2D) + NO -> NO+ + O + 4.37 eV 
              ! -----------

              rr = rr_op2d_p_no__nop_p_o_p_4p37ev ! 1.2e-15   ! R40

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Neutrals(iNO_)

              NeutralLosses(iNO_)  = NeutralLosses(iNO_)  + Reaction
              IonSources(iNOP_) = IonSources(iNOP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction
              NeutralSources(iO_3P_)  = NeutralSources(iO_3P_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.85
				  
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.52
              
              ! -----------
              ! O+(2P) -> O+(2D) + 7320A
              ! -----------

              rr = rr_op2p__op2d_p_7320a ! 0.171  ! R41

              Reaction = &
                   rr * &
                   Ions(iO_2PP_)

              IonSources(iO_2DP_) = IonSources(iO_2DP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              Emission(iE7320_) = Emission(iE7320_) + Reaction

              ! -----------
              ! O+(2D) -> O+(4S) + 3726A
              ! -----------

              rr = rr_op2d__op4s_p_3726a ! 7.7e-5  ! R42
              
              Reaction = &
                   rr * &
                   Ions(iO_2DP_)

              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction

              Emission(iE3726_) = Emission(iE3726_) + Reaction

              ! ----------------------------------------------------------
              ! O(2P)+
              ! ----------------------------------------------------------

              ! Solar EUV

              !Reaction = EuvIonRateS(iLon,iLat,iAlt,iO_2PP_,iBlock) * &
              !     Neutrals(iO_3P_)
!              rr=EuvIonRateS(iLon,iLat,iAlt,iO_2PP_,iBlock) 

              rr = EuvIonRateS(iLon,iLat,iAlt,iO_2PP_,iBlock)
              Reaction = rr

              IonSources(iO_2PP_) = IonSources(iO_2PP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! ----------------------------
              ! Atomic N Photoionization
              ! ----------------------------
              ! ----------------------------------------------------------
              ! N(4S) + hv -> N+
              ! ----------------------------------------------------------
              Reaction = EuvIonRateS(iLon,iLat,iAlt,iNP_,iBlock)

              IonSources(iNP_) = IonSources(iNP_) + Reaction
              NeutralLosses(iN_4S_)  = NeutralLosses(iN_4S_)  + Reaction

              ! ----------------------------
              ! Atomic He Photoionization
              ! ----------------------------
              ! ----------------------------------------------------------
              ! He + hv --> He+  + e-
              ! ----------------------------------------------------------
              Reaction = EuvIonRateS(iLon,iLat,iAlt,iHeP_,iBlock)

              NeutralLosses(iHe_) = NeutralLosses(iHe_) + Reaction
              IonSources(iHeP_) = IonSources(iHeP_) + Reaction

! ----------------------------
! NO Photoionization

!              IonSources(iO_2PP_) = IonSources(iO_2PP_) + Reaction
!              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

!!! Temp change to stop crash
!               ! -----------
!               ! O+(2P) + N2 -> N+ + NO + 0.70 eV
!               ! -----------
!
!               rr = 1.0e-16
!
!               Reaction = &
!                    rr * &
!                    Ions(iO_2PP_) * &
!                    Neutrals(iN2_)
!
!               NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction
!               IonSources(iNP_)     = IonSources(iNP_)     + Reaction
!               NeutralLosses(iN2_)  = NeutralLosses(iN2_)  + Reaction
!               IonLosses(iO_2PP_)   = IonLosses(iO_2PP_)   + Reaction
!
!               ChemicalHeatingSub = &
!                    ChemicalHeatingSub + &
!                    Reaction * 0.70
!
!               ChemicalHeatingS(iop2p_n2) =  &
!                    ChemicalHeatingS(iop2p_n2) + &
!                    Reaction * 0.70

               ! ----------------------------------------------------------
               ! N+
               ! ----------------------------------------------------------
               
              !! Temp change to stop crash
              ! ! -----------
              ! ! O2+ + N(2D) -> N+ + O2 + 0.0 eV
              ! ! -----------
              !
              ! rr = 8.65e-17
              !
              ! Reaction = &
              !      rr * &
              !      Ions(iO2P_) * &
              !      Neutrals(iN_2D_)
              !
              ! NeutralSources(iO2_)  = NeutralSources(iO2_)  + Reaction
              ! IonSources(iNP_)      = IonSources(iNP_)      + Reaction
              ! NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              ! IonLosses(iO2P_)      = IonLosses(iO2P_)      + Reaction
              ! ChemicalHeatingSub = &
              !      ChemicalHeatingSub + &
              !      Reaction * 0.0

              ! -----------
              ! Shunk and Nagy R29
              ! He+ + N2 -> N+ + N + He + 0.28 eV
              ! -----------

              rr = rr_hep_p_n2__np_p_n_p_he_p_0p28ev/1.0e6  ! 7.8e-10 ! R43 

              Reaction = &
                   rr * &
                   Ions(iHeP_) * &
                   Neutrals(iN2_)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralSources(iHe_)   = NeutralSources(iHe_)   + Reaction
              IonSources(iNP_)       = IonSources(iNP_)       + Reaction
              NeutralLosses(iN2_)    = NeutralLosses(iN2_)    + Reaction
              IonLosses(iHeP_)       = IonLosses(iHeP_)       + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.28

              ! -----------
              ! Shunk and Nagy R30
              ! He+ + N2 -> N2+ + He ( + ??? eV)
              ! -----------

              rr = rr_hep_p_n2__n2p_p_he/1.0e6  ! 5.2e-10 ! R44

              Reaction = &
                   rr * &
                   Ions(iHeP_) * &
                   Neutrals(iN2_)

              NeutralSources(iHe_)   = NeutralSources(iHe_)   + Reaction
              IonSources(iN2P_)      = IonSources(iN2P_)      + Reaction
              NeutralLosses(iN2_)    = NeutralLosses(iN2_)    + Reaction
              IonLosses(iHeP_)       = IonLosses(iHeP_)       + Reaction

              !ChemicalHeatingSub = &
              !     ChemicalHeatingSub + &
              !     Reaction * 0.28

              ! -----------
              ! Shunk and Nagy R31
              ! He+ + O2 -> O+ O + He ( + ??? eV)
              ! -----------

              rr = rr_hep_p_o2__o_p_o_p_he/1.0e6  ! 9.7e-10 ! R45

              Reaction = &
                   rr * &
                   Ions(iHeP_) * &
                   Neutrals(iO2_)

              NeutralSources(iHe_) = NeutralSources(iHe_) + Reaction
              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
              IonLosses(iHeP_) = IonLosses(iHeP_) + Reaction

              ! -----------
              ! Shunk and Nagy Radiative Recombination
              ! He+ + e- -> He ( + ??? eV)
              ! -----------

              rr = rr_hep_p_em__he/1.0e6 * te07  ! 4.8e-12 ! R46

              Reaction = &
                   rr * &
                   Ions(iHeP_) * &
                   Ions(ie_)

              NeutralSources(iHe_) = NeutralSources(iHe_) + Reaction
              IonLosses(iHeP_) = IonLosses(iHeP_) + Reaction

              !ChemicalHeatingSub = &
              !     ChemicalHeatingSub + &
              !     Reaction * 0.28

              ! -----------
               ! O+(2P) + N -> N+ + O + 2.7 eV
               ! -----------

               rr = rr_op2p_p_n__np_p_o_p_2p7ev ! 1.0e-16  ! R47

               Reaction = &
                    rr * &
                    Ions(iO_2PP_) * &
                    Neutrals(iN_4S_)

               NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
               IonSources(iNP_)      = IonSources(iNP_)      + Reaction
               NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
               IonLosses(iO_2PP_)    = IonLosses(iO_2PP_)    + Reaction

               ChemicalHeatingSub = &
                    ChemicalHeatingSub + &
                    Reaction * 2.7

              ! -----------
              ! N+ + NO --> N2+ + O + 2.2 eV
              ! -----------

               rr = rr_np_p_no__n2p_p_o_p_2p2ev * ti3m024  ! 8.33e-17 !  R48

               Reaction = &
                    rr * &
                    Ions(iNP_) * &
                    Neutrals(iNO_)

               NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
               IonSources(iN2P_)      = IonSources(iN2P_)      + Reaction
               NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction
               IonLosses(iNP_)    = IonLosses(iNP_)    + Reaction

               ChemicalHeatingSub = &
                    ChemicalHeatingSub + &
                    Reaction * 1.4
				   
               ChemicalHeatingSubI = &
                    ChemicalHeatingSubI + Reaction * 0.8

               ! -----------
               ! N+ + NO --> NO+ + N(4S) + 3.4 eV 
               ! -----------

               rr = rr_np_p_no__nop_p_n4s_p_3p4ev * ti3m024  ! 4.72e-16 ! R49

               Reaction = &
                    rr * &
                    Ions(iNP_) * &
                    Neutrals(iNO_)

               NeutralSources(iN_4S_)   = NeutralSources(iN_4S_)   + Reaction
               IonSources(iNOP_)      = IonSources(iNOP_)      + Reaction
               NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction
               IonLosses(iNP_)    = IonLosses(iNP_)    + Reaction

               ChemicalHeatingSub = &
                    ChemicalHeatingSub + &
                    Reaction * 2.318
				   
               ChemicalHeatingSubI = &
                    ChemicalHeatingSubI + Reaction * 1.082

               ! -----------
               ! N+ + O2 --> NO+ + O(3P) + 6.67 eV  
               ! -----------
              
               if (ti<=1000) then
                  rr = rr_np_p_o2__nop_p_o3p_p_6p67ev_lt * ti3m045  ! 0.495e-16 ! R50
               else 
                  rr = rr_np_p_o2__nop_p_o3p_p_6p67ev_gt            ! 0.855e-16 ! R51
               endif

               Reaction = &
                    rr * &
                    Ions(iNP_) * &
                    Neutrals(iO2_)

               NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
               IonSources(iNOP_)      = IonSources(iNOP_)      + Reaction
               NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
               IonLosses(iNP_)    = IonLosses(iNP_)    + Reaction

               ChemicalHeatingSub = &
                    ChemicalHeatingSub + &
                    Reaction * 4.35
				   
               ChemicalHeatingSubI = &
                    ChemicalHeatingSubI + Reaction * 2.32

              ! -----------
              ! O+(2D) + N -> N+ + O + 1.0 eV
              ! -----------

              rr = rr_op2d_p_n__np_p_o_p_1p0ev ! 1.5e-16  ! R52

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Neutrals(iN_4S_)

              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              IonSources(iNP_)      = IonSources(iNP_)      + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              IonLosses(iO_2DP_)    = IonLosses(iO_2DP_)    + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.467
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 0.533

              ! -----------
              ! N+ + O2 -> NO+ + O(1D) + 4.71 eV
              ! -----------

              if (ti<=1000.0) then
                 rr = rr_np_p_o2__nop_p_o1d_p_4p71ev_lt * ti3m045  ! 1.98e-16 ! R53
              else 
                 rr = rr_np_p_o2__nop_p_o1d_p_4p71ev_gt            ! 3.42e-16 ! R54
              endif

              Reaction = &
                   rr * &
                   Ions(iNP_) * &
                   Neutrals(iO2_)

              if (UseNeutralConstituent(iO_1D_)) then
                 NeutralSources(iO_1D_) = NeutralSources(iO_1D_) + Reaction
              else
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              endif

              IonSources(iNOP_)   = IonSources(iNOP_)   + Reaction
              IonLosses(iNP_)     = IonLosses(iNP_)     + Reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.072
              
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 1.638

              ChemicalHeatingS(inp_o2) =  &
                   ChemicalHeatingS(inp_o2) + &
                   Reaction * 6.67

              ! -----------
              ! N+ + O -> O+ + N + 0.93 eV
              ! -----------

              rr = rr_np_p_o__op_p_n_p_0p93ev ! 2.2e-18  ! R55

              Reaction = &
                   rr * &
                   Ions(iNP_) * &
                   Neutrals(iO_3P_)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              IonSources(iO_4SP_)    = IonSources(iO_4SP_)    + Reaction
              NeutralLosses(iO_3P_)     = NeutralLosses(iO_3P_)     + Reaction
              IonLosses(iNP_)        = IonLosses(iNP_)        + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.496
				   
              ChemicalHeatingSubI = &
                   ChemicalHeatingSubI + Reaction * 0.434

              ChemicalHeatingS(inp_o) =  &
                   ChemicalHeatingS(inp_o) + &
                   Reaction * 0.93

!              ! -----------
!              ! N+ + H -> H+ + N + 0.90 eV
!              ! -----------
!
!              rr = 3.6e-12/1.0e6 
!
!              Reaction = &
!                   rr * &
!                   Ions(iNP_) * &
!                   Neutrals(iH_)
!
!              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
!              IonSources(iHP_)       = IonSources(iHP_)       + Reaction
!              NeutralLosses(iH_)     = NeutralLosses(iH_)     + Reaction
!              IonLosses(iNP_)        = IonLosses(iNP_)        + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 0.9

              ! ----------------------------------------------------------
              ! NO+
              ! ----------------------------------------------------------

              ! -----------
              ! NO+ + e -> O + N(2D) + 0.38 eV 
              ! -----------

              rr = rr_nop_p_e__o_p_n2d_p_0p38ev * te3m085  ! 3.4e-13 ! R56

              Reaction = &
                   rr * &
                   Ions(iNOP_) * &
                   Ions(ie_)

              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction
              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              IonLosses(iNOP_)       = IonLosses(iNOP_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.38

              ChemicalHeatingS(inop_e) =  &
                   ChemicalHeatingS(inop_e) + &
                   Reaction * 0.38
                   
              ! -----------
              ! NO+ + e -> O + N(4S) + 2.77 eV
              ! -----------

              rr = rr_nop_p_e__o_p_n4s_p_2p77ev * te3m085  ! 0.6e-13 ! R57

              Reaction = &
                   rr * &
                   Ions(iNOP_) * &
                   Ions(ie_)


              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              IonLosses(iNOP_)       = IonLosses(iNOP_)       + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.77

               ChemicalHeatingS(inop_e) =  &
                   ChemicalHeatingS(inop_e) + &
                   Reaction * 2.77

              ! ----------------------------------------------------------
              ! N(4S)
              ! ----------------------------------------------------------

              ! -----------
              ! N(2D) + e -> N(4S) + e + 2.38 eV
              ! -----------

              rr = rr_n2d_p_e__n4s_p_e_p_2p38ev * te3m081 ! 3.86e-16 ! R58

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_) * &
                   Ions(ie_)

              ! We create and loose the same amount of e
              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralLosses(iN_2D_)  = NeutralLosses(iN_2D_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.38
				   
              ! -----------
              ! N(2D) + O -> N(4S) + O(3P) + 2.38 eV
              ! N(2D) + O -> N(4S) + O(1D) + 0.42 eV
              ! -----------

              rr = rr_n2d_p_o__n4s_p_o3p_p_2p38ev__n4s_p_o1d_p_0p42ev ! 6.9e-19  ! R59

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_) * &
                   Neutrals(iO_3P_)

              if (UseNeutralConstituent(iO_1D_)) then
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 0.9*Reaction
                 NeutralSources(iO_1D_) = NeutralSources(iO_1D_) + 0.1*Reaction
                 NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction
                 
                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      0.9 * Reaction * 2.38 + &
                      0.1 * Reaction * 0.42
                 
                 ChemicalHeatingS(in2d_o) =  &
                      ChemicalHeatingS(in2d_o) + &
                      0.9 * Reaction * 2.38 + &
                      0.1 * Reaction * 0.42

              else
                 
                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 2.38

              endif

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralLosses(iN_2D_)  = NeutralLosses(iN_2D_)  + Reaction
              
              ! -----------
              ! N(2D) -> N(4S) + 5200A
              ! -----------

              rr = rr_n2d__n4s_p_5200a ! 1.06e-5  ! R60

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralLosses(iN_2D_)  = NeutralLosses(iN_2D_)  + Reaction
              Emission(iE5200_) = Emission(iE5200_) + Reaction

              ! -----------
              ! NO -> N(4S) + O
              ! -----------

              rr = rr_no__n4s_p_o *exp(-1.e-8*(Neutrals(iO2_)*1.e-6)**0.38) ! 4.5e-6 ! R61

              Reaction = &
                   rr * &
                   Neutrals(iNO_)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              NeutralLosses(iNO_)    = NeutralLosses(iNO_)    + Reaction

              ! -----------
              ! N(4S) + O2 -> NO + O + 1.385 eV
              ! -----------
              
              rr = rr_n4s_p_o2__no_p_o_p_1p385ev * tn * exp(-3270/tn) ! 1.5e-20 ! R62

              Reaction = &
                   rr * &
                   Neutrals(iN_4S_) * &
                   Neutrals(iO2_)

              NeutralSources(iNO_)  = NeutralSources(iNO_)  + Reaction
              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              NeutralLosses(iO2_)   = NeutralLosses(iO2_)   + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.385

              ChemicalHeatingS(in_o2) =  &
                   ChemicalHeatingS(in_o2) + &
                   Reaction * 1.385

              ! -----------
              ! N(4S) + NO -> N2 + O + 3.25 eV
              ! -----------
              rr = rr_n4s_p_no__n2_p_o_p_3p25ev ! 3.4e-17   ! R63

              Reaction = &
                   rr * &
                   Neutrals(iN_4S_) * &
                   Neutrals(iNO_)

              NeutralSources(iN2_)  = NeutralSources(iN2_)  + Reaction
              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              NeutralLosses(iNO_)   = NeutralLosses(iNO_)   + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.25
              
               ChemicalHeatingS(ino_n) =  &
                   ChemicalHeatingS(ino_n) + &
                    Reaction * 3.25


              ! -----------
              ! N(2P) -> N(2D) + 10400A
              ! -----------

              rr = rr_n2p__n2d_p_10400a ! 7.9e-2  ! R64

              Reaction = &
                   rr * &
                   Neutrals(iN_2P_)

              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction
              NeutralLosses(iN_2P_)  = NeutralLosses(iN_2P_)  + Reaction

              Emission(iE10400_) = Emission(iE10400_) + Reaction

              ! ----------------------------------------------------------
              ! N(2D)
              ! ----------------------------------------------------------

              ! -----------
              ! N(2D) + O2 -> NO + O(3P) + 3.76 eV
              ! N(2D) + O2 -> NO + O(1D) + 1.80 eV
              ! -----------

              rr = rr_n2d_p_o2__no_p_o3p_p_3p76ev__no_p_o1d_p_1p80ev * exp(-185/tn) ! 9.7e-18 ! R65

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_) * &
                   Neutrals(iO2_)

              if (UseNeutralConstituent(iO_1D_)) then
                 
                 NeutralSources(iO_3P_) = &
                      NeutralSources(iO_3P_) + 0.9 * Reaction
                 NeutralSources(iO_1D_) = &
                      NeutralSources(iO_1D_) + 0.1 * Reaction
                 
                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 3.76 * 0.9 + &
                      Reaction * 1.80 * 0.1
                 
                 ChemicalHeatingS(in2d_o2) =  &
                      ChemicalHeatingS(in2d_o2) + &
                      Reaction * 3.76 * 0.9 + &
                      Reaction * 1.80 * 0.1

              else
                 
                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 3.76
                 
              endif
              NeutralSources(iNO_)   = NeutralSources(iNO_)   + Reaction
              NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              NeutralLosses(iO2_)   = NeutralLosses(iO2_)   + Reaction

              ! -----------
              ! N(2D) + NO -> N2 + O + 5.63 eV
              ! -----------

              rr = rr_n2d_p_no__n2_p_o_p_5p63ev ! 6.7e-17  ! R66
              Reaction = &
                   rr * &
                   Neutrals(iN_2D_) * &
                   Neutrals(iNO_)

              NeutralSources(iN2_)  = NeutralSources(iN2_)  + Reaction
              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              NeutralLosses(iNO_)   = NeutralLosses(iNO_)   + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 5.63

!!! Temp change to stop crash
!              ! ----------------------------------------------------------
!              ! N(2P)
!              ! ----------------------------------------------------------
!			  
!              ! -----------
!              ! N(2P) + e -> N(2D) + e + 1.19 eV
!              ! -----------
!
!              rr = 9.5e-15
!
!              Reaction = &
!                   rr * &
!                   Neutrals(iN_2P_) * &
!                   Ions(ie_)
!
!              ! We create and loose the same amount of e
!              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction
!              NeutralLosses(iN_2P_)  = NeutralLosses(iN_2P_)  + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 1.19
!				   
!              ! -----------
!              ! N(2P) + e -> N(4S) + e + 3.57 eV
!              ! -----------
!
!              rr = 2.04e-16 * (te3m085**(-1))
!
!              Reaction = &
!                   rr * &
!                   Neutrals(iN_2P_) * &
!                   Ions(ie_)
!
!              ! We create and loose the same amount of e
!              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
!              NeutralLosses(iN_2P_)  = NeutralLosses(iN_2P_)  + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 3.57
!				   
!              ! -----------
!              ! N(2P) + NO -> N(4S) + NO + 3.44 eV
!              ! -----------
!
!              rr = 1.8e-16
!
!              Reaction = &
!                   rr * &
!                   Neutrals(iN_2P_) * &
!                   Neutrals(iNO_)
!
!              NeutralSources(iN_4S_)  = NeutralSources(iN_4S_)  + Reaction
!              NeutralLosses(iN_2P_) = NeutralLosses(iN_2P_) + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 3.44
!				   
!              ! -----------
!              ! N(2P) + O(3P) -> N(2D) + O(3P) + 1.19 eV
!              ! -----------
!
!              rr = 1.7e-17
!
!              Reaction = &
!                   rr * &
!                   Neutrals(iN_2P_) * &
!                   Neutrals(iO_3P_)
!
!              NeutralSources(iN_2D_)  = NeutralSources(iN_2D_)  + Reaction
!              NeutralLosses(iN_2P_) = NeutralLosses(iN_2P_) + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 1.19
!				   
!              ! -----------
!              ! N(2P) + O2 -> NO + O(3P) + 4.95 eV
!              ! -----------
!
!              rr = 3.09e-18 * exp(-60/Tn)
!
!              Reaction = &
!                   rr * &
!                   Neutrals(iN_2P_) * &
!                   Neutrals(iO2_)
!
!              NeutralSources(iNO_)  = NeutralSources(iNO_)  + Reaction
!              NeutralSources(iO_3P_)  = NeutralSources(iO_3P_)  + Reaction
!              NeutralLosses(iN_2P_) = NeutralLosses(iN_2P_) + Reaction
!              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 4.95

              ! ----------------------------------------------------------
              ! O(1D)
              ! ----------------------------------------------------------

              if (UseNeutralConstituent(iO_1D_)) then
                 ! ------------
                 ! O(1D) -> O(3P) + 6300A
                 ! ------------
                 
                 rr = rr_o1d__o3p_p_6300a ! 0.0071  ! R67
                 
                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_)
                 
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_)  = NeutralLosses(iO_1D_)  + Reaction
                 
                 Emission(iE6300_) = Emission(iE6300_) + Reaction              
                 
                 
                 ! ------------                                           
                 ! O(1D) -> O(3P) + 6364A                                               
                 ! ------------                                           
                 
                 rr = rr_o1d__o3p_p_6364a ! 0.0022 ! R68
                 
                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_)
                 
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_)  = NeutralLosses(iO_1D_)  + Reaction
                 
                 Emission(iE6364_) = Emission(iE6364_) + Reaction
                 
                 ! ------------                                          
                 ! O(1D) + e -> O(3P) + e + 1.96 eV
                 ! ------------                                             
                 
                 rr = rr_o1d_p_e__o3p_p_e_p_1p96ev * te3m091  ! 2.87e-16 ! R69
                 
                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_) * Ions(ie_)
                 
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_)  = NeutralLosses(iO_1D_)  + Reaction
                 
                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 1.96
                 
                 ChemicalHeatingS(io1d_e) = &
                      ChemicalHeatingS(io1d_e) + &
                      Reaction * 1.96
                 
                 ! ------------
                 ! O(1D) + N2 -> O(3P) + N2 + 1.96 eV
                 ! ------------
                 
                 rr = rr_o1d_p_n2__o3p_p_n2_p_1p96ev * exp(107/Tn)  ! 1.8e-17 ! R70
                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_) * &
                      Neutrals(iN2_)
                 
                 ! We create and loose the same amount of N2
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_) = NeutralLosses(iO_1D_) + Reaction
                 
                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 1.96
                 
                 ChemicalHeatingS(io1d_n2) = &
                      ChemicalHeatingS(io1d_n2) + &
                      Reaction * 1.96

                 ! ------------
                 ! O(1D) + O2 -> O(3P) + O2 + 1.96 eV
                 ! ------------
                 
                 rr = rr_o1d_p_o2__o3p_p_o2_p_1p96ev * exp(67/Tn) ! 3.2e-17 ! R71
                 
                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_) * &
                      Neutrals(iO2_)
                 
                 ! We create and loose the same amount of O2
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_) = NeutralLosses(iO_1D_) + Reaction
                 
                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 1.96
                 
                 ChemicalHeatingS(io1d_o2) = &
                      ChemicalHeatingS(io1d_o2) + &
                      Reaction * 1.96
                 
                 ! ------------
                 ! O(1D) + O(3P) -> O(3P) + O(3P) + 1.96 eV
                 ! ------------
                 
                 rr = rr_o1d_p_o3p__o3p_p_o3p_p_1p96ev * ((Tn/300)**0.14) ! 6.47e-18 ! R72
                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_) * &
                      Neutrals(iO_3P_)
                 
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_) = NeutralLosses(iO_1D_) + Reaction
                 
                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 1.96
                 
                 ChemicalHeatingS(io1d_o) = &
                      ChemicalHeatingS(io1d_o) + &
                      Reaction * 1.96
                 
              endif
              
              ! ----------------------------------------------------------
              ! NO
              ! ----------------------------------------------------------
              ! -----------
              ! NO -> NO+ + e
              ! -----------

!              rr = 6.0e-7

              rr = rr_no__nop_p_e *(1+0.2*(f107-65)/100)*exp(-2.115e-18* &
                   (Neutrals(iO2_)*1.e-6)**0.8855)*szap                        ! 5.88e-7 ! R73
              Reaction = &
                   rr * &
                   Neutrals(iNO_)

              IonSources(iNOP_)   = IonSources(iNOP_)   + Reaction
              NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction
              
              !---- Ions

              if (.not. UseIonChemistry) then
                 IonSources = 0.0
                 IonLosses = 0.0
              else
                 do iIon = 1, nIons-1
                    if (.not.UseIonConstituent(iIon)) then
                       IonSources(iIon) = 0.0
                       IonLosses(iIon) = 0.0
                    endif
                 enddo
              endif

              if (.not. UseNeutralChemistry) then
                 NeutralSources = 0.0
                 NeutralLosses = 0.0
              else
                 do iNeutral = 1, nSpeciesTotal
                    if (.not.UseNeutralConstituent(iNeutral)) then
                       NeutralSources(iNeutral) = 0.0
                       NeutralLosses(iNeutral) = 0.0
                    endif
                 enddo
              endif
              
              ! Take Implicit time step
              Ions(ie_) = 0.0
              do iIon = 1, nIons-1
                 ionso = IonSources(iIon)
                 ionlo = IonLosses(iIon)/(Ions(iIon)+1.0e-6)
                 Ions(iIon) = (Ions(iIon) + ionso * DtSub) / &
                      (1 + DtSub * ionlo)
                 ! sum for e-
                 Ions(ie_) = Ions(ie_) + Ions(iIon)
              enddo

              do iNeutral = 1, nSpeciesTotal

                 neuso = NeutralSources(iNeutral)
                 neulo = NeutralLosses(iNeutral) / (Neutrals(iNeutral)+0.1)

                 !!!
                 if (Neutrals(iNeutral) == 0) &
                      write(*,*) "Neutral is zero : ", iLon,iLat,iAlt,iNeutral

                 Neutrals(iNeutral)=(Neutrals(iNeutral) + neuso * DtSub) / &
                      (1 + DtSub * neulo)

                 NeutralSourcesTotal(ialt,iNeutral) = &
                      NeutralSourcesTotal(ialt,iNeutral) + &
                      NeutralSources(iNeutral) * DtSub

                 NeutralLossesTotal(ialt,iNeutral) = &
                      NeutralLossesTotal(ialt,iNeutral) + &
                      NeutralLosses(iNeutral) * DtSub
                 
              enddo

              ChemicalHeatingRate(iLon,iLat,iAlt) = &
                   ChemicalHeatingRate(iLon,iLat,iAlt) + &
                   ChemicalHeatingSub * DtSub + &
                   ChemicalHeatingSubI * DtSub
				   
              ChemicalHeatingRateIon(iLon,iLat,iAlt) = &
                   0.0
!                   ChemicalHeatingRateIon(iLon,iLat,iAlt) + &
!                   ChemicalHeatingSubI * DtSub

              ChemicalHeatingRateEle(iLon,iLat,iAlt) = &
                   ChemicalHeatingRateEle(iLon,iLat,iAlt) + &
                   ChemicalHeatingSubE * DtSub

              ChemicalHeatingSpecies(iLon,iLat,iAlt,:) = &
                   ChemicalHeatingSpecies(iLon,iLat,iAlt,:) + &
                   ChemicalHeatingS * DtSub

              EmissionTotal = EmissionTotal + Emission(:)*DtSub
              
              DtTotal = DtTotal + DtSub

              if (DtSub < DtMin) DtMin = DtSub

              if (DtSub < 1.0e-9 .and. abs(DtTotal-Dt) > DtSub) then
                 write(*,*) "Chemistry is too fast!!", DtSub

                 ! Check Ions
                 do iIon = 1, nIons
                    write(*,*) "Ion Source/Loss : ", &
                         iIon, IonSources(iIon), IonLosses(iIon)
                 enddo
                 do iNeutral = 1, nSpeciesTotal
                    write(*,*) "Neutral Source/Loss : ", iAlt, &
                         iNeutral, NeutralSources(iNeutral), &
                         NeutralLosses(iNeutral), Neutrals(iNeutral)
                 enddo

                 call stop_gitm("Chemistry is too fast!!")
              endif

              nIters = nIters + 1

           enddo

           IDensityS(iLon,iLat,iAlt,:,iBlock) = Ions
           NDensityS(iLon,iLat,iAlt,:,iBlock) = Neutrals

           Emissions(iLon, iLat, iAlt, :, iBlock) =  &
                Emissions(iLon, iLat, iAlt, :, iBlock) + EmissionTotal
           
           if (DoCheckForNans) then
              do iNeutral=1,nSpeciesTotal
                 if (isnan(Neutrals(iNeutral))) then 
                    write(*,*) "chemistry : Neutral is nan", iLon,iLat,iAlt,iNeutral
                    call stop_gitm("Must stop now.")
                 endif
              enddo
           endif
           
           ChemicalHeating2d(iLon, iLat) =  &
                ChemicalHeating2d(iLon, iLat) + &
                ChemicalHeatingRate(iLon, iLat, iAlt) * &
                Element_Charge * & 
                dAlt_GB(iLon, iLat, iAlt, iBlock)

        enddo ! Alt
     enddo ! Lat
  enddo ! Lon

  ChemicalHeatingRate(:,:,:) = &
       ChemicalHeatingRate(:,:,:) * Element_Charge / &
       TempUnit(1:nLons,1:nLats,1:nAlts) / cp(1:nLons,1:nLats,1:nAlts,iBlock)/&
       rho(1:nLons,1:nLats,1:nAlts,iBlock)
	   
  ChemicalHeatingRateIon(:,:,:) = &
       ChemicalHeatingRateIon(:,:,:) * Element_Charge

  ChemicalHeatingSpecies = ChemicalHeatingSpecies * Element_Charge
  
  if (iDebugLevel > 3) then
     do iIon = 1, nIons
        write(*,*) "====> calc_chemistry: Max Ion Density: ", iIon, &
             maxval(IDensityS(1:nLons,1:nLats,(nAlts*4)/5,iIon,iBlock))
     enddo
  endif

  if (iDebugLevel > 2) &
       write(*,*) "===> calc_chemistry: Average Dt for this timestep : ", &
       (Dt*nLats*nLons*nAlts)/nIters

  call end_timing("calc_chemistry")

end subroutine calc_chemistry
