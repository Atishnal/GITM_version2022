! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

module ModRates

  use ModGITM, only: nLons, nLats, nAlts

  implicit none

!    Temperature Independent (i.e. Static) Reaction Rates

  integer, parameter  :: iRrK4_    =  1
  integer, parameter  :: iRrK5_    =  2
  integer, parameter  :: iRrK6_    =  3
  integer, parameter  :: iRrK7_    =  4
  integer, parameter  :: iRrK8_    =  5
  integer, parameter  :: iRrK9_    =  6
  integer, parameter  :: iRrK10_   =  7
  integer, parameter  :: iRrK16_   =  8
  integer, parameter  :: iRrK17_   =  9
  integer, parameter  :: iRrK18_   = 10
  integer, parameter  :: iRrK21_   = 11
  integer, parameter  :: iRrK22_   = 12
  integer, parameter  :: iRrK23_   = 13
  integer, parameter  :: iRrK24_   = 14
  integer, parameter  :: iRrK26_   = 15
  integer, parameter  :: iRrK27_   = 16
  integer, parameter  :: iRrK32_   = 17
  integer, parameter  :: iRrBeta2_ = 18
  integer, parameter  :: iRrBeta4_ = 19
  integer, parameter  :: iRrBeta6_ = 20
  integer, parameter  :: iRrBeta7_ = 21
  integer, parameter  :: iRrJ1_    = 22
  integer, parameter  :: iRrG9_    = 23
  integer, parameter  :: iRrG12_   = 24
  integer, parameter  :: iRrG13_   = 25
  integer, parameter  :: iRrG21_   = 26
  integer, parameter  :: iRrG22_   = 27
  integer, parameter  :: iRrG24_   = 28
  integer, parameter  :: iRrG25_   = 29
  integer, parameter  :: iRrG26_   = 30
  integer, parameter  :: iRrG27_   = 31
  integer, parameter  :: iRrG30_   = 32
  integer, parameter  :: iRrEA3_   = 33
  integer, parameter  :: iRrEA4_   = 34
  integer, parameter  :: iRrEA5_   = 35
  integer, parameter  :: iRrAlpha3_ = 36

  integer, parameter  :: nRrTempInd = 36

  real :: RrTempInd(nRrTempInd)

!    Temperature Depenedent

  integer, parameter  :: iRrK1_     =  1
  integer, parameter  :: iRrK2_     =  2
  integer, parameter  :: iRrK3_     =  3
  integer, parameter  :: iRrK12_    =  4
  integer, parameter  :: iRrK13_    =  5
  integer, parameter  :: iRrK19_    =  6
  integer, parameter  :: iRrK20_    =  7
  integer, parameter  :: iRrK25_    =  8
  integer, parameter  :: iRrKM12_   =  9
  integer, parameter  :: iRrA1_     = 10
  integer, parameter  :: iRrA2_     = 11
  integer, parameter  :: iRrA3_     = 12
  integer, parameter  :: iRrBeta1_  = 13
  integer, parameter  :: iRrBeta3_  = 14
  integer, parameter  :: iRrBeta5_  = 15
  integer, parameter  :: iRrBeta8_  = 16
  integer, parameter  :: iRrBeta9_  = 17
  integer, parameter  :: iRrBeta9N_ = 18
  integer, parameter  :: iRrBeta17_ = 19
  integer, parameter  :: iRrAlpha5_ = 20
  integer, parameter  :: iRrAlpha7_ = 21
  integer, parameter  :: iRrAlpha8_ = 22
  integer, parameter  :: iRrG19_    = 23

  integer, parameter  :: nRrTempDep = 23

  real :: RrTempDep(nLons, nLats, nAlts, nRrTempInd)

!    Collision Frequencies

  integer, parameter  :: iCfIN_ = 1
  integer, parameter  :: iCfEN_ = 2
  integer, parameter  :: iCfEI_ = 3

  integer, parameter  :: nCf = 3

  real :: Cf(nLons, nLats, nAlts, nRrTempInd)
  
  
  !! Atishnal: declared variables for perturbed reaction rates

  ! R1 ! O+O+M -> O2+M + 5.12 eV  rr=4.7e-33 *((300./tn)**(2))
  real(kind=8) :: rr_o_p_o_p_m__o2_p_m_5p12ev = 4.7E-033

  ! R2 ! O+(2D) + N2 -> N2+ + O + 1.35 eV rr = 1.5e-16 * ti3m055 *
  real(kind=8) :: rr_op2d_p_n2__n2p_p_1p35ev = 1.5000000000000000E-016

  ! R3 ! O+(2P) + N2 -> N2+ + O + 3.05 eV  rr = 2.0e-16 * ti3m055 *
  real(kind=8) :: rr_op2p_p_n2__n2p_p_o_p_3p05ev = 2.0000000000000000E-016

  ! R4 ! N2+ + O2 -> O2+ + N2 + 3.53 eV if (ti<=1000.0) then rr = 5.1e-17 * ti3m116 *
  real(kind=8) :: rr_n2p_p_o2__o2p_p_n2_p_3p53ev_lt =  5.1000000000000003E-017

  ! R5 ! N2+ + O2 -> O2+ + N2 + 3.53 eV if (ti=>1000.0) then rr = 1.26e-17 * ti10m067 *
  real(kind=8) :: rr_n2p_p_o2__o2p_p_n2_p_3p53ev_gt = 1.2600000000000001E-017

  ! R6 ! N2+ + O -> NO+ + N(2D) + 0.70 eV *
  ! -> NO+ + N(4S) + 3.08 eV
  ! if (ti<=1500.0) then rr = 1.33e-16 * ti3m044
  real(kind=8) :: rr_n2p_p_o__nop_p_n2d_p_0p70ev__nop_p_n4s_p_3p08ev_lt = 1.3299999999999999E-016

  ! R7 ! N2+ + O -> NO+ + N(2D) + 0.70 eV *
  ! -> NO+ + N(4S) + 3.08 eV
  ! if (ti=>1500.0) then rr = 6.55e-17 * ti15m023
  real(kind=8) :: rr_n2p_p_o__nop_p_n2d_p_0p70ev__nop_p_n4s_p_3p08ev_gt = 6.5500000000000004E-017

  ! R8 ! N2+ + e -> 2 N(2D) + 1.04 eV ! *
  ! N2+ + e -> 2 N(4S) + 5.77 eV  rr = 2.2e-13 * te3m039
  real(kind=8) :: rr_n2p_p_e__2n2d_p_1p04ev__2n4s_p_5p77ev =  2.2000000000000002E-013

  ! R9 ! N2+ + N(4S) -> N2 + N+ + 2.48 eV  rr = 1.0e-17 *
  real(kind=8) :: rr_n2p_p_n4s__n2_p_np_p_2p48ev = 1.0000000000000001E-017

  ! R10 ! N2+ + O -> O+(4S) + N2 + 1.96 eV rr = 7.0e-18 * ti3m023 *
  real(kind=8) :: rr_n2p_p_o__op4s_p_n2_p_1p96ev = 6.9999999999999997E-018

  ! R11 N2+ + NO -> NO+ + N2 + 6.33 eV rr = 3.6e-16  ! Richards *
  real(kind=8) :: rr_n2p_p_no__nop_p_p_n2_p_6p33ev = 3.5999999999999998E-016

  ! R12 ! O+(4S) + O2 -> O2+ + O + 1.55 eV if (ti<=900.0) then rr = 1.6e-17 * ti3m052
  real(kind=8) :: rr_op4s_p_o2__o2p_p_o_p_1p55ev_lt = 1.6E-017

  ! R13 ! O+(4S) + O2 -> O2+ + O + 1.55 eV if (ti=>900.0) then rr = 9e-18 * ti9m092
  real(kind=8) :: rr_op4s_p_o2__o2p_p_o_p_1p55ev_gt = 9E-018

  ! R14 ! O+(2D) + O2 -> O2+ + O + 4.865 eV  rr = 7.0e-16
  real(kind=8) :: rr_op2d_p_o2__o2p_p_o_p_4p865ev = 7.0e-16

  ! R15 ! O+(2P) + O2 -> O2+ + O + 6.54 eV  rr = 1.3e-16
  real(kind=8) :: rr_op2p_p_o2__o2p_p_o_p_6p54ev = 1.3E-016

  ! R16 O+(2P) + O2 -> O+(4S) + O2 + 5.016 eV rr = 1.3e-16
  real(kind=8) :: rr_op2p_p_o2__op4s_p_o2_p_5p016ev = 1.3E-016

  ! R17 ! N+ + O2 -> O2+ + N(4S) + 2.5 eV if (ti<=1000.0) then rr = 1.925e-16 * ti3m045 *
  real(kind=8) :: rr_np_p_o2__o2p_p_n4s_p_2p5ev_lt = 1.9250000000000000E-016

  ! R18  ! N+ + O2 -> O2+ + N(4S) + 2.5 eV if (ti=>1000.0) then rr = 3.325e-16 *
  real(kind=8) :: rr_np_p_o2__o2p_p_n4s_p_2p5ev_gt = 3.3249999999999998E-016

  ! R19 ! N+ + O2 -> O2+ + N(2D) + 0.1 eV  if (ti<=1000.0) then rr =  0.825e-16 * ti3m045 *
  real(kind=8) :: rr_np_p_o2__o2p_p_n2d_p_0p1ev_lt = 8.2499999999999997E-017

  ! R20 ! N+ + O2 -> O2+ + N(2D) + 0.1 eV  if (ti=>1000.0) then  rr = 1.425e-16 *
  real(kind=8) :: rr_np_p_o2__o2p_p_n2d_p_0p1ev_gt = 1.4250000000000000E-016

  ! R21
  ! O2+ + e -> O(1D) + O(1D) + 3.06 eV if (ti<=1200.0) then rr = 1.95e-13 * te3m07
  ! O2+ + e -> O(3P) + O(1D) + 5.02 eV
  ! O2+ + e -> O(3P) + O(3P) + 6.99 eV
  real(kind=8) :: rr_o2p_p_e__o1d_p_o1d_p_3p06ev_lt = 1.95E-013

  ! R22
  ! O2+ + e -> O(1D) + O(1D) + 3.06 eV if (ti<=1200.0) then rr = 7.39e-14 * te12m056
  ! O2+ + e -> O(3P) + O(1D) + 5.02 eV
  ! O2+ + e -> O(3P) + O(3P) + 6.99 eV
  real(kind=8) :: rr_o2p_p_e__o1d_p_o1d_p_3p06ev_gt = 7.39E-014

  ! R23 ! O2+ + N(4S) -> NO+ + O + 4.21 eV rr = 1.0e-16 ! Richards *
  real(kind=8) :: rr_o2p_p_n4s__nop_p_o_4p21ev = 9.9999999999999998E-017

  ! R24 ! O2+ + N(2D) -> NO+ + O + 6.519 eV rr = 1.8e-16 ! Richards *
   real(kind=8) :: rr_o2p_p_n2d__nop_p_o_p_6p519ev = 1.7999999999999999E-016

  ! R25 ! O2+ + N(2P) -> O2+ + N(4S) + 3.565 eV  rr = 2.2e-17 ! Richards *
  real(kind=8) :: rr_o2p_p_n2p__o2p_p_n4s_p_3p565ev = 2.2000000000000000E-017

  ! R26 ! O2+ + NO -> NO+ + O2 + 2.813 eV rr = 4.5e-16 ! schunk and nagy *
  real(kind=8) :: rr_o2p_p_no_nop_p_o2_p_2p813ev = 4.5000000000000002E-016

  ! R27 ! O+(2D) + O -> O+(4S) + O(3P) + 3.31 eV ! O+(2D) + O -> O+(4S) + O(1D) + 1.35 eV 
  ! rr = 1.0e-17
  real(kind=8) :: rr_op2d_p_o_op4s_p_o3p_p_3p31ev__op4s_p_o1d_p_1p35ev = 1.0E-017

  ! R28 ! O+(2D) + e -> O+(4S) + e + 3.31 eV rr = 6.03e-14 * te3m05
  real(kind=8) :: rr_op2d_p_e__op4s_p_e_p_3p31ev =  6.03E-014

  ! R29 ! O+(2P) + O -> O+(4S) + O + 5.0 eV rr = 4.0e-16
  real(kind=8) :: rr_op2p_p_o__op4s_p_o_p_5p0ev = 4.0E-016

  ! R30 ! O+(2P) + e -> O+(4S) + e + 5.0 eV rr = 3.03e-14 * te3m05
  real(kind=8) :: rr_op2p_p_e__op4s_p_e_p_5p0ev = 3.03E-014

  ! R31 ! O+(2P) -> O+(4S) + 2470A   rr = 0.047
  real(kind=8) :: rr_op2p__op4s_p_2470a = 0.047

  ! R32 ! N+ + O2 -> O+(4S) + NO + 2.31 eV if (ti<=1000.0) then rr = 0.275e-16 * ti3m045 *
  real(kind=8) :: rr_np_p_o2__op4s_p_no_p_2p31ev_lt = 2.7500000000000000E-017

  ! R33 ! N+ + O2 -> O+(4S) + NO + 2.31 eV if (ti=>1000.0) then rr = 0.475e-16 *
  real(kind=8) :: rr_np_p_o2__op4s_p_no_p_2p31ev_gt = 4.7500000000000001E-017

  ! R34 ! O+(4S) + N2 -> NO+ + N(4S) + 1.10 eV if (ti<=1000.0) then rr = 1.2e-18 * ti3m045 *
  real(kind=8) :: rr_op4s_p_n2__nop_p_n4S_p_1p10ev_lt = 1.2000000000000000E-018

  ! R35 ! O+(4S) + N2 -> NO+ + N(4S) + 1.10 eV if (ti=>1000.0) then rr = 7.0e-19 * ti10m212 *
  real(kind=8) :: rr_op4s_p_n2__nop_p_n4S_p_1p10ev_gt = 7.0000000000000003E-019

  ! R36 ! O+(4S) + NO -> NO+ + O + 4.36 eV  rr = 7.0e-19 * ti3m087 *
  real(kind=8) :: rr_op4s_p_no__nop_p_o_p_4p36ev = 7.0000000000000003E-019

  ! R37 ! O+(4S) + N(2D) -> N+ + O + 1.45 eV  rr = 1.3e-16 *
  real(kind=8) :: rr_op4s_p_n2d__np_p_o_p_1p45ev = 1.2999999999999999E-016

  ! R38 ! O+(2P) + e -> O+(2D) + e + 1.69 eV  rr = 1.84e-13 * te3m05 
  real(kind=8) :: rr_op2p_p_e__op2d_p_e_p_1p69ev = 1.84E-013

  ! R39 ! O+(2D) + N2 -> NO+ + N + 4.41 eV  rr = 2.5e-17 *
  real(kind=8) :: rr_op2d_p_n2__nop_p_n_p_4p41ev = 2.4999999999999999E-017

  ! R40 ! O+(2D) + NO -> NO+ + O + 4.37 eV  rr = 1.2e-15 *
  real(kind=8) :: rr_op2d_p_no__nop_p_o_p_4p37ev = 1.2000000000000000E-015

  ! R41 ! O+(2P) -> O+(2D) + 7320A rr = 0.171
  real(kind=8) :: rr_op2p__op2d_p_7320a = 0.171

  ! R42  ! O+(2D) -> O+(4S) + 3726A  rr = 7.7e-5
  real(kind=8) :: rr_op2d__op4s_p_3726a = 7.7E-5

  ! R43  ! He+ + N2 -> N+ + N + He + 0.28 eV rr = 7.8e-10/1.0e6 *
  ! ! Shunk and Nagy R29
  real(kind=8) :: rr_hep_p_n2__np_p_n_p_he_p_0p28ev = 7.7999999999999999E-010

  ! R44  ! He+ + N2 -> N2+ + He ( + ??? eV)  rr = 5.2e-10/1.0e6 *
  ! ! Shunk and Nagy R30
  real(kind=8) :: rr_hep_p_n2__n2p_p_he = 5.1999999999999996E-010

  ! R45 ! He+ + O2 -> O+ O + He ( + ??? eV)  rr = 9.7e-10/1.0e6 
  ! ! Shunk and Nagy R31
  real(kind=8) :: rr_hep_p_o2__o_p_o_p_he = 9.7E-010

  ! R46 ! He+ + e- -> He ( + ??? eV)  rr = 4.8e-12/1.0e6 * te07
  real(kind=8) :: rr_hep_p_em__he =  4.8E-012

  ! R47 ! O+(2P) + N -> N+ + O + 2.7 eV rr = 1.0e-16 *
  real(kind=8) :: rr_op2p_p_n__np_p_o_p_2p7ev = 9.9999999999999998E-017

  ! R48 ! N+ + NO --> N2+ + O + 2.2 eV   rr = 8.33e-17 * ti3m024 *
  real(kind=8) :: rr_np_p_no__n2p_p_o_p_2p2ev = 8.3300000000000017E-017

  ! R49 ! N+ + NO --> NO+ + N(4S) + 3.4 eV  rr = 4.72e-16 * ti3m024 *
  real(kind=8) :: rr_np_p_no__nop_p_n4s_p_3p4ev = 4.7199999999999998E-016
 
  ! R50 ! N+ + O2 --> NO+ + O(3P) + 6.67 eV if (ti<=1000) then rr = 0.495e-16 * ti3m045 *
  real(kind=8) :: rr_np_p_o2__nop_p_o3p_p_6p67ev_lt = 4.9499999999999997E-017

  ! R51 ! N+ + O2 --> NO+ + O(3P) + 6.67 eV if (ti>=1000) then rr = 0.855e-16 line 1617 *
  real(kind=8) :: rr_np_p_o2__nop_p_o3p_p_6p67ev_gt =  8.5499999999999996E-017

  ! R52  ! O+(2D) + N -> N+ + O + 1.0 eV   rr = 1.5e-16*
  real(kind=8) :: rr_op2d_p_n__np_p_o_p_1p0ev = 1.5000000000000000E-016 

  ! R53   ! N+ + O2 -> NO+ + O(1D) + 4.71 eV  if (ti<=1000.0) then  rr = 1.98e-16 * ti3m045 *
  real(kind=8) :: rr_np_p_o2__nop_p_o1d_p_4p71ev_lt = 1.9799999999999999E-016

  ! R54 ! N+ + O2 -> NO+ + O(1D) + 4.71 eV  if (ti=>1000.0) then rr = 3.42e-16 line 1673 *
  real(kind=8) :: rr_np_p_o2__nop_p_o1d_p_4p71ev_gt = 3.4199999999999999E-016

  ! R55 ! N+ + O -> O+ + N + 0.93 eV    rr = 2.2e-18 *
  real(kind=8) :: rr_np_p_o__op_p_n_p_0p93ev = 2.1999999999999998E-018

  ! R56 ! NO+ + e -> O + N(2D) + 0.38 eV  rr = 3.4e-13 * te3m085 *
  real(kind=8) :: rr_nop_p_e__o_p_n2d_p_0p38ev = 3.4000000000000002E-013

  ! R57 ! NO+ + e -> O + N(4S) + 2.77 eV  rr = 0.6e-13 * te3m085 *
  real(kind=8) :: rr_nop_p_e__o_p_n4s_p_2p77ev = 5.9999999999999997E-014

  ! R58  ! N(2D) + e -> N(4S) + e + 2.38 eV  rr = 3.86e-16 * te3m081 *
  real(kind=8) :: rr_n2d_p_e__n4s_p_e_p_2p38ev =  3.8600000000000001E-016

  ! R59 ! N(2D) + O -> N(4S) + O(3P) + 2.38 eV ! N(2D) + O -> N(4S) + O(1D) + 0.42 eV   rr = 6.9e-19 *
  real(kind=8) :: rr_n2d_p_o__n4s_p_o3p_p_2p38ev__n4s_p_o1d_p_0p42ev =  6.8999999999999999E-019

  ! R60  ! N(2D) -> N(4S) + 5200A  rr = 1.06e-5 *
   real(kind=8) :: rr_n2d__n4s_p_5200a = 1.0600000000000000E-005

  ! R61 ! NO -> N(4S) + O   rr=4.5e-6*exp(-1.e-8*(Neutrals(iO2_)*1.e-6)**0.38) *
  real(kind=8) :: rr_no__n4s_p_o = 4.5000000000000001E-006

  ! R62  ! N(4S) + O2 -> NO + O + 1.385 eV   rr = 1.5e-20 * tn * exp(-3270/tn) *
  real(kind=8) :: rr_n4s_p_o2__no_p_o_p_1p385ev = 1.5000000000000001E-020

  ! R63  ! N(4S) + NO -> N2 + O + 3.25 eV   rr = 3.4e-17 *
  real(kind=8) :: rr_n4s_p_no__n2_p_o_p_3p25ev = 3.3999999999999998E-017

  ! R64 ! N(2P) -> N(2D) + 10400A   rr = 7.9e-2 *
  real(kind=8) :: rr_n2p__n2d_p_10400a = 7.9000000000000015E-002

  ! R65 ! N(2D) + O2 -> NO + O(3P) + 3.76 eV *
  ! ! N(2D) + O2 -> NO + O(1D) + 1.80 eV  rr = 9.7e-18 * exp(-185/tn)
  real(kind=8) :: rr_n2d_p_o2__no_p_o3p_p_3p76ev__no_p_o1d_p_1p80ev =  9.6999999999999992E-018

  ! R66  ! N(2D) + NO -> N2 + O + 5.63 eV   rr = 6.7e-17 *
  real(kind=8) :: rr_n2d_p_no__n2_p_o_p_5p63ev = 6.7000000000000004E-017

  ! R67  ! O(1D) -> O(3P) + 6300A  rr = 0.0071
  real(kind=8) :: rr_o1d__o3p_p_6300a = 0.0071

  ! R68  ! O(1D) -> O(3P) + 6364A   rr = 0.0022
  real(kind=8) :: rr_o1d__o3p_p_6364a = 0.0022

  ! R69 ! O(1D) + e -> O(3P) + e + 1.96 eV  rr = 2.87e-16 * te3m091
  real(kind=8) :: rr_o1d_p_e__o3p_p_e_p_1p96ev = 2.87E-016

  ! R70 ! O(1D) + N2 -> O(3P) + N2 + 1.96 eV  rr = 1.8e-17 * exp(107/Tn)
  real(kind=8) :: rr_o1d_p_n2__o3p_p_n2_p_1p96ev = 1.8000000000000000E-017

  ! R71 ! O(1D) + O2 -> O(3P) + O2 + 1.96 eV   rr = 3.2e-17 * exp(67/Tn)
  real(kind=8) :: rr_o1d_p_o2__o3p_p_o2_p_1p96ev = 3.2E-017

  ! R72 ! O(1D) + O(3P) -> O(3P) + O(3P) + 1.96 eV   rr = 6.47e-18 * ((Tn/300)**0.14)
  real(kind=8) :: rr_o1d_p_o3p__o3p_p_o3p_p_1p96ev = 6.47E-018

  ! R73  ! NO -> NO+ + e  rr=5.88e-7*(1+0.2*(f107-65)/100)*exp(-2.115e-18* & (Neutrals(iO2_)*1.e-6)**0.8855)*szap
  real(kind=8) :: rr_no__nop_p_e = 5.8800000000000002E-007
  

contains

  subroutine set_RrTempInd

    RrTempInd(iRrK4_)    = 1.0e-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK5_)    = 4.4E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK6_)    = 4.0E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK7_)    = 2.0E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK8_)    = 1.0E-12 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK9_)    = 6.0E-11 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK10_)   = 1.3E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK16_)   = 4.8E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK17_)   = 1.0E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK18_)   = 4.0E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK21_)   = 0.047 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK22_)   = 0.171 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK23_)   = 8.E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK24_)   = 5.0E-12/ 1.0e6 ! cm-3 to m-3 
    RrTempInd(iRrK26_)   = 7.E-10 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrK27_)   = 7.7E-5 / 1.0e6 ! cm-3 to m-3

    ! Torr et al 79 ( O2+ sink) :
    RrTempInd(iRrK32_)   = 1.8e-10 / 1.0e6 ! cm-3 to m-3

    RrTempInd(iRrBETA2_) = 5.0E-12 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrBETA4_) = 5.0E-13 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrBETA6_) = 7.0E-11 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrBETA7_) = 1.06E-5 / 1.0e6 ! cm-3 to m-3
    RrTempInd(iRrJ1_)    = 5.0e-9 !/ 1.0e6 ! cm-3 to m-3
    !  RrTempInd(J1)    = 0.0

    ! These are all from Reese:
    RrTempInd(iRrG9_)  = 4.8e-10 / 1.0e6
    RrTempInd(iRrG12_) = 6.0e-10 / 1.0e6
    RrTempInd(iRrG13_) = 1.0e-11 / 1.0e6
    RrTempInd(iRrG21_) = 8.0e-13 / 1.0e6
    RrTempInd(iRrG22_) = 7.0e-10 / 1.0e6
    RrTempInd(iRrG24_) = 8.0e-10 / 1.0e6
    RrTempInd(iRrG25_) = 5.2e-11 / 1.0e6
    RrTempInd(iRrG26_) = 1.3e-10 / 1.0e6
    RrTempInd(iRrG27_) = 3.0e-11 / 1.0e6
    RrTempInd(iRrG30_) = 1.0e-10 / 1.0e6
    RrTempInd(iRrEA3_) = 0.047
    RrTempInd(iRrEA4_) = 0.171
    RrTempInd(iRrEA5_) = 7.7e-5
    RrTempInd(iRrAlpha3_) = 1.0e-10 / 1.0e6

  end subroutine set_RrTempInd


end module ModRates
