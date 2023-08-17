! A simple driver routine for testing POP module; could be replaced by a call to POPStep 
! from within host land surface model.
!
! Vanessa Haverd Jan 17, 2014

PROGRAM POP_standalone
USE TypeDef
USE POPModule
USE POP_Types
USE POP_constants


implicit none

integer(i4b), parameter :: nSteps = 450,  NGridCell=1  ! simulation for 450 years and 1 grid-cell
integer(i4b) :: it, k_output, ip
real(dp), DIMENSION(NGridCell,NLayer) :: StemNPP, factor
integer(i4b), DIMENSION(NGridCell,2) ::  disturbance_interval
real(dp), DIMENSION(NGridCell,2) :: disturbance_intensity
integer(i4b) :: out_dens, out_cmass, out_height,out_dens_patch, out_cmass_patch, out_height_patch, out_basalarea, i
integer(i4b) :: out_age_patch
integer(i4b) :: out_biomasslossstress , out_growth, out_diags, out_selfthin, out_biomasslosscrowding, out_biomassloss
CHARACTER(len=1000) :: string1, string2
TYPE ( POP_TYPE )     :: POP
real(dp), DIMENSION(NGridCell) :: LAI 
integer(i4b), PARAMETER :: write_normal_output = 1
integer(i4b) :: count
real(dp) :: meanx, meany
integer(i4b), PARAMETER :: ipMAX = 1


! LAI is only required as in put to POPStep for the purpose of evaluating the projective cover diagnostic
LAI = 2.0
! if NDISTURB = 1 (set in POP_Constants), only complete disturbances are considered, 
! with mean return interval disturbance_interval(1,1). Disturbance intensitites are not used.
! For this test case we use NDISTURB = 1.
! if NDISTURB = 2 (set in POP_Constants), both partial and complete disturbances are considered, 
! with respective mean return intervals disturbance_interval(1,1), disturbance_interval(1,2).
! Partial disturbance intensity is disturbance_intensity(1,1). disturbance_intensity(1,2) is not used. 
disturbance_interval(1,1) = 100. ! mean (catastrophic) disturbance return interval (y)
disturbance_interval(1,2) = 20.  ! mean (partial) disturbance return interval (y)
disturbance_intensity(1,1) = 5. ! Wm-K-1
disturbance_intensity(1,2) = 1000 ! (not used: catastrophic == complete disturbance)

! allocate POP structure for one grid-cell
CALL alloc_POP(POP,1)
! Initialise POP time-step counter
POP%it_pop = 0


! output files
out_dens = 330
out_cmass = 331
out_height = 332
out_biomasslossstress = 328
out_growth = 327
out_basalarea = 326
out_biomasslosscrowding = 325
out_biomassloss = 324

 open (unit=out_growth,file="out_growth.out",status="replace",position="rewind")
 open (unit=out_basalarea,file="basalarea.out",status="replace",position="rewind")
 open (unit=out_biomasslossstress,file="biomasslossstress.out",status="replace",position="rewind")
 open (unit=out_biomasslosscrowding,file="biomasslosscrowding.out",status="replace",position="rewind")
 open (unit=out_biomassloss,file="biomassloss.out",status="replace",position="rewind")
 open (unit=out_dens,file="density.out",status="replace",position="rewind")
 open (unit=out_cmass,file="cmass.out",status="replace",position="rewind")
 open (unit=out_height,file="height.out",status="replace",position="rewind")

 out_dens_patch = 333
 out_cmass_patch = 334
 out_height_patch = 335
 out_diags = 336
 out_age_patch = 337

 open (unit=out_dens_patch,file="density_patch.out",status="replace",position="rewind")
 open (unit=out_cmass_patch,file="cmass_patch.out",status="replace",position="rewind")
 open (unit=out_height_patch,file="height_patch.out",status="replace",position="rewind")
 open (unit=out_diags,file="diags_patch.out",status="replace",position="rewind")
 open (unit=out_age_patch,file="age_patch.out",status="replace",position="rewind")

do ip = 1,ipMAX  ! loop over pixels
	 POP%it_pop = 0
     ! initialise POP variables
     CALL ZeroPOP(pop)
     ! Initialise vector of patches with contrasting prescribed disturbance intervals
     CALL InitPOP1D_Poisson(pop,disturbance_interval)
	 count = 1 ! counter for points to be included on self-thinning curve

	do it = 1,nSteps ! loop over years

	   StemNPP = 0.2 ! annual grid-cell stem biomass increment (kg C m-2 y-1)
	                 ! this could be read in as an annual time series

	   CALL POPStep(POP, StemNPP, disturbance_interval, disturbance_intensity,LAI)

	  if(write_normal_output.eq.1) then
	   if (it.eq.1) then ! output file headers
	   
		  string1 =  '% Year '//trim(pop%pop_grid(1)%bin_labels(1))//' '//trim(pop%pop_grid(1)%bin_labels(2))//' '// &
		  trim(pop%pop_grid(1)%bin_labels(3))//' '//trim(pop%pop_grid(1)%bin_labels(4))//' '// &
		  trim(pop%pop_grid(1)%bin_labels(5))//' '//trim(pop%pop_grid(1)%bin_labels(6))//' '// &
		  trim(pop%pop_grid(1)%bin_labels(7))//' '//trim(pop%pop_grid(1)%bin_labels(8))//' '// &
		  trim(pop%pop_grid(1)%bin_labels(9))//' '//trim(pop%pop_grid(1)%bin_labels(10))//' '// &
		  trim(pop%pop_grid(1)%bin_labels(11))//' '//trim(pop%pop_grid(1)%bin_labels(12)) 

   		  WRITE(out_dens,"(a)") trim(string1)
		  WRITE(out_cmass,"(a)") trim(string1)
		  write(out_height,"(a)") '% Year   Mean_Height     Max_height'

		  string2 =  '% Year '//'Cohort1'//' '//'Cohort2'//' '// &
		  'Cohort3'//' '//'Cohort4'//' '// &
		  'Cohort5'//' '//'Cohort6'//' '// &
		  'Cohort7'//' '//'Cohort8'//' '// &
		  'Cohort9'//' '//'Cohort10'//' '// &
		  'rest'

		  WRITE(out_dens_patch,"(a)") trim(string2)
		  WRITE(out_cmass_patch,"(a)") trim(string2)
		  write(out_height_patch,"(a)") '% Year   trim(string1) '


	   endif ! output file headers

		! write out grid-cell variables (grid-cell # 1)
		WRITE(out_growth,"(i5, 100e16.6)") it, pop%pop_grid(1)%growth
		WRITE(out_basalarea,"(i5, 100e16.6)") it, pop%pop_grid(1)%basal_area
		WRITE(out_dens,"(i5, 100e16.6)") it, pop%pop_grid(1)%densindiv_bin
		WRITE(out_cmass,"(i5, 100e16.6)") it, pop%pop_grid(1)%cmass_stem_bin
		WRITE(out_height,"(i5, 100e16.6)") it, pop%pop_grid(1)%hmean, pop%pop_grid(1)%hmax

		
		! write out patch variables (patch # 6)
		k_output = 6
		WRITE(out_dens_patch,"(i5, 100e16.6)") it, pop%pop_grid(1)%patch(k_output)%Layer(1)%cohort(1:NCOHORT_MAX)%density
		WRITE(out_cmass_patch,"(i5, 100e16.6)") it, pop%pop_grid(1)%patch(k_output)%Layer(1)%cohort(1:NCOHORT_MAX)%biomass
		WRITE(out_age_patch,"(i5, 100e16.6)") it, real(pop%pop_grid(1)%patch(k_output)%Layer(1)%cohort(1:NCOHORT_MAX)%age)
		WRITE(out_biomasslossstress,"(i5, 100e16.6)") it, pop%pop_grid(1)%patch(k_output)%stress_mortality
		!sd edited below line (added &)
		WRITE(out_biomasslosscrowding,"(i5, 100e16.6)") it, pop%pop_grid(1)%patch(k_output)%crowding_mortality, &
				 pop%pop_grid(1)%patch(k_output)%cpc
		WRITE(out_biomassloss,"(i5, 100e16.6)") it, pop%pop_grid(1)%patch(k_output)%mortality

		WRITE(out_diags,"(2i5, 100e16.6)") pop%it_pop, &
				 pop%pop_grid(1)%patch(k_output)%age(1), &
				 pop%pop_grid(1)%patch(k_output)%Layer(1)%density, &
				 pop%pop_grid(1)%patch(k_output)%Layer(1)%biomass, &
				 STEMNPP(1,1)

       endif !write_normal_output
	enddo ! end loop over years
enddo ! end loop over stemNPP

END PROGRAM POP_standalone

