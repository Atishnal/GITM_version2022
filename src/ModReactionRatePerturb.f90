! Declare a module for reaction rate perturbation
module ModReactionRatePerturb
  
  use ModInputs
  
  implicit none
  
  integer, parameter :: nRatesMax = 50
  integer :: i, ierror
  
 ! Declare variables to store the respective reaction rates
  real(kind=8), dimension(nRatesMax) :: PerturbedRates
  character (len=100), dimension(nRatesMax) :: Reactions
  
! Declare a subroutine to read the csv data file
  contains
    subroutine read_csvfile()
       implicit none

       ! Open the data file
       open(10, file = cReactionRateFile, status = 'old', Iostat = ierror)
         
       ! Check if there was an error opening the csv file
       if (ierror /= 0) Then
          print*, "Error opening file reaction_rates.csv"
          stop
       end if

       ! Read in the date from the csv file
       do i = 1, nRatesMax
        read(10,*) Reactions(i), PerturbedRates(i)
        if (ierror /= 0) then
          print*, "Error reading data from file reaction_rates.csv"
          stop
        end if
       end do  
       close(10)
          
     end subroutine read_csvfile


     subroutine get_reaction_rate(reaction, rate)
       character (len =*), intent(in) :: reaction
       character (len =100)  :: tempstring
       real(kind=8), intent(out) :: rate
       logical :: IsFound = .false.

       call read_csvfile

       ! Compare the input string to the read in strings one-by-one
       do i=1, nRatesMax
          ! if they match, then set the rate
          ! (figure out how to compare strings in fortran!!!):
          tempstring = Reactions(i)
          if (trim(reaction) == trim(tempstring)) then
             rate = PerturbedRates(i)
             IsFound = .true.
             exit
          endif
       enddo

       if (.not.IsFound) then
          rate = -1e32
          write(*,*) "I could not find reaction rate : ", reaction
       endif
       
     end subroutine get_reaction_rate  
end module ModReactionRatePerturb







    
