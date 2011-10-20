!     Doug S. Szumski  <d.s.szumski@gmail.com>  01-08-2011
!     Fortran module for processing and analysis software
! 
!     This program is free software; you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation; either version 2 of the License, or
!     (at your option) any later version.
! 
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.


SUBROUTINE s_dev(plat_len_max,pplat,cut_grad,b_toler,toler,max_i,incr,p_dat,p_avg,p_crd,i_in,points)
IMPLICIT NONE

!EXTERNAL
INTEGER, INTENT(IN)    :: points,plat_len_max, pplat
REAL, INTENT(IN)       ::cut_grad, b_toler, toler, max_i, incr
REAL, DIMENSION(points), INTENT(IN):: i_in
REAL, DIMENSION(points), INTENT(OUT):: p_dat, p_avg
INTEGER, DIMENSION(points), INTENT(OUT):: p_crd

!INTERNAL
INTEGER                     :: i,y,plat_len,x, n_plat, counter, counter2
REAL, DIMENSION(3)          :: stats
REAL                        :: dummy, v_pos
LOGICAL                     :: check, fail
LOGICAL, DIMENSION(points)  :: trace_check

!Reset plateau array 
p_dat = 0.0
p_avg = 0.0
p_crd = 0.0
trace_check = .TRUE.  !Array of dim(points)
n_plat = 0
counter = 0
counter2 = 0

!Attempt to fit plateaus of decreasing length to current trace
DO plat_len = plat_len_max, pplat, -1
	!Loop through all data points in current trace
    DO x = 1, points-plat_len_max
		!Check for plateau, start position j, length i: return t/f
		!Calculate coordinates
        y = x + plat_len - 1
        !Check to see if plateau already fitted
        check = .TRUE.
        DO i = x,y
           IF (trace_check(i) .EQV. .FALSE.) THEN 
              check = .FALSE.
           END IF
        END DO
	
        !If no plateaus alreadys fitted to data range then test 
        IF (check .EQV. .TRUE.) THEN
           !Calculate 'plateau' average
           dummy = 0.0
           DO i = x,y
              dummy = dummy + i_in(i)
           END DO
           stats(1) = dummy/plat_len 
           !Calculate 'plateau' post average: Potential bug with upper limit    (check)
           dummy = 0.0
           DO i = y+1, y+pplat
              dummy = dummy + i_in(i)
           END DO
           stats(3) = dummy/pplat
           !Calculate standard deviation of 'plateau'
           dummy = 0.0
           DO i = x,y
              dummy  = dummy + ( i_in(i) - stats(1) )**2
           END DO
           stats(2) = SQRT(dummy/plat_len)
		           
           !Additional plateau criteria
           fail = .FALSE.
           IF (stats(2) < toler*stats(1)) THEN
              !Test plateau gradient 
              dummy = 0.0
              DO i = x, y-1
                 dummy = dummy + ( i_in(i+1) - i_in(i) ) / incr
              END DO
              IF (ABS(dummy) > cut_grad) THEN
	            fail =.TRUE.
                !PRINT*, "FAILED gradient", dummy
              END IF    
              !Get rid of plateaus within b_toler standard deviations of background offset data
              dummy = stats(1) - b_toler*stats(2) 
              IF ( dummy .LE. 0.0 ) THEN
                 fail = .TRUE.
	         !PRINT*, "FAILED plateau within ",b_toler," standard deviations of zero level."
              END IF
              !IF ( 5*stats(2) > stats(1)) THEN 
               !  fail =.TRUE.
              ! END IF
              !Test plateau drop
              !v_pos = stats(1) - 10*stats(2) !Needs variable adding in lieu of '8'
              !IF ( v_pos < stats(3) ) THEN
                 !fail = .TRUE.
	         !PRINT*, "FAILED drop test", v_pos
              !END IF
             !Ignore plateaus with currents above max_i
              IF ( stats(1) > max_i ) THEN
	            fail = .TRUE.
              	!PRINT*, "FAILED max_i", stats(1)
              END IF
           END IF
           !Report plateau pass/fail
           IF ( (stats(2) < toler*stats(1)) .AND. (fail .NEQV. .TRUE.) &
                .AND. (stats(1) < max_i) .AND. (stats(1) > 0.1)) THEN
              !locplat = .TRUE.
              !Copy plateau data points and mark these to prevent re-analysis
              trace_check(x:y) = .FALSE. 
              !Copy plateau data points
              p_dat(x:y) = i_in(x:y)
              !Store plateau average
              p_avg(x:y) = stats(1)
              !Increment global plateau counter
              n_plat = n_plat + 1
              !Store plateau coordinates
              p_crd(x:y) = n_plat
           ELSE
           END IF
        END IF
    END DO
END DO

END SUBROUTINE
