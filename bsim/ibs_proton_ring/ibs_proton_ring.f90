!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! program ibs_proton_ring
!
! Simple program to calculate IBS growth rates for 1 turn in a ring.
!
! Usage: ibs_proton_ring (looks for ibs_proton_ring.in file)
!        ibs_proton_ring file.in
!
!-

program ibs_proton_ring

use beam_mod
use twiss_and_track_mod
use ibs_mod

implicit none

type (lat_struct), target :: lat
type (coord_struct), allocatable :: orbit(:)
type (ibs_sim_param_struct) :: ibs_sim_params
type (ibs_struct) rates1turn

integer :: ix, status
integer :: namelist_file
integer i

character(100) :: lat_name, in_file
character(30), parameter :: r_name = 'ibs_proton_ring'
character(4) :: ibs_formula

real(rp) emit_x, emit_y, sig_z, sig_p, granularity
real(rp) Npart

namelist / ibs_params / &
    lat_name, ibs_formula, Npart, emit_x, emit_y, sig_z, sig_p, ibs_formula, granularity

!------------------------------------------
!Defaults for namelist
lat_name = 'lat.bmad'
ibs_formula = 'bjmt'   ! See ibs_mod

call get_command_argument(1, in_file)

namelist_file = lunget()
print *, 'Opening: ', trim(in_file)
open (namelist_file, file = in_file, status = "old")
read (namelist_file, nml = ibs_params)
close (namelist_file)

! Set params
ibs_sim_params%clog_to_use = 1
ibs_sim_params%formula = ibs_formula

!Parse Lattice
call bmad_parser (lat_name, lat)

lat%param%n_part = Npart

call twiss_and_track(lat, orbit, status)
if (status /= ok$) then
  print *, 'problem with twiss_and_track'
  stop
endif

do i=0,lat%n_ele_track
  lat%ele(i)%a%emit = emit_x
  lat%ele(i)%b%emit = emit_y
  lat%ele(i)%z%sigma_p = sig_p
  lat%ele(i)%z%sigma = sig_z
enddo

call ibs_rates1turn(lat, ibs_sim_params, rates1turn, granularity)

write(*,*) "Ta = ", 1.0d0/rates1turn%inv_Ta/60.0, " (m)"
write(*,*) "Tb = ", 1.0d0/rates1turn%inv_Tb/60.0, " (m)"
write(*,*) "Tz = ", 1.0d0/rates1turn%inv_Tz/60.0, " (m)"

end program
