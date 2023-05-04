! change_continuum.F90 - CKDMIP tool for changing the water vapour continuum
! Author: Robin Hogan <r.j.hogan@ecmwf.int>
! Copyright (C) 2023- ECMWF
!
! This program reads line-by-line absorption spectra in the form of
! layer optical depths, and adds or subtracts a water vapour continuum
! contribution

program change_continuum

  use parkind1,         only : jprb, jpib
  use easy_netcdf,      only : netcdf_file
  use caviar_continuum, only : calc_caviar_continuum
  
  implicit none

#include "ckdmip_version.h"

  real(jprb), parameter :: AccelDueToGravity  = 9.80665_jprb ! m s-2
  real(jprb), parameter :: MolarMassAir       = 28.970_jprb  ! g mol-1

  type(netcdf_file)  :: in_file, out_file
  character(len=512) :: in_file_name, out_file_name, argument_str
  character(len=512) :: constituent_id, profile_type, experiment, experiment_id
  character(len=512) :: new_constituent_id = ' '
  character(len=512) :: add_name = ' ', subtract_name = ' '
  character(len=4000) :: history_str

  ! Number and index of command-line arguments, and number of bands
  integer :: narg, iarg, iarg_in, iarg_out, istatus, ncol, jcol, jlev, nwav, nlev

  ! Layer optical depth of a single column, dimensioned (nwav,nlev)
  real(jprb), allocatable :: layer_od(:,:)

  ! Water vapour continuum absorption coefficient (m2/mol)
  real(jprb), allocatable :: continuum(:,:)

  ! Water vapour column per unit area in a layer (mol/m2)
  real(jprb), allocatable :: h2o_mole_per_layer(:,:)
  
  ! Chunk sizes in compression of "something else"
  integer :: chunksizes(3)

  real(jprb), allocatable :: pressure_hl(:,:), mole_fraction_fl(:,:)
  real(jprb), allocatable :: pressure_fl(:,:), temperature_fl(:,:)
  real(jprb), allocatable :: wavenumber_cm1(:)

  ! Molar mass of gas (g mol-1)
  real(jprb) :: molar_mass

  ! Convert optical depth to (1) mass extinction coefficient or (2)
  ! molar extinction coefficient
  integer :: convert_type

  ! Do we copy existing lines/continuum data or not?
  logical :: no_copy = .false.

  ! Loop over arguments
  narg = command_argument_count()

  if (narg < 2) then
    call print_usage_and_exit()
  end if

  iarg = 1

  write(*,'(a)') '***** CHANGING WATER VAPOUR CONTINUUM *****'
 
  do while (iarg <= narg)
  
    call get_command_argument(iarg, argument_str)
    if (trim(argument_str) == '--subtract') then
      iarg = iarg+1
      call get_command_argument(iarg, subtract_name)
      write(*,'(a,a,a)') 'Subtracting water-vapour continuum model "', trim(add_name), '"'
    else if (trim(argument_str) == '--add') then
      iarg = iarg+1
      call get_command_argument(iarg, add_name)
      write(*,'(a,a,a)') 'Adding water-vapour continuum model "', trim(add_name), '"'
    else if (trim(argument_str) == '--no-copy') then
      write(*,'(a)') 'Not copying existing absorption spectrum'
      no_copy = .true.
    else
      exit
    end if
    iarg = iarg+1
    
  end do

  write(*,'(a)') '*******************************************'
  
  iarg_in = iarg
  iarg_out = iarg+1

  call get_command_argument(iarg_in, in_file_name, status=istatus)
  if (istatus /= 0) then
    error stop '*** Error: failed to read name of input file as string of length < 512'
  end if
  
  call get_command_argument(iarg_out, out_file_name, status=istatus)
  if (istatus /= 0) then
    error stop '*** Error: failed to read name of output file as string of length < 512'
  end if

  ! Open input file
  call in_file%open(trim(in_file_name), iverbose=5)

  call in_file%get_global_attribute('constituent_id', constituent_id)

  if (no_copy) then
    if (len_trim(subtract_name) > 0) then
      error stop '*** Error: cannot subtract continuum if --no-copy specified as well'
    end if
    new_constituent_id = 'h2o-continuum'
  else if (len_trim(subtract_name) > 0) then
    if (len_trim(add_name) > 0) then
      ! We are subtracting then adding a continuum
      if (trim(constituent_id) /= 'h2o') then
        error stop '*** Error: constituent_id must be h2o when subtracting and adding a continuum'
      end if
      new_constituent_id = 'h2o'
    else
      ! We are subtracting only
      if (trim(constituent_id) /= 'h2o') then
        error stop '*** Error: constituent_id must be h2o when subtracting a continuum'
      end if
      new_constituent_id = 'h2o-no-continuum'
    end if
  else
    if (len_trim(add_name) > 0) then
      ! We are adding a continuum only
      if (trim(constituent_id) /= 'h2o-no-continuum') then
        error stop '*** Error: constituent_id must be h2o-no-continuum when adding a continuum'
      end if
      new_constituent_id = 'h2o'
    else
      ! We are neither adding nor subtracting a continuum
      new_constituent_id = constituent_id
    end if
  end if

  ! Get global attributes
  call in_file%get_global_attribute('history', history_str)
  call in_file%get_global_attribute('profile_type', profile_type)
  call in_file%get_global_attribute('experiment', experiment)
  call in_file%get_global_attribute('experiment_id', experiment_id)

  ! Create output file
  call out_file%create(trim(out_file_name), is_hdf5_file=.true., iverbose=5)
  call out_file%copy_dimensions(in_file)
  call out_file%copy_variable_definition(in_file, 'level')
  call out_file%copy_variable_definition(in_file, 'half_level')
  call out_file%copy_variable_definition(in_file, 'pressure_hl')
  call out_file%copy_variable_definition(in_file, 'temperature_hl')
  call out_file%copy_variable_definition(in_file, 'mole_fraction_hl')
  call out_file%copy_variable_definition(in_file, 'pressure_fl')
  call out_file%copy_variable_definition(in_file, 'temperature_fl')
  call out_file%copy_variable_definition(in_file, 'mole_fraction_fl')
  !call out_file%copy_variable_definition(in_file, 'reference_surface_mole_fraction')
  call out_file%copy_variable_definition(in_file, 'wavenumber')

  ! Put global attributes
  call out_file%put_global_attributes( &
       &  title_str = "Spectral absorption profiles of " // trim(to_upper(new_constituent_id)), &
       &  source_str= "Line-By-Line Radiative Transfer Model (LBLRTM)", &
       &  project_str="CKDMIP", &
       &  comment_str="LBLRTM computes layer optical depth from half-level pressure, temperature and mole fraction." &
       &  // new_line('a') // "Layer-mean temperature and mole fractions are computed from the " &
       &  //"half-level values assuming a linear variation with pressure.", &
       &  conventions_str="CF-1.7", prior_history_str=history_str)

  call out_file%put_global_attribute('constituent_id', trim(new_constituent_id))
  call out_file%put_global_attribute('profile_type', trim(profile_type))
  call out_file%put_global_attribute('experiment', trim(experiment))
  call out_file%put_global_attribute('experiment_id', trim(experiment_id))

  ! Choose chunk sizes for compression, allowing for the limit on the
  ! chunk size being 2**32 (which is apparently less)
  nwav = in_file%get_outer_dimension('wavenumber')
  nlev  = in_file%get_outer_dimension('level')
  chunksizes = [nwav,nlev,1]
  if (int(nwav,jpib)*int(nlev,jpib) >= int(2,jpib)**30) then
    chunksizes(2) = 1
  end if

  call out_file%define_variable('optical_depth', &
       &  dim3_name='column', dim2_name='level', dim1_name='wavenumber', &
       &  long_name='Layer optical depth', &
       &  deflate_level=2, shuffle=.true., chunksizes=chunksizes)

  call out_file%copy_variable(in_file, 'level')
  call out_file%copy_variable(in_file, 'half_level')
  call out_file%copy_variable(in_file, 'pressure_hl')
  call out_file%copy_variable(in_file, 'temperature_hl')
  call out_file%copy_variable(in_file, 'mole_fraction_hl')
  call out_file%copy_variable(in_file, 'pressure_fl')
  call out_file%copy_variable(in_file, 'temperature_fl')
  call out_file%copy_variable(in_file, 'mole_fraction_fl')
  !call out_file%copy_variable(in_file, 'reference_surface_mole_fraction')
  call out_file%copy_variable(in_file, 'wavenumber')

  call in_file%get('temperature_fl', temperature_fl)
  call in_file%get('pressure_hl', pressure_hl)
  call in_file%get('pressure_fl', pressure_fl)
  call in_file%get('mole_fraction_fl', mole_fraction_fl)
  call in_file%get('wavenumber', wavenumber_cm1)

  ncol = in_file%get_outer_dimension('optical_depth')

  ! Set continuum initially to zero
  allocate(continuum(nwav,nlev))
  continuum = 0.0_jprb

  ! Water vapour in a layer (mol m-2)
  allocate(h2o_mole_per_layer(nlev,ncol))
  ! The factor of 0.001 converts molar mass from g mol-1 to kg mol-1
  h2o_mole_per_layer = mole_fraction_fl * (pressure_hl(2:nlev+1,:) - pressure_hl(1:nlev,:)) &
       &  / (AccelDueToGravity * MolarMassAir * 0.001_jprb)

  ! Loop over atmospheric column in file
  do jcol = 1,ncol
    write(*,'(a,i0)') 'Column ', jcol
    if (.not. no_copy) then
      ! Load one column of optical depth
      call in_file%get('optical_depth', layer_od, jcol)
    else
      if (.not. allocated(layer_od)) then
        allocate(layer_od(nwav, nlev))
      end if
      layer_od = 0.0_jprb
    end if

    ! Compute continuum absorption to SUBTRACT
    if (len_trim(subtract_name) > 0) then
      write(*,'(a,a,a)') '  Subtracting continuum model "', trim(subtract_name), '"'
      if (trim(subtract_name) == 'CAVIAR') then
        call calc_caviar_continuum(nlev, nwav, pressure_fl(:,jcol), temperature_fl(:,jcol), &
             &                     mole_fraction_fl(:,jcol), wavenumber_cm1, continuum)
      else
        write(*,'(a,a,a)') '*** Error: Continuum model "', trim(subtract_name), '" not understood'
        error stop
      end if

      ! Continuum here has units of m2 mol-1: convert to optical depth
      do jlev = 1,nlev
        continuum(:,jlev) = continuum(:,jlev) * h2o_mole_per_layer(jlev,jcol)
      end do

      ! Subtract continuum optical depth
      layer_od = layer_od - continuum
    end if

    ! Compute continuum absorption to ADD
    if (len_trim(add_name) > 0) then
      write(*,'(a,a,a)') '  Adding continuum model "', trim(add_name), '"'
      if (trim(add_name) == 'CAVIAR') then
        call calc_caviar_continuum(nlev, nwav, pressure_fl(:,jcol), temperature_fl(:,jcol), &
             &                     mole_fraction_fl(:,jcol), wavenumber_cm1, continuum)
      else
        write(*,'(a,a,a)') '*** Error: Continuum model "', trim(add_name), '" not understood'
        error stop
      end if

      ! Continuum here has units of m2 mol-1: convert to optical depth
      do jlev = 1,nlev
        continuum(:,jlev) = continuum(:,jlev) * h2o_mole_per_layer(jlev,jcol)
      end do

      ! Add continuum optical depth
      layer_od = layer_od + continuum
    end if

    ! Write column of modified optical depth
    call out_file%put('optical_depth', layer_od, jcol)
  end do

  call in_file%close()
  call out_file%close()

contains

  subroutine print_usage_and_exit()

    write(*,'(a)') 'Usage: change_continuum [--no-copy] [--subtract <model>] [--add <model>] <input.h5> <output.h5>'
    write(*,'(a)') '  Adds and/or subtracts a water-vapour continuum contribution to the single-gas layer input file,'
    write(*,'(a)') '  storing in the output file. Continuum models:'
    write(*,'(a)') '    CAVIAR'
    write(*,'(a)') '    MT_CKD3.5 (not yet available)'
    stop

  end subroutine print_usage_and_exit

  ! Convert lower case string to upper case.
  ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
  ! Original author: Clive Page
  function to_upper(str_in) result(str_out)
    
    implicit none
    
    character(len=*), intent(in) :: str_in
    character(len=len(str_in)) :: str_out
    integer :: i,j
    
    do i = 1, len(str_in)
      j = iachar(str_in(i:i))
      if (j>= iachar("a") .and. j<=iachar("z") ) then
        str_out(i:i) = achar(iachar(str_in(i:i))-32)
      else
        str_out(i:i) = str_in(i:i)
      end if
    end do
    
  end function to_upper
  
end program change_continuum
