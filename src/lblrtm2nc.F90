! lblrtm2nc - convert LBLRTM files to NetCDF4
!
! Copyright (C) 2019- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author: Robin Hogan <r.j.hogan@ecmwf.int>
!
! This program, lblrtm2nc, converts profiles of layerwise spectral
! optical depth from the line-by-line model "LBLRTM" into NetCDF4
! format. It is designed to produce datafiles for the "CKDMIP"
! project, in which spectral optical depths are provided separately
! for individual gases. It takes eight command-line arguments:
!
! 1. MOLEFRAC: The path to a NetCDF file containing the mole fraction
! (volume mixing ratio) profiles of the various gases.
!
! 2. BASENAME: The base name of LBLRTM output files, which should
! include the full path to the files, along with any prefix that the
! files have, e.g. "/path/to/lblrtm/files/lblrtm-spectra". See below
! for how this is used.
!
! 3. MOLECULE: The name of the molecule to be processed (e.g. "h2o").
! The MOLEFRAC file will then be searched for a variable
! MOLECULE_SCENARIO_mole_fraction or MOLECULE_mole_fraction.
!
! 4. SCENARIO: The name of the scenario to use, used for constructing
! LBLRTM filenames to read.
!
! 5. STARTPROF: The number (1 based) of the first profile to process
! of those in the MOLEFRAC file.
!
! 6. ENDPROF: The number (1 based) of the last profile to process of
! those in the MOLEFRAC file.
!
! 7. NBLOCK: Number of blocks into which each spectrum has been
! computed by LBLRTM.
!
! 8. OUTFILE: Output NetCDF4/HDF5 file.
!
! LBLRTM produces a large number of binary files, each one
! corresponding to a single block of the spectrum, a single layer and
! a single atmospheric column. Given the command-line arguments above,
! lblrtm2nc attempts to read in files with a filename of the form:
!
! BASENAME_MOLECULE_SCENARIO_profilePROFILE_layerLAYER_blockBLOCK.dat
!
! PROFILE is an integer from STARTPROF to ENDPROF, if necessary with
! leading zeros in order that it can express numbers from 000 to
! 999. LAYER is the 1-based layer index counting down from
! top-of-atmosphere, also with leading zeros so that numbers from 000
! to 999 can be accommodated.  The number of layers to expect is
! determined from the MOLEFRAC file. BLOCK is a capital letter (A, B,
! C etc) and up to NBLOCK blocks are searched for.
!
! Typically, usage of of lblrtm2nc would be preceded by the creation
! of symbolic links with filenames that conform to the format
! specified above, linking to the original LBLRTM files.

program lblrtm2nc

  use parkind1,    only : jprb, jpib
  use easy_netcdf, only : netcdf_file

  implicit none

  ! Use easy_netcdf.F90 for reading and writing NetCDF files
  type(netcdf_file) :: prof_file, in_file, out_file

  ! Input NetCDF file name
  character(len=512) :: in_file_name

  ! Prefix (including directory) of LBLRTM input files
  character(len=512) :: in_prefix_name

  ! Output file name
  character(len=512) :: out_file_name

  ! Gas name in lower case, scenario name (e.g. "present")
  character(len=32) :: gas_name, scenario_name

  ! Gas name in upper case
  character(len=32) :: gas_upper_name

  ! Temporary string for reading command line arguments to be
  ! converted to integers
  character(len=32) :: tmp_str

  ! Name of full-level and half-level variables
  character(len=32) :: fl_str, hl_str

  ! Most gases have a CF compliant standard name
  character(len=64) :: mole_fraction_standard_name
  logical :: do_mole_fraction_standard_name

  ! Experiment name & ID, and profile type for this gas
  character(len=128) :: experiment
  character(len=32) :: experiment_id, profile_type

  ! String representing spectral block, A-Z
  character :: block_str

  ! Dimension sizes
  integer :: ncol, nlev, nblock, nspec, ispec, ispec_last, ilevrev

  ! Chunk sizes in compression of optical depth
  integer :: chunksizes(3)

  ! Loop indices
  integer :: jcol, jlev, jblock, jwav

  ! Columns to process
  integer :: istartcol, iendcol

  integer :: istatus

  ! Wavenumber and layerwise optical depth
  real(jprb), allocatable :: wavenumber(:) ! cm-1
  real(jprb), allocatable :: optical_depth(:,:)
  ! Pressure on half and full levels (Pa)
  real(jprb), allocatable :: pressure_hl(:,:)
  real(jprb), allocatable :: pressure_fl(:,:)
  ! Temperature on half and full levels (K)
  real(jprb), allocatable :: temperature_hl(:,:)
  real(jprb), allocatable :: temperature_fl(:,:)
  ! Mole fraction (volume mixing ratio) on half and full levels
  real(jprb), allocatable :: mole_fraction_hl(:,:)
  real(jprb), allocatable :: mole_fraction_fl(:,:)
  ! Well-mixed gases are provided with a reference surface
  ! concentration, to allow for subsequent scaling
  real(jprb) :: reference_surface_mole_fraction
  ! (Half-)level number, starting at 1
  real(jprb), allocatable :: half_level(:)

  logical :: is_allocated

  is_allocated = .false.

  ! Check program called with correct number of arguments
  if (command_argument_count() < 8) then
    write(*,'(a,a)') 'Usage: lblrtm2nc <profile_file> <input_prefix> <gas_name> <scenario_code> ', &
         & '<start_profile> <end_profile> <n_blocks> <output_file>'
    write(*,'(a)') ' e.g.: lblrtm2nc profiles.nc /my/data/lblrtm-od co2 present 1 10 4 spectra_co2_present.nc'
    stop
  end if

  ! Read command-line arguments
  call get_command_argument(1, in_file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read profile file name as string of length < 512'
  end if
  write(*,'(a,a)') 'Input file of pressure, temperature and concentrations: ', &
       &           trim(in_file_name)

  call get_command_argument(2, in_prefix_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read input directory name as string of length < 512'
  end if
  write(*,'(a,a)') 'Prefix of LBLRTM input files: ', trim(in_prefix_name)

  call get_command_argument(3, gas_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read gas name as lowercase string of length < 32'
  end if
  write(*,'(a,a)') 'Molecule: ', trim(gas_name)
  gas_upper_name = to_upper(trim(base_name(gas_name)))

  call get_command_argument(4, scenario_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read scenario as a string of length < 32'
  end if
  write(*,'(a,a)') 'Scenario string: ', trim(scenario_name)

  call get_command_argument(5, tmp_str, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read start profile as a string of length < 32'
  end if
  read(tmp_str,*) istartcol
  call get_command_argument(6, tmp_str, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read end profile as a string of length < 32'
  end if
  read(tmp_str,*) iendcol
  write(*,'(a,i0,a,i0)') 'Processing columns from ', istartcol, ' to ', iendcol

  call get_command_argument(7, tmp_str, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read number of spectral blocks as string of length < 32'
  end if
  read(tmp_str,*) nblock
  write(*,'(a,i0)') 'Number of input spectral blocks: ', nblock

  call get_command_argument(8, out_file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Failed to read output file name as string of length < 512'
  end if
  write(*,'(a,a)') 'Output file name: ', trim(out_file_name)

  ! Read file of concentration profiles
  call prof_file%open(trim(in_file_name), iverbose=3)

  ! Read pressures and temperatures at half levels (layer interfaces)
  ! and full levels (layers)
  call prof_file%get('pressure_hl', pressure_hl)
  call prof_file%get('temperature_hl', temperature_hl)
  call prof_file%get('pressure_fl', pressure_fl)
  call prof_file%get('temperature_fl', temperature_fl)

  ! Load mole fraction (volume mixing ratio) at full and half levels,
  ! first looking for a variable with the scenario in its name
  fl_str = trim(base_name(gas_name)) // '_' // trim(scenario_name) // '_mole_fraction_fl'
  if (.not. prof_file%exists(trim(fl_str))) then
    ! If not present, try without the scenario in its name
    fl_str = trim(base_name(gas_name)) // '_mole_fraction_fl'
    hl_str = trim(base_name(gas_name)) // '_mole_fraction_hl'
  else
    hl_str = trim(base_name(gas_name)) // '_' // trim(scenario_name) // '_mole_fraction_hl'
  end if
  call prof_file%get(trim(fl_str), mole_fraction_fl)
  call prof_file%get(trim(hl_str), mole_fraction_hl)
  if (prof_file%attribute_exists(trim(fl_str), "standard_name", 64)) then
    call prof_file%get_attribute(trim(fl_str), "standard_name", mole_fraction_standard_name)
    do_mole_fraction_standard_name = .true.
  else
    do_mole_fraction_standard_name = .false.
  end if
  call prof_file%get_attribute(trim(fl_str), "profile_type", profile_type)

  ! Only well-mixed gases have a reference surface value
  if (prof_file%attribute_exists(trim(fl_str), 'reference_surface_value')) then
    call prof_file%get_attribute(trim(fl_str), 'reference_surface_value', &
         &                       reference_surface_mole_fraction)
    write(*,'(a,e14.5)') 'Reference surface mole fraction = ', reference_surface_mole_fraction
  else
    reference_surface_mole_fraction = -1.0_jprb
    write(*,'(a)') 'No reference surface mole fraction present'
  end if

  call prof_file%get_global_attribute("experiment", experiment)
  call prof_file%get_global_attribute("experiment_id", experiment_id)

  ! Close file of concentrations
  call prof_file%close()

  ! Prepare to read LBLRTM files

  ! Number of columns and levels
  ncol   = iendcol-istartcol+1
  nlev   = size(pressure_hl,1)-1

  allocate(half_level(nlev+1))
  do jlev = 1,nlev+1
    half_level(jlev) = jlev
  end do

  ! Allocate maximum number of wavenumbers
  nspec = 100000000
  allocate(wavenumber(nspec))
  wavenumber = 0.0_jprb

  ! Open output file
  call out_file%create(trim(out_file_name), is_hdf5_file=.true., iverbose=4)

  ! Define dimensions; defer defining wavenumber dimension until its
  ! length is known
  call out_file%define_dimension("column", ncol)
  call out_file%define_dimension("level", nlev)
  call out_file%define_dimension("half_level", nlev+1)

  ! Global attributes
  call out_file%put_global_attributes( &
       &  title_str = "Spectral optical depth profiles of " // trim(title_str(gas_name)), &
       &  inst_str  = "European Centre for Medium-Range Weather Forecasts (ECMWF)", &
       &  source_str= "Line-By-Line Radiative Transfer Model (LBLRTM)", &
       &  creator_name="Robin Hogan", &
       &  creator_email_str="r.j.hogan@ecmwf.int", &
       &  contributor_name="Marco Matricardi", &
       &  project_str="CKDMIP", &
       &  comment_str="LBLRTM computes layer optical depth from half-level pressure, temperature and mole fraction." &
       &  // new_line('a') // "Layer-mean temperature and mole fractions are computed from the " &
       &  //"half-level values assuming a linear variation with pressure.", &
       &  conventions_str="CF-1.7")
  call out_file%put_global_attribute("constituent_id", gas_name)
  call out_file%put_global_attribute("profile_type", profile_type)
  call out_file%put_global_attribute("experiment", experiment)
  call out_file%put_global_attribute("experiment_id", experiment_id)

  ! Define variables
  call out_file%define_variable("level", dim1_name="level", &
       &   long_name="Full-level (or layer) number", data_type_name='short')
  call out_file%put_attribute("level", "positive", "down")

  call out_file%define_variable("half_level", dim1_name="half_level", &
       &   long_name="Half-level (or layer interface) number", data_type_name='short')
  call out_file%put_attribute("half_level", "positive", "down")

  call out_file%define_variable("pressure_hl", &
       &   dim2_name="column", dim1_name="half_level", &
       &   units_str="Pa", long_name="Pressure at half levels", &
       &   standard_name="air_pressure")

  call out_file%define_variable("temperature_hl", &
       &   dim2_name="column", dim1_name="half_level", &
       &   units_str="K", long_name="Temperature at half levels", &
       &   standard_name="air_temperature")

  if (do_mole_fraction_standard_name) then
    call out_file%define_variable("mole_fraction_hl", &
         &   dim2_name="column", dim1_name="half_level", &
         &   units_str="1", long_name=trim(gas_upper_name)//" mole fraction at half levels", &
         &   standard_name=trim(mole_fraction_standard_name))
  else
    call out_file%define_variable("mole_fraction_hl", &
         &   dim2_name="column", dim1_name="half_level", &
         &   units_str="1", long_name=trim(gas_upper_name)//" mole fraction at half levels")
  end if

  call out_file%define_variable("pressure_fl", &
       &   dim2_name="column", dim1_name="level", &
       &   units_str="Pa", long_name="Layer-mean pressure", &
       &   standard_name="air_pressure")

  call out_file%define_variable("temperature_fl", &
       &   dim2_name="column", dim1_name="level", &
       &   units_str="K", long_name="Layer-mean temperature", &
       &   standard_name="air_temperature")

  if (do_mole_fraction_standard_name) then
    call out_file%define_variable("mole_fraction_fl", &
         &   dim2_name="column", dim1_name="level", &
         &   units_str="1", long_name="Layer-mean "//trim(gas_upper_name)//" mole fraction", &
         &   standard_name=trim(mole_fraction_standard_name))
  else
    call out_file%define_variable("mole_fraction_fl", &
         &   dim2_name="column", dim1_name="level", &
         &   units_str="1", long_name="Layer-mean "//trim(gas_upper_name)//" mole fraction")
  end if

  if (reference_surface_mole_fraction > 0.0_jprb) then
    call out_file%define_variable("reference_surface_mole_fraction", &
         &   units_str="1", long_name="Reference surface "//trim(gas_upper_name)//" mole fraction")
    if (do_mole_fraction_standard_name) then
      call out_file%put_attribute("reference_surface_mole_fraction", "standard_name", &
           &  trim(mole_fraction_standard_name))
    end if
  end if

  ! Put temperature, pressure and mole fraction etc
  call out_file%put("level", half_level(1:nlev))
  call out_file%put("half_level", half_level)
  call out_file%put("pressure_hl", pressure_hl(:,istartcol:iendcol))
  call out_file%put("temperature_hl", temperature_hl(:,istartcol:iendcol))
  call out_file%put("mole_fraction_hl", mole_fraction_hl(:,istartcol:iendcol))
  call out_file%put("pressure_fl", pressure_fl(:,istartcol:iendcol))
  call out_file%put("temperature_fl", temperature_fl(:,istartcol:iendcol))
  call out_file%put("mole_fraction_fl", mole_fraction_fl(:,istartcol:iendcol))
  if (reference_surface_mole_fraction > 0.0_jprb) then
    call out_file%put("reference_surface_mole_fraction", reference_surface_mole_fraction)
  end if

  ! Loop over input files
  
  do jcol = istartcol,iendcol
    do jlev = 1,nlev

      if (jlev == 1 .and. jcol == istartcol) then
        ! Pre-read one set of blocks so we know how much space to
        ! allocate for optical_depth

        ispec = 0
        wavenumber = 0.0
        do jblock = 1,nblock
          ! Create capital letter where 1=A, 2=B ... 26=Z, 27=a, 28=b etc
          if (jblock <= 26) then
            block_str = achar(iachar('A')-1+jblock)
          else
            block_str = achar(iachar('a')-1-26+jblock)
          end if
          write(in_file_name,'(a,a,a,a,a,a,i3.3,a,i3.3,a,a,a)') &
               &  trim(in_prefix_name), '_', &
               &  trim(gas_name), '_', trim(scenario_name),  &
               &  '_profile', jcol, '_layer', jlev, '_block', block_str, '.dat'
          write(*,'(a,a)') 'Pre-reading ', trim(in_file_name)
          ilevrev = nlev+1-jlev
          ispec_last = ispec
          call read_lblrtm_file(trim(in_file_name), ispec, wavenumber(:))
          write(*,'(a,i0,a,f16.8,a,f16.8,a)') '  ', ispec-ispec_last, ' values in the range ', &
               &  wavenumber(ispec_last+1), '-', wavenumber(ispec), ' cm-1'
        end do
        nspec = ispec
        write(*,'(a,i0)') 'Total number of spectral values: ', nspec

        ! Define wavenumber-dependent variables on the first call
        call out_file%define_dimension("wavenumber", nspec)

        call out_file%define_variable("wavenumber", &
             &   dim1_name="wavenumber", units_str="cm-1", long_name="Wavenumber", &
             &   is_double=.true., deflate_level=2, shuffle=.true.)
        
        ! Choose chunk sizes for compression, allowing for the limit
        ! on the chunk size being 2**32 (which is apparently less)
        chunksizes = [nspec,nlev,1]
        if (int(nspec,jpib)*int(nlev,jpib) >= int(2,jpib)**30) then
          chunksizes(2) = 1
        end if

        print *, "chunksizes = ", chunksizes

        call out_file%define_variable("optical_depth", &
             &   dim3_name="column", dim2_name="level", dim1_name="wavenumber", &
             &   long_name="Layer optical depth", &
             &   deflate_level=2, shuffle=.true., chunksizes=chunksizes)

        allocate(optical_depth(nspec,nlev))
        optical_depth = 0.0_jprb

        call out_file%put("wavenumber", wavenumber(1:nspec))
      end if

      ispec = 0
      wavenumber = 0.0
      do jblock = 1,nblock
        ! Create capital letter where 1=A, 2=B ... 26=Z, 27=a, 28=b etc
        if (jblock <= 26) then
          block_str = achar(iachar('A')-1+jblock)
        else
          block_str = achar(iachar('a')-1-26+jblock)
        end if
        write(in_file_name,'(a,a,a,a,a,a,i3.3,a,i3.3,a,a,a)') &
             &  trim(in_prefix_name), '_', &
             &  trim(gas_name), '_', trim(scenario_name),  &
             &  '_profile', jcol, '_layer', jlev, '_block', block_str, '.dat'
        write(*,'(a,a)') 'Reading ', trim(in_file_name)
        ilevrev = nlev+1-jlev
        call read_lblrtm_file(trim(in_file_name), ispec, wavenumber(:), optical_depth(:,ilevrev))
      end do

    end do
    call out_file%put("optical_depth", optical_depth(1:nspec,:), jcol-istartcol+1)
  end do

  call out_file%close()

contains

  ! Read a single LBLRTM file (containing one block of one level of
  ! one column)
  subroutine read_lblrtm_file(file_name, ispec, wavenumber, optical_depth)
    
    use parkind1,    only : jprb, jpib
    
    implicit none
    
    character(*), intent(in)  :: file_name
    integer,    intent(inout) :: ispec
    real(jprb), intent(inout) :: wavenumber(:)
    real(jprb), intent(inout), optional :: optical_depth(:)
    
    character*80 :: header
    real(jprb)   :: temp1(2500)
    
    real(jprb)   :: pv1, pv2, pdv
    integer      :: np_arr, ip, ifreq, ipan, nskip
    
    real(jprb)   :: last_wavenumber, first_wavenumber
    
    open(1,file=file_name,form='unformatted',status='old')
    read(1) header
    ifreq  = 0
    
    first_wavenumber = 0.0_jprb
    
    do ipan = 1,6000
      np_arr = 0
      read(1) pv1,pv2,pdv,np_arr
      if (ifreq == 0 .and. ispec > 0) then
        ! Not the first block: find best place to start
        first_wavenumber= pv1
        !print *, 'wav: ', last_wavenumber, first_wavenumber
        ! Include tolerance of 0.1*dv
        ifreq = ispec
        do while (wavenumber(ifreq) >= first_wavenumber-pdv*0.1_jprb)
          ifreq = ifreq - 1
        end do
      end if
      if (np_arr > 0) then
        read(1,end=991) temp1(1:np_arr)
        do ip = 1,np_arr
          ! Don't store zero wavenumber
          if (ip > 1 .or. pv1 > 0.0_jprb) then
            ifreq = ifreq + 1
            wavenumber(ifreq) = pv1+pdv*(ip-1)
            if (present(optical_depth)) then
!!! FIX ME !!! REMOVE 1.0e6
!!! optical_depth(ifreq) = 1.0e6_jprb*max(0.0_jprb,temp1(ip))
              optical_depth(ifreq) = max(0.0_jprb,temp1(ip))
            end if
          end if
        end do
      else 
        exit
      end if
    end do
    
991 continue
    
    ispec = ifreq
    
    close(1)
    
  end subroutine read_lblrtm_file
  
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
      if (j >= iachar("a") .and. j <= iachar("z") ) then
        str_out(i:i) = achar(iachar(str_in(i:i))-32)
      else
        str_out(i:i) = str_in(i:i)
      end if
    end do
    
  end function to_upper
  
  ! awk -F- '{print $1}', i.e. for input "h2o-no-continuum" return
  ! "h2o"
  function base_name(str_in) result(str_out)
      
    implicit none

    character(len=*), intent(in) :: str_in
    character(len=len(str_in)) :: str_out
    integer :: i,j
  
    str_out = " "

    do i = 1, len(str_in)
      j = iachar(str_in(i:i))
      if (j == iachar("-")) then
        exit
      else
        str_out(i:i) = str_in(i:i)
      end if
    end do

  end function base_name

  ! For input "h2o-no-continuum" return "H2O (no continuum)"
  function title_str(str_in) result(str_out)

    implicit none

    character(len=*), intent(in) :: str_in
    character(len=len(trim(str_in))+2) :: str_out
    integer :: iin, iout, j

    str_out = to_upper(trim(base_name(str_in)))

    if (len(trim(str_out)) < len(trim(str_in))-1) then
      iout = len(trim(str_out))+1
      str_out(iout:iout+1) = " ("
      iout = iout+2
      do iin = iout-1,len(trim(str_in))
        j = iachar(str_in(iin:iin))
        if (j == iachar("-")) then
          str_out(iout:iout) = " "
        else
          str_out(iout:iout) = str_in(iin:iin)
        end if
        iout = iout+1
      end do
      str_out(iout:iout) = ")"
    end if

  end function title_str

end program lblrtm2nc
