program load_analysis
  implicit none
  integer, parameter :: nproc = 128*5 !Please input the number of processores to be used during the actual calculation of matrix elements. 
  integer i, j, dm1, dm2
  integer nmb_file, nmb_non_zero, counter, bgn_mij, end_mij
  integer*8 total_cycles, chk_cycles, load_indv, dm3
  integer, allocatable :: mij_bgn(:), mij_end(:), nonzero_indv(:), read_data(:,:)
  character (len=28) :: dm5
  character (len=9) :: dm7
  character (len=15) :: dm8
  character (len=100) :: bk_temp_1
  
  open(1, file = "MATRIX_FILE_INDEX.DAT", status = 'old', action = 'read')
    read(1,'(a15,i10)') dm8, nmb_file
    read(1,*)
  allocate(mij_bgn(nmb_file), mij_end(nmb_file))
	do i = 1, nmb_file
	  read(1,*) dm1, mij_bgn(i), mij_end(i)
	  if(dm1 /= i) stop "Error in reading MATRIX_FILE_INDEX.DAT"
	end do
  close(1)
  open(1, file = "MATRIX_NONZERO_INDEX.DAT", status = 'old', action = 'read')
  allocate(nonzero_indv(nmb_file))
	read(1,'(a9,i16)') dm7, dm1
	read(1,'(a28,i16)') dm5, nmb_non_zero
	read(1,*)
	read(1,*)
	if(dm1 /= nmb_file) stop "Error in correct number of files"
	do i = 1, nmb_file
	  read(1,*) dm1, nonzero_indv(i)
	  if(dm1 /= i) stop "Error in reading MATRIX_NONZERO_INDEX.DAT"
	end do
  close(1)
  allocate(read_data(2,nmb_non_zero))
  
  counter = 0
  total_cycles = 0
  open(11, file = "Cycles_before_distrubition_Per_Proc.out")
  do i = 1, nmb_file
    write(bk_temp_1, '(a,i0,a,i0,a)') "Info_", mij_bgn(i), "_", mij_end(i), ".DAT"
    bk_temp_1 = trim(bk_temp_1)
	open(1, file = trim(bk_temp_1), status = 'old', action = 'read')
    dm3 = 0
	do j = 1, nonzero_indv(i)
	  counter = counter + 1
	  read(1,*) read_data(1,counter), read_data(2,counter)
	  total_cycles = total_cycles + read_data(2,counter)
	  dm3 = dm3 + read_data(2,counter)
	end do
	close(1)
	write(11,*)i, dm3
	print*, "Done reading file# = ", i, "out of total #files = ", nmb_file
  end do
  close(11)
  if(counter /= nmb_non_zero) then
    print*, "Error in total non-zero elements", counter, nmb_non_zero
	stop
  end if
  Print*, "Total #cycles = ", total_cycles
  load_indv = total_cycles/nproc
  print*, "Total #procs = ", nproc-1, "Load per processor = ", load_indv
  
  j = 0
  dm2 = 0
  dm3 = 0
  open(2, file = "Load.dat")
  open(3, file = "MATRIX_FILE_INDEX_NEW.DAT")
  write(3,'(a15,i10)')"Total #files = ", nproc
  write(3,'(3(a16,2x))') "            #Mij", "      #Mij_begin", "        #Mij_end"
  do i = 1, nmb_non_zero
    write(bk_temp_1, '(a,i0,a)') "Proc_", j, ".DAT"
    bk_temp_1 = trim(bk_temp_1)
	open(1, file = trim(bk_temp_1), action = 'write')
	write(1,*) read_data(1,i)
	if(dm2 == 0) bgn_mij = read_data(1,i)
	dm2 = dm2 + 1
	dm3 = dm3 + read_data(2,i)
	if(dm3 > load_indv .or. i == nmb_non_zero) then
	end_mij = read_data(1,i)
	write(2,*) dm2
	write(3,'(3(i16,2x))') j+1, bgn_mij, end_mij
	print*, "Done with proc# = ", j, read_data(1,i), dm3, dm2
	dm2 = 0
	dm3 = 0
	j = j + 1
	close(1)
	end if
  end do
  close(2)
  close(3)

end program
