	Program MQCT_formatted_matrix_reader
!---------------------------------------------------------------------------------------------------!
! By Dulat Bostan and Dmitri Babikov, Jan. 2026, Cite: MQCT_2026                                    !
! This program is for reading and using matrix elements saved in formated MQCT files.               !
!                                                                                                   !
! Appropriate values for all variables can be found in the file MATRIX_INFO.DAT.                    !
! initial, final - state number, labels all states of the system, including rotational degeneracies !
! ch_ini, ch_fin - channel number, labels non-degenerate states of the system (MQCT channels)       !
! j12_ini, j12_fin - labels degenerate angular momentum states j12 within MQCT channels             !
! m12_ini, m12_fin - labels projections m12 of angular momentum j12 on molecule-molecule axis       !
!                                                                                                   !
! To print one desired matrix element make changes to variables 'initial' and 'final'               !
! or 'ch_ini', 'j12_ini', 'm12_ini', and 'ch_fin', 'j12_fin', 'm12_fin', all after line 115         !
!---------------------------------------------------------------------------------------------------!
	integer i, j, k, initial, final, ini, fin
	integer n_files, n_nonzero, dummy, file_length, MIJ_begin, MIJ_end, n_R, n_ch, n_st, n_mtrx
	integer st_ini, ch_ini, j12_ini, m12_ini, j1_ini, ka1_ini, kc1_ini, j2_ini, ka2_ini, kc2_ini, m12_ind
	integer st_fin, ch_fin, j12_fin, m12_fin, j1_fin, ka1_fin, kc1_fin, j2_fin, ka2_fin, kc2_fin
	integer, allocatable:: st(:), ch(:), j12(:), m12(:), j1(:), ka1(:), kc1(:), j2(:), ka2(:), kc2(:)
	real*8 a, b, c
	real*8, allocatable:: R_com(:), matrix_full(:,:,:,:,:,:)
	character(len=100) dum_char, filename1, filename2

!---------------------------------------------------------------------------------------------------!
! Do not make any changes between two horizontal lines (25-120)
! This is the preparatory part which reads information about quantum states and matrix files

	open(1, file='./MATRIX_TRUNCATED/MATRIX_INFO.DAT', status='old', action='read')
	read(1,*)																						! to skip the line
	read(1,*) dum_char, n_ch																		! Number of channels
	read(1,*) dum_char, n_st																		! Number of states
	read(1,*) dum_char, n_mtrx																		! Number of matrix elements
	read(1,*) dum_char, n_R																			! R grid size
	read(1,*)																						! to skip the header

	j12_max = 0
	m12_max = 0
! Reading and storing state indexes into arrays
	allocate(st(n_st), ch(n_st), j12(n_st), m12(n_st))												! Allocating arrays for molecule state indexes
	allocate(j1(n_st), ka1(n_st), kc1(n_st))														! Allocating arrays for molecule state indexes
	allocate(j2(n_st), ka2(n_st), kc2(n_st))														! Allocating arrays for molecule state indexes
	do i = 1, n_st																					! Do loop over number of j12m12 states
	read(1,*) st(i), ch(i), j12(i), m12(i), j1(i), ka1(i), kc1(i), j2(i), ka2(i), kc2(i)			! Storing molecule state indexes into arrays
	if(j12(i).gt.j12_max) j12_max = j12(i)															! Finding maximum j12 value
	if(m12(i).gt.m12_max) m12_max = m12(i)															! Finding maximum m12 value
	enddo

	print*, 'Salutations from Marquette University!'
	print*, ''
	print*, 'Number of states =', n_st
	print*, 'Number of matrix elements =', n_mtrx
	print*, 'j12_max=', j12_max, ', m12_max=', m12_max

! Allocating full matrix
	allocate(matrix_full(n_ch,n_ch,2*m12_max+1,j12_max+1,j12_max+1, n_R))
	matrix_full = 0.d0
	
! This loop is to skip the transition indexes block
! as this information will be taken from 'Ind_***_***.DAT' files
	read(1,*) 																						! to skip the header
	do i = 1, n_mtrx
	read(1,*)
	enddo

! This loop is to read and store R values of the grid
	read(1,*)																						! to skip the header
	allocate(R_com(n_R))
	do i = 1, n_R																					! Do loop over number of R
	read(1,*) dummy, R_com(i)																		! Storing R center of mass distance
	enddo

	open(2, file='./MATRIX_TRUNCATED/MATRIX_NONZERO_INDEX.DAT', status='old', action='read')
	read(2,*) dum_char, dum_char, n_files 															! Number of files
	read(2,*) dum_char, dum_char, dum_char, dum_char, n_nonzero										! Number of non zero matrix elements
	read(2,*) 																						! to skip the line
	read(2,*) 																						! to skip the header

	open(3, file='./MATRIX_TRUNCATED/MATRIX_FILE_INDEX.DAT', status='old', action='read')
	read(3,*) 																						! to skip the line
	read(3,*)																						! to skip the line

! This is the main part of the code which reads and stores matrix elements into 6D array
! Do loop over number of files. Note: this do loop can be modified to use multiple processors
	do i = 1, 1 !n_files
	read(2,*) dummy, file_length																	! Number of lines in each file
	read(3,*) dummy, MIJ_begin, MIJ_end																! Range of matrix elements in the file
	write(filename1,'(a,i0,a,i0,a)') './MATRIX_TRUNCATED/Ind_', MIJ_begin, '_', MIJ_end, '.DAT'		! Generating the file names for indexes file
	write(filename2,'(a,i0,a,i0,a)') './MATRIX_TRUNCATED/MIJ_', MIJ_begin, '_', MIJ_end, '.DAT'		! Generating the file names for matrix files
! 	write(*,*) 'Reading the file ', filename1
	write(*,*) 'Reading the file ', filename2

	open(4, file=filename1,status='old')															! Opening the file 'Ind_***_***.DAT'
	open(5, file=filename2,status='old')															! Opening the file 'MIJ_***_***.DAT'

	do j = 1, file_length																			! Do loop over the non-zero matrix elements inside the each file
		read(4,*) ini, fin																			! initial and final states
		! Initial state indexes
		ch_ini  = ch(ini)
		j12_ini = j12(ini)
		m12_ini = m12(ini)

		! Final state indexes
		ch_fin 	= ch(fin)
		j12_fin = j12(fin)
		m12_fin = m12(fin)
		m12_ind = m12_ini + m12_max + 1																! Adjusted m12 indexes

		do k = 1, n_R																				! Do loop over number of R
		read(5,*) dummy, dummy, matrix_full(ch_ini,ch_fin,m12_ind,j12_ini+1,j12_fin+1,k)			! Reading matrix into 6D array
		enddo
	enddo
	close(4)
	close(5)
	enddo
    close(1)
	close(2)
	close(3)
	write(*,*) 'Matrix reading is finished.'
	write(*,*)
! Make no changes above this line
!---------------------------------------------------------------------------------------------------!

! The part below can be changed as needed for SACM model
! This is the part for printing one selected matrix element along R grid
	initial = 4
	final = 2
	
	if(initial.lt.final) stop 'Value of Initial state should be higher than Final state'
	if(initial.gt.n_st) stop 'Specified initial state is out of the range'
	if(final.gt.n_st) stop 'Specified final state is out of the range'

! m12 block
	if(m12(initial).ne.m12(final)) stop 'Potential coupling is zero between these two states'
	m12_block = m12(initial)
	m12_ind = m12_block + m12_max + 1

! Initial state indexes
	ch_ini  =  ch(initial)
	j12_ini = j12(initial)
	j1_ini  =  j1(initial)
	ka1_ini = ka1(initial)
	kc1_ini = kc1(initial)
	j2_ini  =  j2(initial)
	ka2_ini = ka2(initial)
	kc2_ini = kc2(initial)

! Final state indexes
	ch_fin 	=  ch(final)
	j12_fin = j12(final)
    j1_fin  =  j1(final)
	ka1_fin = ka1(final)
	kc1_fin = kc1(final)
	j2_fin  =  j2(final)
	ka2_fin = ka2(final)
	kc2_fin = kc2(final)

	open(1999,file='Example.dat',action='write')
	write(1999,'(a7,i5)') 'm12 =', m12_block
	write(1999,'(a15,7(a7))') 'mol-mol_state','j12','j1','ka1','kc1','j2','ka2','kc2'
	write(1999,'(a9,1x,i5,2x,i5)',advance='no') 'initial', ch_ini, j12_ini							! Printing initial states j12 m12 values
	write(1999,'(6(2x,i5))') j1_ini, ka1_ini, kc1_ini, j2_ini, ka2_ini, kc2_ini						! Printing initial states j,ka,kc values
	write(1999,'(a7,3x,i5,2x,i5)',advance='no') 'final', ch_fin, j12_fin							! Printing final states j12 m12 values
	write(1999,'(6(2x,i5))') j1_fin, ka1_fin, kc1_fin, j2_fin, ka2_fin, kc2_fin						! Printing final states j,ka,kc values
	
	write(1999,'(2(2x,a20))') 'R_com(Bohr)','MIJ(Hartree)'

	do k = 1, n_R
	write(1999,'(2x,e20.12)', advance='no') R_com(K)												! Printing R value
	write(1999,'(2x,e20.12)') matrix_full(ch_ini,ch_fin,m12_ind,j12_ini+1,j12_fin+1,k)				! Printing matrix value
	enddo
	write(*,*) 'Writing the output is finished.'
    write(*,*)
	close(1999)
	end program
