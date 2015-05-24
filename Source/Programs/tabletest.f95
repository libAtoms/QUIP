! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! H0 X
! H0 X   libAtoms+QUIP: atomistic simulation library
! H0 X
! H0 X   Portions of this code were written by
! H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! H0 X
! H0 X   Copyright 2006-2010.
! H0 X
! H0 X   These portions of the source code are released under the GNU General
! H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
! H0 X
! H0 X   If you would like to license the source code under different terms,
! H0 X   please contact Gabor Csanyi, gabor@csanyi.net
! H0 X
! H0 X   Portions of this code were written by Noam Bernstein as part of
! H0 X   his employment for the U.S. Government, and are not subject
! H0 X   to copyright in the USA.
! H0 X
! H0 X
! H0 X   When using this software, please cite the following reference:
! H0 X
! H0 X   http://www.libatoms.org
! H0 X
! H0 X  Additional contributions by
! H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
! H0 X
! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

program tabletest

  use libAtoms_module

  implicit none

  type(Table) :: t, t1, t2
  type(Inoutput) :: f
  integer :: i
  character(TABLE_STRING_LENGTH), allocatable, dimension(:) :: tmp_s1
  character(TABLE_STRING_LENGTH), allocatable, dimension(:,:) :: tmp_s2
  logical, allocatable, dimension(:,:) :: tmp_log

  call system_initialise

  call print('empty table')
  call print(t)
  call print('')
  call finalise(t)

  call print('single int table')
  call append(t, 1)
  call append(t, 2)
  call append(t, 3)
  call append(t, (/4,5,6/))
  call print(t)
  call print('int_part(t)')
  call print(int_part(t))
  call finalise(t)
  call print('')

  call print('single real table')
  call append(t, 1.0_dp)
  call append(t, 2.0_dp)
  call append(t, 3.0_dp)
  call append(t, realpart=(/4.0_dp,5.0_dp,6.0_dp/))
  call print(t)
  call print('real_part(t)')
  call print(real_part(t))
  call finalise(t)
  call print('')

  call print('single int and single real table')
  call append(t, 1, 1.0_dp)
  call append(t, 2, 2.0_dp)
  call append(t, 3, 3.0_dp)
  call append(t, intpart=(/4,5,6/), realpart=(/4.0_dp,5.0_dp,6.0_dp/))
  call append(t, intpart=1, realpart=2.0_dp)
  call print(t)
  call print('int_part(t)')
  call print(int_part(t))
  call print('real_part(t)')
  call print(real_part(t))
  call finalise(t)
  call print('')

  call print('single string table')
  call append(t, 'str1      ')
  call append(t, 'str2      ')
  call append(t, 'str3      ')
  call append(t, 'str4      ')
  call print(t)
  call print('str_part(t)')
  allocate(tmp_s2(1,t%n))
  tmp_s2 = str_part(t)
  do i=1,t%N
     call print(tmp_s2(1,i))
  end do
  deallocate(tmp_s2)
  call finalise(t)
  call print('')

  call print('single logical table')
  call append(t, .false.)
  call append(t, .false.)
  call append(t, .true.)
  call append(t, .false.)
  call append(t, .false.)
  call append(t, .false.)
  call print(t)
  call finalise(t)
  call print('')

  call print('multiple int and single string table')
  call append(t, intpart=(/1,2,3/),strpart=(/'0123456789'/))
  call append(t, intpart=(/4,5,6/),strpart=(/'abcdefghij'/))
  call append(t, intpart_2d=reshape((/7,8,9,10,11,12/),(/3,2/)), &
       strpart_2d=reshape((/'hello     ','world     '/),(/1,2/)))
  call print(t)
  call print('int_part(t)')
  call print((int_part(t)))
  call print('int_part(t,1)')
  call print((int_part(t,1)))
  call print('str_part(t)')
  allocate(tmp_s2(1,t%N))
  tmp_s2 = str_part(t)
  do i=1, t%N
     call print(tmp_s2(1,i))
  end do
  deallocate(tmp_s2)
  call print('str_part(t,1)')
  allocate(tmp_s1(t%N))
  tmp_s1 = str_part(t,1)
  do i=1,t%N
     call print(tmp_s1(i))
  end do
  deallocate(tmp_s1)
  call finalise(t)
  call print('')

  call print('multiple real and multiple logical table')
  call append(t, realpart=(/1.0_dp,2.0_dp,3.0_dp/),logicalpart=(/.false.,.true.,.false./))
  call append(t, realpart_2d=reshape((/4.0_dp,5.0_dp,6.0_dp, 7.0_dp,8.0_dp,9.0_dp/),(/3,2/)),&
       logicalpart_2d=reshape((/.false.,.false.,.false.,.true.,.true.,.true./),(/3,2/)))

  call print(t)
  call print('real_part(t)')
  call print((real_part(t)))
  call print('real_part(t,(/1,2/))')
  call print((real_part(t,(/1,2/))))
  call print('logical_part(t)')
  allocate(tmp_log(3,3))
  tmp_log = logical_part(t)
  do i=1,t%N
     call print(tmp_log(:,i))
  end do
  deallocate(tmp_log)

  call print((logical_part(t,2)))
  call finalise(t)
  call print('')

  call print('appending tables to tables')
  call allocate(t1, 2, 1, 0, 0)
  call append(t1, (/1,2/), 3.0_dp)
  call append(t1, (/4,5/), 6.0_dp)
  call append(t1, (/7,8/), 9.0_dp)
  call print('t1=')
  call print(t1)
  call print('t2=')
  call append(t2, t1)
  call append(t2, t1)
  call print(t2)
  call finalise(t1)
  call finalise(t2)

  call allocate(t1, 0, 0, 1, 0)
  call append(t1, 'hello     ')
  call append(t1, 'world     ')
  call print('t1=')
  call print(t1)
  call print('t2=')
  call append(t2, t1)
  call append(t2, t1)
  call print(t2)
  call finalise(t1)
  call finalise(t2)


  call allocate(t1, 0, 0, 1, 2)
  call append(t1, strpart=(/'hello     '/),logicalpart=(/.true.,.false./))
  call append(t1, strpart=(/'world     '/),logicalpart=(/.false.,.true./))
  call print('t1=')
  call print(t1)
  call print('t2=')
  call append(t2, t1)
  call append(t2, t1)
  call print(t2)
  call finalise(t1)

  call print('testing append_column')
  call append_column(t2, 0.0_dp)
  call print(t2)
  call append_column(t2, .true., n_cols=2)
  call print(t2)
  call append_column(t2, '0123456789', n_cols=3)
  call print(t2)
  call append_column(t2, (/.true.,.true.,.true.,.false./))
  call print(t2)
  call append_column(t2, (/1,2,3,4/))
  call print(t2)
  call print('')

  call print('testing subtable')
  t1 = subtable(t2, (/(i,i=1,t2%N-1)/))
  call print(t1)
  call print('')

  call print('testing remove_columns')
  call remove_columns(t2,str_col_min=1)
  call print(t2)
  call remove_columns(t2,logical_col_min=1,logical_col_max=5)
  call print(t2)
  call remove_columns(t2,str_col_min=2)
  call print(t2)
  call remove_columns(t2,real_col_min=1)
  call print(t2)
  call print('')

  call print('testing insert')
  call insert(t2, 3, intpart=(/5/), strpart=(/'alpha     ','beta      '/))
  call print(t2)

  call print('testing find')
  i = find(t2, 5)
  call print('found 5 at position '//i)
  i = find(t2, (/5/))
  call print('found (/5/) at position '//i)
  call print('')

  call print('testing sort')
  call sort(t2)
  call print(t2)
  call print('')

  call print('testing binary search')
  i = search(t2, (/3/))
  call print('found (/5/) at position '//i)

  t1 = t2 
  call print('testing delete')
  call print('deleted first row:')
  call delete(t2, 1, keep_order=.false.) ! by index
  call print('deleted row matching (/5/)')
  call delete(t2, (/5/), keep_order=.false.) ! by value
  call print(t2)
  call print('deleted rows 1 and 2')
  call delete_multiple(t2, (/1,2/))
  call print(t2)
  call wipe(t2)
  call print('wiped:')
  call print(t2)
  call print('')

  call print('copy taken before deletions:')
  call print(t1)

  call finalise(t1)

  call system_finalise
  

end program tabletest
