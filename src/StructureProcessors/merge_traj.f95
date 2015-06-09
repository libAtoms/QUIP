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

#include "error.inc"

! Syntax: merge_traj <outfile> <infiles[@frame]>

program merge_traj
  use libAtoms_module

  integer, parameter  :: FN_MAX = 100

  ! ---

  type(Atoms)      :: at
  type(CInOutput)  :: infile, outfile

  ! ---

  character(FN_MAX)               :: outfn, infn
  character(FN_MAX), allocatable  :: infndescr(:)  

  integer  :: i, j, frame, n
  integer  :: error

  ! ---

  call system_initialise

  n = cmd_arg_count()
  allocate(infndescr(n-1))

  call get_cmd_arg(1, outfn)
  do i = 1, n-1
     call get_cmd_arg(i+1, infndescr(i))
  enddo

  call print("Creating merged trajectory file " // outfn)
  call initialise(outfile, outfn, action=OUTPUT, error=error)
  HANDLE_ERROR(error)
  frame = 0
  do i = 1, n-1
     j = index(infndescr(i), "@")
     if (j == 0) then
        call print("..." // infndescr(i))
        call initialise(infile, infndescr(i), action=INPUT, error=error)
        HANDLE_ERROR(error)
        do j = 1, infile%n_frame
           call read(infile, at, frame=j-1, error=error)
           HANDLE_ERROR(error)
           call write(outfile, at, frame=frame, error=error)
           HANDLE_ERROR(error)
           frame = frame + 1
        enddo
        call finalise(infile)
     else
        infn = infndescr(i)(:j-1)
        read (infndescr(i)(j+1:), *)  j
        call print("..." // trim(infn) // ", frame " // j)
        call initialise(infile, infn, action=INPUT, error=error)
        call read(infile, at, frame=j, error=error)
        HANDLE_ERROR(error)
        call write(outfile, at, frame=frame, error=error)
        HANDLE_ERROR(error)
        frame = frame + 1
        call finalise(infile)
     endif
  enddo
  call finalise(outfile)

  deallocate(infndescr)

  call system_finalise

endprogram merge_traj
