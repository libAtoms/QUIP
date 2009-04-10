program t
implicit none
      integer stat, n_read
      character(len=1024) line

      open (unit=12, file="test_iostat_file", status="UNKNOWN")
      write(unit=12, fmt=*) "BOB"
      close(unit=12)

      open (unit=12, file="test_iostat_file", status="OLD")
      read (unit=12, fmt='(A)', iostat=stat, advance='no', size=n_read) line
      print *, stat
      read (unit=12, fmt='(A)', iostat=stat, advance='no', size=n_read) line
      print *, stat
      close(unit=12)

end program
