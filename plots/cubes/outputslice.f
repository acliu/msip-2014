	! Takes in a binary data cube from 21cmFAST and outputs a single slice
	! in ASCII format for visualization
	! Written by Adrian Liu, 14th July 2010.
	program outputslice
	implicit none
	integer boxlength,slicenum
	character*60 boxfname,outputfname

	open(2,file='args.dat',status='old')
	read(2,*,end=777,err=777) boxfname,boxlength,slicenum,outputfname
	close(2)
		
	call go(boxfname,boxlength,slicenum,outputfname)
	
	stop
777	call usage
	end
	
	subroutine usage
	print *, 'Outputs (in ASCII) a single slice from a binary data cube'
	print *, 'USAGE: call outputslice.x <box fname> <box length> <slice number> <outputfname>'
	print *, 'EXAMPLE: call outputslice.x box.dat 200 50 slice.dat'
	return
	end	
	
	subroutine go(boxfname,boxlength,slicenum,outputfname)
	implicit none
	integer i,j,k,slicenum,boxlength
	real*4 field(boxlength,boxlength,boxlength)
	character*60 boxfname,outputfname
	
	! Check to see that slice is in the box
	if (slicenum .GT. boxlength) then
		print *, 'Umm...your slice number is beyond the box'
		stop
	endif
	
	! Initialize the field array
	do i=1,boxlength
		do j=1,boxlength
			do k=1,boxlength
				field(i,j,k)=0.0
			enddo
		enddo
	enddo

	! Read in the data cube
	open(72,file=boxfname,form='unformatted',access='direct',recl=4)
	do k=1,boxlength
		do j=1,boxlength
			do i=1,boxlength
				read(72,rec=(k-1)*boxlength*boxlength+(j-1)*boxlength+i) field(i,j,k)
			enddo
		enddo
	enddo
	close(72)
	
	! Output the slice
	open(29,file=outputfname)
	do i=1,boxlength
		write(29,*) (field(i,j,slicenum), j=1,boxlength)
	enddo
	close(29)
	return
	end