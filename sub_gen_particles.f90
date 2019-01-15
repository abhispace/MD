! This subroutine generates part_num number of particles at random locations 
! and gives them random velocities in a simulation square of side l
! Author=Abhinav

subroutine gen_particles(part_num,l,k,a,m,temp,v_norm,dt,steps,en,u)

! Setting constants

      character*12 file_name
      character*25 fname,fname_pos
      CHARACTER(*), PARAMETER :: pos_file_loc = "/abhinav/abhinav/Molecular Dynamics/positions/"
      CHARACTER(*), PARAMETER :: fileplace_pos = "/abhinav/abhinav/Molecular Dynamics/output_data/"
      integer part_num,seed,i,steps,xi,j
      real, parameter :: pi = 22.0/7.0
      real, parameter :: c = 299792458						! Speed of light
      real v_mag,l,pos_xarr(part_num),pos_yarr(part_num),v_x(part_num)
      real v_y(part_num),temp,m,k,v_theta(part_num),a,dx,dy,dr,t
!      real dummy_datax,dummy_datay,pos

      if(part_num.le.99) then
         write(fname,'("pos",I2,".txt")')part_num
      endif

      if((part_num.lt.1000).and.(part_num.ge.100)) then
         write(fname,'("pos",I3,".txt")')part_num
      endif

      if(part_num.ge.1000) then
         write(fname,'("pos",I4,".txt")')part_num
      endif

      open (unit=10,file=pos_file_loc//fname,action="read")

      do i=1,part_num
         read(10,'(F12.6,T15,F12.6,T30,F12.6,T45,F12.6,T60,F12.6)')pos_xarr(i),pos_yarr(i),v_x(i),v_y(i),v_mag
!         if(i.eq.1) then
!            write(*,*)':::',pos_xarr(i)
!         endif
!         write(*,*)pos_xarr(i),pos_yarr(i)
      enddo
      close (10)

! Opening and closing energy.txt because it needs to be appended in the next subroutine
! and thus must be flushed at the beginning of the program
      open (unit=13,file="energy.txt",action="write",status="replace")
      close (13)         

!      pos_xarr=pos_xarr*l
!      pos_yarr=pos_yarr*l

!      call random_number(v_theta)
!      v_theta=v_theta*2*pi 

! Setting same velocity magnitude for all particles
!      v_mag=sqrt(2*temp/m)*v_norm
!      v_x=sin(v_theta)*v_mag
!      v_y=cos(v_theta)*v_mag

      write (*,*) 'Velocity magnitude (same for all particles): ', v_mag

      t=0
      do i=1,steps
         call force(i,part_num,l,k,pos_xarr,pos_yarr,v_x,v_y,temp,m,dt,en,u,t,j,steps,a)
!         t=t+steps*dt
         t=t+steps
         write(*,*)'t::::::', t,' i:::::: ',i
         if(mod(i,100000).eq.0) then
            write(*,*)'step number=',i
         endif
      enddo

!      do i=1,steps
!         call force(part_num,l,kmin,kamx,dk,pos_xarr,pos_yarr,v_x,v_y,temp,m,dt)
!         if(mod(i,1000).eq.0) then
!            write(*,*) 'Progress: ',i,'/',steps
!         endif
!      enddo


!---------Following code prints all positions and finds the max y and x axis position values-----------
!
!      do i=1,part_num
!         if(pos_xarr(i).ge.dummy_datax) then
!            dummy_datax=pos_xarr(i)
!         end if
!         if(pos_yarr(i).ge.dummy_datay) then
!            dummy_datay=pos_yarr(i)
!         end if
!         write (*,*) 'x[',i,'] =',pos_xarr(i)
!         write (*,*) 'y[',i,'] =',pos_yarr(i)
!      end do
!      write (*,*) 'x max =',dummy_datax,'y max =',dummy_datay

return
end
