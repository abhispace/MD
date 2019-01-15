program rand_pos
implicit none

      integer i,j,nmin,nmax,n_int,part_steps,counter,part_num
      integer counter1,counter2
      integer, parameter :: seed=1000
!      integer, parameter :: part_num=100
      real, parameter :: pi = 22.0/7.0

!      real, parameter :: m = 9.9765e-26				! Mass of SiO2 (dust)
      real, parameter :: m = 1.66e-25					! Mass of Al2O3 (dust)

      real, parameter :: epsilon_0 = 8.854187817e-12
      real, parameter :: ev_to_joule=1.6e-19
      real, parameter :: e_charge=1.6e-19

      real a,l,v_norm,e0,temp,q

      nmin=0
      nmax=1000
      n_int=100
      part_num=0
      counter1=nmin
      counter2=nmax
      part_steps=(nmax-nmin)/n_int
      temp=2.52                                                 ! in eV (probably kb_t)
      temp=temp*ev_to_joule 

      q=(2e4)*e_charge
      
      do counter=1,part_steps
         part_num=part_num+n_int
        
         write(*,*)'Number of particles: ',part_num
        
         l=4*sqrt(real(part_num))
         a=1/sqrt((l**2)/(pi*part_num))					! Wigner-Seitz radius
	 write(*,*)'l: ',l, ' a:',a

         e0=(q*q)/(4*pi*epsilon_0*a)
         v_norm=1/(sqrt(e0/m))
         call create(part_num,a,l,v_norm,temp,m)
      end do

end program rand_pos

subroutine create(part_num,a,l,v_norm,temp,m)
      integer part_num
      CHARACTER(*), PARAMETER :: pos_file_loc = "/abhinav/abhinav/Molecular Dynamics/positions/"
      real, parameter :: pi = 22.0/7.0
      real pos_xarr(part_num), pos_yarr(part_num),dr,dy,dx,a,l,v_theta,v_mag,v_x(part_num),v_y(part_num)
      real temp, m, v_norm
      character*12 :: fname

         do i=1,part_num

            if(part_num.le.99) then
               write(fname,'("pos",I2,".txt")')part_num
            endif

            if((part_num.lt.1000).and.(part_num.ge.100)) then
               write(fname,'("pos",I3,".txt")')part_num
            endif

            if(part_num.ge.1000) then
               write(fname,'("pos",I4,".txt")')part_num
            endif

            open (unit=10,file=pos_file_loc//fname,action="write")

!            write(*,*)'Filename: ',fname
!            if (mod(i,10).eq.0) then
!               write(*,*)'i: ',i
!            endif

            pos_xarr(i)=rand()
            pos_yarr(i)=rand()
      
            pos_xarr(i)=pos_xarr(i)*l
            pos_yarr(i)=pos_yarr(i)*l

            call random_number(v_theta)
            v_theta=v_theta*2*pi 

! Setting same velocity magnitude for all particles
      	    v_mag=sqrt(2*temp/m)*v_norm
!            write(*,*) v_mag

      	    v_x(i)=sin(v_theta)*v_mag
      	    v_y(i)=cos(v_theta)*v_mag

         enddo

  20       do i=1,part_num
             do j=1,part_num
               if(i.ne.j) then
                  dx=pos_xarr(i)-pos_xarr(j)
                  dy=pos_yarr(i)-pos_yarr(j)
                  dr=sqrt(dx*dx+dy*dy)
                  counter1=counter1+1
                  counter2=counter2+1
               do while(dr.le.2*a)
                  counter1=counter1+1
                  counter2=counter2+1
                  pos_xarr(i)=rand()*l
                  pos_yarr(i)=rand()*l
                  dx=pos_xarr(i)-pos_xarr(j)
                  dy=pos_yarr(i)-pos_yarr(j)
                  dr=sqrt(dx*dx+dy*dy)
!                  write(*,*)i,j
                  goto 20

!                  if ((i.eq.25).and.(j.eq.18)) then
!                     write(*,*)'Before: ',pos_xarr(i),pos_yarr(i) 
 
!                     write(*,*)'pos_xarr',j,pos_xarr(j),pos_yarr(j),'dr',dr
!                 write(*,*)'After: ',pos_xarr(i),pos_yarr(i) 
!                 end if

               enddo
               endif
            end do
         end do
!      pos_xarr=pos_xarr*l
!      pos_yarr=pos_yarr*l

      do i=1,part_num
         write (10,'(F12.6,T15,F12.6,T30,F12.6,T45,F12.6,T60,F12.6)') pos_xarr(i),pos_yarr(i),v_x(i),v_y(i),v_mag
      end do
      close (10)
return
end
