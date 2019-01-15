! This subroutine calculates the pairwise potential energy
! and total energies of each particle. It also calculates the force 
! acting on each particle. It then calls integrator function. 
! Author=Abhinav

subroutine force(count_t,part_num,l,k,pos_xarr,pos_yarr,v_x,v_y,kb_t,m,dt,en,u,t,xi,steps,a)

! Declaring and Setting constants

      integer part_num,i,steps,j,xi
      CHARACTER(*), PARAMETER :: fileplace = "/abhinav/abhinav/Molecular Dynamics/output_data/"
      real v_mag,dt,l,pos_xarr(part_num),pos_yarr(part_num)
      real x_prev(part_num),y_prev(part_num),v_x(part_num),v_y(part_num),kb_t,m,a
      real fx(part_num),fy(part_num),xr,yr,r2,r,maxr,maxxr,maxyr,k,u,ke,en,v(part_num),f_x,f_y,t

      open (unit=11,file=fileplace//"init_force.txt",action="write",status="replace")
!      open (unit=13,file="energy.txt",access="append",status="unknown")

      en=0.0
      u=0.0
!      ke=0.0

      maxr=0
      maxxr=0
      maxyr=0

! Setting forces to zero
!      write(*,*)"##",pos_xarr(1)
      do i=1,part_num
         fx(i)=0
         fy(i)=0
         v(i)=sqrt((v_x(i)**2)+(v_y(i)**2))
         x_prev(i)=pos_xarr(i)-v_x(i)*dt
!         if(i.eq.1) then
!            write(*,*)'>>>>',steps,i,x_prev(i),pos_xarr(i),v_x(i),v_y(i)
!         endif
         y_prev(i)=pos_yarr(i)-v_y(i)*dt
!         write(*,*) ' v(',i,')=',v(i), ' vx=',v_x(i), ' vy=',v_y(i)
      enddo

! Periodic boundary conditions
      do i=1,part_num-1
         do j=i+1,part_num
            xr=pos_xarr(i)-pos_xarr(j)
            yr=pos_yarr(i)-pos_yarr(j)

! Applying periodic boundary condition             
            if(xr.ge.l/2) then
            xr=xr-l
            else if(xr.le.(-l/2)) then
            xr=l+xr  
            end if

            if(yr.ge.l/2) then
            yr=yr-l
            else if(yr.le.(-l/2)) then
            yr=l+yr  
            end if

            r2=(xr*xr+yr*yr)
	    r=sqrt(r2)

!*****************************************************************************************
!Removed the following which was present in last good working back up created on 22 Aug. Although the system had very low U per particle in that backup.
!*****************************************************************************************
! Making sure that the minimum distance of appcoarch is twice the Wigner-Seitz radius,
!            if((r.lt.(2*a)).and.(r.gt.0.0)) then
!            write(*,*)'r=',r
!            xr=2*a*(xr/r)
!            yr=2*a*(yr/r)
!            r2=(xr*xr+yr*yr)
!	    r=sqrt(r2)
!            write(*,*)'r changed to',r
!            end if
   


! Finding the pair wise potential and force (using Yukawa potential for now)
! Force on ith particle is opposite in sign and equal in magnitude to that on jth particle 
            u=u+(1/r2)*exp(-k*r)
!            if(i.eq.4) then
!               write(*,*)'**',xr,pos_xarr(i),pos_xarr(j) 
!               write(*,*) '****',i,u,r,(1/r2)*exp(-k*r),(1/r2),exp(-k*r)
!            endif
!            f_ij=((-exp(-k*r))/(r**2))*((2/r)+k)
            f_x=(xr*exp(-k*r))*((1+k)/(r**3))
            f_y=(yr*exp(-k*r))*((1+k)/(r**3))
            fx(i)=fx(i)+f_x
            fx(j)=fx(j)-f_x
            fy(i)=fy(i)+f_y
            fy(j)=fy(j)-f_y
            
! Finding the largest x and y distances for sanity check
!            if(xr.ge.maxxr) then
!            maxxr=xr
!            end if
!            if(yr.ge.maxyr) then
!            maxyr=yr
!            end if
!            if(r.ge.maxr) then
!            maxr=r
!            end if

        enddo
      enddo

! Finding the kinetic energy
!      do i=1,part_num
!         ke=ke+(v(i)**2)/2
!      end do

! Total energy is sum of kinetic and potential energy      
!      en=ke+u

! Printing the largest x and y distances for sanity check
!      write (*,*) ' max xr', maxxr,' max yr', maxyr,' max r2', maxr

!      write (13,'(F12.9,T15,F12.9)') u/part_num, en/part_num

! Wrinting the intial forces to file - for diagnostic reasons only
!      do i=1,part_num
!         write (11,'(F10.3)') fx(i)		
!      end do
      close (11)
!      close (13)
!      t=0
!      do i=1,steps
      call integrator(count_t,fx,fy,en,pos_xarr,pos_yarr,x_prev,y_prev,v,v_x,v_y,dt,m,part_num,u,k,t,xi,steps)
!         t=t+steps*dt
!      enddo
return 
end
