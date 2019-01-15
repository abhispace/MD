! This subroutine calculates the new positions and velocities of each particle
! Author=Abhinav

subroutine integrator(count_t,fx,fy,en,pos_xarr,pos_yarr,x_prev,y_prev,v,v_x,v_y,dt,m,part_num,u,k,t,xi,steps)

! Declaring and Setting constants

      integer part_num,i,j,steps,xi
      character*25 fname,fname_pos
      character*40 en_fname
      CHARACTER(*), PARAMETER :: fileplace = "/abhinav/abhinav/Molecular Dynamics/output_data/pos/"
      CHARACTER(*), PARAMETER :: en_file = "/abhinav/abhinav/Molecular Dynamics/output_data/energy/"
      real fx(part_num),fy(part_num),en,pos_xarr(part_num),pos_yarr(part_num),k,compare
      real x_prev(part_num),y_prev(part_num),v(part_num),v_x(part_num)
      real temp_x,temp_y,v_y(part_num),dt,m,u,ke,t,x_arr(part_num,steps)

!      if(part_num.le.99) then
!         write(fname,'("k_",F3.0,"_n_",I2,"_t_",F6.1,".txt")')k*100,part_num,t*1000
!      endif

!      if((part_num.lt.1000).and.(part_num.ge.100)) then
!         write(fname,'("k_",F3.0,"_n_",I3,"_t_",F6.1,".txt")')k*100,part_num,t*1000
!      endif

!      if(part_num.ge.1000) then
!         write(fname,'("k_",F3.0,"_n_",I4,"_t_",F6.1,".txt")')k*100,part_num,t*1000
!      endif

      if(part_num.le.99) then
         write(en_fname,'("energy_k_",F3.2,"_n_",I2,"_steps_",I5,"_.txt")')k,part_num,steps
      endif

      if((part_num.lt.1000).and.(part_num.ge.100)) then
         write(en_fname,'("energy_k_",F3.2,"_n_",I3,"_steps_",I5,"_.txt")')k,part_num,steps
      endif

      if(part_num.ge.1000) then
         write(en_fname,'("energy_k_",F3.2,"_n_",I4,"_steps_",I5,"_.txt")')k,part_num,steps
      endif

!      open (unit=12,file=fileplace//fname,action="write",status="replace")
      open (unit=13,file=en_file//en_fname,access="append",status="unknown")
      open (unit=14,file="pos_vel_recordx.txt",action="write",status="replace")
      open (unit=15,file="pos_vel_recordy.txt",action="write",status="replace")

      ke=0

! Finding new positions using Verlet algorithm
      do i=1,part_num
         temp_x=2*pos_xarr(i)-x_prev(i)+(dt**2)*fx(i)
!         if(i.eq.1) then
!            write(*,*)'++++',i,temp_x,'++++++++++++++++++++++++++++++'
!         endif
         temp_y=2*pos_yarr(i)-y_prev(i)+(dt**2)*fy(i)

         v_x(i)=(temp_x-x_prev(i))/(2*dt)
         v_y(i)=(temp_y-y_prev(i))/(2*dt)

         x_prev(i)=pos_xarr(i)
         y_prev(i)=pos_yarr(i)
         

! Applying periodic boundary condition             
         if(temp_x.ge.l) then
         temp_x=temp_x-l
         else if(temp_x.le.(0)) then
         temp_x=l+temp_x  
         end if

         if(temp_y.ge.l) then
         temp_y=temp_y-l
         else if(temp_y.le.(0)) then
         temp_y=l+temp_y  
         end if

!         if((i.eq.4).or.(i.eq.5)) then
!            write(*,*)'----',i,x_prev(i),pos_xarr(i),y_prev(i),pos_yarr(i)
!         endif

         if(temp_x.ge.(2*l)) then
!         int_l=int(2*l)
!         temp_x=int(temp_x)

! Add correction later:::::::::::::::::
!         temp_x=mod(temp_x,int_l)
         else if(temp_x.le.(-2*l)) then
!         int_l=int(2*l)
!         temp_x=int(temp_x)

!         temp_x=mod(temp_x,int_l)
         end if

         if(temp_y.ge.l) then
         temp_y=temp_y-l
         else if(temp_y.le.(0)) then
         temp_y=l+temp_y  
         end if
         pos_xarr(i)=temp_x
         pos_yarr(i)=temp_y



         ke=ke+(v(i)**2)/2
    
!         write (*,*) f(i)*(dt**2)
      enddo
      
!      write(*,*)'Number of particles: ',part_num
!      write(*,*)'m: ',m,' dt: ',dt


! Periodic boundary conditions
!      do i=1,part_num-1
!         do j=i+1,part_num

!        enddo
!      enddo

! Finding the kinetic energy
!      do i=1,part_num
!         ke=ke+(v(i)**2)/2
!      end do

! Total energy is sum of kinetic and potential energy      
      en=ke+u

! Printing the largest x and y distances for sanity check
!      write (*,*) ' max xr', maxxr,' max yr', maxyr,' max r2', maxr

! Writing the next positions and velocities for diagnostics

      do i=1,part_num
!         write (12,'(F9.3,T10,F9.3,T21,F9.6,T32,F9.6)') pos_xarr(i),pos_yarr(i),v_x(i),v_y(i)
!         x_arr(i,xi)=pos_xarr(i)
         if(i.le.99) then
            write(fname_pos,'(I2,"_pos_i_",I2,"_k_",F3.0,".txt")')part_num,i,k*100
         endif

         if((i.lt.1000).and.(part_num.ge.100)) then
            write(fname_pos,'(I3,"_pos_i_",I3,"_k_",F3.0,".txt")')part_num,i,k*100
         endif

         if(i.ge.1000) then
            write(fname_pos,'(I4,"_pos_i_",I4,"_k_",F3.0,".txt")')part_num,i,k*100
         endif
   
        if(count_t.le.steps) then
           write(*,*)' current step::::::: ',count_t
           open (unit=16,file=fileplace//fname_pos,access="append",status="unknown")
           write (16,'(2(F9.3,3X),4(F12.6,3X))')pos_xarr(i),pos_yarr(i),v_x(i),v_y(i),fx(i),fy(i)
           close(16)
        endif

      if(v_x(i).ge.100) then
         write(*,*)v_x(i),v_y(i),i,k,part_num
      endif

!      compare=100.0
!      if((pos_xarr(i).ge.compare).or.(pos_xarr(i).le.(-1*compare))) then
!         write(*,*)'i:',i,'steps:',xi,' and x=',pos_xarr(i), 'time=',t         
!      endif

      end do
!      close (12)
      

      write (13,'(F12.9,T15,F12.9)') u/part_num, en/part_num
!      write (14,'(20(F9.3))') (pos_xarr(i), i=1,part_num)
!      write (15,'(20(F9.3))') (pos_yarr(i), i=1,part_num)

      close (13)
      close (14)
      close (15)
return 
end
