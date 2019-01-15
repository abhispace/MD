program main
implicit none
real, external :: rand_pos

      character*25 file_name
      character*25 fname,fname_pos
      CHARACTER(*), PARAMETER :: fileplace = "/abhinav/abhinav/Molecular Dynamics/output_data/"
      CHARACTER(*), PARAMETER :: fileplace_pos = "/abhinav/abhinav/Molecular Dynamics/output_data/"
      integer part_num,steps,i,part_steps,n_max,n_min,n_int,j,k_step  
      real l,kmin,kmax,k,dk,a,dt,t_norm,v_norm,r_norm,u_norm,e_norm,temp_norm
      real e0,q,temp,en,u
      real, parameter :: pi = 22.0/7.0
      real, parameter :: kb = 1.38064852e-23

!      real, parameter :: m = 9.9765e-26				! Mass of SiO2 (dust)
      real, parameter :: m = 1.66e-25					! Mass of Al2O3 (dust)

      real, parameter :: m_e = 9.11e-31				! Mass of electron
      real, parameter :: epsilon_0 = 8.854187817e-12
      real, parameter :: ev_to_joule=1.6e-19
      real, parameter :: e_charge=1.6e-19
   
      n_min=0
      n_max=1000
      n_int=100
      part_steps=(n_max-n_min)/n_int

!      l=40							! Side length of 2D square; constant for now - change later

      kmin=0.0
      kmax=0.40
      k=kmin
      dk=0.1							! Screening parameter: set to constant for now - change later
      k_step=(kmax-kmin)/dk

      dt=0.0005							! Time step 
      steps=1e+1  						! Total number of steps
      temp=2.52                                                 ! in eV (probably kb_t)
      q=(2e4)*e_charge

      temp=temp*ev_to_joule      
      


      write (*,*) 'Beginning the code!'      
      write (*,*) 'k_step:',k_step

      do j=1,k_step
         k=k+dk
         part_num=0
         write(file_name,'("steps_",I5,"_data",F3.0,"txt")')steps,k*100
         open (unit=20,file=fileplace//file_name,action="write",status="unknown")
         close (20)
      enddo

      k=kmin
      do j=1,k_step
         k=k+dk
         part_num=0
         write(file_name,'("steps_",I5,"_data",F3.0,"txt")')steps,k*100
         write(*,*)'k: ',k,'dk:',dk
         open (unit=20,file=fileplace//file_name,access="append",status="unknown")

         do i=1,part_steps

            part_num=part_num+n_int
            l=4*sqrt(real(part_num))

            a=1/sqrt(pi*part_num/l**2)					! Wigner-Seitz radius
            e0=(q*q)/(4*pi*epsilon_0*a)

! Normalization constants

           t_norm=1/(a*sqrt(m/e0))
           v_norm=1/(sqrt(e0/m))
           r_norm=1/a
           u_norm=2/e0
           e_norm=1/e0
           temp_norm=1/(e0/kb)

           write (*,*) 'Number of particles: ',part_num,' Square side length: ',l
           write (*,*) 'Time step dt: ',dt,' Number of steps: ',steps
!           write (*,*) 'k: ',k,' temp: ',temp
!           write (*,*) 'a: ',a,' e0: ',e0,' q: ',q
!           write (*,*) 'Norm(time): ',t_norm,' Norm(velocity): ',v_norm
!           write (*,*) 'Norm(r): ',r_norm,' Norm(potential): ',u_norm
!           write (*,*) 'Norm(energy): ',e_norm,' Norm(temp): ',temp_norm

           call gen_particles(part_num,l,k,a,m,temp,v_norm,dt,steps,en,u)
           write(*,*)'particle num:',part_num,'en/n:',en/part_num,'u/n:',u/part_num
           write(20,'(I5,T10,F9.4,T25,F9.4)')part_num,en/part_num,u/part_num
           write(*,*)'.........................................................'
         enddo
         write(*,*)'**********************************************************'
         write(*,*)'Writing to file: ',file_name
         close (20)
         write(*,*)'**********************************************************'
      enddo


      write (*,*) 'Done!'
end program main
