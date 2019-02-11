! implementation of the Monte Carlo algorithm for partciles coagulation with the representative particles approach
! based on Zsom & Dullemond (2008) by Joanna Drążkowska, email: asiadrazkowska@gmail.com

! a representative particle nri represents a swarm swarm(nri)%npar of identical particles of mass swarm(nri)%mass each
! mass of the swarms mswrm is constant over time, but as the particles grow their number in a swarm changes
! they all are located in a cell of volume vol and we assume their homogeneous distribution within the cell
! the collision rate is proportional to the sum of the particles masses (the linear coagulation kernel)

! compile with gfortran -fdefault-real-8

program MonteCarlo

   implicit none

   type :: rp
      real           :: npar     ! number of physical particles represented by the swarm
      real           :: mass     ! mass of physical particles represented by the swarm
      ! you can easily add another properties here
   end type rp

   !  ------- set parameters here --------------------------------------------------------------------------------------
   integer, parameter             :: ntot = 200      ! total number of representative particles (swarms)
   real, parameter                :: dens = 1.       ! dust volume mass density
   real, parameter                :: mswrm = 10.e20  ! mass of the swarm
   real, parameter                :: m0 = 1.         ! mass of the monomer (initial mass of the particles)
   real, parameter                :: tend = 20.      ! maximum time of the simulation
   real, parameter                :: dtout = 4.      ! time step for the output
   integer, parameter             :: nbins = 100     ! number of bins for output histograms
   ! -------------------------------------------------------------------------------------------------------------------

   type(rp), dimension(ntot)      :: swarm             ! list of the swarms
   integer                        :: nri               ! index of representative partcile to undergo the next collision
   integer                        :: nrk               ! index of physical partcile to undergo the next collision
   real                           :: vol               ! volume of the cell
   real                           :: time = 0.0        ! time
   real                           :: dtime             ! time step between the subsequent collisions
   real                           :: tout = dtout      ! time of the next output
   real, dimension(ntot,ntot)     :: colrates          ! matrix of the collision rates between the particles
   real, dimension(ntot)          :: colrates_rp       ! probability of choosing given representative particles
   real                           :: totrate           ! total collision rate
   real                           :: rand              ! random number
   real                           :: fin               ! temporary value
   character(len=11)              :: fname             ! output file name
   real, dimension(ntot)          :: colrates_old_col  ! for optimization reasons
   real, dimension(nbins)         :: m2fm              ! mass distribution function
   real, dimension(nbins+1)       :: mgrid             ! mass grid
   integer                        :: i, j
   integer                        :: k = 1
   real                           :: ll, lll, nord

   ! calculate the cell volume
   vol = ntot * mswrm / dens

   ! initialize mass bins for histograms
   mgrid(1) = m0
   nord = log10(mswrm * ntot / m0)  ! how many orders of magnitude in mass should the histogram go through?
   ll = real(nord) / real(nbins)
   lll = ll
   do i = 2, nbins+1
      mgrid(i) = m0 * 10.0 ** lll
      lll = ll * i
   enddo

   ! initialize swarms
   swarm(:)%npar = mswrm / m0
   swarm(:)%mass = m0

   ! calculate the collision rates: this is the linear kernel
   ! and not a physical collision kernel which would be proportional to impact velocities
   do i = 1, ntot
      colrates(i,:) = swarm(:)%npar * 0.5 * (swarm(i)%mass + swarm(:)%mass) / vol
   enddo
   colrates_rp = sum(colrates, dim=2)

   !--------------- MAIN LOOP ------------------------------------------------------------------------------------------
   do while (time < tend)

      ! calculate a total rate for collisions
      totrate = sum(colrates_rp)

      ! determine the time step for the next collision
      call random_number(rand)
      dtime = - 1. / totrate * log(rand)  ! based on the Poisson distribution

      ! update time
      time = time + dtime

      ! select representative particle nri to undergo the collision
      call random_number(rand)
      rand = rand * totrate
      j = 1

      fin = colrates_rp(1)
      do while (rand > fin)
         fin = fin + colrates_rp(j+1)
         j = j + 1
      enddo
      nri = j

      ! select physical particle nrk to undergo the collision
      call random_number(rand)
      rand = rand * colrates_rp(nri)
      j = 1

      fin = colrates(nri,1)
      do while (rand > fin)
         fin = fin + colrates(nri,j+1)
         j = j + 1
      enddo
      nrk = j

      ! perform the collision and update the representative particle nri (here: only sticking)
      swarm(nri)%mass = swarm(nri)%mass + swarm(nrk)%mass
      swarm(nri)%npar = mswrm / swarm(nri)%mass

      ! remember the old column (optimization)
      colrates_old_col(:) = colrates(:,nri)

      ! update the collision rates matrix
      colrates(nri,:) = swarm(:)%npar * 0.5 * (swarm(nri)%mass + swarm(:)%mass) / vol
      colrates(:,nri) = swarm(nri)%npar * 0.5 * (swarm(nri)%mass + swarm(:)%mass) / vol
      colrates_rp(:) = colrates_rp(:) + (colrates(:,nri) - colrates_old_col(:))           ! one could just write colrates_rp = sum(colrates, dim=2) here
      colrates_rp(nri) = sum(colrates(nri,:))                                             ! but this way is much faster

      ! checking for the output
      if (time > tout) then
         write(*,*) 'time: ',time,'producing output ',k

         ! make histogram
         m2fm(:) = 0.0
         do i = 1, ntot
            j = 1
            do while (swarm(i)%mass >= mgrid(j))
               j = j + 1
            enddo
            m2fm(j) = m2fm(j) + ((swarm(i)%npar) * swarm(i)%mass**2) / ((mgrid(j+1) - mgrid(j)) * mswrm * ntot)
         enddo

         write(fname,'(a4,i3.3,a4)') 'out-',k,'.dat'
         open(k,file=fname)
            do j = 1, nbins
               write(k, *) sqrt(mgrid(j+1)*mgrid(j)), m2fm(j)
            enddo
         close(k)
         k = k + 1
         tout = tout + dtout
      endif

   enddo
   !------------------- END OF THE MAIN LOOP ---------------------------------------------------------------------------

   write(*,*) 'tend exceeded, finishing simulation'

end
