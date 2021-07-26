!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!       BINDING ENERGY CALCULATIONS FOR A GIVEN
!!!       LIGAND - MAKING ANGLES -pi/12*i and 
!!!       pi/12*i WITH THE HORIZONTAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!       Both protein and ligand have 3 binding sites

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!             ALGORITHM

!!! pc_cart (lc_cart) cartesian coodinates of protein (ligand)
!!! odd entries x, even entries y
!!! pc_pol (lc_pol) polar coordinates of protein (ligand)
!!! odd entries r, even entries theta
!!! First cartesian coordinates converted to polar
!!! with origin at tip of ligand
!!! polar angle rotated by def (ligand angle)
!!! to realign x-axis along ligand side
!!! energy is minimized in polar coordinates
!!! axis is rotated back and converted back to 
!!! cartesian coordinates


program energyprogram
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! the routines from the paper are imported
  !! I changed ddot etc to make it better inlining
  !! since  the increase is always 1.
  !! one could also use the mkl routines but I didnt do this
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use imported
  implicit none

  character(len=80) :: fmt0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Important integers 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  integer (kind=4) :: namino  !! the number of amino acids
  integer (kind=4) :: nproteins !! the number of protein sequences
  integer (kind=4) :: plinks !! links within protein
  integer (kind=4) :: nligand !! the number of ligands
  integer (kind=4) :: nbind  !! the number of ligand binding sites

  integer (kind=4) :: llinks !! links between prot an ligand
  integer (kind=4) :: naa_color !! tot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! the coupling strengths
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real (kind=8) :: ks
  real (kind=8) :: kw
  real (kind=8) :: kww

  integer (kind=4) :: kstep
  real (kind=8), dimension(:), allocatable :: kmod
  
  real (kind=8) :: els
  real (kind=8) :: elw
  real (kind=8) :: elm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! local structural variation constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real (kind=8) :: dr11, dr22, dr12, dr_table(2,2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Program mode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical :: lstruct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! More fixed parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer (kind=4), parameter :: maxneighbors=12
  real (kind=8) , parameter :: epsilon=0.000001d0
  real (kind=8) , parameter :: pi = 4d00*atan(1d00)

  real (kind=8) , parameter :: sigma = 3d-1   !! Lengthscale of chemical interaction
  real (kind=8) , parameter :: sigma2 = sigma**2
  real (kind=8) , parameter :: sigma3 = sigma2*2.0d00**(0.333333)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Coordinates arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(8), dimension(:), allocatable :: x, y, pc_cart, pc_pol, pc_orig, pc_init
  real(8), dimension(:), allocatable :: lx, ly, lc_cart, lc_pol


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! system parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(4), dimension(:), allocatable ::  plinks1, plinks2, llinks1, llinks2
  real(8), dimension(:), allocatable ::  strength, lstrength, length

  integer(4), dimension(:,:), allocatable ::  amino_list
  integer(4), dimension(:), allocatable ::  amino_strength


  integer*4 i, k, curr_prot, j1, j2, curr_lig, aa, counter  !! loop variables

  real(8), dimension(:), allocatable ::  angles
  real(8), dimension(:,:), allocatable ::  lposy, lposx  !! positions of the ligands
  integer(4), dimension(:), allocatable ::  lposxnumber  !! the x-position of the ligand in relation to the protein
!! this is the index, not the coordinate of the corresponding amino acid
  integer(4), dimension(:,:), allocatable ::  lamino_strength  !! this is 1 or 2

  
 
  integer(4) ::  i1, i2, j, jj
  integer(4), dimension(:,:), allocatable ::  partners   ! for each i the list of partners of i
  integer(4), dimension(:), allocatable ::  npartners
  integer(4), dimension(:,:), allocatable ::  partners2linknumber
  real(8), dimension(:,:), allocatable ::  energyold

  real(8) :: deltaenergy, ligands, newligands
  real(8) :: aatable(2, 2), laatable(2, 2)
  real(8) :: def_energy, chem_energy, energy, f
  real(8) :: x0,y0 ! Origin coords


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! parameters for the gradient programs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer (kind=4) , parameter :: iprint =-1  ! see src.f for explanations

  character*60     task, csave
  integer (kind=4) ::  dimn
  integer (kind=4) , parameter ::  m=5  !! increaing m maight make the minimum
  !! slightly preciswrite (filename1, '( "data/seed_list_", I2, ".dat" )' )  cpuer (not worth it)

  !! derived constants
  integer (kind=4) , parameter :: mmax=m
  integer (kind=4) :: nmax
  integer (kind=4) :: lenwa  
  integer*4 isave(44)
  real*8 dsave(29)
  logical lsave(4)
  !     We specify the tolerances in the stopping criteria.
  !     factr= 0.001d+00/2.2d-16
  real(8), parameter :: factr = 1.0d+7
  real(8), parameter ::  pgtol = 1.0d-5

  integer(4), dimension(:), allocatable :: nbd, iwa
  real(8), dimension(:), allocatable :: upper, lower, wa, derivative

! real :: start, finish

  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read parameters from "parameters.txt"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  open(unit=11, file = 'inputs/parameters.txt')

  read(11, *)  
  read(11, *) nproteins
  read(11, *)  
  read(11, *) namino
  read(11, *)  
  read(11, *) plinks
  read(11, *)  
  read(11, *) nligand
  read(11, *)  
  read(11, *) nbind

  allocate( lposxnumber(nbind) )
  read(11, *)  
  read(11, *) lposxnumber(:)

  read(11, *)  
  read(11, *) ks
  read(11, *)  
  read(11, *) kw
  read(11, *)  
  read(11, *) kww
  read(11, *)  
  read(11, *) kstep

  allocate( kmod(kstep) )

  read(11, *)  
  read(11, *) kmod
  read(11, *)  
  read(11, *) els
  read(11, *)  
  read(11, *) elw
  read(11, *)  
  read(11, *) lstruct
  read(11, *)  
  read(11, *) dr11
  read(11, *)  
  read(11, *) dr12
  read(11, *)  
  read(11, *) dr22
  close(11)

  elm = (els + elw) / 2.0d0 

  dr_table(1,1) = dr11
  dr_table(1,2) = dr12
  dr_table(2,1) = dr12
  dr_table(2,2) = dr22
  
  llinks = nbind
  naa_color = namino + nbind
  dimn = 2 * namino
  nmax = dimn
  lenwa = 2*mmax*nmax +  4*nmax+ 11*mmax*mmax + 8*mmax



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Allocate arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate( pc_cart(dimn), pc_pol(dimn), pc_orig(dimn), pc_init(dimn) )
  allocate( lc_cart(2*nbind), lc_pol(2*nbind) )

  allocate( x(namino), y(namino) )
  allocate( lx(nbind), ly(nbind) )


  allocate( plinks1(plinks), plinks2(plinks), strength(plinks), length(plinks) )
  allocate( llinks1(llinks), llinks2(llinks), lstrength(llinks) )

  allocate( amino_list(nproteins,naa_color), amino_strength(naa_color) )


  allocate( angles(nligand), lposx(nligand,nbind), lposy(nligand,nbind) )
  allocate( lamino_strength(nligand,nbind) )

  
  allocate( partners(namino,maxneighbors), npartners(namino), partners2linknumber(namino,namino) )
  allocate( energyold(namino,maxneighbors) )

  allocate( upper(dimn),lower(dimn), nbd(dimn), iwa(3*dimn), wa(lenwa), derivative(nmax) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! interaction tables

  !!  amino acid type interactions

  laatable(1, 1) = elw 
  laatable(1, 2) = elm
  laatable(2, 1) = elm
  laatable(2, 2) = els


  ! ligand - protein interaction indices

! lposxnumber(1)=2
! lposxnumber(2)=7
! lposxnumber(3)=3
  
  !! Create lists with protein-ligand link indices 
  do i = 1, nbind ! the 3 ligand binding sites
     llinks1(i) = lposxnumber(i)
     llinks2(i) = i
  enddo
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read protein positions

  !! reading the original positions of the aminoacids
  open(unit = 11, file = 'inputs/prot_xy.dat')
  do i = 1, namino
    read(11, *) pc_orig(2*i-1:2*i)
    x(i) = pc_orig(2*i-1)
    y(i) = pc_orig(2*i)
  enddo
  close(11)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read ligand types

  open(unit = 13, file = 'inputs/ligands.dat')
  do i = 1, nligand
!    read(13, *) lamino_strength(i,:), angles(i)  !! read the ligands
     read(13, *) lamino_strength(i,:), angles(i), lposx(i,:), lposy(i,:)  !! read the ligands
     do j = 1, nbind
        if(lamino_strength(i,j) .ne. 1 .and. lamino_strength(i,j) .ne. 2) then
           print *, "inconsistent ligand data", i, j
           stop
        endif
      enddo
  enddo
  close(13)


  !! Convert ligand angles to positions

  !! shape of ligand with respect to horizontal 
  !! ligand angles def - ligand1 - 0, ligand2 - pi/12, ligand3 - 2pi/12
  !! ligand4 - 3pi/12, ligand5 - 4pi/12, ligand6 - 5pi/12

! do i = 1, nligand
!    lposx(i,1) = -1d0*cos(angles(i))
!    lposx(i,2) = 0
!    lposx(i,3) = 1d0*cos(angles(i))
!    lposy(i,1) = -1d0*sin(angles(i))
!    lposy(i,2)= 0
!    lposy(i,3)=-1d0*sin(angles(i))
! enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read protein adjacency file
  !! Reading and creating link tables for fast gradient calculation
  !! for protein internal elastic bonds
  do i=1, namino
     npartners(i) = 0
  enddo

  open(unit=20, file = 'inputs/prot_adjacency.dat')
  do k = 1, plinks
     read(20,*) plinks1(k), plinks2(k)

     if (plinks1(k) > namino .or. plinks2(k) > namino &
         .or. plinks1(k) <= 0 .or. plinks2(k) <=0) then
        print*, "inconsistent links", plinks1(k), plinks2(k)
        stop
     endif

     if((x(plinks1(k)) - x(plinks2(k)))**2 &
         + (y(plinks1(k)) - y(plinks2(k)))**2 - 1 > 0.001) then
        print*, "long distance", k, plinks1(k), plinks2(k)
!       stop
     endif

     i1 = plinks1(k)
     i2 = plinks2(k)
     npartners(i1) = npartners(i1) + 1
     npartners(i2) = npartners(i2) + 1

     if (npartners(i1) > maxneighbors .or. npartners(i2) > maxneighbors) then
        print *, "too many neighbors"
        stop 
     endif

     partners(i1,npartners(i1)) = i2
     partners(i2,npartners(i2)) = i1
     partners2linknumber(i1,i2) = k
     partners2linknumber(i2,i1) = k
  enddo
  close(20)

  !! Set the equilibrium length based on the original coordinates
  do i=1, plinks
     length(i) = sqrt((x(plinks1(i)) - x(plinks2(i)))**2 + &
     & (y(plinks1(i)) - y(plinks2(i)))**2)
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read all protein sequences 

  open(unit=23, file = 'inputs/prot_seq.dat') !! these are the generations
  do i = 1, nproteins    !! the loop over all proteins
     read(23,*) amino_list(i,:)  !! read the mutation
     do j = 1, naa_color
        if(amino_list(i,j) .ne.1 .and. amino_list(i,j) .ne. 2) then
            print *, "inconsistent generations", curr_prot, j
            stop
        endif
     enddo
  enddo
  close(23)


  if (lstruct) then
    open(unit=50, file='config_local.dat') !! File for saving local configurations
  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! the output file(s)
! call cpu_time(start)

  open(unit=45, file="config.out")
  open(unit=46, file="energy.out")
  open(unit=47, file="config_init.dat") !! File for saving local configurations


  !! looping over all proteins inside
  !! loop aver all ligands
  do curr_prot = 1, nproteins   !! loop over proteins

     do aa = 1, naa_color
        amino_strength(aa) = amino_list(curr_prot, aa)
     enddo

     if (lstruct) then

         !! Update structure according to dr_table
         call update_structure(dimn, namino, plinks, pc_orig, pc_init, amino_list(curr_prot,:), plinks1, plinks2, dr_table)

         do i = 1, namino
            x(i) = pc_init(i*2-1)
            y(i) = pc_init(i*2)
         end do

         !! Set the equilibrium length based on the updated coordinates
         do i=1, plinks
            length(i) = sqrt((x(plinks1(i)) - x(plinks2(i)))**2 + &
            & (y(plinks1(i)) - y(plinks2(i)))**2)
         enddo

         write(50,'(*(f6.3,x))')  pc_init

     else

        pc_init(:) = pc_orig(:)

     end if

    

     do curr_lig = 1, nligand  !! loop over ligands

        do i = 1, nbind
           ! Update bond strength
           j1 = lamino_strength(curr_lig, llinks2(i))
           j2 = amino_strength(namino + i)
           lstrength(i) = laatable(j1, j2)

           ! Update ligand position
           lc_cart(2*i-1) = lposx(curr_lig,i)!+ x(lposxnumber(2))
           lc_cart(2*i) = lposy(curr_lig,i)!+ y(lposxnumber(1))
        enddo
        

        !!! Only run on a single ligand for testing!
!       if ((curr_lig /= 6) .or. (curr_lig /= 5)) cycle

        !! tip of the ligand - chosen as origin
        x0 = lc_cart(3)
        y0 = lc_cart(4)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!           ENERGY MINIMIZATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !! initialize protein coordinates
        pc_cart(:) = pc_init(:)

        !! establish bounds on protein coordinates
        call get_bounds(dimn, nbind, angles(curr_lig), lposxnumber, nbd, upper, lower)


        !! Translate protein so that it almost overlaps with the ligand
        call translation_step(dimn, pc_cart, lposxnumber, x0, y0, angles(curr_lig), lower, upper)
        write(47,'(*(f6.3,x))')  pc_cart

        !! Loop over spring constant steps
        do counter = 1, kstep

          !!  amino acid type interactions

          aatable(1, 1) = kmod(counter) * kww 
          aatable(1, 2) = kmod(counter) * kw
          aatable(2, 1) = kmod(counter) * kw
          aatable(2, 2) = kmod(counter) * ks
  
          ! compute linkstrength
          do i = 1, plinks
             j1 = amino_strength(plinks1(i))
             j2 = amino_strength(plinks2(i))
             strength(i) = aatable(j1, j2)   
          enddo


          !! convert protein to polar coordinates with origin at tip of the ligand
          !! pc_pol are polar coordinates - odd entries r, even are theta
          call aaxy2polar(pc_cart, x0, y0, angles(curr_lig), pc_pol)

          !! convert ligand to polar coordinates
          !! lc_pol are polar coordinates - odd entries r, even are theta
          call lxy2polar(lc_cart, x0, y0, angles(curr_lig), lc_pol)



          !!    comments from the gradient people
          !     We wish to have output at every iteration.

          !     We specify the tolerances in the stopping criteria.


!         !     We specify the dimension n of the sample problem and the number
          !     m of limited memory corrections stored.  (n and m should not
          !     exceed the limits nmax and mmax respectively.)


          !     We now provide nbd which defines the bounds on the variables:
          !     l   specifies the lower bounds,
          !     u   specifies the upper bounds. 

          !     First set bounds on the odd-numbered variables.
          !! 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!     minimization angles constrained between 0 and pi + 2 def
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! only restriction on angles
          !! "nbd = 1" means only a lower bound is used,
          !! "nbd = 2" means both lower and upper bounds are used

!         do  i = 1, dimn
!            if (mod(i,2) .eq. 0) then
!               nbd(i) = 2
!               if (int(i/2) == lposxnumber(1)) then
!                 lower(i)   = 0.5 * pi + angles(curr_lig)
!                 upper(i)   = pi + 2d0 * angles(curr_lig)
!               elseif (int(i/2) == lposxnumber(3)) then
!                 lower(i)   = 0
!                 upper(i)   = 0.5 * pi + angles(curr_lig)
!               else
!                 lower(i)   = 0
!                 upper(i)   = pi + 2d0 * angles(curr_lig)
!               endif
!            else
!               nbd(i) = 1
!               lower(i)   = 0d0
!               upper(i)   = 10d0
!            endif
!         enddo


          !! !     We now define the starting point.
          !!  
          !! !     We now write the heading of the output.

          task='START'
  111     continue
          call setulb(dimn, m, pc_pol, lower, upper, nbd, f, derivative, factr, pgtol, &
                      wa, iwa, task, iprint, csave, lsave, isave, dsave)


          if (task(1:2) .eq. 'FG') then

             call totalenergy(pc_pol, f)
             call  ligandenergy(pc_pol, ligands)

             do i = 1, namino
                do k = 1, npartners(i)
                   j = partners(i,k)
                   !! |r1-r2| = ri^2 + r2^2 - 2r1*r2*cos(t1-t2)
                   energyold(i,k) = (sqrt(pc_pol(2*i-1)**2 + pc_pol(2*j-1)**2 - &
                                  & pc_pol(2*i-1)*pc_pol(2*j-1)*2.0*cos(pc_pol(2*i) - pc_pol(2*j))) &
                                  & - length(partners2linknumber(i,j)))**2 * strength(partners2linknumber(i,j))

                enddo
             enddo


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !! this is where the huge time gain comes from
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


             !! we only compute the energy change for the affected links
             !! in "direction" j
             do j=1 , dimn
                i = int((j + 1) / 2)
                pc_pol(j) = pc_pol(j) + epsilon
                deltaenergy = 0
                do k = 1, npartners(i)
                   jj = partners(i,k)
                   !! |r1-r2| = ri^2 + r2^2 - 2r1*r2*cos(t1-t2)
                   deltaenergy = deltaenergy &
                        +(sqrt(pc_pol(2*i-1)**2 + pc_pol(2*jj-1)**2 - pc_pol(2*i-1)*pc_pol(2*jj-1)*&
                        &2.0*cos(pc_pol(2*i) - pc_pol(2*jj))) &
                        &-length(partners2linknumber(i,jj)))**2 * &
                        strength(partners2linknumber(i,jj)) - energyold(i,k)
                enddo
                call ligandenergy(pc_pol, newligands)
                derivative(j) = (deltaenergy + newligands - ligands)/epsilon
                pc_pol(j) = pc_pol(j) - epsilon
             enddo


             goto 111
          elseif (task(1:5) .eq. 'NEW_X') then
             goto 111
             stop
          else

             !        We terminate execution when task is neither FG nor NEW_X.
             !        We print the information contained in the string task
             !        if the default output is not used and the execution is
             !        not stopped intentionally by the user. 

             if (iprint .le. -1 .and. task(1:4) .ne. 'STOP') write(6,*) task

          endif



          call totalenergy(pc_pol, energy)
          call proteinenergy(pc_pol, def_energy)
          call ligandenergy(pc_pol, chem_energy)


          !! convert polar back to cartesian coordinate
          call aapolar2xy(pc_pol, x0, y0, angles(curr_lig), pc_cart)

          !! Calculate energy from initial guess
          call aaxy2polar(pc_init, x0, y0, angles(curr_lig), pc_pol)
          call totalenergy(pc_pol, f)

          !! **redundant** return to the final pc_pol, just in case I forget about this
          call aaxy2polar(pc_cart, x0, y0, angles(curr_lig), pc_pol)


!          write(fmt, 
           write(45,'(*(f6.3,x))')  pc_cart
           write(46,'(4(e18.10,x))') energy, def_energy, chem_energy, f


          !! Finished one protein-ligand combination

        end do
!       write(45,*)

!       call cpu_time(finish)
!       print '("Time = ",f20.1," seconds.")',finish-start

     enddo ! curr_lig=1,nligand
  enddo ! curr_prot=1,nproteins

  close(45)
  close(46)
  close(47)
  close(50)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

  !! plinks1(k) = number of amino 1 in link
  !! plinks2(k) = number of amino 2 in link
  !! strength(k) = strength of that link
  !! length(k) = natural length of link k (for the moment, I take it equal to 1)
  !!   so its not used

  !! llinks1(l) = number of ligandamino 1 in link
  !! llinks2(l) = number of amino  in link
  !! lstrength(l) = strength of that link
  !! llength(l) = natural length of link l
  !! this IS used


  !DIR$ ATTRIBUTES FORCEINLINE :: energy

  !!  Total energy = spring + chemical lj energies
  subroutine totalenergy (xp,result)

    implicit none
    real*8 xp(2*namino)
    real*8 result, dr, rs
    integer*4 k,l
    result=0
    do k= 1,plinks
       result=result+ &
            (sqrt(xp(2*plinks1(k)-1)**2+xp(2*plinks2(k)-1)**2 &
            -xp(2*plinks1(k)-1)*xp(2*plinks2(k)-1)*2.0* &
            &cos(xp(2*plinks1(k))-xp(2*plinks2(k)))) &
            &-length(k))**2*strength(k)
    enddo
    !        print *,"---",result      

    do l=1,llinks
       dr = xp(2*llinks1(l)-1)**2+lc_pol(2*llinks2(l)-1)**2 &
            -xp(2*llinks1(l)-1)*lc_pol(2*llinks2(l)-1)*&
            &2.0*cos(xp(2*llinks1(l))-lc_pol(2*llinks2(l)))

       rs = dr/sigma2
       result=result - lstrength(l)*exp(-rs)

    enddo
    !        print *,"+++",result

  end subroutine totalenergy


!!!       only spring energies of the protein 

  subroutine proteinenergy (xp,result)

    implicit none
    real*8 xp(2*namino)
    real*8 result
    integer*4 k
    result=0
    do k= 1,plinks
       result=result+ &
            (sqrt(xp(2*plinks1(k)-1)**2+xp(2*plinks2(k)-1)**2 &
            -xp(2*plinks1(k)-1)*xp(2*plinks2(k)-1)*2.0* &
            &cos(xp(2*plinks1(k))-xp(2*plinks2(k)))) &
            &-length(k))**2*strength(k)
    enddo

  end subroutine proteinenergy

!!!   only lj chemical energies of the protein-ligand interaction

  subroutine ligandenergy (xp,result)

    implicit none
    real*8 xp(2*namino), dr, rs
    real*8 result
    integer*4 l
    result=0

    do l=1,llinks
       dr = xp(2*llinks1(l)-1)**2+lc_pol(2*llinks2(l)-1)**2 &
            -xp(2*llinks1(l)-1)*lc_pol(2*llinks2(l)-1)*&
            &2.0*cos(xp(2*llinks1(l))-lc_pol(2*llinks2(l)))

       rs = dr/sigma2
       result=result - lstrength(l)*exp(-rs)

    enddo
  end subroutine ligandenergy



  !! convert cartesian to polar coordinates for proteins
  !! rotate by angle theta to align x-axis along the ligand
  subroutine aaxy2polar(pc_cart, x0, y0, theta, pc_pol)

    implicit none
    real*8 pc_cart(2*namino), x0, y0, pc_pol(2*namino), theta
    integer*4 aa

    do aa = 1,namino
        pc_pol(2*aa-1) = sqrt((pc_cart(2*aa-1)-x0)**2 + (pc_cart(2*aa)-y0)**2)
        pc_pol(2*aa) = atan2(pc_cart(2*aa) - y0, pc_cart(2*aa-1) - x0) + theta

        ! Make sure that angle is non-negative, since this is required
        ! by the constraint keeping the protein above the ligand
        if (pc_pol(2*aa) < 0.0d0) then
            pc_pol(2*aa) = pc_pol(2*aa) + 2.0d0 * pi
        end if
    enddo
  end subroutine aaxy2polar

  !! convert cartesian to polar coordinates for ligands
  subroutine lxy2polar(lc_cart, x0, y0, theta, lc_pol)

    implicit none
    real*8 lc_cart(2*nbind), x0, y0, lc_pol(2*nbind), theta
    integer*4 ligand

    do ligand = 1,nbind
       lc_pol(2*ligand-1) = sqrt((lc_cart(2*ligand-1)-x0)**2+(lc_cart(2*ligand)-y0)**2)
       lc_pol(2*ligand) = atan2(lc_cart(2*ligand)-y0,lc_cart(2*ligand-1)-x0)+theta
    enddo
  end subroutine lxy2polar

  !! convert polar to cartesian coordinates
  subroutine aapolar2xy(pc_pol, x0, y0, theta, pc_cart)

    implicit none
    real*8 pc_pol(2*namino), x0, y0, pc_cart(2*namino), theta
    integer*4 aa

    
    do aa = 1, namino
       pc_cart(2*aa-1) = pc_pol(2*aa-1)*cos(pc_pol(2*aa)-theta) + x0
       pc_cart(2*aa) = pc_pol(2*aa-1)*sin(pc_pol(2*aa)-theta) + y0
    enddo

  end subroutine aapolar2xy

  !! calculate boundaries
  subroutine get_bounds(dimn, nbind, angle, lposxnumber, nbd, upper, lower)

    implicit none
    integer (kind=4), intent(in)::  dimn, nbind
    real(8), intent(in) :: angle
    integer(4), dimension(nbind), intent(in) ::  lposxnumber(:) 

    real(8), dimension(dimn), intent(inout) :: upper, lower
    integer(4), dimension(dimn), intent(inout) :: nbd

          do  i = 1, dimn
             if (mod(i,2) .eq. 0) then
                nbd(i) = 2
                if (int(i/2) == lposxnumber(1)) then
                  lower(i)   = 0.5 * pi + angle
                  upper(i)   = pi + 2d0 * angle
                elseif (int(i/2) == lposxnumber(3)) then
                  lower(i)   = 0
                  upper(i)   = 0.5 * pi + angle
                else
                  lower(i)   = 0
                  upper(i)   = pi + 2d0 * angle
                endif
             else
                nbd(i) = 1
                lower(i)   = 0d0
                upper(i)   = 10d0
             endif
          enddo

  end subroutine get_bounds

  ! Optimize initial configuration by translating the entire protein
  ! so that it almost overlaps with the ligand.
  ! This is achieved using the same set of constraints that are needed
  ! for the deformation optimization step, plus an additional constraint
  ! which is needed in case theta = pi/2, to stop the ligand from going down indefinitely
  subroutine translation_step(dimn, pc_cart, lposxnumber, x0, y0, theta, lower, upper)

    implicit none

    integer(4), intent(in) :: dimn 
    real(8), intent(in):: x0, y0, theta
    real(8), dimension(dimn), intent(in) :: upper, lower
    real(8), dimension(dimn), intent(inout) :: pc_cart
    integer(4), dimension(nbind), intent(in) ::  lposxnumber
    
    real(8) :: delta, pc_save(dimn), pc_pol(dimn)
    integer(4) ::  i, j, c, reject, cnt
    logical ::  accept
    real(8), parameter :: delta0 = 0.2
    integer(4), parameter :: steps = 5
    integer(4), parameter :: maxcnt = 50

    delta = delta0
    c = 1
    reject = 0

    do i = 1, steps

        ! Need to limit the number of times this can repeat,
        ! as it is possible to reach a steady-state which endlessly
        ! oscillates between left-to-right states

        cnt = 0

        do while ( (c < 4) .and. (cnt < maxcnt))

            cnt = cnt + 1

            pc_save(:) = pc_cart(:)

            select case(c)
            
              ! Move downwards (y)
              case(1)
                do j = 1, int(dimn / 2)
                    pc_cart(j*2) = pc_cart(j*2) - delta
                end do

              ! Move left (x)
              case(2)
                do j = 1, int(dimn / 2)
                    pc_cart(j*2-1) = pc_cart(j*2-1) - delta
                end do
                ! After moving left, try downwards again
                c = 1
            
              ! Move right (x)
              case(3)
                do j = 1, int(dimn / 2)
                    pc_cart(j*2-1) = pc_cart(j*2-1) + delta
                end do
                ! After moving right, try downwards again
                c = 1
            
            end select

            ! Check constraints
            accept = .true.
            ! Angle constraints
            call aaxy2polar(pc_cart, x0, y0, theta, pc_pol)
            do j = 1, int(dimn / 2)
                if (pc_pol(j*2) < lower(j*2)) then
                    accept = .false.
                elseif (pc_pol(j*2) > upper(j*2)) then
                    accept = .false.
                end if
            end do
            ! y-axis constraint
            ! Central protein binding site cannot go below central ligand site
            if (pc_cart(lposxnumber(2)*2) < y0) accept = .false.

            if (accept) then
                reject = 0
            else
                pc_cart(:) = pc_save(:)
                reject = reject + 1
                ! If "downwards" failed, then try "left"
                if (reject == 1) then
                    c = 2
                ! If "left" also failed, then try "right"
                elseif (reject == 2) then
                    c = 3
                ! If "right" also failed, then end loop
                elseif (reject == 3) then
                    c = 4
                end if
            end if

        end do

        delta = delta * 0.5d0

    end do


  end subroutine translation_step


  subroutine update_structure(dimn, namino, plinks, pc_orig, pc_cart, seq, plinks1, plinks2, table)

    implicit none

    integer(4), intent(in) :: dimn, namino, plinks
    integer(4), dimension(namino), intent(in) :: seq
    integer(4), dimension(plinks), intent(in) :: plinks1, plinks2
    real(8), dimension(dimn), intent(in) :: pc_orig
    real(8), dimension(dimn), intent(inout) :: pc_cart
    real(8), dimension(2,2), intent(in) :: table

    integer(4) ::  i, j, k
    real(8) :: vec(dimn), vec0(2)

    vec(:) = 0.0d0
    vec0(:) = 0.0d0

    do i = 1, plinks
        j = plinks1(i)
        k = plinks2(i)
        vec0 = (pc_orig(k*2-1:k*2) - pc_orig(j*2-1:j*2)) * table(seq(j), seq(k))
        vec(j*2-1:j*2) = vec(j*2-1:j*2) + vec0
        vec(k*2-1:k*2) = vec(k*2-1:k*2) - vec0
    end do

    pc_cart(:) = pc_orig(:) + vec(:)


  end subroutine update_structure

end program energyprogram

