!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!       BINDING ENERGY CALCULATIONS FOR A GIVEN
!!!       LIGAND - MAKING ANGLES -pi/12*i and 
!!!       pi/12*i WITH THE HORIZONTAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!       Both protein and ligand have 3 binding sites

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!             ALGORITHM

!!! xp0 (lxf) cartesian coodinates of protein (ligand)
!!! odd entries x, even entries y
!!! xp (lxp) polar coordinates
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
  character (len=90) :: dirname
  character (len=90) :: cmd
  character (len=90) :: outname
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! the parameters
  !! links and llinks are parameters here, because
  !!I assume for now that we do not change the models too often
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  integer (kind=4), parameter :: naa=13  !! the number of amino acids
  integer (kind=4), parameter :: nlinks_lig_aa=3  !! the number of ligand aa
  integer (kind=4)  , parameter :: links=25 !! links within protein
  integer (kind=4)  , parameter :: llinks = nlinks_lig_aa !! links between prot an ligand
  integer (kind=4)  , parameter :: naa_color = naa+nlinks_lig_aa
  integer (kind=4) , parameter :: maxneighbors=12
  real (kind=8) , parameter :: epsilon=0.000001d0
  real (kind=8) , parameter :: pi = 4d00*atan(1d00)
  !! the dx for the derivative calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!nproteins is fixed for these tests,
  !! but will be dynamically changing when we do mutations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer (kind=4)  , parameter :: nproteins=2**naa_color !! the number of generations of mutations (for this test)
  integer (kind=4)  , parameter :: nligands = 6  !! because of symmetry
  integer (kind=4) , parameter :: nangles = 180 !! the niumber of ligand angles we consider
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! the coupling strengths
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real (kind=8) :: ks != 9.9d0
  real (kind=8) :: kw != 6.9d0
  real (kind=8) :: kww != 3.9d0
  
  real (kind=8) , parameter :: l0 = 1.0d0

  real (kind=8) :: els != 0.20d0

  real (kind=8) :: elw != 0.06d0
  real (kind=8) :: elm != (els+elw)/2

  real (kind=8) , parameter :: sigma = 3d-1
  real (kind=8) , parameter :: sigma2 = sigma**2
  real (kind=8) , parameter :: sigma3 = sigma2*2.0d00**(0.333333)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! note : nlinks_lig_aa is the number ofligands
  !!        lbonds is the number of bonds from the ligands
  !!        currently :: lbonds=nlinks_lig_aa  (1 per ligand)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer (kind=4) , parameter :: lbonds=nlinks_lig_aa
  real(kind=8), parameter :: downshift=0!sigma
  !! how much the ligands are pushed down
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! system parameters
  integer (kind=4) , parameter :: iprint =-1  ! see src.f for explanations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real*8 x(naa),y(naa),xp(2*naa),xp0(2*naa), xf(2*naa), xorig(2*naa)
  real*8 lxp(2*nlinks_lig_aa),lx(nlinks_lig_aa),ly(nlinks_lig_aa),lxp0(2*nlinks_lig_aa)
  integer*4 links1(links),links2(links)
  integer*4 llinks1(llinks), llinks2(llinks)
  real*8 strength(links),lstrength(llinks)
  real*8 length(links)
  integer*4 amino_list(nproteins,naa_color)
  integer*4 amino_strength(naa_color)
  integer*4 temp(naa_color), ltemp(nlinks_lig_aa)
  integer*4 i,k,currentprotein,j1,j2, currentligand,aa,cnt  !! loop variables
  real (kind=8) :: angles(nangles)
  real*8 e, energies(2*nligands), etemp, xtemp(2*naa)
  real*8 lposy(nangles,nlinks_lig_aa),lposx(nangles,nlinks_lig_aa)  !! the vertical positions of the ligands
  integer*4 lposxnumber(nlinks_lig_aa), colors(nlinks_lig_aa) !! the x-position of the ligand
  !! this is the index, not the coordinate of the corresponding amino
  integer*4 lamino_strength(nligands,nlinks_lig_aa)  !! this is 1 or 2
  integer*4 partners(naa,maxneighbors)  ! for each i the list of partners of i
  integer*4 npartners(naa),i1,i2,j,jj
  integer*4 partners2linknumber(naa,naa)
  real*8 energyold(naa,maxneighbors),deltaenergy,ligands,newligands
  real*8 aatable(2, 2), laatable(2, 2), energy
  real*8 def_energy, chem_energy
  real*8 x0,y0
  integer*4 currentangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! parameters for the gradient programs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  character*60     task, csave
  real*8 upper(2*naa),lower(2*naa)
  integer (kind=4) , parameter ::  dimn=2*naa
  integer (kind=4) , parameter ::  m=5  !! increaing m maight make the minimum
  !! slightly preciswrite (filename1, '( "data/seed_list_", I2, ".dat" )' )  cpuer (not worth it)

  !! derived constants
  integer (kind=4) , parameter :: mmax=m
  integer (kind=4) , parameter :: nmax=dimn
  integer (kind=4) , parameter :: lenwa = 2*mmax*nmax +  4*nmax+ 11*mmax*mmax + 8*mmax  
  integer*4 isave(44), nbd(dimn), iwa(3*dimn)
  real*8 dsave(29),factr,pgtol,wa(lenwa),derivative(nmax),f,fold
  logical lsave(4)
  real :: start, finish
  integer*4 counter
  !!

  !     We specify the tolerances in the stopping criteria.
  factr  = 1.0d+7
  !     factr= 0.001d+00/2.2d-16
  pgtol  = 1.0d-5



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read the interaction energies

  open(unit = 31, file = 'inputwedge/energy_const.dat')
  read(31, *)  
  read(31, *) ks
  read(31, *)  
  read(31, *) kw
  read(31, *)  
  read(31, *) kww
  read(31, *)  
  read(31, *) els
  read(31, *)  
  read(31, *) elw
  close(31)

  elm = (els + elw) / 2.0d0 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! geoetry and interaction tables

  !!  amino acid type interactions

  aatable(1, 1) = kww 
  aatable(1, 2) = kw
  aatable(2, 1) = kw
  aatable(2, 2) = ks
  
  laatable(1, 1) = elw 
  laatable(1, 2) = elm
  laatable(2, 1) = elm
  laatable(2, 2) = els


  !! ligand bindng site positions lposy
  !!   and shapes angles

  lposxnumber(1)=2
  lposxnumber(2)=7
  lposxnumber(3)=3
  
  !!! colors at the 3 binding sites
    
  colors(1)=naa+1
  colors(2)=naa+2
  colors(3)=naa+3
  

  !! shape of ligand with respect to horizontal 
  !! ligand angles def - ligand1 - 0, ligand2 - pi/12, ligand3 - 2pi/12
  !! ligand4 - 3pi/12, ligand5 - 4pi/12, ligand6 - 5pi/12

  do currentangle  = 1,nangles
     angles(currentangle)= (currentangle-1)*pi/(2d00*nangles)
     lposy(currentangle,1)=-1d0*sin(angles(currentangle))-downshift
     lposy(currentangle,2)= 0 - downshift
     lposy(currentangle,3)=-1d0*sin(angles(currentangle)) - downshift
     lposx(currentangle,1) = -1d0*cos(angles(currentangle))
     lposx(currentangle,2) = 0
     lposx(currentangle,3) = 1d0*cos(angles(currentangle))
  enddo
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! reading input data from files

  open(unit = 18, file = 'inputwedge/x13.dat')
  open(unit = 19, file = 'inputwedge/y13.dat')
  open(unit = 20, file = 'inputwedge/adj_list13.dat')
  open(unit = 22, file = 'inputwedge/ligand3.dat')
  open(unit = 23, file = 'inputwedge/amino16.dat') !! these are the generations
  !! I consider that all links are of length 1.



  !! reading the positions of the aminoacids
  read(18, *) x
  read(19, *) y
  close(18)
  close(19)
  do i=1,naa
     xp(2*i-1)=x(i)
     xp(2*i)=y(i)
  enddo
  xorig=xp  !! saving xp (for the test loops)

  !! read the 6 ligand types
  do currentligand=1,nligands
     read(22,*) ltemp  !! read the ligands
     do k=1,nlinks_lig_aa
        lamino_strength(currentligand,k)=ltemp(k)
        if(ltemp(k).ne.1 .and. ltemp(k).ne.2) then
           print *, "inconsistent ligand data", currentligand,k
           stop
        endif
      enddo
  enddo
  close(22)


!!!!!!!!!!!!!!!!!
  do k = 1, nlinks_lig_aa ! the 3 ligand binding sites
     llinks1(k) = lposxnumber(k)
     llinks2(k) = k
  enddo


  !! reading and creating link tables for fast gradient calculation
  do i=1,naa
     npartners(i)=0
  enddo
  do k=1,links
     read(20,*)links1(k),links2(k)
     if(links1(k)>naa .or. links2(k) >naa &
          .or.links1(k)<=0 .or. links2(k)<=0)then
        print*, "inconsistent links",links1(k),links2(k)
        stop
     endif
     if((x(links1(k))-x(links2(k)))**2 &
          +(y(links1(k))-y(links2(k)))**2-1>0.001)then
        print*, "long distance",k,links1(k),links2(k)
!        stop
     endif

     i1=links1(k)
     i2=links2(k)
     npartners(i1)=npartners(i1)+1
     npartners(i2)=npartners(i2)+1
     if(npartners(i1)>maxneighbors .or. npartners(i2)>maxneighbors) then
        print *, "too many neighbors"
        stop 
     endif
     partners(i1,npartners(i1))=i2
     partners(i2,npartners(i2))=i1
     partners2linknumber(i1,i2)=k
     partners2linknumber(i2,i1)=k
  enddo

  !! strength the links amino - ligands
  do i=1,links
     length(i)=sqrt((x(links1(i))-x(links2(i)))**2+&
     &(y(links1(i))-y(links2(i)))**2)
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! the main input of mutations
  !! we read here all possible mutations

  do currentprotein=1,nproteins    !! the loop over all proteins
     read(23,*)temp  !! read the mutation
     do k=1,naa_color
        amino_list(currentprotein,k)=temp(k)
        if(temp(k).ne.1 .and. temp(k).ne.2) then
           print *, "inconsistent generations", currentprotein,k
           stop
        endif
     enddo
  enddo
  close(23)





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! the output file(s)
  call cpu_time(start)
! write (dirname, '( "out", F06.2, "_", F06.2, "_",F06.2,"_",F06.2,"_",F06.2,"_",F06.2 )' )  els,elm,elw,ks,kw,kww
! 
! write (cmd,'( "mkdir -p out",F06.2,"_",F06.2,"_",F06.2,"_",F06.2,"_",F06.2,"_",F06.2)') els,elm,elw,ks,kw,kww
! print *,cmd
! call system(cmd)
! print *,dirname
! write (outname, '( "out", F06.2, "_",F06.2, "_",F06.2,"_",F06.2,"_",F06.2,"_",F06.2,"/output.dat" )' )  els,elm,elw,ks,kw,kww
! print *,outname

 
!  open(unit=41, file='outputifort0.60.3/xf13.dat')
!  open(unit=42, file='outputifort0.60.3/energy13.dat')
!  open(unit=43, file='outputifort0.60.3/deformation_energy.dat')
!  open(unit=44, file='outputifort0.60.3/chemical_energy.dat')
! open(unit=45, file=outname)
  open(unit=45, file="output.dat")


  write(45,'(a155)') &
       'protein, ligand, angle, lig coordinates(6(f6.3)),&
       protein no(i5),amino sequence(14(i1)),&
       total,deformation,chemical(3(f6.3)),protein coordinates(28(f6.3))'
  write(45,'("els=",F06.2," elm=",F06.2," elw=",F06.2," ks=",F06.2," kw=",F06.2," kww=",F06.2)')  els,elm,elw,ks,kw,kww


  xf = 0

!!! looping over all proteins inside
  !! loop aver all ligands
  !! loop over all angles
  counter=0
  do currentprotein=1,1    !! loop over proteins
! do currentprotein=1,nproteins   !! loop over proteins
     counter=counter+1;



     do aa =1,naa_color
        amino_strength(aa)=amino_list(currentprotein,aa)
     enddo


     ! compute linkstrength
     do k=1,links
        j1=amino_strength(links1(k))
        j2=amino_strength(links2(k))
        strength(k)= aatable(j1, j2)   
     enddo

     energies = 0
     
     cnt = 1

     do currentligand = 1, nligands  !! loop over ligands
        do k=1,nlinks_lig_aa
           j1=lamino_strength(currentligand,llinks2(k))
           j2=amino_strength(colors(k))
           lstrength(k)= laatable(j1, j2)
        enddo
        

        do currentangle = 1, nangles
           !! lxf - ligand coordinates
           !! creating the links amino - ligands
           do i=1,lbonds
              lx(i)=lposx(currentangle,i)+x(lposxnumber(2))
              ly(i)=lposy(currentangle,i)+y(lposxnumber(1))
              lxp0(2*i-1)=lx(i)
              lxp0(2*i)=ly(i)

           enddo

           if (currentligand .ne. 6) cycle

           !! tip of the ligand - chosen as origin
           x0 = lxp0(3)
           y0 = lxp0(4)

           etemp = 0

           xtemp = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!           ENERGY MINIMIZATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


           !! initialize protein coordinates			
           xp0 = xorig

           !! convert protein to polar coordinates with origin at tip of the ligand
           !! xp are polar coordinates - odd entries r, even are theta
           call aaxy2polar(xp0, x0, y0, angles(currentangle), xp)

           !! convert ligand to polar coordinates
           !! lxp are polar coordinates - odd entries r, even are theta
           call lxy2polar(lxp0, x0, y0, angles(currentangle), lxp)


           call totalenergy(xp,e)

           fold=e



           !!    comments from the gradient people
           !     We wish to have output at every iteration.

           !     We specify the tolerances in the stopping criteria.


           !     We specify the dimension n of the sample problem and the number
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
           do  i = 1, dimn
              if (mod(i,2) .eq. 0) then
                 nbd(i) = 2
                 lower(i)   = 0
                 upper(i)   = pi+2d0*angles(currentangle)
              else
                 nbd(i) = 1
                 lower(i)   = 0d0
                 upper(i)   = 10d0
              endif
           enddo


           !! !     We now define the starting point.
           !!  
           !! !     We now write the heading of the output.

           task='START'
111        continue
           call setulb(dimn,m,xp,lower,upper,nbd,f,derivative,factr,pgtol,wa,iwa,task,iprint, &
                csave,lsave,isave,dsave)


           if (task(1:2) .eq. 'FG') then

              call totalenergy(xp,f)


              fold=f

              call  ligandenergy(xp,ligands)

              do i=1,naa
                 do k=1,npartners(i)
                    j=partners(i,k)
                    !! |r1-r2| = ri^2 + r2^2 - 2r1*r2*cos(t1-t2)
                    energyold(i,k)=(sqrt(xp(2*i-1)**2+xp(2*j-1)**2 -xp(2*i-1)*xp(2*j-1)*2.0*cos(xp(2*i)-xp(2*j))) &
                         &-length(partners2linknumber(i,j)))**2*strength(partners2linknumber(i,j))

                 enddo
              enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !! this is where the huge time gain comes from
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


              !! we only compute the energy change for the affected links
              !! in "direction" j
              do j=1,dimn
                 i=int((j+1)/2)
                 xp(j)=xp(j)+epsilon
                 deltaenergy=0
                 do k=1,npartners(i)
                    jj=partners(i,k)
                    !! |r1-r2| = ri^2 + r2^2 - 2r1*r2*cos(t1-t2)
                    deltaenergy=deltaenergy &
                         +(sqrt(xp(2*i-1)**2+xp(2*jj-1)**2 -xp(2*i-1)*xp(2*jj-1)*&
                         &2.0*cos(xp(2*i)-xp(2*jj))) &
                         &-length(partners2linknumber(i,jj)))**2* &
                         strength(partners2linknumber(i,jj))-energyold(i,k)
                 enddo
                 call ligandenergy(xp,newligands)
                 derivative(j)=(deltaenergy+newligands-ligands)/epsilon
                 xp(j)=xp(j)-epsilon
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



           energy = f

           call  proteinenergy(xp,f)

           def_energy = f

           call  ligandenergy(xp,f)

           chem_energy = f

           !! convert polar back to cartesian coordinate
           call aapolar2xy(xp, x0, y0, angles(currentangle), xf)

            write(45,'(i6,i4,i4,x,6(f6.3,x),x,16(i1),3(e18.10,x),26(f6.3,x))') &
                counter,currentligand,currentangle,lxp0,amino_strength,&
                &energy,def_energy,&
                &chem_energy,xf
           !!,xf


        enddo !! currentangle=1,nangles
!       write(45,*)
        call cpu_time(finish)

!        print '("Time = ",f20.1," seconds.")',finish-start


     enddo ! currentligand=1,nligands
  enddo ! currentprotein=1,nproteins


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

  !! links1(k) = number of amino 1 in link
  !! links2(k) = number of amino 2 in link
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
    real*8 xp(2*naa)
    real*8 result, dr, rs
    integer*4 k,l
    result=0
    do k= 1,links
       result=result+ &
            (sqrt(xp(2*links1(k)-1)**2+xp(2*links2(k)-1)**2 &
            -xp(2*links1(k)-1)*xp(2*links2(k)-1)*2.0* &
            &cos(xp(2*links1(k))-xp(2*links2(k)))) &
            &-length(k))**2*strength(k)
    enddo
    !        print *,"---",result      

    do l=1,llinks
       dr = xp(2*llinks1(l)-1)**2+lxp(2*llinks2(l)-1)**2 &
            -xp(2*llinks1(l)-1)*lxp(2*llinks2(l)-1)*&
            &2.0*cos(xp(2*llinks1(l))-lxp(2*llinks2(l)))

       rs = dr/sigma2
       result=result - lstrength(l)*exp(-rs)

    enddo
    !        print *,"+++",result

  end subroutine totalenergy


!!!       only spring energies of the protein 

  subroutine proteinenergy (xp,result)

    implicit none
    real*8 xp(2*naa)
    real*8 result
    integer*4 k
    result=0
    do k= 1,links
       result=result+ &
            (sqrt(xp(2*links1(k)-1)**2+xp(2*links2(k)-1)**2 &
            -xp(2*links1(k)-1)*xp(2*links2(k)-1)*2.0* &
            &cos(xp(2*links1(k))-xp(2*links2(k)))) &
            &-length(k))**2*strength(k)
    enddo

  end subroutine proteinenergy

!!!   only lj chemical energies of the protein-ligand interaction

  subroutine ligandenergy (xp,result)

    implicit none
    real*8 xp(2*naa), dr, rs
    real*8 result
    integer*4 l
    result=0

    do l=1,llinks
       dr = xp(2*llinks1(l)-1)**2+lxp(2*llinks2(l)-1)**2 &
            -xp(2*llinks1(l)-1)*lxp(2*llinks2(l)-1)*&
            &2.0*cos(xp(2*llinks1(l))-lxp(2*llinks2(l)))

       rs = dr/sigma2
       result=result - lstrength(l)*exp(-rs)

    enddo
  end subroutine ligandenergy



  !! convert cartesian to polar coordinates for proteins
  !! rotate by angle theta to align x-axis along the ligand
  subroutine aaxy2polar(xp0, x0, y0, theta, xp)

    implicit none
    real*8 xp0(2*naa), x0, y0, xp(2*naa), theta
    integer*4 aa

    do aa = 1,naa
       xp(2*aa-1) = sqrt((xp0(2*aa-1)-x0)**2+(xp0(2*aa)-y0)**2)
       xp(2*aa) = atan2(xp0(2*aa)-y0,xp0(2*aa-1)-x0)+theta
    enddo
  end subroutine aaxy2polar

  !! convert cartesian to polar coordinates for ligands
  subroutine lxy2polar(lxp0, x0, y0, theta, lxp)

    implicit none
    real*8 lxp0(2*nlinks_lig_aa), x0, y0, lxp(2*nlinks_lig_aa), theta
    integer*4 ligand

    do ligand = 1,nlinks_lig_aa
       lxp(2*ligand-1) = sqrt((lxp0(2*ligand-1)-x0)**2+(lxp0(2*ligand)-y0)**2)
       lxp(2*ligand) = atan2(lxp0(2*ligand)-y0,lxp0(2*ligand-1)-x0)+theta
    enddo
  end subroutine lxy2polar

  !! convert polar to cartesian coordinates
  subroutine aapolar2xy(xp, x0, y0, theta, xf)

    implicit none
    real*8 xp(2*naa), x0, y0, xf(2*naa), theta
    integer*4 aa

    
    do aa = 1, naa
       xf(2*aa-1) = xp(2*aa-1)*cos(xp(2*aa)-theta) + x0
       xf(2*aa) = xp(2*aa-1)*sin(xp(2*aa)-theta) + y0
    enddo

  end subroutine aapolar2xy

end program energyprogram

