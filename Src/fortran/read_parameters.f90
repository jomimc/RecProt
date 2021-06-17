subroutine read_parameters(nproteins, namino, plinks, nligand, ks, kw, kww, kstep, kmod, els, elw, elm)

  implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! integer parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer (kind=4) :: namino  !! the number of amino acids
  integer (kind=4) :: nproteins !! the number of protein sequences
  integer (kind=4) :: plinks !! links within protein
  integer (kind=4) :: nligand !! the number of ligands
  integer (kind=4) :: nbind  !! the number of ligand binding sites



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
  close(11)

  elm = (els + elw) / 2.0d0 

end subroutine read_parameters
