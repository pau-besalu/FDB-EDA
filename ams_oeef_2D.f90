program ams_construct_2d
!This program was coded by Pau Besalú-Sala, 2022. PhD candidate from the
!University of Girona.
!It is linked with the script eef.sh which aims to Taylor-expansion the
!ASM model for the EF catalysis
!This software computes the expansions as well as the numerical
!derivatives for the 2D analysis
implicit none
integer :: i,j,numF,numpos
!double precision, parameter :: kcal=627.5096080d0 !1 a.u. = 627.5 kcal/mol
double precision, parameter :: kcal=627.50 !1 a.u. = 627.5 kcal/mol
double precision, parameter :: debye=0.393456 !1 debye = 0.39 a.u.
double precision :: be1, be2, be3, be4 !bonding energies
double precision, dimension (3) :: deriv1,deriv2,err1,err2 !RR 1st and 2nd derivatives, and associated numerical errors
double precision, dimension (3,9) :: s1,s2,stot 
double precision, dimension (9) :: int1,int2,int3,inttot
double precision, dimension (9) :: f
double precision, dimension (4,3) :: d !dipole moments
double precision, dimension (4,3,3) :: a ! read polarizabilities, polarizabilities
double precision, allocatable, dimension (:,:,:) :: iex,iey,iez !interaction energies for 1d 
double precision, dimension (2) :: ie2d !reading the energies from summary file
double precision, allocatable, dimension (:,:) :: pauli2De,elstat2De,oi2De !2D energies (numF,identifyer,value)
double precision :: E0 !energy zero for 2D derivatives
double precision, dimension (3) :: first,strf,paulif,elstatf,oif !first derivative for 2D analysis
double precision, dimension (3,3) :: second,strs,paulis,elstats,ois ! second derivative for 2D analysis
double precision, allocatable, dimension (:) :: f2d !field for 2D calculations 
double precision :: minf,step
double precision :: stre0,paulie0,elstate0,oie0
integer :: grid2d
character*80 :: junk1,junk2,junk3,junk4 

!this is read from the screen
!those will be the nine outputed fields
write(*,*) "Please, enter nine fields"
do i=1,9
  read(*,*) f(i)
end do
read(*,*) grid2d
write(*,'("Fields displayed will be: ",9es11.2)') f(:)
write(*,'("Number of points in 2D map will be: ",I8)') grid2d**2

open(unit=1,file='summary.dat',status="old")
 !Reading basic info
 read(1,*) junk1
 read(1,*) numpos 
 numF=numpos*2+1
 !3 Types of Interaction energy (Elstat, Pauli,OI), 2 columns (field and energy
 !value) and numF number of rows
 allocate(pauli2De(numF,13),elstat2De(numF,13),oi2De(numF,13))
 allocate(iex(3,numF,2))
 allocate(iey(3,numF,2))
 allocate(iez(3,numF,2))

 read(1,*) junk1
 read(1,*) junk1
 !reading bonding energies of the distorted/optimized subsystems
 read(1,*) junk1,junk2,junk3,junk4,be1
 read(1,*) junk1,junk2,junk3,junk4,be2
 read(1,*) junk1,junk2,junk3,junk4,be3
 read(1,*) junk1,junk2,junk3,junk4,be4
 !Reading the dipole moments of the four compounds
 do i=1,4 !for each system
   read(1,*) junk1
   read(1,*) junk1,d(i,1)
   read(1,*) junk1,d(i,2)
   read(1,*) junk1,d(i,3)
 end do
 do i=1,4 !for each system read the polarizability 
   read(1,*) junk1
   read(1,*) junk1
   read(1,*) a(i,1,1),a(i,1,2),a(i,1,3)
   read(1,*) a(i,2,1),a(i,2,2),a(i,2,3)
   read(1,*) a(i,3,1),a(i,3,2),a(i,3,3)
 end do
 do i=1,18 !for each energy component
   read(1,*) junk1
   read(1,*) junk1
   read(1,*) junk1
   do j=1,numF
      read(1,*) ie2d(1),ie2d(2)

      !distribute the data
      !the first 10 if's are 1D data
      if (i==1) then !is X elstat
        iex(1,j,1)=ie2d(1); iex(1,j,2)=ie2d(2)
        if (j==1) then
          elstat2De(1,1)=ie2d(2)
          !write(*,*) j,"1",elstat2De(j,1)
          !write(*,*) "enter",j,1
        else if (j.gt.1.and.j.le.(numpos+1)) then
          elstat2De(j,2)=ie2d(2)
          !write(*,*) "enter",j,1
          !write(*,*) j,"2",elstat2De(j,2)
        else
          elstat2De(j,3)=ie2d(2)
          !write(*,*) "enter",j,1
          !write(*,*) j,"3",elstat2De(j,3)
        end if

      else if (i==2) then !is X pauli
        iex(2,j,1)=ie2d(1); iex(2,j,2)=ie2d(2)
        if (j==1) then
          pauli2De(1,1)=ie2d(2)
          !write(*,*) j,"1",pauli2De(j,1)
        else if (j.gt.1.and.j.le.(numpos+1)) then
          pauli2De(j,2)=ie2d(2)
          !write(*,*) j,"2",pauli2De(j,2)
        else
          pauli2De(j,3)=ie2d(2)
          !write(*,*) j,"3",pauli2De(j,3)          
        end if

      else if (i==3) then ! is X oi
        iex(3,j,1)=ie2d(1); iex(3,j,2)=ie2d(2)
        if (j==1) then
          oi2De(1,1)=ie2d(2)
          !write(*,*) j,"1",oi2De(j,1)          
        else if (j.gt.1.and.j.le.(numpos+1)) then
          oi2De(j,2)=ie2d(2)
          !write(*,*) j,"2",oi2De(j,2)          
        else
          oi2De(j,3)=ie2d(2)
          !write(*,*) j,"3",oi2De(j,3)          
        end if

      else if (i==4) then ! is Y elstat
        iey(1,j,1)=ie2d(1); iey(1,j,2)=ie2d(2)
        if (j==1) then
          elstat2De(1,1)=ie2d(2)
          !write(*,*) j,"1",elstat2De(j,1)
        else if (j.gt.1.and.j.le.(numpos+1)) then
          elstat2De(j,4)=ie2d(2)
          !write(*,*) j,"4",elstat2De(j,4)
        else
          elstat2De(j,5)=ie2d(2)
          !write(*,*) j,"5",elstat2De(j,5)
        end if
        
      else if (i==5) then ! is Y pauli
        iey(2,j,1)=ie2d(1); iey(2,j,2)=ie2d(2)
        if (j==1) then
          pauli2De(1,1)=ie2d(2)
          !write(*,*) j,"1",pauli2De(j,1)
        else if (j.gt.1.and.j.le.(numpos+1)) then
          pauli2De(j,4)=ie2d(2)
          !write(*,*) j,"4",pauli2De(j,4)
        else
          pauli2De(j,5)=ie2d(2)
          !write(*,*) j,"5",pauli2De(j,5)
        end if
        
      else if (i==6) then ! is Y oi
        iey(3,j,1)=ie2d(1); iey(3,j,2)=ie2d(2)
        if (j==1) then
          oi2De(1,1)=ie2d(2)
          !write(*,*) j,"1",oi2De(j,1)
        else if (j.gt.1.and.j.le.(numpos+1)) then
          oi2De(j,4)=ie2d(2)
          !write(*,*) j,"4",oi2De(j,4)
        else
          oi2De(j,5)=ie2d(2)
          !write(*,*) j,"5",oi2De(j,5)
        end if

      else if (i==7) then ! is Z elstat
        iez(1,j,1)=ie2d(1); iez(1,j,2)=ie2d(2)
        if (j==1) then
          elstat2De(1,1)=ie2d(2)
          !write(*,*) j,"1",elstat2De(j,1)
        else if (j.gt.1.and.j.le.(numpos+1)) then
          elstat2De(j,6)=ie2d(2)
          !write(*,*) j,"6",elstat2De(j,6)
        else
          elstat2De(j,7)=ie2d(2)
          !write(*,*) j,"7",elstat2De(j,7)
        end if

      else if (i==8) then ! is Z pauli
        iez(2,j,1)=ie2d(1); iez(2,j,2)=ie2d(2)
        if (j==1) then
          pauli2De(1,1)=ie2d(2)
          !write(*,*) j,"1",pauli2De(j,1)
        else if (j.gt.1.and.j.le.(numpos+1)) then
          pauli2De(j,6)=ie2d(2)
          !write(*,*) j,"6",pauli2De(j,6)
        else
          pauli2De(j,7)=ie2d(2)
          !write(*,*) j,"7",pauli2De(j,7)
        end if
        
      else if (i==9) then ! is Z oi
        iez(3,j,1)=ie2d(1); iez(3,j,2)=ie2d(2)
        if (j==1) then
          oi2De(1,1)=ie2d(2)
        else if (j.gt.1.and.j.le.(numpos+1)) then
          oi2De(j,6)=ie2d(2)
        else
          oi2De(j,7)=ie2d(2)
        end if

      ! from here is 2D data
      else if (i==10) then !is XY elstat
        if (j==1) then
          elstat2De(1,1)=ie2d(2)
        else if (j.gt.1.and.j.le.(numpos+1)) then
          elstat2De(j,8)=ie2d(2)
        else
          elstat2De(j,9)=ie2d(2)
        end if

      else if (i==11) then !is XY pauli
        if (j==1) then
          pauli2De(1,1)=ie2d(2)
        else if (j.gt.1.and.j.le.(numpos+1)) then
          pauli2De(j,8)=ie2d(2)
        else
          pauli2De(j,9)=ie2d(2)
        end if

      else if (i==12) then !is XY oi
        if (j==1) then
          oi2De(1,1)=ie2d(2)
        else if (j.gt.1.and.j.le.(numpos+1)) then
          oi2De(j,8)=ie2d(2)
        else
          oi2De(j,9)=ie2d(2)
        end if

      else if (i==13) then !is XZ elstat
        if (j==1) then
          elstat2De(1,1)=ie2d(2)
        else if (j.gt.1.and.j.le.(numpos+1)) then
          elstat2De(j,10)=ie2d(2)
        else
          elstat2De(j,11)=ie2d(2)
        end if
              
      else if (i==14) then !is XZ pauli
        if (j==1) then
          pauli2De(1,1)=ie2d(2)
        else if (j.gt.1.and.j.le.(numpos+1)) then
          pauli2De(j,10)=ie2d(2)
        else
          pauli2De(j,11)=ie2d(2)
        end if

      else if (i==15) then !is XZ oi
        if (j==1) then
          oi2De(1,1)=ie2d(2)
        else if (j.gt.1.and.j.le.(numpos+1)) then
          oi2De(j,10)=ie2d(2)
        else
          oi2De(j,11)=ie2d(2)
        end if

      else if (i==16) then !is YZ elstat
        if (j==1) then
          elstat2De(1,1)=ie2d(2)
        else if (j.gt.1.and.j.le.(numpos+1)) then
          elstat2De(j,12)=ie2d(2)
        else
          elstat2De(j,13)=ie2d(2)
        end if

      else if (i==17) then !is YZ pauli
        if (j==1) then
          pauli2De(1,1)=ie2d(2)
          !write(*,*) j,"1",pauli2De(j,1)
        else if (j.gt.1.and.j.le.(numpos+1)) then
          pauli2De(j,12)=ie2d(2)
          !write(*,*) j,"12",pauli2De(j,12)
        else
          pauli2De(j,13)=ie2d(2)
          !write(*,*) j,"13",pauli2De(j,13)
        end if

      else if (i==18) then !is YZ oi
        if (j==1) then
          oi2De(1,1)=ie2d(2)
        else if (j.gt.1.and.j.le.(numpos+1)) then
          oi2De(j,12)=ie2d(2)
        else
          oi2De(j,13)=ie2d(2)
        end if

      else
        write(*,*)"Error on data collection"
      endif
   enddo
 end do

 do i=1,3 ; do j=1,numF !the fields (column 1) are scaled for the 1D analysis
   iex(i,j,1)=iex(i,j,1)*1.0d-4
   iey(i,j,1)=iey(i,j,1)*1.0d-4 
   iez(i,j,1)=iez(i,j,1)*1.0d-4 
 end do; end do
close(1)
! End of summary.dat reading

!if you want, just uncomment and use it
!call print_summary_read(elstat2De,pauli2De,oi2De,numF)

!Preparing output
write(*,*)
write(*,*) "   ASM analysis of OEEF catalysis  "
write(*,*) "   ------------------------------  "
write(*,*)
write(*,*) "   Strain energies at F=0 for each reactant (kcal/mol)"
write(*,'("                Reactant 1: ",es18.9)') (be1-be3)*kcal
write(*,'("                Reactant 2: ",es18.9)') (be2-be4)*kcal
write(*,*)
write(*,*) "   The 1st derivative of the strain wrt the field is (kcal/mol)" 
write(*,'("                Reactant 1 ")')
write(*,'("                wrt Fx: ", es18.9)') (d(1,1)-d(3,1))*debye*kcal
write(*,'("                wrt Fy: ", es18.9)') (d(1,2)-d(3,2))*debye*kcal
write(*,'("                wrt Fz: ", es18.9)') (d(1,3)-d(3,3))*debye*kcal
write(*,*)
write(*,'("                Reactant 2 ")')
write(*,'("                wrt Fx: ", es18.9)') (d(2,1)-d(4,1))*debye*kcal
write(*,'("                wrt Fy: ", es18.9)') (d(2,2)-d(4,2))*debye*kcal
write(*,'("                wrt Fz: ", es18.9)') (d(2,3)-d(4,3))*debye*kcal
write(*,*) 
write(*,*)
write(*,*) "   The 2nd derivative of the strain wrt the field is (kcal/mol)"
write(*,*) "   Please, note that the crossed-derivatives are not printed nor coded"
write(*,'("                Reactant 1 ")')
write(*,'("                wrt Fx^2: ", es18.9)') (a(1,1,1)-a(3,1,1))*kcal
write(*,'("                wrt Fy^2: ", es18.9)') (a(1,2,2)-a(3,2,2))*kcal
write(*,'("                wrt Fz^2: ", es18.9)') (a(1,3,3)-a(3,3,3))*kcal
write(*,*)
write(*,'("                Reactant 2 ")')
write(*,'("                wrt Fx^2: ", es18.9)') (a(2,1,1)-a(4,1,1))*kcal
write(*,'("                wrt Fy^2: ", es18.9)') (a(2,2,2)-a(4,2,2))*kcal
write(*,'("                wrt Fz^2: ", es18.9)') (a(2,3,3)-a(4,3,3))*kcal
write(*,*)

!Compute the strain energy derivatives (up to polar).
! E_i(F)=E_i(0)-mu_i*F-(1/2)*a_i*F**2 where _i=TS,R
! dE(F)=E_TS(F)-E_R(F) 

do i=1,9 !all tested fields
  do j=1,3 !all directions
    s1(j,i)=((be1-be3)-(d(1,j)-d(3,j))*debye*f(i)-0.5d0*(a(1,j,j)-a(3,j,j))*f(i)**2)*kcal
    s2(j,i)=((be2-be4)-(d(2,j)-d(4,j))*debye*f(i)-0.5d0*(a(2,j,j)-a(4,j,j))*f(i)**2)*kcal
    stot(j,i)=s1(j,i)+s2(j,i)
  enddo
end do

write(*,*)
write(*,*) "   With this information, we can compute the strain(F) "
write(*,'("    Fx (a.u.)      Strain(F) R1       Strain(F) R2       Strain(F)")')
write(*,'("    ---------      ------------       ------------       ---------")')
do i=1,9
write(*,'("  ",es11.2,3x,es18.9,x,es18.9,x,es18.9)') f(i),s1(1,i),s2(1,i),stot(1,i)
enddo
write(*,*) 
write(*,'("    Fy (a.u.)      Strain(F) R1       Strain(F) R2       Strain(F)")')
write(*,'("    ---------      ------------       ------------       ---------")')
do i=1,9
write(*,'("  ",es11.2,3x,es18.9,x,es18.9,x,es18.9)') f(i),s1(2,i),s2(2,i),stot(2,i)
enddo
write(*,*)
write(*,'("    Fz (a.u.)    Strain(F) R1         Strain(F) R2       Strain(F)")')
write(*,'("    ---------    ------------         ------------       ---------")')
do i=1,9
write(*,'("  ",es11.2,3x,es18.9,x,es18.9,x,es18.9)') f(i),s1(3,i),s2(3,i),stot(3,i)
enddo
write(*,*)
write(*,*)
write(*,*) "The interaction energy in function of the EF is the following:"
write(*,*)

! Compute the Taylor expansion for the interaction energy
deriv1=0.0d0; deriv2=0.0d0 ; err1=0.0d0 ; err2=0.0d0
call RR(numF,iex,deriv1,deriv2,err1,err2,kcal)
!initialize all values
elstate0=0.0d0;paulie0=0.0d0;oie0=0.0d0
elstatf=1.0d9;paulif=1.0d9;oif=1.0d9
elstats=1.0d9;paulis=1.0d9;ois=1.0d9
!save E0's and 1st derivs along X for later
elstate0=iex(1,1,2)
paulie0=iex(2,1,2)
oie0=iex(3,1,2)
elstatf(1)=deriv1(1)/kcal
paulif(1)=deriv1(2)/kcal
oif(1)=deriv1(3)/kcal
elstats(1,1)=deriv2(1)/kcal
paulis(1,1)=deriv2(2)/kcal
ois(1,1)=deriv2(3)/kcal

!Compute the interaction energy (contributions and total) for the direction of
!the field applied
do i=1,9 !all tested fields
  int1(i)=iex(1,1,2)*kcal+deriv1(1)*f(i)+0.5d0*deriv2(1)*f(i)**2
  int2(i)=iex(2,1,2)*kcal+deriv1(2)*f(i)+0.5d0*deriv2(2)*f(i)**2
  int3(i)=iex(3,1,2)*kcal+deriv1(3)*f(i)+0.5d0*deriv2(3)*f(i)**2
  inttot(i)=int1(i)+int2(i)+int3(i)
end do

call print_IntEnDecomp(deriv1,deriv2,err1,err2,f,inttot,int1,int2,int3,"X")

! Compute the Taylor expansion for the interaction energy
deriv1=0.0d0; deriv2=0.0d0 ; err1=0.0d0 ; err2=0.0d0
call RR(numF,iey,deriv1,deriv2,err1,err2,kcal)
!save 1st and diagonal second derivs along Y for later
elstatf(2)=deriv1(1)/kcal
elstats(2,2)=deriv2(1)/kcal
paulif(2)=deriv1(2)/kcal
paulis(2,2)=deriv2(2)/kcal
oif(2)=deriv1(3)/kcal
ois(2,2)=deriv2(3)/kcal

!Compute the interaction energy (contributions and total) for the direction of
!the field applied
do i=1,9 !all tested fields
  int1(i)=iey(1,1,2)*kcal+deriv1(1)*f(i)+0.5d0*deriv2(1)*f(i)**2
  int2(i)=iey(2,1,2)*kcal+deriv1(2)*f(i)+0.5d0*deriv2(2)*f(i)**2
  int3(i)=iey(3,1,2)*kcal+deriv1(3)*f(i)+0.5d0*deriv2(3)*f(i)**2
  inttot(i)=int1(i)+int2(i)+int3(i)
end do

call print_IntEnDecomp(deriv1,deriv2,err1,err2,f,inttot,int1,int2,int3,"Y")

! Compute the Taylor expansion for the interaction energy
deriv1=0.0d0; deriv2=0.0d0 ; err1=0.0d0 ; err2=0.0d0
call RR(numF,iez,deriv1,deriv2,err1,err2,kcal)
!save 1st derivs along Z for later
elstatf(3)=deriv1(1)/kcal
elstats(3,3)=deriv2(1)/kcal
paulif(3)=deriv1(2)/kcal
paulis(3,3)=deriv2(2)/kcal
oif(3)=deriv1(3)/kcal
ois(3,3)=deriv2(3)/kcal

!Compute the interaction energy (contributions and total) for the direction of
!the field applied
do i=1,9 !all tested fields
  int1(i)=iez(1,1,2)*kcal+deriv1(1)*f(i)+0.5d0*deriv2(1)*f(i)**2
  int2(i)=iez(2,1,2)*kcal+deriv1(2)*f(i)+0.5d0*deriv2(2)*f(i)**2
  int3(i)=iez(3,1,2)*kcal+deriv1(3)*f(i)+0.5d0*deriv2(3)*f(i)**2
  inttot(i)=int1(i)+int2(i)+int3(i)
end do

call print_IntEnDecomp(deriv1,deriv2,err1,err2,f,inttot,int1,int2,int3,"Z")
! end of 1D analysis

!!!!! 2D analysis and priting
! We require the field and the 1st, 2nd derivatives including crossed terms.
! Then we can compute whatever we want with this input

!preparation of the field grid
minf=-abs(f(1))
step=2*abs(minf)/(grid2d-1)
allocate(f2d(grid2d))
f2d=0.0d0; f2d(1)=minf
do i=2,grid2d; f2d(i)=f2d(i-1)+step; enddo
!write(*,*) "Fields in 2D will be :" ; do i=1,grid2d ;  write(*,*) f2d(i) ; end do

! compute 2D strain of R1 in the planes XY XZ and YZ
E0=be1-be3
first(1)=-(d(1,1)-d(3,1))*debye
first(2)=-(d(1,2)-d(3,2))*debye
first(3)=-(d(1,3)-d(3,3))*debye
second(1,1)=-(a(1,1,1)-a(3,1,1))
second(1,2)=-(a(1,1,2)-a(3,1,2))
second(1,3)=-(a(1,1,3)-a(3,1,3))
second(2,1)=second(1,2)
second(2,2)=-(a(1,2,2)-a(3,2,2))
second(2,3)=-(a(1,2,3)-a(3,2,3))
second(3,1)=second(1,3)
second(3,2)=second(2,3)
second(3,3)=-(a(1,3,3)-a(3,3,3))
call compute_2d(E0,first,second,f2d,kcal,grid2d)
call system("mv EXY.dat StrR1_XY.dat")
call system("mv EXZ.dat StrR1_XZ.dat")
call system("mv EYZ.dat StrR1_YZ.dat")

! compute 2D strain of R2 in the planes XY XZ and YZ
E0=be2-be4
first(1)=-(d(2,1)-d(4,1))*debye
first(2)=-(d(2,2)-d(4,2))*debye
first(3)=-(d(2,3)-d(4,3))*debye
second(1,1)=-(a(2,1,1)-a(4,1,1))
second(1,2)=-(a(2,1,2)-a(4,1,2))
second(1,3)=-(a(2,1,3)-a(4,1,3))
second(2,1)=second(1,2)
second(2,2)=-(a(2,2,2)-a(4,2,2))
second(2,3)=-(a(2,2,3)-a(4,2,3))
second(3,1)=second(1,3)
second(3,2)=second(2,3)
second(3,3)=-(a(2,3,3)-a(4,3,3))
call compute_2d(E0,first,second,f2d,kcal,grid2d)
call system("mv EXY.dat StrR2_XY.dat")
call system("mv EXZ.dat StrR2_XZ.dat")
call system("mv EYZ.dat StrR2_YZ.dat")

! compute 2D strain of R2 in the planes XY XZ and YZ
stre0=(be1-be3)+(be2-be4)
strf(1)=-((d(2,1)-d(4,1))+(d(1,1)-d(3,1)))*debye
strf(2)=-((d(2,2)-d(4,2))+(d(1,2)-d(3,2)))*debye
strf(3)=-((d(2,3)-d(4,3))+(d(1,3)-d(3,3)))*debye
strs(1,1)=-(a(2,1,1)-a(4,1,1)+a(1,1,1)-a(3,1,1))
strs(1,2)=-(a(2,1,2)-a(4,1,2)+a(1,1,2)-a(3,1,2))
strs(1,3)=-(a(2,1,3)-a(4,1,3)+a(1,1,3)-a(3,1,3))
strs(2,1)=second(1,2)
strs(2,2)=-(a(2,2,2)-a(4,2,2)+a(1,2,2)-a(3,2,2))
strs(2,3)=-(a(2,2,3)-a(4,2,3)+a(1,2,3)-a(3,2,3))
strs(3,1)=second(1,3)
strs(3,2)=second(2,3)
strs(3,3)=-(a(2,3,3)-a(4,3,3)+a(1,3,3)-a(3,3,3))
call compute_2d(stre0,strf,strs,f2d,kcal,grid2d)
call system("mv EXY.dat Str_XY.dat")
call system("mv EXZ.dat Str_XZ.dat")
call system("mv EYZ.dat Str_YZ.dat")

!print grid for isocontour lines
!call printgrid(grid2d,f2d)

! compute 2D Pauli in the planes XY XZ and YZ ! compute first the 2nd crossed derivatives
! now we compute the crossed terms
err2=0.0d0
call RR_2D(numF,pauli2De,paulis)
!do i=1,3; write(*,*) (paulis(i,j),j=1,3); end do

! now we do the expansion
call compute_2d(paulie0,paulif,paulis,f2d,kcal,grid2d)
call system("mv EXY.dat pau_XY.dat")
call system("mv EXZ.dat pau_XZ.dat")
call system("mv EYZ.dat pau_YZ.dat")

! compute 2D elstat in the planes XY XZ and YZ ! compute first the 2nd crossed derivatives
! now we compute the crossed terms
err2=0.0d0
call RR_2D(numF,elstat2De,elstats)
!do i=1,3; write(*,*) (elstats(i,j),j=1,3); end do

! now we do the expansion
call compute_2d(elstate0,elstatf,elstats,f2d,kcal,grid2d)
call system("mv EXY.dat elst_XY.dat")
call system("mv EXZ.dat elst_XZ.dat")
call system("mv EYZ.dat elst_YZ.dat")

! compute 2D OI in the planes XY XZ and YZ ! compute first the 2nd crossed derivatives
! now we compute the crossed terms
err2=0.0d0
call RR_2D(numF,oi2De,ois)
!do i=1,3; write(*,*) (ois(i,j),j=1,3); end do

! now we do the expansion
call compute_2d(oie0,oif,ois,f2d,kcal,grid2d)
call system("mv EXY.dat oi_XY.dat")
call system("mv EXZ.dat oi_XZ.dat")
call system("mv EYZ.dat oi_YZ.dat")

!use OI variable as the one for total interaction
oie0=oie0+elstate0+paulie0
do i=1,3; oif(i)=oif(i)+elstatf(i)+paulif(i); enddo
do i=1,3;do j=1,3; ois(i,j)=ois(i,j)+elstats(i,j)+paulis(i,j); enddo;enddo

! compute 2D OI in the planes XY XZ and YZ ! compute first the 2nd crossed derivatives
! now we compute the crossed terms
err2=0.0d0
call RR_2D(numF,oi2De,ois)
!do i=1,3; write(*,*) (ois(i,j),j=1,3); end do

! now we do the expansion
call compute_2d(oie0,oif,ois,f2d,kcal,grid2d)
call system("mv EXY.dat inten_XY.dat")
call system("mv EXZ.dat inten_XZ.dat")
call system("mv EYZ.dat inten_YZ.dat")

!use OI variable as the one for total energy (it was already the total int. energy)
oie0=oie0+stre0
do i=1,3; oif(i)=oif(i)+strf(i); enddo
do i=1,3;do j=1,3; ois(i,j)=ois(i,j)+strs(i,j); enddo;enddo

! compute 2D OI in the planes XY XZ and YZ ! compute first the 2nd crossed derivatives
! now we compute the crossed terms
err2=0.0d0
call RR_2D(numF,oi2De,ois)
!do i=1,3; write(*,*) (ois(i,j),j=1,3); end do

! now we do the expansion
call compute_2d(oie0,oif,ois,f2d,kcal,grid2d)
call system("mv EXY.dat toten_XY.dat")
call system("mv EXZ.dat toten_XZ.dat")
call system("mv EYZ.dat toten_YZ.dat")

deallocate(f2d)
deallocate(iex,iey,iez)
deallocate(pauli2De,elstat2De,oi2De)

end

!print the output for the 1D analysis
subroutine print_IntEnDecomp(deriv1,deriv2,err1,err2,f,inttot,int1,int2,int3,direction)
implicit none
integer :: i
double precision, dimension (3) :: deriv1, deriv2, err1, err2 !output values
double precision, dimension (9) :: int1,int2,int3,inttot,f
character*1 :: direction

write(*,'( "FIELD DIRECTION IS :",a)') direction
write(*,*) "Derivatives of the Electrostatic Interaction (kcal/mol)"
write(*,*) "       Derivative          Num. Error"
write(*,'( " 1st ",es18.9,x,es18.9)') deriv1(1), err1(1)
write(*,'( " 2nd ",es18.9,x,es18.9)') deriv2(1), err2(1)
write(*,*)
write(*,*) "Derivatives of the Pauli Repulsion (kcal/mol)"
write(*,*) "       Derivative          Num. Error"
write(*,'( " 1st ",es18.9,x,es18.9)') deriv1(2), err1(2)
write(*,'( " 2nd ",es18.9,x,es18.9)') deriv2(2), err2(2)
write(*,*)
write(*,*) "Derivatives of the Orbital Interaction (kcal/mol)"
write(*,*) "       Derivative          Num. Error"
write(*,'( " 1st ",es18.9,x,es18.9)') deriv1(3), err1(3)
write(*,'( " 2nd ",es18.9,x,es18.9)') deriv2(3), err2(3)
write(*,*)
write(*,'("    F  (a.u.)        Elstat(F)            Pauli(F)            OI(F)      ")')
write(*,'("    ---------      ------------         ------------       ---------     ")')
do i=1,9
write(*,'("  ",es11.2,3x,es18.9,x,es18.9,x,es18.9)') f(i),int1(i),int2(i),int3(i)
enddo
write(*,*)
write(*,*) " Total interaction energy(F) "
do i=1,9
write(*,'("  ",es11.2,3x,es18.9,x,es18.9,x,es18.9)') f(i),inttot(i)
end do
write(*,*)

end

!Romberg-Rutishauser method for the numerical derivatives.
subroutine RR(numF,ie,deriv1,deriv2,err1,err2,kcal)
implicit none
integer :: numF,numpos !number of EF, number of positive fields
integer :: i,j,k
double precision :: kcal
double precision, dimension (3,numF,2) :: ie !interaction energy tensor
double precision, dimension (3) :: deriv1, deriv2, err1, err2 !output values
integer, dimension (3) :: savej,savek
double precision, dimension (3,((numF-1)/2)-1,((numF-1)/2)-1) :: r,er1,er2
double precision, dimension (3,(numF-1)/2) :: ff
double precision, dimension (3,(numF-1)/2,(numF-1)/2) :: disp1,disp2

!Here we start the numerical derivatives for the three contributions within the
!interaction energy (Pauli repulsion, orbital interaction and electrostatic
!interaction)
!ff,r,err are the FF vector, the RR triangle and the error triangle
!respectively (for each energy component).

!initialization
numpos=(numF-1)/2
ff=1.0d8 ; r=1.0d8 

!compute first derivative, FF
do i=1,3 !for each component
  do j=1,numpos
    ff(i,j)=(ie(i,j+1,2)-ie(i,j+1+numpos,2))/(2*ie(i,j+1,1))
  end do
end do

!Compute the RR correction to FF vector
do i=1,3 !for each energy component

  !first RR iteration is done respect to the FF vector
  j=1
  do while (j<=(numpos-1))
     r(i,1,j)=(4*ff(i,j)-ff(i,j+1))/3.0d0
     j=j+1     
  end do
 
  ! from 2nd RR iteration onwards is is computed from the r tensor 
  do k=2,numpos-1 !each RR iteration
     j=1 ! used for counting
     do while (j<=(numpos-k))  !Only upper side triangle
       r(i,k,j)=(4**k*r(i,k-1,j)-r(i,k-1,j+1))/(4**k-1)
       j=j+1
     end do
  end do
end do

!save the triangles to a new array for displaying purposes (compact array with
!FF and RR matrixes)
disp1=1.0d8
do i=1,3 ! for each component
  do k=1,numpos !each RR iteration
    do j=1,numpos
      if (k==1) then
        disp1(i,k,j)=ff(i,j)
      else if (j<=(numpos-1)) then !the sizes of the tensors are different
        disp1(i,k,j)=r(i,k-1,j)
      end if
    end do
  end do
end do

!Compute error triangle from the display triangle
er1=1.0d8 ; err1=1.0d8
do i=1,3 !for each energy component
 do k=1,numpos-1 !each iteration
   do j=1,numpos-1 !each field
     er1(i,k,j)=abs(disp1(i,k,j+1)-disp1(i,k,j))

     ! this is to avoid having zeros in the lower part of the triangle (to catch the minimal error this must be avoided)
     if (er1(i,k,j)==0.0d0) er1(i,k,j)=1.0d8 

     !save the smallest error
     if (er1(i,k,j)<err1(i)) then
       savek(i)=k
       savej(i)=j
       err1(i)=er1(i,k,j)
     end if 
   end do
 end do
end do

!save final results for first derivative
do i=1,3 ; deriv1(i)=disp1(i,savek(i),savej(i)+1) ; end do

!DISPLAY FOR THE USER TO SEE
call display_triangles(1,numpos,disp1,er1,deriv1,err1,savek,savej)

!Convert to kcal/mol
do i=1,3
  deriv1(i)=deriv1(i)*kcal
  err1(i)=err1(i)*kcal
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! End of first derivative. Starts second derivative calculation!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!initialization
ff=1.0d8 ; r=1.0d8

!compute second derivative, FF
do i=1,3 !for each component
  do j=1,numpos
    ff(i,j)=(ie(i,j+1,2)+ie(i,j+1+numpos,2)-2.0d0*ie(i,1,2))/(ie(i,j+1,1)**2)
  end do
end do

!Compute the RR correction to FF vector
do i=1,3 !for each energy component

  !first RR iteration is done respect to the FF vector
  j=1
  do while (j<=(numpos-1))
     r(i,1,j)=(4*ff(i,j)-ff(i,j+1))/3.0d0
     j=j+1
  end do

  ! from 2nd RR iteration onwards is is computed from the r tensor
  do k=2,numpos-1 !each RR iteration
     j=1 ! used for counting
     do while (j<=(numpos-k))  !Only upper side triangle
       r(i,k,j)=(4**k*r(i,k-1,j)-r(i,k-1,j+1))/(4**k-1)
       j=j+1
     end do
  end do
end do

!save the triangles to a new array for displaying purposes (compact array with
!FF and RR matrixes)
disp2=1.0d8
do i=1,3 ! for each component
  do k=1,numpos !each RR iteration
    do j=1,numpos
      if (k==1) then
        disp2(i,k,j)=ff(i,j)
      else if (j<=(numpos-1)) then !the sizes of the tensors are different
        disp2(i,k,j)=r(i,k-1,j)
      end if
    end do
  end do
end do

!Compute error triangle from the display triangle
er2=1.0d8 ; err2=1.0d8
do i=1,3 !for each energy component
 do k=1,numpos-1 !each iteration
   do j=1,numpos-1 !each field
     er2(i,k,j)=abs(disp2(i,k,j+1)-disp2(i,k,j))

     ! this is to avoid having zeros in the lower part of the triangle (to catch
     ! the minimal error this must be avoided)
     if (er2(i,k,j)==0.0d0) er2(i,k,j)=1.0d8

     !save the smallest error
     if (er2(i,k,j)<err2(i)) then
       savek(i)=k
       savej(i)=j
       err2(i)=er2(i,k,j)
     end if
   end do
 end do
end do

!save final results for first derivative
do i=1,3 ; deriv2(i)=disp2(i,savek(i),savej(i)+1) ; end do

!DISPLAY FOR THE USER TO SEE
call display_triangles(2,numpos,disp2,er2,deriv2,err2,savek,savej)

!Convert to kcal/mol
do i=1,3
  deriv2(i)=deriv2(i)*kcal
  err2(i)=err2(i)*kcal
end do

end

subroutine display_triangles(n,s,d,e,save_deriv,save_er,savek,savej)
implicit none
integer :: s ! size of the deriv matrix
integer :: i,n, j, k,l,print_size,position
integer, dimension (3) :: savej,savek
double precision, dimension (3) :: save_er,save_deriv !output values
double precision, dimension (3,s,s) :: d,p !derivatives matrix
double precision, dimension (3,s-1,s-1) :: e !error matrix
character(len=100) :: format,nchar
character(len=103) :: ch
logical :: filexist

if (s>12) then
  write(*,*) "NOT POSSIBLE TO PRINT. TOO BIG MATRIX"
  return
else
 
!remove possible old files
write(nchar,'("deriv",I1,".dat")') n
ch="rm "//nchar
filexist=.false.
inquire(file=trim(nchar),exist=filexist)
if (filexist.eqv..true.) then ; call system(ch) ;end if

open(unit=1,file=trim(nchar),status='new')
  write(1,*) "Units are a.u."
  do k=1,3 !each energy component
    write(1,*)
    write(1,'("Energy component ",i1)') k
    write(1,'(I1," Derivative, error and location (RR iter, field) are:")') n
    write(1,*) save_deriv(k),save_er(k),savek(k),savej(k)+1
    write(1,*)
    write(1,*) "Triangle is"

    !conditional printing of the proper matrixes
    !Derivatives matrix
    print_size=s
    do i=1,3
      do j=1,s
        do l=1,s ;p(i,l,j)=d(i,j,l) ; end do
      end do
    end do
    write(1,"(/1X,A)") "Derivative matrix:"
    position = 2
    do i = 1, print_size
      write(format,"(A,I4.4,A)")  "(12ES15.6,T", Position,")"
      write(1,'(9ES15.6)')  (p(k,i,j), j = 1, print_size-i+1)
      position = position + 15 
    end do

    !Error matrix
    print_size=s-1
    do i=1,3
      do j=1,s-1
        do l=1,s-1 ;p(i,l,j)=e(i,j,l) ; end do
      end do
    end do
    write(1,"(/1X,A)") "Error matrix:"
    position = 2
    do i = 1, print_size
      write(format,"(A,I4.4,A)")  "(12ES15.6,T", Position,")"
      write(1,'(9ES15.6)')  (p(k,i,j), j = 1, print_size-i+1)
      position = position + 15
    end do

  end do
close(1)

end if

end

!Romberg-Rutishauser method for the numerical derivatives
!2nd derivatives matrix (2D, with crossed terms).
subroutine RR_2D(numF,e,d)
implicit none
integer :: numF,numpos !number of EF in the numerical procedure,number of positive fields
integer :: i,j,savei,savej
double precision, dimension (numF,13) :: e !grid of energies
double precision, dimension (3,3) :: d !2nd derivative (matrix)
double precision, dimension (((numF-1)/2)-1,((numF-1)/2)-1) :: r,er2
double precision, dimension ((numF-1)/2,(numF-1)/2) :: disp2
double precision, dimension ((numF-1)/2) :: ffxy,ffxz,ffyz
double precision :: h,err2!stepsize for each iteration

!CALCULATION OF THE 2ND DERIVATIVE MATRIX ONLY

!The expected entered energies matrix e is (numF,7) which means that for each field
!we have 7 possible energy values. The values are the following:
!
! origin(1), x(2), -x(3), y(4), -y(5), z(6),-z(7),xy(8),-xy(9),xz(10),-xz(11),yz(12),-yz(13)
!
!
!               Y axis  
!               |
!               x       x                      4       8
!               |
!               |
!-------x-------x-------x------Xaxis  3        1       2
!               |
!               |
!       x       x                     9        5
!               |
!

!
!               Z axis
!               |
!               x       x                      6       10
!               |
!               |
!-------x-------x-------x------Xaxis  3        1       2
!               |
!               |
!       x       x                     11       7
!               |
!

!
!               Z axis
!               |
!               x       x                      6       12
!               |
!               |
!-------x-------x-------x------Yaxis  5        1       4
!               |
!               |
!       x       x                     13       7
!               |
!

!initialization
numpos=(numF-1)/2
ffxy=1.0d8; ffxz=1.0d8 ;ffyz=1.0d8 ; r=1.0d8

!compute second derivative, FFxy
do i=2,numpos+1
   h=2**(i-2)*1.0d-4 !stepsize (value of the EF)
   ffxy(i-1)=-1/(2*h**2)*(e(i,2)+e(i+numpos,3)+e(i,4)+e(i+numpos,5)-2*e(1,1)-e(i,8)-e(i+numpos,9))
end do

!compute second derivative, FFxz
do i=2,numpos+1
   h=2**(i-2)*1.0d-4 !stepsize (value of the EF)
   ffxz(i-1)=-1/(2*h**2)*(e(i,2)+e(i+numpos,3)+e(i,6)+e(i+numpos,7)-2*e(1,1)-e(i,10)-e(i+numpos,11))
end do

!compute second derivative, FFyz
do i=2,numpos+1
   h=2**(i-2)*1.0d-4 !stepsize (value of the EF)
   ffyz(i-1)=-1/(2*h**2)*(e(i,4)+e(i+numpos,5)+e(i,6)+e(i+numpos,7)-2*e(1,1)-e(i,12)-e(i+numpos,13))
end do

!Compute the RR correction to FF_XY vector

!first RR iteration is done respect to the FF vector
j=1
do while (j<=(numpos-1))
  r(1,j)=(4*ffxy(j)-ffxy(j+1))/3.0d0
  j=j+1
end do

! from 2nd RR iteration onwards is is computed from the r tensor
do i=2,numpos-1 !each RR iteration
  j=1 ! used for counting
  do while (j<=(numpos-i))  !Only upper side triangle
    r(i,j)=(4**i*r(i-1,j)-r(i-1,j+1))/(4**i-1)
    j=j+1
 end do
end do

!save the triangles to a new array for displaying purposes (compact array with
!FF and RR matrixes)
disp2=1.0d8
do i=1,numpos !each RR iteration
  do j=1,numpos
    if (i==1) then
      disp2(i,j)=ffxy(j)
    else if (j<=(numpos-1)) then !the sizes of the tensors are different
      disp2(i,j)=r(i-1,j)
    end if
  end do
end do

!Compute error triangle from the display triangle
er2=1.0d8 ; err2=1.0d8
do i=1,numpos-1 !each iteration
  do j=1,numpos-1 !each field
    er2(i,j)=abs(disp2(i,j+1)-disp2(i,j))

   ! this is to avoid having zeros in the lower part of the triangle (to catch
   ! the minimal error this must be avoided)
    if (er2(i,j)==0.0d0) er2(i,j)=1.0d8

   !save the smallest error
    if (er2(i,j)<err2) then
      savei=i
      savej=j
      err2=er2(i,j)
    end if
  end do
end do

!save final results for second derivative FXY
d(1,2)=disp2(savei,savej+1)
if (abs(d(1,2))>1E-2.and.(err2/d(1,2)) > 5.0d-2) then !if the value is significant + has >5% error
   write(*,*) "WARNING IN 2D PRINTING"
   write(*,*) "The 2D derivative xy has an error higher than 5%"
   write(*,*) " % error         value of derivative      romberg error"
   write(*,*) err2/d(1,2)*100.0d0, "%",d(1,2),err2
else
   !write(*,*) "nowarn",d(1,2),err2
end if

!Compute the RR correction to FF_XZ vector
!first RR iteration is done respect to the FF_XZ vector
r=1.0d8
j=1
do while (j<=(numpos-1))
  r(1,j)=(4*ffxz(j)-ffxz(j+1))/3.0d0
  j=j+1
end do

! from 2nd RR iteration onwards is is computed from the r tensor
do i=2,numpos-1 !each RR iteration
  j=1 ! used for counting
  do while (j<=(numpos-i))  !Only upper side triangle
    r(i,j)=(4**i*r(i-1,j)-r(i-1,j+1))/(4**i-1)
    j=j+1
 end do
end do

!save the triangles to a new array for displaying purposes (compact array with
!FF and RR matrixes)
disp2=1.0d8
do i=1,numpos !each RR iteration
  do j=1,numpos
    if (i==1) then
      disp2(i,j)=ffxz(j)
    else if (j<=(numpos-1)) then !the sizes of the tensors are different
      disp2(i,j)=r(i-1,j)
    end if
  end do
end do

!Compute error triangle from the display triangle
er2=1.0d8 ; err2=1.0d8
do i=1,numpos-1 !each iteration
  do j=1,numpos-1 !each field
    er2(i,j)=abs(disp2(i,j+1)-disp2(i,j))

   ! this is to avoid having zeros in the lower part of the triangle (to catch
   ! the minimal error this must be avoided)
    if (er2(i,j)==0.0d0) er2(i,j)=1.0d8

   !save the smallest error
    if (er2(i,j)<err2) then
      savei=i
      savej=j
      err2=er2(i,j)
    end if
  end do
end do

!save final results for second derivative FXZ
d(1,3)=disp2(savei,savej+1)
if (abs(d(1,3))>1E-2.and.(err2/d(1,3)) > 5.0d-2) then !if the value is significant + has >5% error
   write(*,*) "WARNING IN 2D PRINTING"
   write(*,*) "The 2D derivative xz has an error higher than 5%"
   write(*,*) " % error         value of derivative      romberg error"
   write(*,*) err2/d(1,3)*100.0d0, "%",d(1,3),err2
else
   !write(*,*) "nowarn",d(1,3),err2
end if


!Compute the RR correction to FF_YZ vector
!first RR iteration is done respect to the FF vector
r=1.0d8
j=1
do while (j<=(numpos-1))
  r(1,j)=(4*ffyz(j)-ffyz(j+1))/3.0d0
  j=j+1
end do

! from 2nd RR iteration onwards is is computed from the r tensor
do i=2,numpos-1 !each RR iteration
  j=1 ! used for counting
  do while (j<=(numpos-i))  !Only upper side triangle
    r(i,j)=(4**i*r(i-1,j)-r(i-1,j+1))/(4**i-1)
    j=j+1
 end do
end do

!save the triangles to a new array for displaying purposes (compact array with
!FF and RR matrixes)
disp2=1.0d8
do i=1,numpos !each RR iteration
  do j=1,numpos
    if (i==1) then
      disp2(i,j)=ffyz(j)
    else if (j<=(numpos-1)) then !the sizes of the tensors are different
      disp2(i,j)=r(i-1,j)
    end if
  end do
end do

!Compute error triangle from the display triangle
er2=1.0d8 ; err2=1.0d8
do i=1,numpos-1 !each iteration
  do j=1,numpos-1 !each field
    er2(i,j)=abs(disp2(i,j+1)-disp2(i,j))

   ! this is to avoid having zeros in the lower part of the triangle (to catch
   ! the minimal error this must be avoided)
    if (er2(i,j)==0.0d0) er2(i,j)=1.0d8

   !save the smallest error
    if (er2(i,j)<err2) then
      savei=i
      savej=j
      err2=er2(i,j)
    end if
  end do
end do

!save final results for second derivative FYZ
d(2,3)=disp2(savei,savej+1)
if (abs(d(2,3))>1E-2.and.(err2/d(2,3)) > 5.0d-2) then !if the value is significant + has >5% error
   write(*,*) "WARNING IN 2D PRINTING"
   write(*,*) "The 2D derivative yz has an error higher than 5%"
   write(*,*) " % error         value of derivative      romberg error"
   write(*,*) err2/d(2,3)*100.0d0, "%",d(2,3),err2
else
   !write(*,*) "nowarn",d(2,3),err2
end if

!saving others by symmetry
!need to fill the xx, yy and zz derivs yet
d(2,1)=d(1,2)
d(3,1)=d(1,3)
d(3,2)=d(2,3)

end

subroutine locate(iunit,string)
character*80 :: var
integer :: ii,iunit
character string*(*)

do
   read(iunit,'(a80)',end=9,err=9) var
end do

9 ii=0

do while(ii==0)
   backspace iunit
   backspace iunit
   read(iunit,'(a80)') var
   if(index(var,string).ne.0) then
      ii=1
!      write(*,*) "Section found"
   end if
   if(index(var,"Entering Gaussian System").ne.0) then
      write(*,*) string, " section not found"
      return
   end if
end do
return
end


subroutine compute_2d(E0,p,s,f,kcal,grid)
implicit none
integer :: grid,i,j
double precision :: E0,E
double precision, dimension(3) :: p
double precision, dimension(3,3) :: s
double precision, dimension(grid) :: f
double precision :: kcal

!write(*,*) "Debugging mode: printing derivatives"
!write(*,*) "p1 p2 s(1,1) s(2,2) s(1,2)"
!write(*,'(7F11.6)') p(1),p(2),s(1,1),s(2,2),s(1,2)
open(unit=43,file="EXY.dat",status="new")
do i=1,grid !field along X
  do j=1,grid !field along Y
    E=E0+(p(1)*f(i)+p(2)*f(j))
    E=E+0.5d0*(s(1,1)*f(i)*f(i)+s(2,2)*f(j)*f(j)+2.0d0*s(1,2)*f(i)*f(j))
    E=E*kcal
    write(43,*) f(i),f(j),E
  end do
end do
close(43)

!write(*,*) "Debugging mode: printing derivatives"
!write(*,*) "p1, p3, s(1,1), s(3,3), s(1,3)"
!write(*,*) p(1),p(3),s(1,1),s(3,3),s(1,3)
open(unit=44,file="EXZ.dat",status="new")
do i=1,grid !field along X
  do j=1,grid !field along Z
    E=E0+(p(1)*f(i)+p(3)*f(j))
    E=E+0.5d0*(s(1,1)*f(i)*f(i)+s(3,3)*f(j)*f(j)+2.0d0*s(1,3)*f(i)*f(j))
    E=E*kcal
    write(44,*) f(i),f(j),E
  end do
end do
close(44)

!write(*,*) "Debugging mode: printing derivatives"
!write(*,*) "p1, p3, s(1,1), s(3,3), s(1,3)"
!write(*,*) p(1),p(3),s(1,1),s(3,3),s(1,3)
open(unit=45,file="EYZ.dat",status="new")
do i=1,grid !field along Y
  do j=1,grid !field along Z
    E=E0+(p(2)*f(i)+p(3)*f(j))
    E=E+0.5d0*(s(2,2)*f(i)*f(i)+s(3,3)*f(j)*f(j)+2.0d0*s(2,3)*f(i)*f(j))
    E=E*kcal
    write(45,*) f(i),f(j),E
  end do
end do
close(45)

end

subroutine print_summary_read(elstat2De,pauli2De,oi2De,numF)
implicit none
integer :: i,j,numF,numpos
double precision, dimension (numF,13) :: elstat2De,pauli2De,oi2De

write(*,*)
write(*,*) "Values read from summary.dat"

numpos=(numF-1)/2

do i=2,12,2
  write(*,*) 1,1,elstat2De(1,1)
  do j=2,numF
     if (j<=(numpos+1)) then
       write(*,*) j,i,elstat2De(j,i)
     else if (j>(numpos+1)) then
       write(*,*) j,i,elstat2De(j,i+1)
     end if
  end do

  write(*,*) 1,1,pauli2De(1,1)
  do j=2,numF
     if (j<=(numpos+1)) then
       write(*,*) j,i,pauli2De(j,i)
     else if (j>(numpos+1)) then
       write(*,*) j,i,pauli2De(j,i+1)
     end if
  end do

  write(*,*) 1,1,oi2De(1,1)
  do j=2,numF
     if (j<=(numpos+1)) then
       write(*,*) j,i,oi2De(j,i)
     else if (j>(numpos+1)) then
       write(*,*) j,i,oi2De(j,i+1)
     end if
  end do

end do

write(*,*)

end

subroutine printgrid(grid,f)
implicit none
integer :: i, grid
double precision, dimension(grid) :: f

open(unit=11,file="grid.plt",status="new")
do i=1,grid
 write(11,*) f(i)
end do
close(11)

end
