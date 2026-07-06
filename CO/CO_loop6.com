
BASIS=aug-cc-pVQZ

geometry={angstrom;
C;
O,C,r}

rgeom=[0.7,0.8,0.9,1.0,1.1,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.78,1.80] ang


field=[-0.005, 0.005]

do i=1,#rgeom

 r=rgeom(i)
 rx(i) = rgeom(i)
 XX(i) = 'xxxx'
 YY(i) = 'yyyy'

 {rhf}
 CCSD(T)
 ener1(i) = energy
 dip,,,field(1)

 {rhf}
 CCSD(T)
 ez1=energy
 dip,,,field(2)

 {rhf}
 CCSD(T)
 ez2=energy
 dipau=(ez2-ez1)/(field(2)-field(1))
 dm1(i)=dipau*TODEBYE

enddo


table, XX,rx,ener1
DIGITS, 0, 5,  12


table, YY,rx,dm1
DIGITS, 0, 5,  12

!table, rx,ener1
!DIGITS, 5,  12
!save,co_loop1.dat,new

---,  
  
