program tridigonal_matrix_test

 implicit none
 integer,parameter :: nm=30,nn=30
 real,parameter :: Di=1.0,rx=1.0,ry=1.0,deltat=1.0,ipu=1.0
 real a(nm),b(nm),c(nm),d(nm,nn),un05(nm,nn),un(nm,nn)
 integer i,j

 do i=1,nm
    a(i)=-Di*rx
    b(i)=1+2*Di*rx
    c(i)=-Di*rx
 enddo

 
 b(1)=0.5*b(1)
 b(nm)=0.5*b(nm)


 do i=1,nm
    do j=1,nn
       un(i,j)=0
    enddo
 enddo

do j=1,nn


 if(j=1 .or. j=nn) then

    d(1,j)=0.5*(un(1,j)+Di*ry*(2*un(1,j)-2*un(1,j))+0.5*deltat*ipu*1)
    d(nm,j)=0.5*(un(nm,j)+Di*ry*(2*un(nm,j)-2*un(nm,j))+0.5*deltat*ipu*1)

     do i=2,nm-1
        d(i,j)=un(i,j)+Di*ry*(2*un(i,j+1)-2*un(i,j))+0.5*deltat*ipu*1
     enddo

 else
    do i=1,nm
       d(i,j)=un(i,j)+Di*ry*(un(i,j+1)-2*un(i,j)+un(i,j-1))+0.5*deltat*ipu*1
    enddo

 endif

 call tridiagonal_matrix(a,b,c,d,nm,un05(:,j))

 do i = 1,nm
    print*,i,j,un05(i,j)
 enddo

 enddo

stop
end program tridigonal_matrix_test

subroutine tridiagonal_matrix(a,b,c,d,n,x)

 implicit none
 real a(n),b(n),c(n),d(n),x(n)
 real G(n),H(n),den
 integer i,n
 G(1) = -c(1)/b(1)
 H(1) = d(1)/b(1)
 
 do i = 2,n
    den = 1/(b(i)+a(i)*G(i-1))
    G(i) = -c(i)*den
    H(i) = (d(i)-a(i)*H(i-1))*den
 enddo

 x(n) = H(n)

 do i = n-1,1,-1
    x(i) = G(i)*x(i+1)+H(i)
 enddo

end subroutine tridiagonal_matrix

