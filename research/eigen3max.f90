program eigenvalue

 implicit none
 integer,parameter :: n=20,m=20
 real a(2,2),x(2,2),eigen(2),maxeigen(0:n)
 real u,v,d,x0(2),Discri,apara(0:n),bpara(0:m)
 real,parameter :: eps=1e-14,ipu=0.03
 integer i,j,k,icount(0:n)


 open(11,file='eigen.txt')

 apara(0)=-10
 bpara(0)=-10

 do j=0,m
 
   write(*,*) "j=",j
 
    do i=0,n
 
       write(*,*) 'i=',i

       Discri=(apara(i)+bpara(j)-1)**2+4*apara(i)*(1+bpara(j)) !ルートの中身が負にならない条件

       if(Discri < 0) then
          write(*,*) 'jimeinakainomi'
          cycle
       endif

       x0(1)=(1-apara(i)-bpara(j)+((apara(i)+bpara(j)-1)**2+4*apara(i)*(1+bpara(j)))**0.5)/2
       x0(2)=(1-apara(i)-bpara(j)+((apara(i)+bpara(j)-1)**2+4*apara(i)*(1+bpara(j)))**0.5)/2

        a(1,1)=1/ipu*(1-2*x0(1)-2*apara(i)*bpara(j)*x0(2)/(x0(1)+apara(i))**2) !線形化行列の値の代入
        a(1,2)=1/ipu*(bpara(j)*(x0(1)-apara(i))/(x0(1)+apara(i)))
        a(2,1)=1
        a(2,2)=-1
 
        u=a(1,1)+a(2,2)
        v=a(1,1)*a(2,2)-a(1,2)*a(2,1)

        d=sqrt(u*u-4*v)
        eigen(1)=(u+d)/2
        eigen(2)=(u-d)/2

        do k=1,2
          if(abs(a(1,1)-eigen(k)) > eps .or. abs(a(1,2)) > eps) then
            x(1,k)=-a(1,2)
            x(2,k)=a(1,1)-eigen(k)
           else
            x(1,k)=a(2,2)-eigen(k)
            x(2,k)=-a(2,1)
           endif
            d=sqrt(x(1,k)*x(1,k)+x(2,k)*x(2,k))
            x(1,k)=x(1,k)/d
            x(2,k)=x(2,k)/d
         enddo
 
         call check_eigenvalue(a,eigen(1),x(1,1),2)
         call check_eigenvalue(a,eigen(2),x(1,2),2)
    
   maxeigen(i)=max(eigen(1),eigen(2))

   apara(i+1)=0.1+apara(i)
  icount(i+1)=i+1

    enddo

    write(11,'(2x,a,(100i12))') 'i=',icount(0:n)
    write(11,'(a,i3,(100e12.4))') "j=",j,maxeigen(0:n)

  bpara(j+1)=0.1+bpara(j)

 enddo

 close(11)

end program eigenvalue

subroutine check_eigenvalue(a,eigen,x,n)

 implicit none
 real a(n,n),eigen,x(n),xe,err
 integer n,i,j
 character(80) form

 err=0
 do i=1,n
    xe=0
    do j=1,n
       xe=xe+a(i,j)*x(j)
     enddo
 
     err=max(err,abs(eigen*x(i)-xe))
 enddo

 print "(' eigenvalue = ',f14.7,'  error =',es11.4)",eigen,err
 form=" (' eigenvector = (',f12.5,@(',',f12.5),')')"

 i=index(form,'@')

 form(i:i)=char(ichar('0')+n-1)
 print form,x(1),(x(1),i=2,n)
 
end subroutine check_eigenvalue


