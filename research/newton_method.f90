program newton_method

 implicit none
 real x1,x2,x,eps,y,dy,xmin,xmax,dx
 external f
 integer ind,n.nmax

 xmin=0
 xmax=6
 nmax=10
 dx=(xmax-xmin)/nmax
 x1=xmin
 eps=1e-10

 do n=1,nmax
    x2=xmin+dex*n
    call newton(f,x1,x2,x,eps,ind)
    if(ind>-0) then
      call f
      print*,'x='x,y,ind
    else if(ind<)
