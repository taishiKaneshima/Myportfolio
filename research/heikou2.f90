program heikou_newton

 implicit none
 real x1(2),x2(2),x(2),y(2),dx(2),df(2,2),xmin(2),xmax(2),eps(2)
 integer ind,m,n,x1nmax,x2nmax
 external function

 xmin(1)=0.0    !��ԕ��ƁA�������̐ݒ�
 xmax(1)=6.0
 xmin(2)=0.0
 xmax(2)=6.0
 x1nmax=2
 x2nmax=2

 dx(1)=(xmax(1)-xmin(1))/x1nmax   !�e���W���ł̍��ݕ��̌v�Z
 dx(2)=(xmax(2)-xmin(2))/x2nmax

 x1=xmin   !�j���[�g���@�̊J�n�_�̑��

 eps(1)=1; eps(2)=1   !��������l
 
 do m=1,x2nmax   !������������Ԃɑ΂��ăj���[�g���@��K�p
    x2(2)=xmin(2)+dx(2)*m
    do n=1,x1nmax
       x2(1)=xmin(1)+dx(1)*n
       call newton(function,x1,x2,x,eps,ind)
       if(ind >= 0) then   !���̎��X�V�l���������a���Ɏ��܂�
         call function(x,y,df)
         print *,'(n,m)=',n,m
         print *,'X(1,2)=',x(1),x(2)
         print *,'Y(1,2)=',y(1),y(2),ind
       else if(ind <-1) then   !���̎��j���[�g���@�ł͔��U�A�������͎������a���Ɏ��܂�Ȃ�
       print*,'No convergence !'
       endif
       x1(1)=x2(1)   !���̏���Ԃ��l���邽�߂̑���
    enddo
    x1(2)=x2(2)
 enddo

end program heikou_newton


subroutine function(x,f,df)   !�֐��Ɠ��֐��̐ݒ�

 real x(2),f(2),df(2,2)
 real,parameter :: ipu=1,a=1,b=1
 
 f(1)=1/ipu*(x(1)*(1-x(1))-b*x(2)*(x(1)-a)/(x(1)+a))
 f(2)=x(1)-x(2)
 df(1,1)=1/ipu*((1-x(1))-x(1)-2*a*b*x(2)/(x(1)+a)**2)
 df(1,2)=1/ipu*(1-2*x(1)-2*a*b*x(2)/(x(1)+a)**2)
 df(2,1)=1
 df(2,2)=-1

end subroutine function


subroutine newton(subr,x10,x20,x,eps,ind)

 implicit none
 real x10(2),x20(2),x2xm(2),xmx2(2),x(2),x1(2),x2(2),y2(2),dy2(2,2),&
&xm(2),ym(2),dym(2,2),y2ym(2),dy2ym(2,2),ymy2(2),dymy2(2,2),&
&xs(2),ys(2),xsx2(2),xsxm(2),xsx2xm(2),xsxmx2(2),ysx2xm(2),ysxmx2(2),eps(2),det,invdy2(2,2)
 external subr
 integer ind,it,ib
 integer,parameter :: itmax=100   !�j���[�g���@�̔�����

 xm=x10   !����Ԃ̂S���̍��W�̌v�Z
 x2=x20
 x2xm(1)=x2(1); x2xm(2)=xm(2)
 xmx2(1)=xm(1); xmx2(2)=x2(2)

 call subr(xm,ym,dym)   !�S���ł̊֐��Ɠ��֐��̒l�̌Ăяo��
 call subr(x2,y2,dy2)
 call subr(x2xm,y2ym,dy2ym)
 call subr(xmx2,ymy2,dymy2)
ind=0   !�j���[�g���@��p���邩�A�񕪖@��p���邩�̃p�����[�^

 if(all(ym == 0)) then   !���̎��̂������ƂȂ�
   x=xm
   return
 else if(all(y2 == 0)) then   !���̎��̂������ƂȂ�
   x=x2
   return
 else if(all(y2ym == 0)) then   !���̎��̂������ƂȂ�
   x=x2xm
   return
 else if(all(ymy2 == 0)) then   !���̎��̂������ƂȂ�
   x=xmx2
   return
 else if(all(ym*y2*y2ym*ymy2 > 0 )) then   !���Ԓl�̒藝���K�p�ł��Ȃ��ꍇ
   ind=-1   !���̎����̋�Ԃɉ������݂���Ƃ͌���Ȃ�
   return
 endif

!4���̂����A�֐��̐�Βl���ŏ��ƂȂ�_�ɑ΂��j���[�g���@��K�p����
 if(max(abs(ym(1))*abs(ym(2)),abs(y2(1))*abs(y2(2)),abs(y2ym(1))*abs(y2ym(2)),&
&abs(ymy2(1)*ymy2(2))) == abs(y2(1))*abs(y2(2))) then
   xs=xm
   xm=x2
   x2=xs
   y2=ym
   dy2=dym
   xsx2xm=x2xm
   x2xm=xmx2
   xmx2=xsx2xm
   else if(max(abs(ym(1))*abs(ym(2)),abs(y2(1))*abs(y2(2)),abs(y2ym(1))*abs(y2ym(2)),&
&abs(ymy2(1)*ymy2(2))) == abs(y2ym(1))*abs(y2ym(2))) then
          xs=x2xm
          xsxm=xm
          xm=x2xm
          xsx2=x2
          x2=xmx2
          y2=ymy2
          dy2=dymy2
          xmx2=xsx2
          x2xm=xsxm
   else if(max(abs(ym(1))*abs(ym(2)),abs(y2(1))*abs(y2(2)),abs(y2ym(1))*abs(y2ym(2)),&
&abs(ymy2(1)*ymy2(2))) == abs(ymy2(1))*abs(ymy2(2))) then
          xs=xmx2
          xsxm=xm
          xm=xmx2
          xsx2=x2
          x2=x2xm
          y2=y2ym
          x2xm=xsx2
          xmx2=xsxm
          dy2=dy2ym
   endif

 xs=xm   !xm�̒l�i�Ԑ��l�̐�Βl����������Ԃ̒[�_�̌��������̓_�j��xs�ɕۑ����Ă���
 x1=xm   !xm�𔽕��@�̊J�n�_�Ƃ���i���ۂɂ͂��̌��������̓_�Ƀj���[�g���@�K�p�j


 do it=1,itmax
   det=dy2(1,1)*dy2(2,2)-dy2(2,1)*dy2(1,2)
   invdy2(1,1)=det*dy2(2,2)
   invdy2(1,2)=-det*dy2(1,2)
   invdy2(2,1)=-det*dy2(2,1)
   invdy2(2,2)=det*(1,1)
   xm=x2-matmul(invdy2,y2)   !�j���[�g���@�̑Q����
   ib=0   !���̎��̓j���[�g���@��K�p
   if(xm(1) < min(x2(1),xs(1)) .or. xm(1) > max(x2(1),xs(1))&
      &.or. xm(2) < min(x2(2),xs(2)) .or. xm(2) > max(x2(2),xs(2))) then   !�X�V�l����ԊO�ɏo��Ƃ�
     xm=(xs+x2)/2   !xm�ɒ��_����
     ib=1   !���̎��͓񕪖@��K�p
   endif

 xsx2xm=x2xm   !��Ԃ��X�V����O�ɍ���4���̒l��ۑ����Ă���
 ysx2xm=y2ym
 x2xm(1)=x2(1); x2xm(2)=xm(2)
 xsxmx2=xmx2
 ysxmx2=ymy2
 xmx2(1)=xm(1); xmx2(2)=x2(2)
 ys=ym

   call subr(xm,ym,dym) !�X�V�lxm�ł̊֐��l�̌Ăяo��
   call subr(x2xm,y2ym,dy2ym)
   call subr(xmx2,ymy2,dymy2)

   if(all(ym == 0)) then   !���̎��̂������ƂȂ�
     x=xm
     exit
   endif

   if(all(y2ym == 0)) then   !���̎��̂������ƂȂ�
     x=x2xm
     exit
   endif

   if(all(ymy2 == 0)) then   !���̎��̂������ƂȂ�
     x=xmx2
   endif

   if(all(y2*ym*y2ym*ymy2 < 0 .or. y2*ym*y2ym < 0)) then
     xs=x2
     x1=x2   !�X�V�l�͌��x�Q�Ƃ���
   else if(ib > 0) then   !���̎��͓񕪖@��K�p
     x1=xs
   else if(all(ymy2*ysx2xm*y2ym*ym < 0 .or. ymy2*ysx2xm*y2ym < 0)) then   !�ȉ��ł͑���3�̋�Ԃɉ�������ꍇ�𒲂ׂ�
     xs=xsx2xm
     x1=xsx2xm   
   else if(all(y2ym*ym*ymy2*ysxmx2 < 0 .or. y2ym*ym*ymy2 < 0)) then
     xs=xsxmx2
     x1=xsxmx2
   else if(all(ys*ymy2*ym*y2ym < 0 .or. ys*ymy2*ym < 0)) then
     x1=x2
   endif

   x2=xm   !�X�V�l�̓Y�������Q�Ƃ���
   y2=ym
   dy2=dym

   if(all(abs(x2-x1) < eps)) then   !�X�V�l�ƒ[�_�̕�����������l���Ɏ��܂邩�̊m�F
     x=xm
     exit
   endif
 enddo

 if(it > itmax) then   !���̎��͌J��Ԃ������Ɏ������a�Ɏ��܂�Ȃ�
   ind=-itmax
 else   !�������a�Ɏ��܂������̔����񐔂�ind�ɑ������
   ind=it
 endif

end subroutine newton
