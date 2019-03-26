program heikou_newton

 implicit none
 real x1(2),x2(2),x(2),y(2),dx(2),df(2,2),xmin(2),xmax(2),eps(2)
 integer ind,m,n,x1nmax,x2nmax
 external function

 xmin(1)=0.0    !区間幅と、分割数の設定
 xmax(1)=6.0
 xmin(2)=0.0
 xmax(2)=6.0
 x1nmax=2
 x2nmax=2

 dx(1)=(xmax(1)-xmin(1))/x1nmax   !各座標軸での刻み幅の計算
 dx(2)=(xmax(2)-xmin(2))/x2nmax

 x1=xmin   !ニュートン法の開始点の代入

 eps(1)=1; eps(2)=1   !収束判定値
 
 do m=1,x2nmax   !分割した小区間に対してニュートン法を適用
    x2(2)=xmin(2)+dx(2)*m
    do n=1,x1nmax
       x2(1)=xmin(1)+dx(1)*n
       call newton(function,x1,x2,x,eps,ind)
       if(ind >= 0) then   !この時更新値が収束半径内に収まる
         call function(x,y,df)
         print *,'(n,m)=',n,m
         print *,'X(1,2)=',x(1),x(2)
         print *,'Y(1,2)=',y(1),y(2),ind
       else if(ind <-1) then   !この時ニュートン法では発散、もしくは収束半径内に収まらない
       print*,'No convergence !'
       endif
       x1(1)=x2(1)   !次の小区間を考えるための操作
    enddo
    x1(2)=x2(2)
 enddo

end program heikou_newton


subroutine function(x,f,df)   !関数と導関数の設定

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
 integer,parameter :: itmax=100   !ニュートン法の反復回数

 xm=x10   !小区間の４隅の座標の計算
 x2=x20
 x2xm(1)=x2(1); x2xm(2)=xm(2)
 xmx2(1)=xm(1); xmx2(2)=x2(2)

 call subr(xm,ym,dym)   !４隅での関数と導関数の値の呼び出し
 call subr(x2,y2,dy2)
 call subr(x2xm,y2ym,dy2ym)
 call subr(xmx2,ymy2,dymy2)
ind=0   !ニュートン法を用いるか、二分法を用いるかのパラメータ

 if(all(ym == 0)) then   !この時のｘが解となる
   x=xm
   return
 else if(all(y2 == 0)) then   !この時のｘが解となる
   x=x2
   return
 else if(all(y2ym == 0)) then   !この時のｘが解となる
   x=x2xm
   return
 else if(all(ymy2 == 0)) then   !この時のｘが解となる
   x=xmx2
   return
 else if(all(ym*y2*y2ym*ymy2 > 0 )) then   !中間値の定理が適用できない場合
   ind=-1   !この時この区間に解が存在するとは限らない
   return
 endif

!4隅のうち、関数の絶対値が最小となる点に対しニュートン法を適用する
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

 xs=xm   !xmの値（間数値の絶対値が小さい区間の端点の向かい側の点）をxsに保存しておく
 x1=xm   !xmを反復法の開始点とする（実際にはその向かい側の点にニュートン法適用）


 do it=1,itmax
   det=dy2(1,1)*dy2(2,2)-dy2(2,1)*dy2(1,2)
   invdy2(1,1)=det*dy2(2,2)
   invdy2(1,2)=-det*dy2(1,2)
   invdy2(2,1)=-det*dy2(2,1)
   invdy2(2,2)=det*(1,1)
   xm=x2-matmul(invdy2,y2)   !ニュートン法の漸化式
   ib=0   !この時はニュートン法を適用
   if(xm(1) < min(x2(1),xs(1)) .or. xm(1) > max(x2(1),xs(1))&
      &.or. xm(2) < min(x2(2),xs(2)) .or. xm(2) > max(x2(2),xs(2))) then   !更新値が区間外に出るとき
     xm=(xs+x2)/2   !xmに中点を代入
     ib=1   !この時は二分法を適用
   endif

 xsx2xm=x2xm   !区間を更新する前に今の4隅の値を保存しておく
 ysx2xm=y2ym
 x2xm(1)=x2(1); x2xm(2)=xm(2)
 xsxmx2=xmx2
 ysxmx2=ymy2
 xmx2(1)=xm(1); xmx2(2)=x2(2)
 ys=ym

   call subr(xm,ym,dym) !更新値xmでの関数値の呼び出し
   call subr(x2xm,y2ym,dy2ym)
   call subr(xmx2,ymy2,dymy2)

   if(all(ym == 0)) then   !この時のｘが解となる
     x=xm
     exit
   endif

   if(all(y2ym == 0)) then   !この時のｘが解となる
     x=x2xm
     exit
   endif

   if(all(ymy2 == 0)) then   !この時のｘが解となる
     x=xmx2
   endif

   if(all(y2*ym*y2ym*ymy2 < 0 .or. y2*ym*y2ym < 0)) then
     xs=x2
     x1=x2   !更新値は後にx２とおく
   else if(ib > 0) then   !この時は二分法を適用
     x1=xs
   else if(all(ymy2*ysx2xm*y2ym*ym < 0 .or. ymy2*ysx2xm*y2ym < 0)) then   !以下では他の3つの区間に解がある場合を調べる
     xs=xsx2xm
     x1=xsx2xm   
   else if(all(y2ym*ym*ymy2*ysxmx2 < 0 .or. y2ym*ym*ymy2 < 0)) then
     xs=xsxmx2
     x1=xsxmx2
   else if(all(ys*ymy2*ym*y2ym < 0 .or. ys*ymy2*ym < 0)) then
     x1=x2
   endif

   x2=xm   !更新値の添え字を２とおく
   y2=ym
   dy2=dym

   if(all(abs(x2-x1) < eps)) then   !更新値と端点の幅が収束判定値内に収まるかの確認
     x=xm
     exit
   endif
 enddo

 if(it > itmax) then   !この時は繰り返し操作後に収束半径に収まらない
   ind=-itmax
 else   !収束半径に収まった時の反復回数をindに代入する
   ind=it
 endif

end subroutine newton
