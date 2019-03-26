program ADI_test

 implicit none
 integer,parameter :: nm=99,nn=99,nl=1000       !空間領域と時間幅の分割数

 real,parameter :: Di=1.0,Di2=1.0,dx=1,dy=1,dt=0.1,ipu=1.0,aa=1,bb=1   !拡散係数、刻み幅、時定数,tysonモデルのパラメータ

 real an(0:nm),bn(0:nm),cn(0:nm),an05(0:nm),bn05(0:nm),cn05(0:nm),&                  !係数行列の成分（uに関して）
      &an2(0:nm),bn2(0:nm),cn2(0:nm),an205(0:nm),bn205(0:nm),cn205(0:nm),rx,ry,&     !係数行列（vに関して）
      &d(0:nm,0:nn),un(0:nm,0:nn,0:nl),un05(0:nm,0:nn,0:nl),vn(0:nm,0:nn,0:nl),vn05(0:nm,0:nn,0:nl),&
      &non1(0:nm,0:nn),non2(0:nm,0:nn) 
!dは定数項ベクトル成分、u,vが濃度変数、nonは非線形成分、n,n05は時間ステップを表す

 integer i,j,k

 rx=dt/2/(dx**2)      !定数をまとめてすっきりさせる
 ry=dt/2/(dy**2)

!初期条件の設定

!濃度uの初期条件の設定
 do i=0,nm
    do j=0,nn
       un(i,j,0)=20
    enddo
 enddo

!濃度vの初期条件の設定
 do i=0,nm
    do j=0,nn
       vn(i,j,0)=20
    enddo
 enddo

!境界条件の設定(濃度uについて)
 do i=0,nm        !ADIスキームによる行列成分（時刻nのとき）
    an(i)=-Di*rx
    bn(i)=1+2*Di*rx       !サブルーチンの計算より0だと発散するので注意
    cn(i)=-Di*rx
 enddo

 bn(0)=0.5*bn(0)  !斉次ノイマン条件
 bn(nm)=0.5*bn(nm)

 do i=0,nm        !ADIスキームの行列成分（時刻n+0.5のとき）
    an05(i)=-Di*ry
    bn05(i)=1+2*Di*ry
    cn05(i)=-Di*ry
 enddo

!境界条件の設定（濃度ｖについて）
 do i=0,nm        !ADIスキームの行列成分（時刻nのとき）
    an2(i)=-Di2*rx
    bn2(i)=1+2*Di2*rx
    cn2(i)=-Di2*rx
 enddo

 bn2(0)=0.5*bn2(0)  !斉次ノイマン条件
 bn2(nm)=0.5*bn2(nm)

 do i=0,nm        !ADIスキームの行列成分（時刻n+0.5のとき）
    an205(i)=-Di2*ry
    bn205(i)=1+2*Di2*ry
    cn205(i)=-Di2*ry
 enddo

 bn205(0)=0.5*bn205(0)  !斉次ノイマン条件
 bn205(nm)=0.5*bn205(nm)

 write(11,*) "t=",0
 write(11,*) "u(0,0,0)"

 do i=0,nm
    write(12,'(100e12.4)') un(i,0:nn,0)
  enddo

write(13,*) "v(0,0,0)"

 do i=0,nm
    write(14,'(100e12.4)') vn(i,0:nn,0)
  enddo


 do k=0,nl

!n - >n+0.5の時間ステップ後の濃度の計算

!非線形項の計算(tysonモデルの場合)（un=aaだと発散するので注意）
 do i=0,nm
    do j=0,nn
       non1(i,j)=ipu*dt*0.5*(un(i,j,k)*(1-un(i,j,k))-bb*vn(i,j,k)*(un(i,j,k)-aa)/(un(i,j,k)+aa))
       non2(i,j)=un(i,j,k)-vn(i,j,k)
    enddo
 enddo


!濃度uの漸化式を解く
!ステップ1（時間がn -> n+0.5に変化した時の濃度変化）
!定数項ベクトルdの計算
    do j=0,nn
       if(j==0 .or. j==nn) then    !境界条件により列ベクトルの両端は別個に扱う
          d(0,j)=0.5*(un(0,j,k)+Di*ry*(2*un(0,j,k)-2*un(0,j,k))+0.5*dt*ipu*non1(i,j))
          d(nm,j)=0.5*(un(nm,j,k)+Di*ry*(2*un(nm,j,k)-2*un(nm,j,k))+0.5*dt*ipu*non1(i,j))
          do i=1,nm-1
             d(i,j)=un(i,j,k)+Di*ry*(2*un(i,j+1,k)-2*un(i,j,k))+0.5*dt*ipu*non1(i,j)
          enddo
       else      !列ベクトルの両端以外の成分のときの計算
          do i=0,nm
             d(i,j)=un(i,j,k)+Di*ry*(un(i,j+1,k)-2*un(i,j,k)+un(i,j-1,k))+0.5*dt*ipu*non1(i,j)
          enddo
       endif

       call tridiagonal_matrix(an,bn,cn,d,nm,un(:,j,k))  !係数行列、定数項ベクトルを求めたので変数列ベクトルを求めるサブルーチンを呼び出す
        
    enddo
 un05=un     !0.5ステップ後のuの濃度が求まる


!濃度vの漸化式を解く
!ステップ1（時間がn -> n+0.5に変化した時の濃度変化）
!定数項ベクトルdの計算
    do j=0,nn
       if(j==0 .or. j==nn) then    !境界条件により列ベクトルの両端は別個に扱う
          d(0,j)=0.5*(vn(0,j,k)+Di2*ry*(2*vn(0,j,k)-2*vn(0,j,k))+0.5*dt*ipu*non2(i,j))
          d(nm,j)=0.5*(vn(nm,j,k)+Di2*ry*(2*vn(nm,j,k)-2*vn(nm,j,k))+0.5*dt*ipu*non2(i,j))
          do i=1,nm-1
             d(i,j)=vn(i,j,k)+Di2*ry*(2*vn(i,j+1,k)-2*vn(i,j,k))+0.5*dt*ipu*non2(i,j)
          enddo
       else      !列ベクトルの両端以外の成分のときの計算
          do i=0,nm
             d(i,j)=vn(i,j,k)+Di2*ry*(vn(i,j+1,k)-2*vn(i,j,k)+vn(i,j-1,k))+0.5*dt*ipu*non2(i,j)
          enddo
       endif

       call tridiagonal_matrix(an2,bn2,cn2,d,nm,vn(:,j,k))  !係数行列、定数項ベクトルを求めたので変数列ベクトルを求めるサブルーチンを呼び出す
        
    enddo
 vn05=vn     !0.5ステップ後のvの濃度が求まる





!0.5 -> 1の時間ステップ後の濃度の計算

!非線形項の計算（un=aaだと発散するので注意）
 do i=0,nm
    do j=0,nn
       non1(i,j)=ipu*dt*0.5*(un05(i,j,k)*(1-un05(i,j,k))-bb*vn05(i,j,k)*(un05(i,j,k)-aa)/(un05(i,j,k)+aa))
       non2(i,j)=un05(i,j,k)-vn05(i,j,k)
    enddo
 enddo


!濃度uの漸化式を解く
!ステップ2（n+0.5->n+1に変化した時の濃度変化）
!定数項ベクトルdの計算
    do i=0,nm
       if(i==0 .or. i==nm) then    !境界条件により行ベクトルの両端は別個に扱う
          d(i,0)=0.5*(un05(i,0,k)+Di*rx*(2*un05(i,0,k)-2*un05(i,0,k))+0.5*dt*ipu*non1(i,j))
          d(i,nn)=0.5*(un05(i,nn,k)+Di*rx*(2*un05(i,nn,k)-2*un05(i,nn,k))+0.5*dt*ipu*non1(i,j))
          do j=1,nn-1
             d(i,j)=un05(i,j,k)+Di*rx*(2*un05(i+1,j,k)-2*un05(i,j,k))+0.5*dt*ipu*non1(i,j)
          enddo
       else      !行ベクトルの両端以外の成分のときの計算
          do j=0,nn
          d(i,j)=un05(i,j,k)+Di*rx*(un05(i+1,j,k)-2*un05(i,j,k)+un05(i-1,j,k))+0.5*dt*ipu*non1(i,j)
          enddo
       endif

       call tridiagonal_matrix(an05,bn05,cn05,d,nn,un05(i,:,k))  !係数行列、定数項ベクトルを求めたので変数行ベクトルを求めるサブルーチンを呼び出す
        
    enddo

 do i=0,nm
    do j=0,nn
       un(i,j,k+1)=un05(i,j,k)
    enddo
 enddo


!濃度vの漸化式を解く
!ステップ2（n+0.5->n+1に変化した時の濃度変化）
!定数項ベクトルdの計算
    do i=0,nm
       if(i==0 .or. i==nm) then    !境界条件により行ベクトルの両端は別個に扱う
          d(i,0)=0.5*(vn05(i,0,k)+Di2*rx*(2*vn05(i,0,k)-2*vn05(i,0,k))+0.5*dt*ipu*non2(i,j))
          d(i,nn)=0.5*(vn05(i,nn,k)+Di2*rx*(2*vn05(i,nn,k)-2*vn05(i,nn,k))+0.5*dt*ipu*non2(i,j))
          do j=1,nn-1
             d(i,j)=vn05(i,j,k)+Di2*rx*(2*vn05(i+1,j,k)-2*vn05(i,j,k))+0.5*dt*ipu*non2(i,j)
          enddo
       else      !行ベクトルの両端以外の成分のときの計算
          do j=0,nn
          d(i,j)=vn05(i,j,k)+Di2*rx*(vn05(i+1,j,k)-2*vn05(i,j,k)+vn05(i-1,j,k))+0.5*dt*ipu*non2(i,j)
          enddo
       endif

       call tridiagonal_matrix(an205,bn205,cn205,d,nn,vn05(i,:,k))  !係数行列、定数項ベクトルを求めたので変数行ベクトルを求めるサブルーチンを呼び出す
        
    enddo
 
    do i=0,nm
       do j=0,nn
       vn(i,j,k+1)=vn05(i,j,k)
       enddo
    enddo

    write(11,*) "t=",k+1

    write(11,*) "u"

    write(12,'(//)')
    write(14,'(//)')

    do i=0,nm
        write(12,'(100e12.4)') un(i,0:nn,k)
     enddo

    write(11,*) "v"
    do i=0,nm
        write(14,'(100e12.4)') vn(i,0:nn,k)
    enddo

 enddo

stop
end program ADI_test









subroutine tridiagonal_matrix(a,b,c,d,n,x)  !直接法（ガウスの消去法）で連立方程式を解くプログラム

 implicit none
 real a(0:n),b(0:n),c(0:n),d(0:n),x(0:n)
 real G(0:n),H(0:n),den
 integer i,n
 G(0) = -c(0)/b(0)
 H(0) = d(0)/b(0)
 
 do i = 1,n
    den = 1/(b(i)+a(i)*G(i-1))
    G(i) = -c(i)*den
    H(i) = (d(i)-a(i)*H(i-1))*den
 enddo

 x(n) = H(n)

 do i = n-1,0,-1
    x(i) = G(i)*x(i+1)+H(i)
 enddo

end subroutine tridiagonal_matrix

