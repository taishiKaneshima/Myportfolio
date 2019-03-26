program ADI_test

 implicit none
 integer,parameter :: nm=99,nn=99,nl=1000       !��ԗ̈�Ǝ��ԕ��̕�����

 real,parameter :: Di=1.0,Di2=1.0,dx=1,dy=1,dt=0.1,ipu=1.0,aa=1,bb=1   !�g�U�W���A���ݕ��A���萔,tyson���f���̃p�����[�^

 real an(0:nm),bn(0:nm),cn(0:nm),an05(0:nm),bn05(0:nm),cn05(0:nm),&                  !�W���s��̐����iu�Ɋւ��āj
      &an2(0:nm),bn2(0:nm),cn2(0:nm),an205(0:nm),bn205(0:nm),cn205(0:nm),rx,ry,&     !�W���s��iv�Ɋւ��āj
      &d(0:nm,0:nn),un(0:nm,0:nn,0:nl),un05(0:nm,0:nn,0:nl),vn(0:nm,0:nn,0:nl),vn05(0:nm,0:nn,0:nl),&
      &non1(0:nm,0:nn),non2(0:nm,0:nn) 
!d�͒萔���x�N�g�������Au,v���Z�x�ϐ��Anon�͔���`�����An,n05�͎��ԃX�e�b�v��\��

 integer i,j,k

 rx=dt/2/(dx**2)      !�萔���܂Ƃ߂Ă������肳����
 ry=dt/2/(dy**2)

!���������̐ݒ�

!�Z�xu�̏��������̐ݒ�
 do i=0,nm
    do j=0,nn
       un(i,j,0)=20
    enddo
 enddo

!�Z�xv�̏��������̐ݒ�
 do i=0,nm
    do j=0,nn
       vn(i,j,0)=20
    enddo
 enddo

!���E�����̐ݒ�(�Z�xu�ɂ���)
 do i=0,nm        !ADI�X�L�[���ɂ��s�񐬕��i����n�̂Ƃ��j
    an(i)=-Di*rx
    bn(i)=1+2*Di*rx       !�T�u���[�`���̌v�Z���0���Ɣ��U����̂Œ���
    cn(i)=-Di*rx
 enddo

 bn(0)=0.5*bn(0)  !�Ď��m�C�}������
 bn(nm)=0.5*bn(nm)

 do i=0,nm        !ADI�X�L�[���̍s�񐬕��i����n+0.5�̂Ƃ��j
    an05(i)=-Di*ry
    bn05(i)=1+2*Di*ry
    cn05(i)=-Di*ry
 enddo

!���E�����̐ݒ�i�Z�x���ɂ��āj
 do i=0,nm        !ADI�X�L�[���̍s�񐬕��i����n�̂Ƃ��j
    an2(i)=-Di2*rx
    bn2(i)=1+2*Di2*rx
    cn2(i)=-Di2*rx
 enddo

 bn2(0)=0.5*bn2(0)  !�Ď��m�C�}������
 bn2(nm)=0.5*bn2(nm)

 do i=0,nm        !ADI�X�L�[���̍s�񐬕��i����n+0.5�̂Ƃ��j
    an205(i)=-Di2*ry
    bn205(i)=1+2*Di2*ry
    cn205(i)=-Di2*ry
 enddo

 bn205(0)=0.5*bn205(0)  !�Ď��m�C�}������
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

!n - >n+0.5�̎��ԃX�e�b�v��̔Z�x�̌v�Z

!����`���̌v�Z(tyson���f���̏ꍇ)�iun=aa���Ɣ��U����̂Œ��Ӂj
 do i=0,nm
    do j=0,nn
       non1(i,j)=ipu*dt*0.5*(un(i,j,k)*(1-un(i,j,k))-bb*vn(i,j,k)*(un(i,j,k)-aa)/(un(i,j,k)+aa))
       non2(i,j)=un(i,j,k)-vn(i,j,k)
    enddo
 enddo


!�Z�xu�̑Q����������
!�X�e�b�v1�i���Ԃ�n -> n+0.5�ɕω��������̔Z�x�ω��j
!�萔���x�N�g��d�̌v�Z
    do j=0,nn
       if(j==0 .or. j==nn) then    !���E�����ɂ���x�N�g���̗��[�͕ʌɈ���
          d(0,j)=0.5*(un(0,j,k)+Di*ry*(2*un(0,j,k)-2*un(0,j,k))+0.5*dt*ipu*non1(i,j))
          d(nm,j)=0.5*(un(nm,j,k)+Di*ry*(2*un(nm,j,k)-2*un(nm,j,k))+0.5*dt*ipu*non1(i,j))
          do i=1,nm-1
             d(i,j)=un(i,j,k)+Di*ry*(2*un(i,j+1,k)-2*un(i,j,k))+0.5*dt*ipu*non1(i,j)
          enddo
       else      !��x�N�g���̗��[�ȊO�̐����̂Ƃ��̌v�Z
          do i=0,nm
             d(i,j)=un(i,j,k)+Di*ry*(un(i,j+1,k)-2*un(i,j,k)+un(i,j-1,k))+0.5*dt*ipu*non1(i,j)
          enddo
       endif

       call tridiagonal_matrix(an,bn,cn,d,nm,un(:,j,k))  !�W���s��A�萔���x�N�g�������߂��̂ŕϐ���x�N�g�������߂�T�u���[�`�����Ăяo��
        
    enddo
 un05=un     !0.5�X�e�b�v���u�̔Z�x�����܂�


!�Z�xv�̑Q����������
!�X�e�b�v1�i���Ԃ�n -> n+0.5�ɕω��������̔Z�x�ω��j
!�萔���x�N�g��d�̌v�Z
    do j=0,nn
       if(j==0 .or. j==nn) then    !���E�����ɂ���x�N�g���̗��[�͕ʌɈ���
          d(0,j)=0.5*(vn(0,j,k)+Di2*ry*(2*vn(0,j,k)-2*vn(0,j,k))+0.5*dt*ipu*non2(i,j))
          d(nm,j)=0.5*(vn(nm,j,k)+Di2*ry*(2*vn(nm,j,k)-2*vn(nm,j,k))+0.5*dt*ipu*non2(i,j))
          do i=1,nm-1
             d(i,j)=vn(i,j,k)+Di2*ry*(2*vn(i,j+1,k)-2*vn(i,j,k))+0.5*dt*ipu*non2(i,j)
          enddo
       else      !��x�N�g���̗��[�ȊO�̐����̂Ƃ��̌v�Z
          do i=0,nm
             d(i,j)=vn(i,j,k)+Di2*ry*(vn(i,j+1,k)-2*vn(i,j,k)+vn(i,j-1,k))+0.5*dt*ipu*non2(i,j)
          enddo
       endif

       call tridiagonal_matrix(an2,bn2,cn2,d,nm,vn(:,j,k))  !�W���s��A�萔���x�N�g�������߂��̂ŕϐ���x�N�g�������߂�T�u���[�`�����Ăяo��
        
    enddo
 vn05=vn     !0.5�X�e�b�v���v�̔Z�x�����܂�





!0.5 -> 1�̎��ԃX�e�b�v��̔Z�x�̌v�Z

!����`���̌v�Z�iun=aa���Ɣ��U����̂Œ��Ӂj
 do i=0,nm
    do j=0,nn
       non1(i,j)=ipu*dt*0.5*(un05(i,j,k)*(1-un05(i,j,k))-bb*vn05(i,j,k)*(un05(i,j,k)-aa)/(un05(i,j,k)+aa))
       non2(i,j)=un05(i,j,k)-vn05(i,j,k)
    enddo
 enddo


!�Z�xu�̑Q����������
!�X�e�b�v2�in+0.5->n+1�ɕω��������̔Z�x�ω��j
!�萔���x�N�g��d�̌v�Z
    do i=0,nm
       if(i==0 .or. i==nm) then    !���E�����ɂ��s�x�N�g���̗��[�͕ʌɈ���
          d(i,0)=0.5*(un05(i,0,k)+Di*rx*(2*un05(i,0,k)-2*un05(i,0,k))+0.5*dt*ipu*non1(i,j))
          d(i,nn)=0.5*(un05(i,nn,k)+Di*rx*(2*un05(i,nn,k)-2*un05(i,nn,k))+0.5*dt*ipu*non1(i,j))
          do j=1,nn-1
             d(i,j)=un05(i,j,k)+Di*rx*(2*un05(i+1,j,k)-2*un05(i,j,k))+0.5*dt*ipu*non1(i,j)
          enddo
       else      !�s�x�N�g���̗��[�ȊO�̐����̂Ƃ��̌v�Z
          do j=0,nn
          d(i,j)=un05(i,j,k)+Di*rx*(un05(i+1,j,k)-2*un05(i,j,k)+un05(i-1,j,k))+0.5*dt*ipu*non1(i,j)
          enddo
       endif

       call tridiagonal_matrix(an05,bn05,cn05,d,nn,un05(i,:,k))  !�W���s��A�萔���x�N�g�������߂��̂ŕϐ��s�x�N�g�������߂�T�u���[�`�����Ăяo��
        
    enddo

 do i=0,nm
    do j=0,nn
       un(i,j,k+1)=un05(i,j,k)
    enddo
 enddo


!�Z�xv�̑Q����������
!�X�e�b�v2�in+0.5->n+1�ɕω��������̔Z�x�ω��j
!�萔���x�N�g��d�̌v�Z
    do i=0,nm
       if(i==0 .or. i==nm) then    !���E�����ɂ��s�x�N�g���̗��[�͕ʌɈ���
          d(i,0)=0.5*(vn05(i,0,k)+Di2*rx*(2*vn05(i,0,k)-2*vn05(i,0,k))+0.5*dt*ipu*non2(i,j))
          d(i,nn)=0.5*(vn05(i,nn,k)+Di2*rx*(2*vn05(i,nn,k)-2*vn05(i,nn,k))+0.5*dt*ipu*non2(i,j))
          do j=1,nn-1
             d(i,j)=vn05(i,j,k)+Di2*rx*(2*vn05(i+1,j,k)-2*vn05(i,j,k))+0.5*dt*ipu*non2(i,j)
          enddo
       else      !�s�x�N�g���̗��[�ȊO�̐����̂Ƃ��̌v�Z
          do j=0,nn
          d(i,j)=vn05(i,j,k)+Di2*rx*(vn05(i+1,j,k)-2*vn05(i,j,k)+vn05(i-1,j,k))+0.5*dt*ipu*non2(i,j)
          enddo
       endif

       call tridiagonal_matrix(an205,bn205,cn205,d,nn,vn05(i,:,k))  !�W���s��A�萔���x�N�g�������߂��̂ŕϐ��s�x�N�g�������߂�T�u���[�`�����Ăяo��
        
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









subroutine tridiagonal_matrix(a,b,c,d,n,x)  !���ږ@�i�K�E�X�̏����@�j�ŘA���������������v���O����

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

