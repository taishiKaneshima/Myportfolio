program kuukanitiyou

 implicit none
 integer,parameter :: n=1000
 real u(0:n),v(0:n)
 integer i
 real,parameter :: ipu=1,apara=2,bpara=20,dt=0.001
 


!‰Šú’l‚Ìİ’è
 u(0)=0.0001
 v(0)=0.0001

 open(14,file="kuukan.text")

 write(14,'(2f12.8)') u(0),v(0)

 do i=0,n-1

    u(i+1)=u(i)+dt*1/ipu*(u(i)*(1-u(i))-bpara*v(i)*(u(i)-apara)/(u(i)+apara))
   v(i+1)=v(i)+dt*(u(i)-v(i))
   
   write(14,'(2f12.8)') u(i+1),v(i+1)

 end do
 
 close(14)

stop

end program kuukanitiyou
