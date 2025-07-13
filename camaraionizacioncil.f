program camaraionizacion
use randomnumber
implicit none
real*8::h,PI,H_elec,rhomed,rho,Vol,hmax,hmin,haux,dt,alfa,mu_pos,mu_neg,k,t,V,recol_neg,recol_pos,a,b,error,er,fsum,fcuadsum,fmed,f
real*8::x
real*8,dimension(:),allocatable::rholocal,r_neg,r_pos,phi_neg,phi_pos,z_neg,z_pos,p
integer::N0,Npos,Nneg,i,j,m,i3,i4,m1

!se calcula el volumen interior del caso cilindrico conociendo la longtiud del mismo, así como los radios de los electrodos exterior(a) e interior(b)
H_elec=0.207

write(6,*)'ELECCIÓN DE RADIO DE ELECTRODOS. Elija un valor de los radio mayor a y menor b. Radio a (entre 0 y 0.5cm):'
read(5,*) a
write(6,*)' Radio b (entre 0 y a):'
read(5,*) b
!a=0.305
!b=0.055

PI=3.141592
Vol=PI*(a**2-b**2)*H_elec

!se puede calcular  la densidad media de iones positivos

write(6,*) 'ELECCIÓN DE DENSIDAD DE IONES.'
write(6,*)'Elija un valor de x entre 1 y 10 siendo la densidad inicial de iones 2.183E8/(xE4) cm**(-3):'
read(5,*) x

!Calculamos la densidad de iones sabiendo que una radiación de 1 mGy genera 2.183E8/cm³ y dividirla por un cierto factor que luego se añadirá a la constante de recombinación
rho=2.183E8/(x*10**4)
!Tomamos como constante de recombinación k =1.6E-6 cm³/s  y multiplicamos por el mismo factor que habiamos dividido la densidad
k=(1.6E-6)*(x*10**4)
N0=int(rho*Vol)
Nneg=N0
Npos=NNeg


write(6,*) 'ELECCIÓN DE VOLTAJE. Elija un valor de V entre 50 y 400 V:'
read(5,*) V
!V=400
mu_pos=1.80
mu_neg=2.10
alfa=100

allocate(rholocal(1:Nneg))
allocate(z_neg(1:Nneg))
allocate(z_pos(1:Npos))
allocate(r_neg(1:Nneg))
allocate(r_pos(1:Npos))
allocate(phi_neg(1:Nneg))
allocate(phi_pos(1:Npos))
allocate(p(1:Nneg))

error=1
er=1
m1=1
fsum=0.0
fcuadsum=0.0
!do while(er.gt.(0.0001))
	!inicializamos de nuevo el numero de iones(que antes ha disminuido).
	Nneg=N0
	Npos=NNeg
	t=0.0
	recol_pos=0
	recol_neg=0
	
	OPEN(unit=1,file='posicionesneg',status='old')
	OPEN(unit=2,file='posicionespos',status='old')
	call dran_ini(m1+52*m1)
	!elijamos una distribucion inicial de iones.
	call distribucionini
	!hallamos un valor apropiado de h para las densidades locales.
	call calcularh
	!calculamos el intervalo temporal y elegimos el campo máximo, el cual será para el radio mínimo->b
	dt=h/(mu_neg*E(b)*alfa)
	!do while (rhomed.gt.0)
	do while (Npos.gt.1)

		!ahora vamos a escribir la posicion de cada ion negativo y positivo con objetivo de hacer gift
		do i=1,N0
			if(i.le.Nneg)then
			write(1,*)r_neg(i)*cos(phi_neg(i)),r_neg(i)*sin(phi_neg(i)),z_neg(i)
			else 
			write(1,*) b/1000.0,10,10
			end if
		end do
		write(1,*)' '
		write(1,*)' '
		do i=1,N0
			if(i.le.Npos)then
			write(2,*)r_pos(i)*cos(phi_pos(i)),r_pos(i)*sin(phi_pos(i)),z_pos(i)
			else 
			write(2,*) b/1000.0,10,10
			end if
		end do
		write(2,*)' '
		write(2,*)' '
		!evolucionamos el sistema
		call pasotemporal
		!write(6,*)'paso temporal'
	end do
	!Acabada una simulación calculamos eficiencia
	f=100.0*(recol_neg+Nneg)/(1.0*N0)
	fsum=fsum+f
	fcuadsum=fcuadsum+f**2
	fmed=fsum/(1.0*m1)
	if(m1.gt.1)then
		error=sqrt(1/(m1*1.0*(m1-1))*(fcuadsum-fsum**2/(1.0*m1)))
		er=error/fmed
	end if
	write(6,*)'La eficiencia de recolección es del ',fmed,' %.'
	!write(6,*)m1,' ',fmed,' ',error,' ',er
	m1=m1+1
!end do
	



stop
contains



!función para calcular la distancia entre un ión negativo k y un ion positivo l
real*8 function dist(o,l)
	integer::o,l
	real*8::auxx,auxy,auxz

auxx=r_pos(l)*cos(phi_pos(l))-r_neg(o)*cos(phi_neg(o))
auxy=r_pos(l)*sin(phi_pos(l))-r_neg(o)*sin(phi_neg(o))
auxz=z_pos(l)-z_neg(o)
dist=sqrt(auxx**2+auxy**2+auxz**2)
end function
			
real*8 function E(r1)
	real*8::r1
	E=V/(r1*log(a/b))
end function
	
	
	



subroutine calcularh
	real*8::aux
	!elegimos un valor inicial del error(aux) mayor que 0.01 para que se produzca el do while
	aux=1
	!implementamos el algoritmo de biseccion estimando un intervalo que puede contener a h, por ejemplo entre el origen la  semialtura del cilindro
	hmax=H_elec/2.0
	hmin=0.0
	!tenemos que realizar el bucle las neces que sea necesario hasta alcanzar un valor apropiado de h
	do while(aux.gt.(0.01))
		!do m=1,3
		h=(hmax+hmin)/2.0
		call densidades
		!veamos el error relativo de pholocal para un h dado
		aux=abs(rhomed-rho)/rho
		
		if(rhomed.gt.rho) then
			!hmax=h
			hmin=h
		else 
			!hmin=h
			hmax=h
		end if
	end do
	!end do		
end subroutine




!tmb necesirtaremos una subroutina que cree una distribución uniforme a partir de una densidad dada. Se ha considerado un espacio cilíndrico homogeneo
subroutine distribucionini
	do m=1,Nneg
		
		r_neg(m)=b+abs((a-b)*sqrt(dran_u()))
		r_pos=r_neg
		phi_neg(m)=abs(dran_u()*2*PI)
		phi_pos(m)=phi_neg(m)
		z_pos(m)=abs(H_elec*dran_u())
		z_neg(m)=z_pos(m)
		
		
	end do
end subroutine

!Esta subrutina permite calcular las densidades(positivas) locales asociadas a cada ion negativo, así como la densidad media usando la aportación de cada densidad local.
subroutine densidades
	rholocal=0
		rhomed=0

			
		
		!primero debemos hallar la densidad local positiva de cada ion  negativo para poder realizar el promedio, para ello sumaremos la contribucion de los iones positivos j que rodean a un ion negativo i

		do i=1,Nneg
			do j=1,Npos
				if(dist(i,j).lt.h)then
					rholocal(i)=rholocal(i)+3.0/(4*PI*h**3) 
					
				end if
			
			end do	
			rhomed=rhomed+rholocal(i)/(Nneg*1.0)
		end do
end subroutine



subroutine pasotemporal
	integer::l
	real*8::aux
	!en primer lugar elegimos el valor del salto temporal dt.Para ello emplearemos la mobilidad iónica negtiva y el h hallado anteriormente. El dt será el tiempo que tarda una partícula en moverse 			una distancia h, estando sometida la particula a un campo electrico uniforme E=V/d en el caso plano-paralelo.
	
	!write(6,*)'El dt es: ',dt
	call densidades
	!write(6,*)rhomed, ' es la densidad promedio de las locales'
	!en segundo lugar tenemos en cuenta la recombinación.
	l=1
	do while (l.le.Nneg)
		!calculamos la probabilidad de recombinación
		!p(l)=k*rholocal(l)*dt
		p(l)=1-exp(-k*rholocal(l)*dt)
		!write(6,*)'La probabilidad de recombinación del ion negativo ',l,' es',p(l)
		aux=dran_u()
		!eliminamos el ion negativo si la probabilidad es mayor que el nº aleatorio obtenido. Tmb eliminamos el ion positivo mas cercano
		if(p(l).gt.aux)then
			call eliminariones(l)
			!reducimos el indice l para que se repita el calculo de probabilidad para la partícula N, ahora en el lugar l al haber elminiado la partícula l 
			l=l-1
		end if
		l=l+1
	end do
	
	!Por ultimo movemos cada partícula 
	do l=1,Nneg
		r_neg(l)=r_neg(l)-mu_neg*E(r_neg(l))*dt
		
	end do
	do l=1,Npos
		r_pos(l)=r_pos(l)+mu_pos*E(r_pos(l))*dt
		
	end do
	call reordenar
	t=t+dt
end subroutine




subroutine eliminariones(i2)
integer::iaux,i2
iaux=1
	!veamos cual es el ion positivo m mas cercano al negativo i2.!empezamos el bucle comparando las distancias a i de lios iones positivos 1 y 2
	do m=2,Npos
		if(dist(i2,m).lt.dist(i2,iaux))then
			iaux=m
			
		end if
	end do
	!Como realizamos la eliminación de iones? Por ejemplo sustituimos el N por el ion i o k en cuestion.
	!write(6,*)' Los iones negativo nº',i2,' y positivo nº ',iaux,' se han recombinado.'
	
	r_neg(i2)=r_neg(Nneg)
	phi_neg(i2)=phi_neg(Nneg)
	z_neg(i2)=z_neg(Nneg)
	rholocal(i2)=rholocal(Nneg)
	Nneg=Nneg-1	
	
	r_pos(iaux)=r_pos(Npos)
	phi_pos(iaux)=phi_pos(Npos)
	z_pos(iaux)=z_pos(Npos)
	Npos=Npos-1
		
end subroutine
!creamos una función que reordene los coeficientes de los iones eliminando aquellos que se hayan recombinado o esten en los electrodos.
	
subroutine reordenar
	i4=1
	!ordenamos iones negativos
	do while(i4.le.Nneg)
		if(r_neg(i4).lt.b)then
			r_neg(i4)=r_neg(Nneg)
			phi_neg(i4)=phi_neg(Nneg)
			z_neg(i4)=z_neg(Nneg)
			rholocal(i4)=rholocal(Nneg)
			Nneg=Nneg-1	
			recol_neg=recol_neg+1
			!write(6,*)'El ion negativo ',i4,' ha desaparecido.'	
		else
		i4=i4+1
		end if
	end do
	i4=1	
	do while(i4.le.Npos)
		if(r_pos(i4).gt.a)then
			r_pos(i4)=r_pos(Npos)
			phi_pos(i4)=phi_pos(Npos)
			z_pos(i4)=z_pos(Npos)
			Npos=Npos-1
			recol_pos=recol_pos+1
			!write(6,*)'El ion positivo ',i4,' ha desaparecido.'		
		else
		i4=i4+1
		end if
		
	end do
end subroutine

end program
