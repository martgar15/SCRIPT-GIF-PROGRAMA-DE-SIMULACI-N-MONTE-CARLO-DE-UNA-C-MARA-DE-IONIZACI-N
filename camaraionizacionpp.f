!Programa de simulación Monte Carlo de una cámara de ionización plano-paralela en reducción a 1D
program camaraionizacion
use randomnumber
implicit none
real*8::h,PI,d,rhomed,rho,Vol,hmax,hmin,haux,rd,dt,alfa,mu_pos,mu_neg,k,t,V,recol_neg,recol_pos,error,er,fsum,fcuadsum,fmed,f,x
real*8,dimension(:),allocatable::rholocal,x_neg,y_neg,z_neg,x_pos,y_pos,z_pos,r,phi,zp,p
integer::N0,Npos,Nneg,i,j,m,i3,i4,m1

!se calcula el volumen activo de la cámara plano-paralela conociendo la distancia entre placas y el radio de la sección circular de las mismas en cm
write(6,*)'ELECCIÓN DE SEPARACIÓN DE ELECTRODOS. Elija un valor de d entre 0.1 y 0.3 cm:'
read(5,*) d
!d=0.2
rd=0.5
PI=3.141592
Vol=PI*(rd**2)*d

write(6,*) 'ELECCIÓN DE DENSIDAD DE IONES.'
write(6,*)'Elija un valor de x entre 1 y 10 siendo la densidad inicial de iones 2.183E8/(xE4) cm**(-3):'
read(5,*) x

!Calculamos la densidad de iones sabiendo que una radiación de 1 mGy genera 2.183E8/cm³ y dividirla por un cierto factor que luego se añadirá a la constante de recombinación
rho=2.183E8/(x*10**4)
!Tomamos como constante de recombinación k =1.6E-6 cm³/s  y multiplicamos por el mismo factor que habiamos dividido la densidad
k=(1.6E-6)*(x*10**4)
!Calculamos el número inicial de iones
N0=int(rho*Vol)
Nneg=N0
Npos=NNeg
!Debemos elegir un voltaje, así como un valor a las movilidades en cm²/(V*s) y el parámetro alfa.

write(6,*) 'ELECCIÓN DE VOLTAJE. Elija un valor de V entre 50 y 400 V:'
read(5,*) V
!V=400
mu_pos=1.80
mu_neg=2.10
alfa=100.0

allocate(rholocal(1:Nneg))
allocate(x_neg(1:Nneg))
allocate(y_neg(1:Nneg))
allocate(z_neg(1:Nneg))
allocate(x_pos(1:Npos))
allocate(y_pos(1:Npos))
allocate(z_pos(1:Npos))
allocate(r(1:Nneg))
allocate(phi(1:Nneg))
allocate(zp(1:Nneg))
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
	
	!hacemos evolucionar el sistema. Los iones se recombinarán, se moverán debido al campo eléctrico generado en los electrodos y aquellos que alcancen los electrodos serán contados.
	OPEN(unit=1,file='posicionesneg',status='old')
	OPEN(unit=2,file='posicionespos',status='old')
	call dran_ini(m1+1412*m1)
	!elijamos una distribucion inicial de iones.
	call distribucionini
	!hallamos un valor apropiado de h para las densidades locales.
	call calcularh
	!calculamos el intervalo temporal y elegimos el campo máximo, el cual será para el radio mínimo->b
	dt=h/(mu_neg*(V/d)*alfa)
	!do while (rhomed.gt.0)
	do while (Npos.gt.1)
		!ahora vamos a escribir la posicion de cada ion negativo y positivo con objetivo de hacer gift
		do i=1,N0
			if(i.le.Nneg)then
			write(1,*)x_neg(i),y_neg(i),z_neg(i)
			else 
			write(1,*) 10,10,10
			end if
			end do
			write(1,*)' '
			write(1,*)' '
		do i=1,N0
			if(i.le.Npos)then
			write(2,*)x_pos(i),y_pos(i),z_pos(i)
			else 
			write(2,*) 10,10,10
			end if
			end do
			write(2,*)' '
			write(2,*)' '
		!evolucionamos el sistema
		call pasotemporal
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

	!write(6,*)m1,' ',fmed,' ',error,' ',er
	write(6,*)'La eficiencia de recolección es del ',fmed,' %.'
	m1=m1+1
!end do

stop
contains

!función para calcular la distancia entre un ión negativo k y un ion positivo l
real*8 function dist(o,l)
	integer::o,l
	dist=sqrt((x_pos(l)-x_neg(o))**2+(y_pos(l)-y_neg(o))**2+(z_pos(l)-z_neg(o))**2)
end function

subroutine calcularh
	real*8::aux
	!elegimos un valor inicial del error(aux) mayor que 0.01 para que se produzca el do while
	aux=1
	!implementamos el algoritmo de biseccion estimando un intervalo que puede contener a h, por ejemplo entre el origen y ek radio de la seccióncilíndrica de la lamina plano-paralela
	hmax=rd
	hmin=0.0
	!tenemos que realizar el bucle las neces que sea necesario hasta alcanzar un valor apropiado de h
	do while(aux.gt.(0.01))
		!do m=1,3
		h=(hmax+hmin)/2.0
		!write(6,*) 'El valor de h es: ',h
		call densidades
		!veamos el error relativo de pholocal para un h dado
		aux=abs(rhomed-rho)/rho
		!write(6,*)'Densidad promedio de densidades locales: ',rhomed
		!write(6,*)'Densidad total: ',rho
		!write(6,*)'error relativo', aux
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

!tmb necesirtaremos una subroutina que cree una distribución uniforme a partir de una densidad dada. Se ha considerado un espacio cilíndrico homogeneo cuyas bases son las 2 placas(electrodos).
subroutine distribucionini
	do m=1,Nneg
		!x(m)=1*dran_u()
		!y(m)=1*dran_u()
		!z(m)=1*dran_u()
		!tambien podemos emplear coordenadas cilíndricas que sería el caso equivalente a la cámara plano paralela teniendo 3 coordenadas r,phi y la distancia al primer electrdo zp que será la 'altura' del cilindro.
		phi(m)=abs(dran_u()*2*PI)
		!r(m)=abs(rd*dran_u())
		!tras crear una imagen vemos que es una distribución en el interior de un cilindro pero no es homogénea :(.
		!el culpable de que la distribución no sea homogénea es que la superficie de un disco crece con r² luego para subsanar este efecto le calcularemos la raiz a la distribucion unitaria aleatoria haciendo que ahora sí la densidad superficial se mantenga constante en todo el disco y no haya acumulación de partículas para radios pequeños.
		r(m)=abs(rd*sqrt(dran_u()))
		zp(m)=abs(d*dran_u())
		!pasamos a coordenadas cartesianas
		x_neg(m)=r(m)*cos(phi(m))
		y_neg(m)=zp(m)
		z_neg(m)=r(m)*sin(phi(m))!
		!En el instante inicial los iones positivos estarán junto a los negativos
		x_pos(m)=x_neg(m)
		y_pos(m)=y_neg(m)
		z_pos(m)=z_neg(m)
		
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
					!write(6,*) rholocal(i)
					!write(6,*)'Densidad local de ion numero ',i,' :',rholocal(i)
					!write(6,*)'El ionnegativo nº ',i,'interactua con el ion positivo ',j,' por proximidad'
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
		y_neg(l)=y_neg(l)-mu_neg*(V/d)*dt
		
		!veamos que iones son recolectados
		!if(y_neg(l).lt.0)then
		!	recol_neg=recol_neg+1
		!	Write(6,*)'El ion negativo ',l,' ha alcanzado el electrodo positivo'
			
		!end if
	end do
	do l=1,Npos
		y_pos(l)=y_pos(l)+mu_pos*(V/d)*dt
		!if (y_pos(l).gt.d)then
		!	recol_pos=recol_pos+1
		!	Write(6,*)'El ion positivo ',l,' ha alcanzado el electrodo negativo'
		!end if
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
	
	x_neg(i2)=x_neg(Nneg)
	y_neg(i2)=y_neg(Nneg)
	z_neg(i2)=z_neg(Nneg)
	rholocal(i2)=rholocal(Nneg)
	Nneg=Nneg-1	
	
	x_pos(iaux)=x_pos(Npos)
	y_pos(iaux)=y_pos(Npos)
	z_pos(iaux)=z_pos(Npos)
	Npos=Npos-1
		
end subroutine
!creamos una función que reordene los coeficientes de los iones eliminando aquellos que se hayan recombinado o esten en los electrodos.
	
subroutine reordenar
	i4=1
	!ordenamos iones negativos
	do while(i4.le.Nneg)
		if(y_neg(i4).lt.0)then
			x_neg(i4)=x_neg(Nneg)
			y_neg(i4)=y_neg(Nneg)
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
		if(y_pos(i4).gt.d)then
			x_pos(i4)=x_pos(Npos)
			y_pos(i4)=y_pos(Npos)
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
