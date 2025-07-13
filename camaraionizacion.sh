#!/bin/bash

echo "Selecciona un tipo de simulación pulsando 1 o 2:"
echo "1 - Cámara plano-paralela"
echo "2 - Clamara cilíndrica"

read -p "Introduce 1 o 2: " opcion

if [ "$opcion" = "1" ]; then

	echo " Compilando el programa Fortran..."

	#Compilar el programa principal junto con el objeto randomnumber.o
	gfortran-11 camaraionizacionpp.f -c -ffree-form -fcheck=all -g -O3
	
	gfortran-11 camaraionizacionpp.o randomnumber.o -o camaraionizacionpp.exe
	
	
	#Verificar si la compilación ha sido exitosa
	if [ $? -ne 0 ]; then
	    echo "Error en la compilación."
	    exit 1
	fi
	
	echo " Ejecutando el programa Fortran..."
	./camaraionizacionpp.exe
	
	#Verificar si se ha ejecutado correctamente
	if [ $? -ne 0 ]; then
	    echo " Error al ejecutar el programa."
	    exit 1
	fi
	
	echo " Generando GIF con Gnuplot..."
	gnuplot animacionpp.plt
	
	#Verificar si Gnuplot ha creado la grafica
	if [ $? -ne 0 ]; then
	    echo "Error al generar la animación con Gnuplot."
	    exit 1
	fi
	
	echo "¡GIF generado exitosamente!"
	
	xdg-open animacionpp.gift
elif [ "$opcion" = "2" ]; then
	
	echo " Compilando el programa Fortran..."

	#Compilar el programa principal junto con el objeto randomnumber.o
	gfortran-11 camaraionizacioncil.f -c -ffree-form -fcheck=all -g -O3
	
	gfortran-11 camaraionizacioncil.o randomnumber.o -o camaraionizacioncil.exe
	
	
	#Verificar si la compilación ha sido exitosa
	if [ $? -ne 0 ]; then
	    echo "Error en la compilación."
	    exit 1
	fi
	
	echo " Ejecutando el programa Fortran..."
	./camaraionizacioncil.exe
	
	#Verificar si se ha ejecutado correctamente
	if [ $? -ne 0 ]; then
	    echo " Error al ejecutar el programa."
	    exit 1
	fi
	
	echo " Generando GIF con Gnuplot..."
	gnuplot animacioncil.plt
	
	#Verificar si Gnuplot ha creado la grafica
	if [ $? -ne 0 ]; then
	    echo "Error al generar la animación con Gnuplot."
	    exit 1
	fi
	
	echo "¡GIF generado exitosamente!"
	
	xdg-open animacioncil.gift
fi
