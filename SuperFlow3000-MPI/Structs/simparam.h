/*
 * simparam.h
 *
 *  Created on: 05.11.2013
 *      Author: felix
 */

#ifndef SIMPARAM_H_
#define SIMPARAM_H_

#include "../Misc/template.h"
struct Simparam {
	RealType xLength;  //Gebietslänge in x-Richtung
	RealType yLength;  //Gebietslänge in y.Richtung
	int iMax;          //Anzahl der inneren Zellin in x-Richtung
	int jMax;          //Anzahl der inneren Zellen in y-Richtung

	RealType tEnd;		//Endzeit
	RealType deltaT;
	RealType tau;   //Sicherheitsfaktor
	RealType deltaVec;   //Zeitabstand für die Ausgabe

	IndexType iterMax;     // Maximale Anzahl an Iterationen
	RealType eps;	//Toleranz für Druckiteration
	RealType omg;	//Relaxationsfaktor
	RealType alpha;	//Upwind-Differencing-Faktor
	RealType RE;			//Reynoldszahl
	RealType GX;			//Äußere Kräfte g_x und g_y
	RealType GY;
	RealType VI;			//Initailwerte für u,v undp
	RealType UI;
	RealType PI;
};


#endif /* SIMPARAM_H_ */
