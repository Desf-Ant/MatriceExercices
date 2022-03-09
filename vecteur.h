#ifndef VECTEUR_H
#define VECTEUR_H

#include <math.h>
#include <vector>
#include <iostream>

class Vecteur
{
    private :
        unsigned int taille;
        std::vector<double> lesCoeffs;

    public :
        Vecteur( unsigned int size );
        Vecteur( const Vecteur& vec );
        ~Vecteur();

        void affiche ( void ) const;
        void setCoeff ( unsigned int index, float valeur );
        float getCoeff ( unsigned int index ) const;
        unsigned int getTaille ( void ) const;
        double getNorme ( void ) const;
        Vecteur normalise ( void );
        double produitScalaireNorme ( const Vecteur& autre );
};

Vecteur operator* ( const double lambda, const Vecteur& V );
bool operator== ( const Vecteur& A, const Vecteur& B );
bool operator!= ( const Vecteur& A, const Vecteur& B );
Vecteur operator+ ( const Vecteur& A, const Vecteur& B );
Vecteur operator- ( const Vecteur& A, const Vecteur& B );
double operator* ( const Vecteur& A, const Vecteur& B );

#endif // VECTEUR_H
