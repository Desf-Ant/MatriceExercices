#ifndef MATRICE_H
#define MATRICE_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "vecteur.h"

class Matrice
{
    private :
        unsigned int nb_lignes;
        unsigned int nb_colonnes;
        std::vector <double> lesCoeffs;
        bool eigen_analysis_effectuee;

        std::vector < double > lesValeursPropres;
        std::vector < Vecteur > lesVecteursPropres;

    public :
        // Constructeurs
        Matrice ( unsigned int size );
        Matrice ( unsigned int L, unsigned int C );
        Matrice ( const Matrice& M );
        //Matrice ( const Eigen::MatrixXd& eigen );
        Matrice ( const Vecteur& diagonale );
        ~Matrice();

        // Accesseurs
        void affiche ( void ) const;
        std::vector<double> getCoeffs ( void ) const;
        unsigned int getNbLignes ( void ) const;
        unsigned int getNbColonnes ( void ) const;
        bool isSquare( void ) const;
        double getCoeff ( unsigned int i, unsigned int j ) const;
        void setCoeff ( unsigned int i, unsigned int j, double val );

        // Méthodes
        Matrice inverse ( void ) const;
        void eigenAnalysis ( void );
        Vecteur powerIteration ( void );
        double calculValeurPropreAssociee ( Vecteur v );
        void ordonneValeurPropre ( void );
        void ordonneValeurPropreAbsolue ( void );

        // Accesseurs après analyse
        Matrice getMatriceValeursPropres ( void );
        Vecteur getVecteurValeursPropres ( void );
        float getValeurPropre ( unsigned int n );
        Matrice getMatriceVecteursPropres ( void );
        Vecteur getVecteurPropre ( unsigned int n );

    public :
        static Matrice identite ( unsigned int size );
        static Matrice tensoriel (const Vecteur& A, const Vecteur& B);
};

Matrice operator* ( const double lambda, const Matrice& M );
Vecteur operator* ( const Matrice& M, const Vecteur& V );
Vecteur operator* ( const Vecteur& V, const Matrice& M );
Matrice operator* ( const Matrice& A, const Matrice& B );
Matrice operator+ ( const Matrice& A, const Matrice& B );
Matrice operator- ( const Matrice& A, const Matrice& B );

#endif // MATRICE_H
