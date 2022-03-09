#include "vecteur.h"

Vecteur::Vecteur( unsigned int size )
{
    this->taille = size;
    this->lesCoeffs = std::vector<double>();
    for ( unsigned int i = 0 ; i < this->taille ; i++ )
    {
        this->lesCoeffs.push_back( 0.0 );
    }
}

Vecteur::Vecteur( const Vecteur& vec )
{
    this->taille  = vec.getTaille();
    this->lesCoeffs = std::vector<double>();
    for ( unsigned int i = 0 ; i < this->taille ; i++ )
    {
        this->lesCoeffs.push_back( vec.getCoeff(i) );
    }
}

Vecteur::~Vecteur()
{
    // Rien
}

void Vecteur::affiche ( void ) const
{
    if ( this->taille == 0 )
    {
        std::cout << "[ ]" <<std::endl;
    }
    else
    {
        std::cout << "[ ";
        for ( unsigned int i = 0 ; i < this->taille-1 ; i++ )
        {
            std::cout << ((int)(100 * this->getCoeff(i)))/100.0 << " ; ";
        }
        std::cout << ((int)(100 * this->getCoeff( this->taille-1 )))/100.0 << " ] " << std::endl;
    }
}

void Vecteur::setCoeff ( unsigned int index, float valeur )
{
    this->lesCoeffs[index] = valeur;
}

float Vecteur::getCoeff ( unsigned int index ) const
{
    return this->lesCoeffs[index];
}

unsigned int Vecteur::getTaille ( void ) const
{
    return this->taille;
}

double Vecteur::getNorme ( void ) const
{
    double norme = 0.0;
    for ( unsigned int i = 0 ; i < this->taille ; i++ )
    {
        norme += this->getCoeff(i) * this->getCoeff(i);
    }
    return sqrt ( norme );
}

Vecteur Vecteur::normalise ( void )
{
    Vecteur res = Vecteur ( this->taille );
    double norme = this->getNorme();
    for ( unsigned int i = 0 ; i < res.getTaille() ; i++ )
    {
        res.setCoeff(i, this->getCoeff(i)/norme );
    }
    return res;
}

double Vecteur::produitScalaireNorme ( const Vecteur& autre )
{
    double scal = (*this) * autre;
    scal /= this->getNorme();
    scal /= autre.getNorme();
    return scal;
}

Vecteur operator* ( const double lambda, const Vecteur& V )
{
    Vecteur res = Vecteur ( V.getTaille() );
    for ( unsigned int i = 0 ; i < res.getTaille() ; i++ )
    {
        res.setCoeff( i, lambda*V.getCoeff(i) );
    }
    return res;
}

bool operator== ( const Vecteur& A, const Vecteur& B )
{
    if ( A.getTaille() != B.getTaille() ) return false;
    double somme_diff = 0.0;
    for ( unsigned int i = 0 ; i < A.getTaille() ; i++ )
    {
        somme_diff += abs ( A.getCoeff(i) - B.getCoeff(i) );
    }
    return somme_diff < 0.00001;
}

bool operator!= ( const Vecteur& A, const Vecteur& B )
{
    return !( A == B );
}

Vecteur operator+ ( const Vecteur& A, const Vecteur& B )
{
    if ( A.getTaille() != B.getTaille() ) return Vecteur(0);
    Vecteur res = Vecteur ( A.getTaille() );
    for ( unsigned int i = 0 ; i < res.getTaille() ; i++ )
    {
        res.setCoeff( i, A.getCoeff(i) + B.getCoeff(i) );
    }
    return res;
}

Vecteur operator- ( const Vecteur& A, const Vecteur& B )
{
    if ( A.getTaille() != B.getTaille() ) return Vecteur(0);
    Vecteur res = Vecteur ( A.getTaille() );
    for ( unsigned int i = 0 ; i < res.getTaille() ; i++ )
    {
        res.setCoeff( i, A.getCoeff(i) - B.getCoeff(i) );
    }
    return res;
}

double operator* ( const Vecteur& A, const Vecteur& B )
{
    if ( A.getTaille() != B.getTaille() ) return 0.0;
    double scal = 0.0;
    for ( unsigned int i = 0 ; i < A.getTaille() ; i++ )
    {
        scal += A.getCoeff(i) * B.getCoeff(i);
    }
    return scal;
}
