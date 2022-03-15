    #include "matrice.h"

Matrice random ( unsigned int L, unsigned int C )
{
    Matrice M = Matrice ( L, C );
    for ( unsigned int i = 0 ; i < L ; i++ )
    {
        for ( unsigned int j = 0 ; j < C ; j++ )
        {
            M.setCoeff( i, j, (rand()%1000)/100.0 );
        }
    }
    return M;
}

int main(int argc, char *argv[])
{
    Matrice M = random ( 5, 5 );
    M.affiche();

    Matrice inv_M = M.inverse();;
    inv_M.affiche();

    Matrice prod = M*inv_M;
    prod.affiche();

    Vecteur vecPropre = M.powerIteration();
    vecPropre.affiche();

    M.calculValeurPropreAssociee(vecPropre);
    Vecteur vecPropreMax = M.getVecteurValeursPropres();
    vecPropreMax.affiche();

    M.eigenAnalysis();
    for ( unsigned int i = 0 ; i < 5 ; i++ )
    {
        Vecteur produit = M.getValeurPropre(i) * M.getVecteurPropre(i);
        Vecteur transfo = M * M.getVecteurPropre(i);

        std::cout << M.getValeurPropre(i) << " -> ";
        produit.affiche();
        std::cout << M.getValeurPropre(i) << " -> ";
        transfo.affiche();
        std::cout << std::endl;
    }

    return 0;
}
