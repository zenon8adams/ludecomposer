#include <vector>
#include <string>
#include <initializer_list>
#include <iterator>
#include <unordered_map>
#include <numeric>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "algebra-basic.hpp"

class LUDecomposer
{
public:

    using MatrixD = std::vector<std::vector<value_t>>;
    using MatrixAL = std::vector<std::vector<AlUnit>>;
    using Matrix2D = std::pair<MatrixD, MatrixD>;

    LUDecomposer( const MatrixD& given)
    : alpha( given)
    {
        resize( ldelta);
        resize( udelta);
        for( size_t i = 0; i < alpha.size(); ++i)
        {
            for( size_t j = 0; j < alpha.size(); ++j)
            {
                size_t idx = (i+1) * 10 + (j+1);
                std::string ld = "L"s + std::to_string( idx),
                            ud = "U"s + std::to_string( idx);
                if( j <= i)
                {
                    ldelta[ i][ j] = Monomial( 1, {ld});
                    answers[ ld] = { 1, i == j};
                }
                
                if( j >= i)
                {
                    udelta[ i][ j] = Monomial( 1, {ud});
                    answers[ ud] = { 0, false};
                }
            }
        }
    }

    Matrix2D soln()
    {
        auto lu = matrixProduct();
        
        for( size_t i = 0; i < alpha.size(); ++i)
        {
            for( size_t j = 0; j < alpha[ i].size(); ++j)
            {
                lu[ i][ j] = lu[ i][ j].eval( answers);
                auto rem = lu[ i][ j].unknown();
                if( rem.size() == 1)
                {
                    auto weight = -(lu[ i][ j] - rem[0] - Monomial(alpha[ i][ j])).value()
                                            / rem[ 0].coefficient();
                    answers[ rem[ 0].hashVal()] = { weight, true};
                }
            }
        }

        return toMatrix2D( lu);
    }

private:

    Matrix2D toMatrix2D( const MatrixAL& lu)
    {
        Matrix2D result;
        auto& [ lm, um] = result;
        resize( lm);
        resize( um);

        for( size_t i = 0; i < alpha.size(); ++i)
        {
            for( size_t j = 0; j < alpha[ i].size(); ++j)
            {
                lm[ i][ j] = answers[ ldelta[ i][ j].hashVal()].first;
                um[ i][ j] = answers[ udelta[ i][ j].hashVal()].first;
            }
        }

        return result;
    }

    MatrixAL matrixProduct()
    {
        MatrixAL result;
        resize( result);
        for( size_t i = 0; i < alpha.size(); ++i)
            for( size_t j = 0; j < alpha.size(); ++j)
                for( size_t k = 0; k < alpha.size(); ++k)
                    result[ i][ j] += ldelta[ i][ k] * udelta[ k][ j];

        return result;
    }


    template <typename Matrix>
    void resize( Matrix& mat)
    {
        mat.resize( alpha.size());
        for( size_t i = 0; i < alpha.size(); ++i)
            mat[ i].resize( alpha[ i].size());
    }


    MatrixD alpha;

    MatrixAL ldelta,
             udelta;
    
    result_t answers;

};

template <typename Container>
void printMatrix( const Container& matrix)
{
    for( size_t i = 0; i < matrix.size(); ++i )
    {
        for( size_t j = 0; j < matrix[ i].size(); ++j)
            std::cout << std::setprecision( 3) << (ZERO(matrix[ i][ j]) ? 0 : matrix[ i][ j]) <<'\t';
        putchar('\n');
    }

}

int main()
{
    LUDecomposer::MatrixD mat({ { 1, -1, 0}, { 2, 2, 3}, { -1, 3, 2}});
    LUDecomposer dcmp( mat);
    const auto& [ lm, um]  = dcmp.soln();

    std::cout << "Given:\n";
    printMatrix( mat);

    std::cout << "Lower Triangular Matrix:\n";
    printMatrix( lm);

    std::cout << "Upper Triangular Matrix:\n";
    printMatrix( um);
}
