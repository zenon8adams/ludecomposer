#include <iostream>
#include <cstdio>
#include <iomanip>
#include <unordered_map>
#include <deque>
#include "macros.hpp"
#include "algebra-basic.hpp"

template <typename Matrix>
Matrix crop(  const Matrix& old_mat,  size_t col)
{
    size_t oldsize = old_mat.size();
    size_t newsize = oldsize - 1;
    Matrix new_mat( newsize);
    for( size_t i = 0; i < newsize; ++i)
        new_mat[ i].resize( newsize);

    for( size_t i = 1, n = 0; i < oldsize; ++i)
    {
        for( size_t j = 0; j < oldsize; ++j)
        {
            if( j != col)
            {
                new_mat[ n / newsize][ n % newsize] = old_mat[ i][ j];
                ++n;
            }
        }
    }

    return new_mat;
}

template <typename Matrix>
typename Matrix::value_type::value_type
        det2X2( const Matrix& matrix)
{
    return matrix[ 0][ 0] * matrix[ 1][ 1]
                          -
           matrix[ 0][ 1] * matrix[ 1][ 0];
}
             
template <typename Matrix>
typename Matrix::value_type::value_type determinantOf( const Matrix& matrix)
{
    if( matrix.size() == 2)
        return det2X2( matrix);
    static size_t level = 0;

    printf( "Level: %zu\n", ++level);
    
    typename Matrix::value_type::value_type sum{};
    for( size_t i = 0, sz = matrix.size(); i < sz; ++i)
        sum += ((((int)i & 1) << 1) - 1) * matrix[ 0][ i] 
                        * 
                determinantOf( crop( matrix, i));

    return sum;
}

template <typename Matrix>
void printMatrix( const Matrix& matrix)
{
    for( size_t i = 0; i < matrix.size(); ++i)
    {
        for( size_t j = 0; j < matrix[ i].size(); ++j)
            std::cout << matrix[ i][ j] << '\t';
        putchar( '\n');
    }
}


template <typename Matrix, 
          typename InnerType = typename Matrix::value_type::value_type>
InnerType determinantOf_( const Matrix& matrix)
{
    std::unordered_map<size_t, InnerType> latch;
    std::deque<size_t> indexes;
    size_t N = matrix.size();
    for( size_t i = 0; i < N; ++i)
    {
        indexes.push_front( i);
        for( size_t j = 0; j < N; ++j)
            latch[ i * N + j] = matrix[ i][ j];
    }
    printf( "Size: %zu\n", indexes.size());
    size_t start, end = N;
    while( !indexes.empty())
    {
        auto current = indexes.back();
        indexes.pop_back();
        size_t i = current / N, j = current % N;
        if( i == j)
        {
            start = i + 1;
            end -= 1;
            printf( "s: %zu, e: %zu\n", start, end);
        }
        for( size_t k = 0, idx = i; k < end; ++k, ++idx)
        {
            if( j == idx)
                ++idx;
            printf( "( %zu, %zu) = %d\n", start, idx, matrix[ start][ idx]);
            indexes.push_front( start * N + idx);
        }
        putchar('\n');
        if( start == end)
            break;
    }
    return {};
}

            


int main()
{
    Variable v1( "a", 2), v2( "b", 3), v3( "c", -2),
             v4( "a", 1);
    Monomial m1( 4, {v1}), m2( -2, {v2}), m3( 1, {v3}),
             m4( 6, {v4}), m5( 0);
    AlUnit a1( { m1, m3, m4}), a2( m2), a3( m3), a4( m4), a5( {});

    std::vector<std::vector<AlUnit>> matrix{
        { m2, m2, m3, m4, m1, m3, m1, m3},
        { m4, m2, m4, m1, m5, m2, m2, m4},
        { m2, m1, m3, m5, m3, m4, m2, m1},
        { m1, m2, m2, m4, m1, m3, m1, m5},
        { m2, m4, m1, m1, m5, m2, m2, m1},
        { m5, m4, m3, m3, m1, m2, m2, m2},
        { m1, m4, m2, m1, m1, m3, m5, m3},
        { m5, m5, m1, m3, m1, m2, m3, m1}
        };

    //std::cout << determinantOf( matrix).toString() <<'\n';
    //
    /*std::vector<std::vector<int>> mat{
    { 1, 4, -1, 5, 6},
    { 3, 0,  1, 4, 2},
    { 2, 7,  4, 1, 3},
    { 5, 9, -2, 8, 1},
    { 2, 7,  3, 2, 9}
    };

    determinantOf_( mat);
    */
    std::cout << "dF("<< a1.toString() << ")/da: "
              << a1.integral( "a").toString() <<'\n';


}
