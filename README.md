# ludecomposer
A matrix decomposition algorithm

## Extras
Has support for multiprecision operation

## Dependency
<ul>
 <li>gmp - Gnu multiprecision library</li>
</ul>

## Example
```LUDecomposer::MatrixD mat({ { 1, -1, 0}, { 2, 2, 3}, { -1, 3, 2}});
    LUDecomposer dcmp( mat);
    const auto& [ lm, um]  = dcmp.soln();

    std::cout << "Given:\n";
    printMatrix( mat);

    std::cout << "Lower Triangular Matrix:\n";
    printMatrix( lm);

    std::cout << "Upper Triangular Matrix:\n";
    printMatrix( um);
```

## Building
 ### Without gmp
```
g++ ludecomposer.cpp -o lu
```
 ### With gmp
```
g++ -D__MP ludecomposer.cpp -lgmp -lgmpxx -o ludecomposer
```
