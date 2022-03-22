#include <iostream>
#include <unordered_map>
#include <gmpxx.h>

int main()
{
    std::unordered_map<char, std::string> charmap = {

        { '0', "\u2070"},
        { '1', "\u00B9"},
        { '2', "\u00B2"},
        { '3', "\u00B3"},
        { '4', "\u2074"},
        { '5', "\u2075"},
        { '6', "\u2076"},
        { '7', "\u2077"},
        { '8', "\u2078"},
        { '9', "\u2079"},
        { '+', "\u207A"},
        { '-', "\u207B"}
    };

    std::string text = "3x" + charmap[ '2'] + "+"
                        "5x" + charmap[ '-'] + charmap[ '3']
                        + charmap[ '('] +
                        charmap[ '8'] + charmap[ '+']
                        + charmap[ '9'] + charmap[ ')'];

    std::cout << text <<'\n';
}
