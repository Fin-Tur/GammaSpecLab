#include <iostream>
#include "driver/GammaMBibSearchTest.h"
#include "driver/SecondDSearchTest.h"

int main() {
    /*GammaMBibSearchTest tester;
    tester.test_gamma_m_bib_search();*/
    SecondDSearchTest::test_seconddsearch();
    std::cout << "Finished";
}
