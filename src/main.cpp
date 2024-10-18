// [[Rcpp::plugins(cpp20)]]

#undef NDEBUG
//needed for Rstudio Rcpp compilation!!!

#include <glm/common.hpp>
#include <glm/fwd.hpp>
#include <glm/geometric.hpp>
#include <iostream>
#include <Rcpp.h>
#include "createfeatures.h"
//#include <RInside.h>
using namespace Rcpp;
using namespace std;
using glm::vec3;
#include <boost/numeric/ublas/storage.hpp>

int main()
{
    cout << "Starting processing!" << endl;

    UGridAndPoints gp;

    gp.CalculateFeatures();

    UGrid g(0.5, {0,0,0},{1,2,3});
    g(vec3{0.5f,0.5f,0.5f})++;
    g(vec3{0.5f,0.5f,0.5f})++;

    cout << "Done processing!" << endl;
    return 0;
}
