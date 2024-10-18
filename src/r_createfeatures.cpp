// [[Rcpp::plugins(cpp20)]]

#undef NDEBUG

#include "createfeatures.h"
#include <Rcpp.h>
using namespace Rcpp;

//' @name UGridAndPoints
//' @title Represents a grid of tree loaded from XYZ file
//' @description Most important method is CalculateFeaturesR
//' @field new Constructor \itemize{
//' \item FileNamePath: path to the XYZ file
//' \item swapYZ: TRUE/FALSE
//' \item cellSize: Size of the cell
//' }
//' @field CalculateFeaturesR Calculates features by using symmetry measure
//' @field ReadFileToGridR Reads the tree from XYZ file
//' @field GetPointsR Gets a list of all points of the XYZ file
struct UGridAndPointsR: public UGridAndPoints{
    void ReadFileToGridR(const Rcpp::String path, bool swapYZ, float cellSize);
    Rcpp::NumericMatrix GetPointsR();
    Rcpp::NumericMatrix GetGridAsPointsR();
    Rcpp::List CalculateFeaturesR();
};

void Clean(UGridAndPointsR *f){
    //delete f; // TODO: error: double free or corruption
}

RCPP_MODULE(UGridAndPointsR) {
    //Rcpp::function("norm",&main);
    Rcpp::class_<UGridAndPointsR>("UGridAndPoints","SAS")
        .constructor("Creates new Grid")
        //.finalizer(&Clean)
        .method("ReadFileToGridR", &UGridAndPointsR::ReadFileToGridR, "Reads the tree from XYZ file")
        .method("GetPointsR", &UGridAndPointsR::GetPointsR, "Gets a list of all points of the XYZ file")
        .method("GetGridAsPointsR", &UGridAndPointsR::GetGridAsPointsR)
        .method("CalculateFeaturesR", &UGridAndPointsR::CalculateFeaturesR);
    //  .field()
}

void UGridAndPointsR::ReadFileToGridR(const String path, bool swapYZ, float cellSize){
    ReadFileToGrid(path.get_cstring(), swapYZ, cellSize);
}

NumericMatrix UGridAndPointsR::GetPointsR(){
    NumericMatrix res(Pts.size(), 3);
    int r=0;
    for(auto p:Pts){
        res(r, 0) = p.x;
        res(r, 1) = p.y;
        res(r, 2) = p.z;
        r++;
    }
    return res;
}

NumericMatrix UGridAndPointsR::GetGridAsPointsR(){
    NumericMatrix res(G.Cells.size(), 4);
    int r=0;
    for(int x = 0; x < G.Cells.extents()[0]; x++){
        for(int y = 0; y < G.Cells.extents()[1]; y++){
            for(int z = 0; z < G.Cells.extents()[2]; z++){
                res(r, 0) = x;
                res(r, 1) = y;
                res(r, 2) = z;
                res(r, 3) = G(x,y,z);
                r++;
            }
        }
    }
    return res;
}



List UGridAndPointsR::CalculateFeaturesR(){
    CalculateFeatures();
    //TODO: remove outliers!
    Rcpp::List l;
    l.push_back(Feat.UnionRevolutionAndOrg, "UnionRevolutionAndOrg");
    l.push_back(Feat.IntersectionRevolutionAndOrg, "IntersectionRevolutionAndOrg");
    l.push_back(Feat.IOU_RevolutionAndOrg,"IOU_RevolutionAndOrg");

    {
        Rcpp::NumericVector v(Feat.NumPointsPerSliceZ.size());
        int r=0;
        for(auto i: Feat.NumPointsPerSliceZ){
            v(r++)=i;
        }
        l.push_back(v,"NumPointsPerSlice");
    }
    {
        Rcpp::NumericVector v(Feat.RadiousPerSliceZ.size());
        int r=0;
        for(auto i: Feat.RadiousPerSliceZ){
            v(r++)=i;
        }
        l.push_back(v,"RadiousPerSliceZ");
    }
    {
        Rcpp::NumericVector v(Feat.AvgDistPerSliceZ.size());
        int r=0;
        for(auto i: Feat.AvgDistPerSliceZ){
            v(r++)=i;
        }
        l.push_back(v,"AvgDistPerSliceZ");
    }

    {
        Rcpp::NumericVector v(Feat.IOUs_RevolutionAndOrgPerSlice.size());
        int r=0;
        for(auto i: Feat.IOUs_RevolutionAndOrgPerSlice){
            v(r++)=i;
        }
        l.push_back(v,"IOUs_RevolutionAndOrgPerSlice");
    }

    {
        Rcpp::NumericVector v(Feat.NumPointsFromVerticalLine.size());
        int r=0;
        for(auto i: Feat.NumPointsFromVerticalLine){
            v(r++)=i;
        }
        l.push_back(v,"NumPointsFromVerticalLine");
    }

    return l;
}

