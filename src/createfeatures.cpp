#undef NDEBUG

#include "createfeatures.h"
//needed for Rstudio Rcpp compilation!!!

#include <glm/common.hpp>
#include <glm/fwd.hpp>
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>
#include <iostream>

#include <filesystem>
//#include <RInside.h>

using namespace std;
using glm::vec3;

#include <boost/numeric/ublas/tensor/expression.hpp>

auto operator< (UGrid::itensor const& lhs, UGrid::itensor const& rhs) {
    return boost::numeric::ublas::detail::make_binary_tensor_expression<UGrid::itensor> (lhs, rhs, [](int const& l, int const& r){ return l < r; });
}
auto operator> (UGrid::itensor const& lhs, UGrid::itensor const& rhs) {
    return boost::numeric::ublas::detail::make_binary_tensor_expression<UGrid::itensor> (lhs, rhs, [](int const& l, int const& r){ return l > r; });
}
auto operator|| (UGrid::itensor const& lhs, UGrid::itensor const& rhs) {
    return boost::numeric::ublas::detail::make_binary_tensor_expression<UGrid::itensor> (lhs, rhs, [](int const& l, int const& r){ return l || r; });
}
auto operator&& (UGrid::itensor const& lhs, UGrid::itensor const& rhs) {
    return boost::numeric::ublas::detail::make_binary_tensor_expression<UGrid::itensor> (lhs, rhs, [](int const& l, int const& r){ return l && r; });
}

#include <glm/ext.hpp>
#include <fstream>
#include <algorithm>


UGridAndPoints::UGridAndPoints(){

}

UGridAndPoints::~UGridAndPoints(){
}

void UGridAndPoints::ReadFileToGrid(const string fpath, bool swapYZ, float cellSize){
    Pts.clear();
    FPath=fpath;
    vec3 bmin(std::numeric_limits<float>::max());
    vec3 bmax(std::numeric_limits<float>::min());
    {
        std::ifstream fileIn(fpath);
        vec3 p;

        std::string line;
        while(std::getline(fileIn,line)){
            if (true){
                std::stringstream ss(line);
                std::string token;
                if(ss >> p.x >> p.y >> p.z){
                    if(swapYZ){
                        swap(p.y,p.z);
                    }
                    Pts.push_back(p);
                    bmin=glm::min(bmin, p);
                    bmax=glm::max(bmax, p);
                }
            }
        }

        /*while(fileIn >> p.x >> p.y >> p.z){
            Pts.push_back(p);
            bmin=glm::min(bmin, p);
            bmax=glm::max(bmax, p);
        }*/
    }
    if(Pts.size()==0){
        throw std::runtime_error("Error reading file " + fpath);
    }
    assert(Pts.size()>0);

    vec3 center = std::accumulate(Pts.begin(), Pts.end(), vec3(0.0f)) / float(Pts.size());//middle voxel // BEST!

    vec3 d=bmax-bmin;
    G=UGrid(cellSize,center-d,center+d);
    G.Allocate(Pts);

    Ones=UGrid::itensor(G.Cells.extents());
    Zeros=UGrid::itensor(G.Cells.extents());
    Ones=1;
    Ones=0;
}



UGrid::itensor UGridAndPoints::CalcUnion(UGrid::itensor X, UGrid::itensor Y)
{
    UGrid::itensor  uni  = (X>Zeros) || (Y>Zeros);
    return uni;
}

UGrid::itensor UGridAndPoints::CalcIntersection(UGrid::itensor X, UGrid::itensor Y)
{
    UGrid::itensor inter  = (X>Zeros) && (Y>Zeros);
    return inter;
}

void UGridAndPoints::CalculateFeatures(){

    G.Allocate(Pts);
    auto cOriginal = G.Cells;

    G.AllocateSphericalZ(Pts);

    Feat.RadiousPerSliceZ.clear();
    Feat.AvgDistPerSliceZ.clear();
    Feat.NumPointsPerSliceZ.clear();
    glm::vec2 vertCenter(G.CenterBBox().x, G.CenterBBox().y);
    for(const auto &slice: G.PointsPerSliceZ){
        float maxD=0;
        for(const auto &p: slice){
            maxD=max(maxD, glm::distance(vertCenter,glm::vec2(p.x, p.y)));
        }
        Feat.RadiousPerSliceZ.push_back(maxD);


        float avgD=std::accumulate(slice.cbegin(),slice.cend(),0.0f,[vertCenter](float lhs, auto p){
            return lhs + glm::distance(vertCenter,glm::vec2(p.x, p.y));
        });

        Feat.AvgDistPerSliceZ.push_back(
            slice.size()>0?
                avgD/slice.size():
                0);

        Feat.NumPointsPerSliceZ.push_back(slice.size());
    }


    auto cRound = G.Cells;

    UGrid::itensor uni = CalcUnion(cRound, cOriginal);

    UGrid::itensor inter = CalcIntersection(cOriginal, cRound);

    Feat.UnionRevolutionAndOrg=std::count_if(uni.begin(), uni.end(),[](int v){return v==1;});
    Feat.IntersectionRevolutionAndOrg=std::count_if(inter.begin(), inter.end(),[](int v){return v==1;});
    Feat.IOU_RevolutionAndOrg=(double)Feat.IntersectionRevolutionAndOrg/Feat.UnionRevolutionAndOrg;

    Feat.IOUs_RevolutionAndOrgPerSlice.clear();
    for(unsigned int z=0;z<uni.size(2);z++){
        int unis=0;
        int inters=0;
        for(unsigned int x=0;x<uni.size(0);x++){
            for(unsigned int y=0;y<uni.size(1);y++){
                if(uni.at(x,y,z))unis++;
                if(inter.at(x,y,z))inters++;
            }
        }
        if(unis==0)
            Feat.IOUs_RevolutionAndOrgPerSlice.push_back(0);
        else
            Feat.IOUs_RevolutionAndOrgPerSlice.push_back((float)inters/unis);
    }

    Feat.NumPointsFromVerticalLine.clear();
    for(float dist=0;dist<glm::distance(G.BBoxMin,G.BBoxMax);dist+=G.CellSize){
        int num=0;
        for(const auto &p:G.Pts){
            if(glm::distance(vec3(vertCenter.x,vertCenter.y,p.z),p)>=dist &&
                glm::distance(vec3(vertCenter.x,vertCenter.y,p.z),p)<(dist+G.CellSize) ){
                num++;
            }
        }
        Feat.NumPointsFromVerticalLine.push_back(num);
    }

    // symmetry per slice: IOU
}

UGrid::UGrid():BBoxMin(glm::vec3(std::numeric_limits<float>::min())),BBoxMax(glm::vec3(std::numeric_limits<float>::max())){

}

UGrid::~UGrid(){

}

UGrid::UGrid(float cellSize, glm::vec3 bboxMin, glm::vec3 bboxMax):
    BBoxMin(bboxMin),BBoxMax(bboxMax)
{
    CellSize = cellSize;
    BBoxMin = bboxMin;
    BBoxMax = bboxMax;

    //center
    auto d=bboxMax-bboxMin;
    auto cells=glm::ivec3(glm::floor(d/CellSize))+1;

    Cells=itensor{unsigned(cells.x), unsigned(cells.y), unsigned(cells.z)};
    for(auto&c:Cells)c=0;
}

vec3 UGrid::CenterBBox(){
    return (BBoxMin+BBoxMax)/2.0f;
}
using glm::ivec3;
ivec3 UGrid::CellCoord(glm::vec3 p){
    return glm::floor((p-BBoxMin)/CellSize);
}

vec3 UGrid::CellCenterPos(glm::ivec3 c){
    return BBoxMin + vec3(c)*CellSize + CellSize/2.0f;
}

int &UGrid::operator()(glm::vec3 p){
    const glm::ivec3 ci=CellCoord(p);
    return Cells.at(unsigned(ci.x), unsigned(ci.y), unsigned(ci.z));
}

int &UGrid::operator()(glm::ivec3 ci){
    return Cells.at(unsigned(ci.x), unsigned(ci.y), unsigned(ci.z));
}

int &UGrid::operator()(int x, int y, int z){
    return Cells.at(unsigned(x), unsigned(y), unsigned(z));
}

void UGrid::LoadSlices(const std::vector<glm::vec3> &points)
{
    Pts = points;
    PointsPerSliceZ=std::vector<std::vector<vec3>>(Cells.extents()[2]);
    for(auto p:points){
        int zi=CellCoord(p).z;
        PointsPerSliceZ[zi].push_back(p);
    }
}

void UGrid::Allocate(const std::vector<glm::vec3> &points){
    for(auto&c:Cells)c=0;
    LoadSlices(points);

    for(auto p:points){
        (*this)(p)++;
    }
}

void UGrid::AllocateSphericalZ(const std::vector<glm::vec3> &points){
    for(auto&c:Cells) c=0;
    LoadSlices(points);

    for(int z = 0; z < Cells.extents()[2];z++){

        for(const vec3 &p:PointsPerSliceZ[z]){
            const glm::vec3 centerPos=(BBoxMin + BBoxMax) / 2.0f;
            const float pDistToCenter = glm::distance(glm::vec2(centerPos.x, centerPos.y), glm::vec2(p.x, p.y));

            for(int x = 0; x < Cells.extents()[0];x++){

                for(int y = 0; y < Cells.extents()[1];y++){
                    ivec3 cell=ivec3(x,y,z);
                    vec3 cc=CellCenterPos(cell);
                    //dist from central voxel

                    float cDistToCenter = glm::distance(glm::vec2(centerPos.x, centerPos.y), glm::vec2(cc.x, cc.y));

                    if(glm::distance(cDistToCenter, pDistToCenter) < CellSize){
                        (*this)(cell)++;
                    }
                }
            }
        }
    }
}

void UGrid::ExportVoxelsToXYZ(const string fName){
    if(std::filesystem::exists(fName)){
        throw std::runtime_error("File already exist");
    }
    std::ofstream f(fName);

    for(int x = 0; x < Cells.extents()[0];x++){
        for(int y = 0; y < Cells.extents()[1];y++){
            for(int z = 0; z < Cells.extents()[2];z++){
                if((*this)(x,y,z)){
                    f << x << " "<< y << " "<< z << std::endl;
                }
            }
        }
    }
}
