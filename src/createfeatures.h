#ifndef CREATEFEATURES_H
#define CREATEFEATURES_H

#include <glm/fwd.hpp>
#include <vector>
#include <boost/numeric/ublas/tensor.hpp>
#include <glm/vec3.hpp>

struct Cell{
    std::vector<glm::vec3> Points;
};

struct UGrid{
    using vv3tensor  = boost::numeric::ublas::tensor<std::vector<glm::vec3>>;
    using itensor  = boost::numeric::ublas::tensor<int>;
    itensor Cells{2,2,2};
    float CellSize=0.5;

    glm::vec3 BBoxMin;
    glm::vec3 BBoxMax;

    std::vector<std::vector<glm::vec3>> PointsPerSliceZ;
    std::vector<glm::vec3> Pts;

    UGrid();
    UGrid(float cellSize, glm::vec3 bboxMin, glm::vec3 bboxMax);

    ~UGrid();

    glm::vec3 CenterBBox();

    glm::ivec3 CellCoord(glm::vec3 p);
    glm::vec3 CellCenterPos(glm::ivec3 c);

    int& operator()(glm::vec3 p);
    int& operator()(glm::ivec3 ci);
    int& operator()(int x, int y, int z);

    void LoadSlices(const std::vector<glm::vec3> &points);

    void Allocate(const std::vector<glm::vec3> &points);

    void AllocateSphericalZ(const std::vector<glm::vec3> &points);

    void ExportVoxelsToXYZ(const std::string fName);
};


struct UGridAndPoints{
    UGrid G;
    std::vector<glm::vec3> Pts;
    std::string FPath;
    UGrid::itensor Ones;
    UGrid::itensor Zeros;

    struct Features{

        int UnionRevolutionAndOrg;
        int IntersectionRevolutionAndOrg;
        double IOU_RevolutionAndOrg;
        std::vector<int> NumPointsPerSliceZ;
        std::vector<float> RadiousPerSliceZ;
        std::vector<float> AvgDistPerSliceZ;
        std::vector<float> IOUs_RevolutionAndOrgPerSlice;
        std::vector<float> NumPointsFromVerticalLine;
    };
    Features Feat;

    UGridAndPoints();

    ~UGridAndPoints();

    void ReadFileToGrid(const std::string fpath, bool swapYZ, float cellSize);

    UGrid::itensor CalcUnion(UGrid::itensor X, UGrid::itensor Y);
    UGrid::itensor CalcIntersection(UGrid::itensor X, UGrid::itensor Y);

    void CalculateFeatures();
};



#endif // CREATEFEATURES_H
