#ifndef MESH_H
#define MESH_H
#include <iostream>
#include <vector>
/*class boundaryField
{
    //follow the arrangement: upper,left,lower,right
    std::vector<boundaryPatch> boundary; 
    Mesh& mesh;

    boundaryField()
    {}
};
class boundaryPatch
{
    std::vector<double>
}
*/
class boundaryPatch 
{
    public:
    std::vector<std::vector<int> > index;
    std::vector<std::vector<int> > neighbour;
    std::vector<double> dist;
    std::vector<double> bcValue;
    int type;
    
    boundaryPatch
    (
        const std::vector<std::vector<int> >& Index,
        const std::vector<std::vector<int> >& Neighbour,
        const std::vector<double>& Value,
        const std::vector<double>& bcDist,
        int type
    ):type(type)
    {
        for(int i=0;i<Index.size();i++)
        {
            index.push_back(std::vector<int>());
            neighbour.push_back(std::vector<int>());
            for(int j = 0; j<Index[i].size();j++)
            {
                index[i].push_back(Index[i][j]);
                neighbour[i].push_back(Neighbour[i][j]);
            }
            dist.push_back(bcDist[i]);
            bcValue.push_back(Value[i]);
        }
    }
    ~boundaryPatch(){}

};

class Mesh
{
    public:

  //  std::vector<double> dx;
  //  std::vector<double> dxInv;
    int nx;
    int ny;
    const std::vector<double>& x;
    const std::vector<double>& y;
    const std::vector<boundaryPatch>& boundary;
    std::vector<std::vector<bool> > onBoundary;
    std::vector<std::vector<int> > bcType;
    std::vector<std::vector<double> > bcValueMap;
    //constructor
    Mesh(const std::vector<double>& x, const std::vector<double>& y, int nx, int ny, std::vector<boundaryPatch>& boundary,std::vector<std::vector<bool> > boundarymap,std::vector<std::vector<int> > BC_typemap,std::vector<std::vector<double> > bcvaluemap)
    :x(x),y(y),nx(nx),ny(ny),boundary(boundary)
    {
        for(int i=0;i<ny;i++)
        {
            onBoundary.push_back(std::vector<bool>());
            bcType.push_back(std::vector<int>());
            bcValueMap.push_back(std::vector<double>());
            for(int j = 0; j< nx;j++)
            {
                onBoundary[i].push_back(boundarymap[i][j]);
                bcType[i].push_back(BC_typemap[i][j]);
                bcValueMap[i].push_back(bcvaluemap[i][j]);
            }
        }
      //  dx.push_back( (xlim.second-xlim.first)/nx);
      //  dx.push_back( (ylim.second-ylim.first)/ny);
      //  dxInv.push_back(1/dx[0]);
      //  dxInv.push_back(1/dx[1]);
    }
    

    ~Mesh(){}
};
#endif