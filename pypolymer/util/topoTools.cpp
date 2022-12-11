// example.cpp

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <cstdlib>
#include <iostream>
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

using namespace Eigen;
using namespace std;

namespace py = pybind11;

float add(float a, float b){
    return a+b;
}

Array3f compute_rcm(const MatrixXd &coords){
    Array3f Rcm;
    Rcm << 0, 0, 0;
    
    float *rcm = Rcm.data();
    const double *crd = coords.data();
    const int ncoords = coords.rows();
    const int *ncd = &ncoords;
    
    for(int i=0; i<ncoords; i++){
        *rcm = *rcm + *(crd+i);
        *(rcm+1) = *(rcm+1) + *(crd+i+*ncd);
        *(rcm+2) = *(rcm+2) + *(crd+i+*ncd+*ncd);
    }
    
    *rcm = *rcm/(*ncd);
    *(rcm+1) = *(rcm+1)/(*ncd);
    *(rcm+2) = *(rcm+2)/(*ncd);
    
    return Rcm;
}

Matrix3d compute_gyration(const MatrixXd &coords_to_rcm){
    
    Matrix3d gyration;
    gyration << 0, 0, 0, 
                0, 0, 0,
                0, 0, 0;
    
    
    double *gyra=gyration.data();
    const double *coord=coords_to_rcm.data();
    
    double a=0, b=0, c=0;
    double *x=&a, *y=&b, *z=&c;
    
    const int ncoords = coords_to_rcm.rows();

    for(int i=0; i<ncoords; i++){
        
        *x = *(coord+i), *y=*(coord+i+ncoords), *z=*(coord+i+2*ncoords);
        
        *gyra += *x*(*x), *(gyra+1) += *y**x, *(gyra+2) += *z**x;
        *(gyra+3) += *y**x, *(gyra+4) += *y**y, *(gyra+5) += *y**z;
        *(gyra+6) += *x**z, *(gyra+7) += *y**z, *(gyra+8) += *z**z;
    }
    
    *gyra = *gyra/ncoords, *(gyra+1) = *(gyra+1)/ncoords, *(gyra+2) = *(gyra+2)/ncoords;
    *(gyra+3) = *(gyra+3)/ncoords, *(gyra+4) = *(gyra+4)/ncoords, *(gyra+5) = *(gyra+5)/ncoords;
    *(gyra+6) = *(gyra+6)/ncoords, *(gyra+7) = *(gyra+7)/ncoords, *(gyra+8) = *(gyra+8)/ncoords;
    
    return gyration;
}

MatrixXd parse_boundary(const MatrixXd &coords, const MatrixXd &bonds, const MatrixXd &box){
    
    int ncoords = coords.rows();
    const int nbonds = bonds.rows();
    
    MatrixXd new_coords(ncoords,3);
    new_coords << coords;
    double *n_cd = new_coords.data();
    
    double p0_[3]={0.0}, p1_[3]={0.0};
    int bond_[2] = {0};
    int pid_[ncoords] = {0};
    double *p0=p0_, *p1=p1_;
    int *bond=bond_, *pid=pid_;
    
    const double *cd=coords.data();
    const double *bd=bonds.data();
    
    double vec_[3] = {0, 0, 0};
    double *vec = vec_;
    
    for(int i=0; i<nbonds; i++){
        
        *bond = (int)*(bd+i);
        *(bond+1) = (int)*(bd+i+nbonds);
        
        int aid1=*bond, aid2=*(bond+1);
        
        *p0 = *(n_cd+aid1+0);
        *(p0+1) = *(n_cd+aid1+ncoords);
        *(p0+2) = *(n_cd+aid1+2*ncoords);

        *p1 = *(cd+aid2);
        *(p1+1) = *(n_cd+aid2+ncoords);
        *(p1+2) = *(n_cd+aid2+2*ncoords);
        
        *(pid+i+1) = *(bond+1);

        *vec = *p1-*p0;
        *(vec+1) = *(p1+1)-*(p0+1);
        *(vec+2) = *(p1+2)-*(p0+2);
        
        bool State = true;
        for(int j=0; j<i; j++){
            if(*(bond+1) == *(pid+j)){
                State = false;
                break;
            }
        }
                 
        if(State){
            for(int j=0; j<3; j++){
                float span = abs(*(vec+j))-box(j)/2;
                if (span>=0){
                    *(p1+j) = -(*(vec+j)/abs(*(vec+j)))*box(j) + *(p1+j);
                }
            }
        }
         
        *(n_cd+aid2) = *p1;
        *(n_cd+aid2+ncoords) = *(p1+1);
        *(n_cd+aid2+2*ncoords) = *(p1+2); 
    }
    
    return new_coords;
}


PYBIND11_MODULE(topo_tools_cpp, m){
    m.doc() = "pybind11 example plugin";
    m.def("parse_boundary", &parse_boundary, "A function which restores polymer's topology which is broken by periodical boundary.");
    m.def("compute_gyration", &compute_gyration, "A function which computes the gyration tensor of polymer.");
    m.def("compute_rcm", &compute_rcm, "A function which computes the center of mass of polymer.");
    m.def("add", &add);
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}




