// geometry_approximate.cpp

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <cstdlib>
#include <iostream>
#include <math.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

using namespace Eigen;
using namespace std;

namespace py = pybind11;

void print_arr(double *p, const int &ncoords){
    for(int i=0; i<ncoords*3; i=i+3)
        cout<<*p<<"\t"<<*(p+1)<<"\t"<<*(p+2)<<endl;
}

MatrixXd geometry_approximate(MatrixXd &coords, const int &ncoords){

    // Center of hypotenuse:
    // Apexes for nex triangle search.
    MatrixXd Apex(ncoords,3);
    double *apex=Apex.data();
    // double Apex[ncoords][3]={0};
    // double *apex=NULL;
    // apex = Apex[0];
    
    Matrix<double, 1, 7> Nvec;
    double *nvec = Nvec.data();
    for(int i=0; i<7; i++)
        *(nvec+i)=0;
    
    //
    double p1[3]={0}, p2[3]={0}, p3[3]={0}, papex[3]={0};
    //
    Vector3d v1(0,0,0), v2(0,0,0), nv(0,0,0);
    double *_v1=v1.data(), *_n = nv.data(), *_v2=v2.data();
    
    Vector3d vcurr(0,0,0), vprev(0,0,0), nr(0,0,0);
    double *_nr = nr.data();
    
    // // Store x,y,z data:
    // const int length=ncoords*5; 
    // double xs[length]={0}, ys[length]={0}, zs[length]={0};
    
    // loop to find all triangles
    int ntags = ncoords/2 + (ncoords%2); // Number of triangles current points have.
    // Number of atoms before each search.
    int _nps = ncoords;
    double *coord=coords.data();
    
    int total_tags = 0;

    for(int nm=0; nm<1000; nm++){
            
        int i=0, k=0, nps = _nps; 
        _nps = 0;
        
        // if(nm==0){
        //     cout<<coords<<endl;   
        // }
        // else{
        //     cout<<"Round "<<nm<<endl;
        //     cout<<Apex<<endl;
        // }
        
        // cout<<"n triangle apexes:\t"<<nps<<endl;
        // Make all points be a closed shape.
        /*      The center of hypotenuse for last two and first 1 points. 
                When num of points is even, the last bond and first bond are seen
            as the two sides of triangle. 
                When # of points is odd, the last bond can be cut to 2 bonds and 
            become the two side of triangle.
        */
        
        // cout<<"\nNumber of apexes: "<<nps<<endl;
        double x0=*(coord), y0=*(coord+ncoords), z0=*(coord+2*ncoords);
        double x1_=*(coord+nps-1), y1_=*(coord+ncoords+nps-1), z1_=*(coord+2*ncoords+nps-1);
        double x2_=*(coord+nps-2), y2_=*(coord+ncoords+nps-2), z2_=*(coord+2*ncoords+nps-2);
        
        /*
            The break out of this for loop.
            When 
        */
        if(nps==3){
            vprev<<x0-x1_, y0-y1_, z0-z1_;
            vcurr<<x2_-x1_, y2_-y1_, z2_-z1_;
            
            nr = vprev.cross(vcurr);
            *(nvec) += *_nr;
            *(nvec+1) += *(_nr+1);
            *(nvec+2) += *(_nr+2);
            *(nvec+3) += abs(*_nr)/2;
            *(nvec+4) += abs(*(_nr+1))/2;
            *(nvec+5) += abs(*(_nr+2))/2;
            *(nvec+6) += nr.norm()/2;
            total_tags+=1;
            // cout<<"Total triangles:"<<total_tags<<endl;
            break;
        }
        
        if(nps==4){
            
            double x3_=*(coord+nps-3), y3_=*(coord+ncoords+nps-3), z3_=*(coord+2*ncoords+nps-3);
            
            vprev<<x0-x3_, y0-y3_, z0-z3_;
            vcurr<<x2_-x3_, y2_-y3_, z2_-z3_;
            nr = vprev.cross(vcurr);
            
            v1<<x0-x1_, y0-y1_, z0-z1_;
            v2<<x2_-x1_, y2_-y1_, z2_-z1_;
            nv = v1.cross(v2);
            
            *(nvec) += (*_n + *_nr);
            *(nvec+1) += (*(_n+1) + *(_nr+1));
            *(nvec+2) += (*(_n+2) + *(_nr+2));
            *(nvec+3) += (abs(*_n)/2 + abs(*_nr)/2);
            *(nvec+4) += (abs(*(_n+1))/2 + abs(*(_nr+1))/2);
            *(nvec+5) += (abs(*(_n+2))/2 + abs(*(_nr+2))/2);
            *(nvec+6) += (nv.norm() + nr.norm())/2;
            
            total_tags+=2;
            // cout<<"Total triangles:"<<total_tags<<endl;
            break;
        } 
        
        if(nps%2==0){
            *(apex+k) = (x2_ + x0)/2;
            *(apex+k+ncoords) = (y2_ + y0)/2;
            *(apex+k+2*ncoords) = (z2_ + z0)/2;
            _nps = _nps+1;
            // Calculate normal vector and area of triangle:
            *_v1 = x2_-x1_, *(_v1+1) = y2_-y1_, *(_v1+2)=z2_-z1_;
            *_v2 = x0 -x1_, *(_v2+1) = y0 -y1_, *(_v2+2)=z0 -z1_;
            nv = v1.cross(v2);
            
            *(nvec) += *_n;
            *(nvec+1) += *(_n+1);
            *(nvec+2) += *(_n+2);
            *(nvec+3) += abs(*_n)/2;
            *(nvec+4) += abs(*_n+1)/2;
            *(nvec+5) += abs(*_n+2)/2;
            *(nvec+6) += nv.norm()/2;
            total_tags+=1;
            
        }else{ // In this case, the angle between two bonds is 0 so that area is 0.
            *(apex+k) = (x1_ + x0)/2;
            *(apex+k+ncoords) = (y1_ + y0)/2;
            *(apex+k+2*ncoords) = (z1_ + z0)/2;
            _nps = _nps+1;
        }
        

        // nps=1 means above add of points.
        // When ncoords is odd, the i+2 in loop can reach npoints-1,
        // When ncoords is even, the i+2 in loop can reach npoints-2.
        while(i<nps-2){
            *p1 = *(coord+i), *(p1+1)=*(coord+i+ncoords), *(p1+2)=*(coord+i+2*ncoords);
            *p2 = *(coord+i+1), *(p2+1)=*(coord+i+ncoords+1), *(p2+2)=*(coord+i+2*ncoords+1);
            *p3 = *(coord+i+2), *(p3+1)=*(coord+i+ncoords+2), *(p3+2)=*(coord+i+2*ncoords+2);
            
            // Calculate center of hypotenuse:
            // The center of hypotenuse as the apex of triangle for next split.
            *(apex+k+1)=(*p1+*p3)/2, 
            *(apex+k+ncoords+1)=(*(p1+1)+*(p3+1))/2, 
            *(apex+k+2*ncoords+1)=(*(p1+2)+*(p3+2))/2;
            _nps = _nps+1;
            
            vprev<<*(apex+k)-*(p1), *(apex+k+ncoords)-*(p1+1), *(apex+2*ncoords+k)-*(p1+2);
            vcurr<<*(apex+k+1)-*(p1), *(apex+1+k+ncoords)-*(p1+1), *(apex+2*ncoords+k+1)-*(p1+2);
            nr = vprev.cross(vcurr); // Normal vector remain.
            
            
            // Calculate normal vector of triangle:
            *_v1 = *p1-*p2, *(_v1+1) = *(p1+1)-*(p2+1), *(_v1+2)=*(p1+2)-*(p2+2);
            *_v2 = *p3-*p2, *(_v2+1) = *(p3+1)-*(p2+1), *(_v2+2)=*(p3+2)-*(p2+2);
            nv = v1.cross(v2);
            
            *(nvec) += (*_n + *_nr);
            *(nvec+1) += (*(_n+1) + *(_nr+1));
            *(nvec+2) += (*(_n+2) + *(_nr+2));
            *(nvec+3) += (abs(*_n)/2 + abs(*_nr)/2);
            *(nvec+4) += (abs(*(_n+1))/2 + abs(*(_nr+1))/2);
            *(nvec+5) += (abs(*(_n+2))/2 + abs(*(_nr+2))/2);
            *(nvec+6) += (nv.norm() + nr.norm())/2;
            total_tags+=2;

            i = i+2;
            k = k+1;
        }
        
        // Count remaining triangles.
        if(nps%2==0){

            // Calculate normal vector and area of triangle:
            vprev<<*(apex+k)-x2_, *(apex+k+ncoords)-y2_, *(apex+2*ncoords+k)-z2_;
            vcurr<<*(apex)-x2_, *(apex+ncoords)-y2_, *(apex+2*ncoords)-z2_;
            
        }else{ 
              
            vprev<<*(apex+k)-x1_, *(apex+k+ncoords)-y1_, *(apex+2*ncoords+k)-z1_;
            vcurr<<*(apex)-x1_, *(apex+ncoords)-y1_, *(apex+2*ncoords)-z1_;
            
        }
        
        nr = vprev.cross(vcurr);
        *(nvec) += *_nr;
        *(nvec+1) += *(_nr+1);
        *(nvec+2) += *(_nr+2);
        *(nvec+3) += abs(*_nr)/2;
        *(nvec+4) += abs(*_nr+1)/2;
        *(nvec+5) += abs(*_nr+2)/2;
        *(nvec+6) += nr.norm()/2;
        total_tags+=1;
        
        coord=apex;
    }
    return Nvec;
}

PYBIND11_MODULE(geometry_approximation_cpp, m){
    m.doc() = "pylmp geometry approximation plugin";
    m.def("geometry_approximate", &geometry_approximate, "Compute the Area of geometry shape for many points.");
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}