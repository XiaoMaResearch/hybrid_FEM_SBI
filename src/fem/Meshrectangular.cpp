//
//  Meshrectangular.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/5/18.
//
//

#include "Meshrectangular.hpp"

void Meshrectangular(double x_min, double x_max ,double y_min, double y_max, double dx ,double dy ,int nx, int ny , MatrixXd &Node ,MatrixXi_rm &Element)
{
    for (int i =0; i<=ny;i++)
    {
        for (int j= 0; j<=nx; j++)
        {
            Node((nx+1)*i+j,0)=x_min+j*dx;
            Node((nx+1)*i+j,1)=y_min+i*dy;
        }
    }
    for (int i =0; i<ny;i++)
    {
        for (int j= 0; j<nx; j++)
        {
            double Element_temp = nx*i+j;
            Element.row(nx*i+j) << Element_temp+i ,Element_temp+i+1 , Element_temp+i+nx+2, Element_temp+i+nx+1;
        }
    }
}
