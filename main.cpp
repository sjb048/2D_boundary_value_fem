/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no
 * SINTEF ICT, Department of Applied Mathematics,
 * P.O. Box 124 Blindern,
 * 0314 Oslo, Norway.
 *
 * This file is part of TTL.
 *
 * TTL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * TTL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with TTL. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using TTL.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the TTL library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT.
 * g++ -std=c++11 -Wno-c++11-extensions \
    -I./TTL/include \
    -I/opt/homebrew/opt/eigen/include/eigen3 \
    main.cpp \
    ./TTL/src/utils/HandleId.cpp \
    ./TTL/src/halfedge/HeTriang.cpp \
    -o main

 */
#include <iostream>
#include <ttl/ttl.h>
#include <ttl/halfedge/HeTriang.h>
#include <ttl/halfedge/HeDart.h>
#include <ttl/halfedge/HeTraits.h>

using namespace hed; // (to avoid using prefix hed::)
#define _USE_MATH_DEFINES


#include <cmath>
#include <fstream>
#include<array>
#include<vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>
using namespace Eigen;
using namespace std;


// double M_PI = 3.14159265358979323846;

class FEMobject{

public:
    vector<Node*>* nodes = new vector<Node*>;
    using Index = int; // or as appropriate

    string problemtype; //specific problem
    Triangulation triang; //mesh
    SparseMatrix<double> A; //stiffness matrix
    SparseMatrix<double> M; //mass matrix
    SparseMatrix<double> R; //Robin matrix
    VectorXd b; //load vector
    VectorXd r; //Robin vector

    int lambda = 81;

    void setproblemtype(string problem){

        problemtype = problem;

    }
void circleMesh(int n, int m, double r){

    // Initialize with three nodes forming a non-collinear triangle
    for (int i = 0; i < 3; i++){
        double angle = (i * 2 * M_PI) / 3.0; // Equally spaced angles
        double radius = (r * 1.0) / n; // Radius for the first circle
        double x = radius * cos(angle);
        double y = radius * sin(angle);
        Node* p = new Node(x, y);
        nodes->push_back(p);
    }

    // Insert remaining nodes on concentric circles
    for (int k = 1; k < n; k++){ // Start from k=1 since k=0 is the first circle
        for (int i = 0; i < (m * (k + 1)); i++ ){
            double radius = (r * (k + 1)) / n;
            double angle = (i * 2 * M_PI) / (m * (k + 1));

            double x = radius * cos(angle);
            double y = radius * sin(angle);

            Node* p = new Node(x, y);
            nodes->push_back(p);
        }
    }

    // Print number of nodes for verification
    std::cout << "Total nodes: " << nodes->size() << std::endl;

    // Create Delaunay triangulation
    triang.createDelaunay(nodes->begin(), nodes->end());
}

void squareMesh(int n, double d, Vector2d Op){

    for (int j = 0;  j<= n; j++){
            for (int i = 0;  i <= n ; i++){

                Node* p = new Node(Op(0) + ((i * d)/n), Op(1) + ((j*d)/n));

                nodes->push_back(p);
            }

        }

    triang.createDelaunay(nodes->begin(), nodes->end());
}
std::array<Vector3d, 2> gradients(Vector3d x, Vector3d y, double area){
        Vector3d b((y(1)-y(2))/(2*area), (y(2)-y(0))/(2*area), (y(0)-y(1))/(2*area));
        Vector3d c((x(2)-x(1))/(2*area), (x(0)-x(2))/(2*area), (x(1)-x(0))/(2*area));
        return {b, c};
    }

double triarea(Node* N1, Node* N2, Node* N3){
    Vector3d a;
    a << N2->x() - N1->x(), N2->y() - N1->y(), 0;
    Vector3d b;
    b << N3->x() - N1->x(), N3->y() - N1->y(), 0;
    double area = 0.5 * (b.cross(a)).norm();
    return area;
}

SparseMatrix<double> stiffMat(list<Edge*> trilist, int np){
    trilist = triang.getLeadingEdges();
    SparseMatrix<double> A(np,np);
    A.setZero();
        // list <Edge*>::iterator K;

    for (auto K = trilist.begin(); K != trilist.end(); K++){

            Edge* edg = *K;
            Node* N1 = edg->getSourceNode();
            Node* N2 = edg->getTargetNode();
            Node* N3 = edg->getNextEdgeInFace()->getTargetNode();
           
           
            // ...and their global indices
            Vector3i loc2glb;
            loc2glb << N1->id(), N2->id(), N3->id();
            

            // extract coordinates of the nodes
            Vector3d x;
            x << N1->x(), N2->x(), N3->x();
            Vector3d y;
            y << N1->y(), N2->y(), N3->y();
            
            double area = triarea(N1,N2,N3); // compute the triangle area
            // Vector2<Vector3d> gradphi = gradients(x, y, area); // compute the gradients
            std::array<Vector3d, 2> gradphi = gradients(x, y, area); // compute the gradients

            // compute AK
            Vector3d b = gradphi[0];
            Vector3d c = gradphi[1];

            Matrix3d AK = (b*b.transpose() + c*c.transpose())*area;

            // insert AK to the global stiffness matrix A
            // Insert AK into global stiffness matrix A
            for (int i = 0; i < 3; i++){
                for (int j = 0; j < 3; j++){
                    A.coeffRef(loc2glb(i), loc2glb(j)) += AK(i, j);
                }
            }
    
        }
    return A;
}

    SparseMatrix<double> massMat(list<Edge*> trilist, int np){
        SparseMatrix<double> M(np, np); // initialization of the mass matrix M
        M.setZero();

        trilist = triang.getLeadingEdges();
        for (auto K = trilist.begin(); K != trilist.end(); K++){ // loop over triangles
            Edge* edg = *K; // extract the edge of the triangle
            // extract nodes...
            Node* N1 = edg->getSourceNode();
            Node* N2 = edg->getTargetNode();
            Node* N3 = edg->getNextEdgeInFace()->getTargetNode();
            // ...and their global indices
            Vector3i loc2glb;
            loc2glb << N1->id(), N2->id(), N3->id();

            // Vector3< Index> loc2glb = {N1->id(), N2->id(), N3->id()};
            // extract coordinates of the nodes
            Vector3d x;
            x << N1->x(), N2->x(), N3->x();
            Vector3d y;
            y << N1->y(), N2->y(), N3->y();
            double area = triarea(N1,N2,N3); // compute the triangle area
            // compute MK
            Matrix3d mat;
            mat << 2, 1, 1,
                    1, 2, 1,
                    1, 1, 2;
            Matrix3d MK = (1.0 / 12.0) * mat * area;

            SparseMatrix<double> N(np, np);
            for (int i = 0; i < 3; i++){
                for (int j = 0; j < 3; j++){
                    M.coeffRef(loc2glb(i), loc2glb(j)) += MK(i, j);
                }
            }
        }
        return M;
    }

    VectorXd loadVect(list<Edge*> trilist, int np){
        VectorXd b = VectorXd::Zero(np); // initialization of the load vector b
       
        trilist = triang.getLeadingEdges();

        for (auto K = trilist.begin(); K != trilist.end(); K++){ // loop over triangles
            Edge* edg = *K; // extract the edge of the triangle
            // extract nodes...
            Node* N1 = edg->getSourceNode();
            Node* N2 = edg->getTargetNode();
            Node* N3 = edg->getNextEdgeInFace()->getTargetNode();
            // ...and their global indices
            Vector3i loc2glb;
            loc2glb << N1->id(), N2->id(), N3->id();

            // extract coordinates of the nodes
            Vector3d x;
            x << N1->x(), N2->x(), N3->x();
            Vector3d y;
            y << N1->y(), N2->y(), N3->y();
            double area = triarea(N1,N2,N3); // compute the triangle area

            // compute bK
            Vector3d F;
            F << f(x(0), y(0)), f(x(1), y(1)), f(x(2), y(2));
            Vector3d bK = (F / 3.0) * area;
            //add bK to the global load vector b
            for (int j = 0; j < 3; j++){
                b(loc2glb(j)) += bK(j);
            }

        }
        return b;
    }

    SparseMatrix<double> robinMat(list<Dart> boundary, int np){
        SparseMatrix<double> R(np, np);
        R.setZero();

        Edge* edge = triang.getBoundaryEdge(); //returns an arbitrary boundary edge
        Dart b_dart(edge); //creates a boundary dart
        //creates a list of boundary darts
        ttl::getBoundary(b_dart, boundary);
       

        // initialization of the boundary matrix R
        for (auto E = boundary.begin(); E != boundary.end(); E++){ // loop over boundary edges
            // extract nodes...
            Node* N1 = E->getNode();
            Node* N2 = E->getOppositeNode();
            // ...and their global indices
            Vector2i loc2glb;
            loc2glb << N1->id(), N2->id();

            // Vector2<Index> loc2glb = {N1->id(), N2->id()};
            // extract coordinates of the nodes
            Vector2d x;
            x << N1->x(), N2->x();
            Vector2d y;
            y << N1->y(), N2->y();
            // compute the length of the edge
            double len = (x - y).norm();
               
            // find edge centroid
            double xc = (x(0) + x(1)) / 2.0;
            double yc = (y(0) + y(1)) / 2.0;
            double k = kappa(xc, yc); // compute the value of κ at centroid
            // compute RE
            Matrix2d mat;
            mat << 2, 1,
                   1, 2;
            Matrix2d RE = (k / 6.0) * mat * len;

            SparseMatrix<double> S(np, np);
            for (int i = 0; i < 2; i++){
                for (int j = 0; j < 2; j++){
                    R.coeffRef(loc2glb(i), loc2glb(j)) += RE(i, j);
                }
            }

        }
        return R; }

    VectorXd RobinVect(list<Dart> boundary, int np){
        VectorXd r =  VectorXd::Zero(np);
        Edge* edge = triang.getBoundaryEdge(); //returns an arbitrary boundary edge
        Dart b_dart(edge); //creates a boundary dart
        //creates a list of boundary darts
        ttl::getBoundary(b_dart, boundary);
        
        // initialization of the boundary vector r
        for ( auto E = boundary.begin(); E != boundary.end(); E++ ){ // loop over boundary edges
            // extract nodes...
            Node* N1 = E->getNode();
            Node* N2 = E->getOppositeNode();
            // ...and their global indices
            Vector2i loc2glb;
            loc2glb << N1->id(), N2->id();

            // Vector2< Index> loc2glb = {N1->id(), N2->id()};
            // extract coordinates of the nodes
            Vector2d x = {N1->x(), N2->x()};
            Vector2d y = {N1->y(), N2->y()};
            // compute the length of the edge
            double len = sqrt((x(0) - x(1))*(x(0) - x(1)) + (y(0) - y(1))*(y(0) - y(1)));
            // find edge centroid
            double xc = (x(0) + x(1))/2.0;
            double yc = (y(0) + y(1))/2.0;
            double tmp = kappa(xc,yc)*gD(xc,yc)+gN(xc,yc); // compute the value of the boundary conditions
            
            Vector2d vec;
            vec << 1.0, 1.0;
            Vector2d rE = (tmp * vec * len) / 2.0;
            //add rE to the global boundary vector r
            for (int j = 0; j < 2; j++){
                r(loc2glb(j)) += rE(j);
            }
        }
        return r;
    }

    double kappa(double x, double y){
        if(this->problemtype == "Laplace") {
                return 1e6;
            } else if(this->problemtype == "poisson") {
                return 1e6;
            }else if(this->problemtype == "HelmHoltz") {
                if (x > 0.0)
                    return 0.0;
                else
                    return 1e6;
            }
            return 0; //default value
    }
    double gN (double x, double y){
        return 0;
    }
    double gD(double x, double y){
         if (this->problemtype == "Laplace" ){
            double phi = atan2(y,x);
            return cos(4*phi);
        } else if (this->problemtype == "poisson" ){
            return (y*y)/2;
        } else if (this->problemtype == "HelmHoltz" ){
            return 0.25;
        } else {
            return 0.0; // Default return value
        }
    }
    double f(double x, double y){
        return 1;
    }
    double uexact(double x, double y){
        if (this->problemtype == "Laplace"){
            double rho = sqrt(x * x + y * y);
            double phi = atan2(y, x);
            return pow(rho, 4) * cos(4 * phi);
        } else if(this->problemtype == "poisson") {
            return (1.0 - pow(x, 2)) / 2.0;
        }
        else if(this->problemtype == "HelmHoltz"){
            return 0.25 * (cos(sqrt(lambda) * x) + tan(sqrt(lambda)) * sin(sqrt(lambda) * x));
        } else{
            return 0.0;
        }
    }
    // 
    double getError(){
        double error = 0.0;
        list<hed::Edge*> trilist = triang.getLeadingEdges();
        
        for (auto K = trilist.begin(); K != trilist.end(); K++){ //loop over triangles
            Edge* edge = *K;
            Node* nodea = edge->getSourceNode();
            Node* nodeb = edge->getTargetNode();
            Node* nodec = edge->getNextEdgeInFace()->getTargetNode();
            Vector3d x;
            x << nodea->x(), nodeb->x(), nodec->x();
            Vector3d y;
            y << nodea->y(), nodeb->y(), nodec->y();
            double area = triarea(nodea, nodeb, nodec);
            double xc = (x(0) + x(1) + x(2)) / 3.0;
            double yc = (y(0) + y(1) + y(2)) / 3.0;
            double ubar = uexact(xc, yc);
            double uhbar = (nodea->z() + nodeb->z() + nodec->z()) / 3.0;
            double eK = pow((ubar - uhbar), 2) * area;
            error += eK;
        }
        return sqrt(error);
    }

    int getDoFs(){
        return this->triang.getNodes()->size(); //returns #DoF
    }

    void solve(){
        
        list<hed::Edge*> trilist = triang.getLeadingEdges();
        Edge* edge = triang.getBoundaryEdge(); //returns an arbitrary boundary edge
        Dart b_dart(edge); //creates a boundary dart
        list<Dart> boundary;
        list<Node*>* nodelist = triang.getNodes(); //returns a list of nodes
        int np = nodelist->size();
        list<Node*>::iterator L;
        // Debug: Print node IDs and positions
        for(auto node : *nodelist){
            // std::cout << "Node ID: " << node->id() << " Position: (" << node->x() << ", " << node->y() << ")\n";
        }

        A = stiffMat( trilist, np); //returns the stiffness matrix A
        R = robinMat(boundary, np); //returns the Robin matrix R
        M = massMat(trilist, np);
        r = RobinVect(boundary, np);
        b = loadVect(trilist, np);
        VectorXd zeta(np);
        SimplicialLDLT<SparseMatrix<double> > solver;
        if (this->problemtype == "Laplace"){

            zeta = solver.compute(A + R).solve(r);

        }else if(this->problemtype == "poisson"){

            zeta = solver.compute(A + R).solve( r + b);

        }else if(this->problemtype == "HelmHoltz"){
             zeta = solver.compute(A + R - 81 * M).solve( r);
        }

        for (L = nodelist->begin(); L != nodelist->end(); L++){ //loop over nodes
            Node* node = *L;
            node->init(node->x(),node->y(), zeta(node->id())); //set the node values
        }


    }
    void visualization(string filename){
        ofstream objfile;
        objfile.open(filename);


        list<hed::Node*>* nodelist = triang.getNodes(); //returns a list of nodes
        
        for (auto L = nodelist->begin(); L != nodelist->end(); L++){ //loop over nodes
            hed::Node* node = *L;
            objfile << "v " << node->x() << " " << node->z()  << " " << node->y() << "\n";
        }

        list<hed::Edge*> trilist = triang.getLeadingEdges(); //returns a list of triangles
        // Generate normals
        for (auto K = trilist.begin(); K != trilist.end(); K++){
            hed::Edge* edge = *K;
            hed::Node* N1 = edge->getSourceNode();
            hed::Node* N2 = edge->getTargetNode();
            hed::Node* N3 = edge->getNextEdgeInFace()->getTargetNode();
            Vector3d a(N2->x() - N1->x(), N2->y() - N1->y(), N2->z() - N1->z());
            Vector3d b(N3->x() - N1->x(), N3->y() - N1->y(), N3->z() - N1->z());
            Vector3d normal = a.cross(b);
            objfile << "vn " << normal.x() << " " << normal.z() << " " << normal.y() << "\n";
        }


        int normind = 1; //counts the index of a normal vector
        
        for (auto T = trilist.begin(); T != trilist.end(); T++){ //loop over triangles
            hed::Edge* edge = *T;
            hed::Node* N1 = edge->getSourceNode();
            hed::Node* N2 = edge->getTargetNode();
            hed::Node* N3 = edge->getNextEdgeInFace()->getTargetNode();

            int np = nodelist->size();
            objfile << "f " << N1->id()-np-3 << "//" << normind << " " <<
            N2->id()-np-3 << "//" << normind << " " <<
            N3->id()-np-3 << "//" << normind << "\n"; // f v1//vn1 v2//vn2 v3//vn3
        }
        objfile.close();
    }
    void writeNodeErrors(const std::string& filename) {
        std::ofstream outfile(filename);
        if (!outfile.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
            return;
        }
        
        outfile << "NodeID,X,Y,Error\n"; // Write the CSV header
        
        for (auto node : *nodes) { // Iterate over all nodes
            if (!node) continue; // Check for null pointers
            double exact = uexact(node->x(), node->y());
            double numerical = node->z(); // Assuming node->z() holds the FEM solution
            double error = exact - numerical;
            
            outfile << node->id() << "," << node->x() << "," << node->y() << "," << error << "\n";
        }
        
        outfile.close();
        // std::cout << "Node errors successfully written to " << filename << std::endl;
    }

};

int main() {

    FEMobject model;
    Vector2d Op(0, -0.5)  ;

     model.squareMesh(40, 1, Op );

    // model.circleMesh(1, 3, 1.0 );
    // model.circleMesh(5, 12, 1.0);


     model.setproblemtype("poisson");
    // model.setproblemtype("HelmHoltz");
    // model.setproblemtype("Laplace");


    model.solve();

    model.visualization("Model.obj");

    auto error = model.getError();

    auto dof = model.getDoFs();
    

    cout<< "\n Error : "  << error;

    cout<< "\n Dof : " << dof;
    // Save node errors to CSV
    // std::ofstream outfile("node_errors_square_Laplace.txt");
    // std::ofstream outfile("node_errors_square_HelmHoltz.txt");
    std::ofstream outfile("node_errors_square_poisson.txt");
    if (outfile.is_open()) {
        outfile << "Error : " << error << "\n";
        outfile << "Dof : " << dof << "\n";
        outfile.close();
    } else {
        std::cerr << "Error: Unable to open results.txt for writing." << std::endl;
    }
    // model.writeNodeErrors("node_errors_circle_Laplace.txt");

    return 0;
}

