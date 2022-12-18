#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <iostream>
#include <ostream>
#include "HalfedgeBuilder.cpp"
#include "utils.cpp"

using namespace Eigen;
using namespace std;

MatrixXd V1; // vertex coordinates of the input mesh
MatrixXi F1; // incidence relations between faces and edges

double lambda = 1.0e+1;
double epsilon = 3.0e-3;
double tau_incr = 2.0;
double tau_decr = 2.0;
double mu = 10;

// MATHEMATIC FUNCTIONS //

// Return the norm 2 of a vector : OK
double norm_2(RowVector3d x){
    return sqrt(pow(x(0), 2) + pow(x(1), 2) + pow(x(2), 2));
}

// Return the norm ||A||_W = Tr(AWA.T)
double norm_matrix(MatrixXd A, MatrixXd W){
    return (A * W * A.transpose()).trace();
}

// One solution of the lasso problem [Boyd et al. 2011; Tibshirani 1996] OK
RowVector3d shrinkage(double kappa, RowVector3d x){
    RowVector3d z;
    for(int i=0; i<3; i++){
        z(i) = (1 - kappa/abs(x(i))) + x(i);
    }
    return z;
}

// Return the cotangent of x
float cotan(float x){
    return 1 / tan(x);
}


// GEOMETRIC FUNCTIONS //

// Return the angle between three points: OK
double get_angle(MatrixXd V, int pt_ref, int pt_1, int pt_2){
    RowVector3d e_1 = V.row(pt_ref) - V.row(pt_1);
    RowVector3d e_2 = V.row(pt_ref) - V.row(pt_2);
    float x = double(e_1.dot(e_2)) / double(norm_2(e_1) * norm_2(e_2));
    return acos(x);
}

// Return the normal of a triangle: OK
MatrixXd get_face_weigh_normal(HalfedgeDS he, MatrixXd V, int e_1, int e_2){

    int pt_0 = he.getTarget(e_1);
    int pt_1 = he.getTarget(he.getOpposite(e_1));
    int pt_2 = he.getTarget(he.getOpposite(e_2));
    std::cout << "points concernés: " << pt_0 << "  " << pt_1 << "   " << pt_2 << std::endl;

    // Caclulate the angle and the normal of the triangle i, pt_1, pt_2 (with respect to i)
    double angle = get_angle(V, pt_0, pt_1, pt_2);
    std::cout << "angle: " << angle << std::endl;

    // Compute the normal with the cross product of two edges of the triangle
    RowVector3d vec_1 = V.row(pt_1) - V.row(pt_0);
    RowVector3d vec_2 = V.row(pt_2) - V.row(pt_0);
    MatrixXd face_normal = vec_2.cross(vec_1);

    face_normal.normalize();

    // Weigh the normal
    face_normal *= angle;

    std::cout<< "normal : " << face_normal << std::endl;
    return face_normal.transpose();
}

// Return the normal of a triangle: OK
RowVector3d get_face_normal(HalfedgeDS he, MatrixXd V, RowVector3d vec_1, RowVector3d vec_2){
    RowVector3d face_normal = vec_2.cross(vec_1);
    face_normal.normalize();
    std::cout<< "normal : " << face_normal << std::endl;
    return face_normal.transpose();
}

// Return the area of a triangle based on points of V
float triangle_area(MatrixXd V, RowVector3i triangle){
    float base = norm_2(V.row(triangle(0)) - V.row(triangle(1)));
    float edge_1 = norm_2(V.row(triangle(2)) - V.row(triangle(1)));
    float edge_2 = norm_2(V.row(triangle(2)) - V.row(triangle(0)));

    float x = (pow(edge_2, 2) - pow(edge_1, 2) + pow(base, 2)) / (2 * base);
    float height = sqrt(pow(edge_2, 2) - pow(x, 2));

    return (base * height) / 2;

}

// Return the barycenter of a triangle
RowVector3d get_face_barycenter(MatrixXd V, RowVector3i triangle){
    RowVector3d barycenter(1, 3);
    for(int i=0; i<3; i++){
        barycenter += V.row(triangle(i));
    }
    barycenter /=3;
    return barycenter;
}

// HALFEDGE FUNCTIONS //

int vertex_Degree(HalfedgeDS he, int v) {
    // TO BE COMPLETED
    int result = 0;
    int e = he.getEdge(v);

    int e_next = he.getNext(e);
    int p_edge = he.getOpposite(e_next);

    while(p_edge != e){
        e_next = he.getNext(p_edge);
        p_edge = he.getOpposite(e_next);
        result += 1;
    }

    return result + 1;
    }

// COMPUTATIONAL FUNCTIONS //

// Retrurn a 3 x |N(i)| matrice of rim/spoke edge vectors: OK
MatrixXd compute_distances(HalfedgeDS he, MatrixXd V, int i, int N_neighbour){
    // Create D_i
    MatrixXd D(3, N_neighbour);
    // Initialize the he data structure course
    int p_edge = he.getEdge(i);
    int e_next;
    RowVector3d dist_vect;
    // For each neighbour
    for (int j=0; j<N_neighbour; j++){
        e_next = he.getNext(p_edge);
        // Calculate the distance vector
        dist_vect = V.row(he.getTarget(e_next)) - V.row(i);
        D.block<3,1>(0, j) = dist_vect;
        p_edge = he.getOpposite(e_next);
    }
    return D;
}



// Return  the unit area-weighted normal vector of a vertex i with the he data structure: OK
MatrixXd compute_area_weighted_normal(HalfedgeDS he, MatrixXd V, int i, int N_neighbour){
    // Create n_i
    MatrixXd normal(3, 1);
    // Initialize the he data structure course
    vector<int> edge_i; // Contain the edge (in the he ds) of each face related to i
    int p_edge = he.getEdge(i);
    int e_next;
    // For each neigbhour,
    for (int i=0; i<N_neighbour; i++){
        e_next = he.getNext(p_edge);
        p_edge = he.getOpposite(e_next);
        edge_i.push_back(p_edge);
    }
    normal += get_face_weigh_normal(he, V, edge_i.at(0), edge_i.at(1));
    normal += get_face_weigh_normal(he, V, edge_i.at(1), edge_i.at(2));
    normal += get_face_weigh_normal(he, V, edge_i.at(2), edge_i.at(0));

    return normal;
}


// Return the diagonal matrix of cotangent weights: OK
MatrixXd compute_cotangent_weight(HalfedgeDS he, MatrixXd V, int i, int N_neighbour){

    // Declaration of the variables
    MatrixXd W(N_neighbour, N_neighbour);
    int p_edge = he.getEdge(i);
    float alpha, beta;
    int e_next, v_1, v_2, v_ref;

    for(int j=0; j<N_neighbour; j++){
        e_next = he.getNext(p_edge);

        v_1 = he.getTarget(p_edge);
        v_2 = he.getTarget(e_next);
        v_ref = he.getTarget(he.getOpposite(p_edge));
        alpha = get_angle(V, v_ref, v_1, v_2);

        p_edge = he.getOpposite(e_next);

        v_1 = he.getTarget(p_edge);
        v_2 = he.getTarget(e_next);
        v_ref = he.getTarget(he.getNext(p_edge));
        beta = get_angle(V, v_ref, v_1, v_2);

        W.row(j)(j) = cotan(alpha) + cotan(beta);

    }
    return W;
}

// Return the mean of the area of the adjacent faces area: OK
double compute_barycenter_area(HalfedgeDS he, MatrixXd V, MatrixXi F, int i, int N_neighbour){
    double a_i;
    int p_edge = he.getEdge(i);
    int e_next;
    RowVector3i triangle;
    int face;
    for(int i=0; i<N_neighbour; i++){
        face = he.getFace(p_edge);
        triangle = F.row(face);
        a_i += triangle_area(V, triangle);
        e_next = he.getNext(p_edge);
        p_edge = he.getOpposite(e_next);
    }
    a_i /= N_neighbour;
    return a_i;
}

// CLASS STRUCTURE //

class vertex_features{
public:
    MatrixXd R_i;
    MatrixXd D_i;
    MatrixXd D_i_tilde;
    MatrixXd normal_i;
    MatrixXd W_i;
    double a_i;

    void initialize(int N_neighbour){
        W_i = MatrixXd::Zero(N_neighbour, N_neighbour);
        D_i = MatrixXd::Zero(3, N_neighbour);
        D_i_tilde = MatrixXd::Zero(3, N_neighbour);
        normal_i = MatrixXd::Zero(3, 1);
        R_i = MatrixXd::Zero(3, 3);
        a_i = 0;
    }
    void compute_parameters(HalfedgeDS he, MatrixXd V, MatrixXd V_tilde, MatrixXi F, int i, int N_neighbour){
        D_i = compute_distances(he, V, i, N_neighbour);
        D_i_tilde = compute_distances(he, V_tilde, i, N_neighbour);
        normal_i = compute_area_weighted_normal(he, V, i, N_neighbour);
        W_i = compute_cotangent_weight(he, V, i, N_neighbour);
        a_i = compute_barycenter_area(he, V, F, i, N_neighbour);
    }


};

class Trf
{
public:
    vector<RowVector3d> z;
    vector<RowVector3d> u;
    vector<MatrixXd> R;
    double energy;
    double penalty = 10e-4;

    void local_step(HalfedgeDS he, MatrixXd V, MatrixXi F, MatrixXd V_tilde, int i, vertex_features param);

    void initialize(MatrixXd V){
        for(int i=0; i<V.rows(); i++){
            z.push_back(RowVector3d(0.0, 0.0, 0.0));
            u.push_back(RowVector3d(0.0, 0.0, 0.0));
            R.push_back(MatrixXd::Zero(3, 3));
        }
    }
};


// CUBE STYLIZATION FUNCTIONS //

double update_penalty(double penalty, RowVector3d u){
    if(norm_2(r) > mu * norm_2(s)){
        penalty = tau_incr * penalty;
    }
    else if(norm_2(s) > mu * norm_2(r)){
        penalty = penalty / tau_decr;
    }
    // else : penalty = penalty : useless to compute
    return penalty;
}



void Trf::local_step(HalfedgeDS he, MatrixXd V, MatrixXi F, MatrixXd V_tilde, int i, vertex_features param){
    // For all vertices
    energy = 0;
    for(int i=0; i<V.rows(); i++){
        int N_neighbour = vertex_Degree(he, i);

        // Compute Di ni Wi ai Di_tilde
        param.compute_parameters(he, V, V_tilde, F, i, N_neighbour);

        //std::cout<< "norme " << norm_matrix(R_i * D_i - D_i_tilde, W_i) << std::endl;
        energy += norm_matrix(param.R_i * param.D_i - param.D_i_tilde, param.W_i);

        MatrixXd U_i(3, N_neighbour + 1);
        U_i.block(0, 0, 3, N_neighbour) = param.D_i;
        U_i.block<3, 1>(0, N_neighbour) = param.normal_i;

        MatrixXd V_i(3, N_neighbour + 1);
        V_i.block(0, 0, 3, N_neighbour) = param.D_i_tilde;
        V_i.block<3, 1>(0, N_neighbour) = (z.at(i) - u.at(i)).transpose();
        //std::cout << "V_i " << V_i << std::endl;
        //std::cout << "U_i" << U_i << std::endl;

        param.R_i = V_i * U_i.transpose();
        //std::cout << "R_i" << R_i << std::endl;

        R.at(i) = param.R_i;

        // Update the parameters z, u, penalty
        z.at(i) = shrinkage(lambda * a_i / penalty, (R_i * normal_i).transpose() + u.at(i)); // as in the lasso problem [Boyd et al. 2011; Tibshirani 1996]
        // We update the penalty according to Sec. 3.4.1 in [Boyd et al.2011] where u needs to be rescaled accordingly after updating the penalty.
        update_penalty(penalty, u)
    }
}

// Transform V_tilde with the transformation R
MatrixXd global_step(vector<MatrixXd> R, MatrixXd V, MatrixXd V_tilde){
    for(int i=0; i<V.rows(); i++){
        std::cout << "R_i " << R.at(i) << std::endl;
        V_tilde.block<1, 3>(i, 0) = (R.at(i) * V.block<1, 3>(i, 0).transpose()).transpose();
    }
    return V_tilde;
}

MatrixXd cube_stylization(HalfedgeDS he, MatrixXd V, MatrixXi F){

    Trf trf;
    MatrixXd V_tilde = V;
    //MatrixXd V_temp; // will be useful to store matrix for a while
    trf.initialize(V);
    vertex_features param;

    while(trf.energy > epsilon){
        trf.energy = 0;
        trf.local_step(he, V, F, V_tilde, i, param);
        V_tilde = global_step(trf.R, V, V_tilde);
    }
    return V_tilde;
}


// OTHERS //

void draw_bounding_box(igl::opengl::glfw::Viewer &viewer, MatrixXd &V) {
  // compute the corners of the bounding box
  Vector3d m = V.colwise().minCoeff();
  Vector3d M = V.colwise().maxCoeff();

  MatrixXd V_box(8,3); // Corners of the bounding box
  MatrixXi E_box(12,2); // edges of the bounding box

  V_box <<
  m(0), m(1), m(2),
  M(0), m(1), m(2),
  M(0), M(1), m(2),
  m(0), M(1), m(2),
  m(0), m(1), M(2),
  M(0), m(1), M(2),
  M(0), M(1), M(2),
  m(0), M(1), M(2);

  E_box <<
  0, 1,
  1, 2,
  2, 3,
  3, 0,
  4, 5,
  5, 6,
  6, 7,
  7, 4,
  0, 4,
  1, 5,
  2, 6,
  7 ,3;

  viewer.append_mesh();
  viewer.data(1).add_points(V_box,RowVector3d(1,0,0));

  for (unsigned i=0;i<E_box.rows(); ++i) // Plot the edges of the bounding box
    viewer.data().add_edges
    (
      V_box.row(E_box(i,0)),
      V_box.row(E_box(i,1)),
      RowVector3d(1,0,0)
    );

}


void draw_face_normals(igl::opengl::glfw::Viewer &viewer, const MatrixXd &V, const MatrixXi &F)
{
    std::cout << "Drawing face normals " << std::endl;
    HalfedgeBuilder* builder=new HalfedgeBuilder();
    HalfedgeDS he=builder->createMesh(V.rows(), F);
    RowVector3d c = RowVector3d(1.0, 0.0, 0.0);
    MatrixXd face_normal;
    //viewer->append_mesh();
    for (int i = 0; i < F.rows(); i++)
    {

      RowVector3d barycenter = get_face_barycenter(V, F.row(i));
      std::cout << "barycenter" << barycenter << std::endl;
      RowVector3d vec_1 = V.row(F.row(i)(0))-V.row(F.row(i)(1));
      RowVector3d vec_2 = V.row(F.row(i)(0))-V.row(F.row(i)(2));
      face_normal = get_face_normal(he, V, vec_1, vec_2);
      std::cout << i << " face normal of  " << F.row(i) << ": " << face_normal << std::endl;
      viewer.data(2).add_edges(
          barycenter,
          barycenter+ face_normal,
          c);
    }
  }

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
  std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;
  if(key == 'T'){

      HalfedgeBuilder* builder=new HalfedgeBuilder();
      HalfedgeDS he=builder->createMesh(V1.rows(), F1);

      V1 = cube_stylization(he, V1, F1);

      viewer.data().clear();
      viewer.data().set_mesh(V1, F1);
  }
  return false;
}

void createTetraedron(MatrixXd &Vertices, MatrixXi &Faces)
{
    Vertices = MatrixXd(4, 3);
    Faces = MatrixXi(4, 3);

    Vertices <<
    0.0, 0., 0.,
    1.0, 0.0, 0.,
    0., 0., 1.,
    0.0, 1.0, 0.0;

    Faces <<
    0, 3, 1,
    0, 2, 3,
    0, 1, 2,
    1, 3, 2;
}

void createOctagon(MatrixXd &Vertices, MatrixXi &Faces)
{
    Vertices = MatrixXd(6, 3);
    Faces = MatrixXi(8, 3);

    Vertices << 0.0, 0.0, 1.0,
        1.000000, 0.000000, 0.000000,
        0.000000, 1.000000, 0.000000,
        -1.000000, 0.000000, 0.000000,
        0.000000, -1.000000, 0.000000,
        0.000000, 0.000000, -1.000000;

    Faces << 0, 1, 2,
        0, 2, 3,
        0, 3, 4,
        0, 4, 1,
        5, 2, 1,
        5, 3, 2,
        5, 4, 3,
        5, 1, 4;
}

int main(int argc, char *argv[]){


    createTetraedron(V1, F1);
    //igl::readOBJ("../meshes/spot.obj", V1, F1); // Load an input mesh in OFF format
    //createOctagon(V1, F1);

    igl::opengl::glfw::Viewer viewer;
    draw_face_normals(viewer, V1, F1);
    viewer.callback_key_down = &key_down;
    std::cout << "key_down ok" << std::endl;
    viewer.data().set_mesh(V1, F1);
    std::cout << "set mesh ok" << std::endl;
    viewer.launch();

}


