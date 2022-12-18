#include <igl/opengl/glfw/Viewer.h>
#include <ostream>

using namespace Eigen;

/**
 * A class for representing linear transformations on 3D points (using homogeneous coordinates)
 * */
class MeshTransformation
{
public:
/*
Initialize the identity transformation
**/
  MeshTransformation()
  {
    MatrixXd m(4, 4);
    m(0, 0) = 1.0; m(1, 0) = 0.0; m(2, 0) = 0.0; m(3, 0) = 0.0;
    m(0, 1) = 0.0; m(1, 1) = 1.0; m(2, 1) = 0.0; m(3, 1) = 0.0;
    m(0, 2) = 0.0; m(1, 2) = 0.0; m(2, 2) = 1.0; m(3, 2) = 0.0;
    m(0, 3) = 0.0; m(1, 3) = 0.0; m(2, 3) = 0.0; m(3, 3) = 1.0;

    M = m;
  }


/*
Initialize a scaling transformation
**/
  MeshTransformation(double s1, double s2, double s3)
  {
      MatrixXd m(4, 4);
      m(0, 0) = s1; m(1, 0) = 0.0; m(2, 0) = 0.0; m(3, 0) = 0.0;
      m(0, 1) = 0.0; m(1, 1) = s2; m(2, 1) = 0.0; m(3, 1) = 0.0;
      m(0, 2) = 0.0; m(1, 2) = 0.0; m(2, 2) = s3; m(3, 2) = 0.0;
      m(0, 3) = 0.0; m(1, 3) = 0.0; m(2, 3) = 0.0; m(3, 3) = 1.0;

      M = m;
  }

/*
Initialize a rotation transformation around a given axis (X, Y or Z) <br><br>

 @param  direction  a value 0, 1 or 2 indicating the direction (X, Y or Z respectively)
**/
  MeshTransformation(double theta, int direction)
  {
      // Create float based on boolean to test the value of direction
      float test_0 = float(direction == 0);
      float test_1 = float(direction == 1);
      float test_2 = float(direction == 2);

      MatrixXd m(4, 4);
      m(0, 0) = test_0 + (test_1 + test_2) * cos(theta); m(1, 0) = test_2 * (-sin(theta)); m(2, 0) = test_1 * (-sin(theta)); m(3, 0) = 0.0;
      m(0, 1) = test_2 * sin(theta); m(1, 1) = test_1 + (test_0 + test_2) * cos(theta) ; m(2, 1) = test_0 * (-sin(theta)); m(3, 1) = 0.0;
      m(0, 2) = test_1 * sin(theta); m(1, 2) = test_0 * sin(theta); m(2, 2) =test_2 + (test_0 + test_1) * cos(theta); m(3, 2) = 0.0;
      m(0, 3) = 0.0; m(1, 3) = 0.0; m(2, 3) = 0.0; m(3, 3) = 1.0;

      M = m;
  }

/*
Initialize a translation
**/
  MeshTransformation(RowVector3d t)
  {
      MatrixXd m(4, 4);
      m(0, 0) = 1.0; m(1, 0) = 0.0; m(2, 0) = 0.0; m(3, 0) = t[0];
      m(0, 1) = 0.0; m(1, 1) = 1.0; m(2, 1) = 0.0; m(3, 1) = t[1];
      m(0, 2) = 0.0; m(1, 2) = 0.0; m(2, 2) = 1.0; m(3, 2) = t[2];
      m(0, 3) = 0.0; m(1, 3) = 0.0; m(2, 3) = 0.0; m(3, 3) = 1.0;

      M = m;
  }

/*
Matrix accessor

@return  the matrix transformation
**/
  MatrixXd get_matrix() {
    return M;
  }

/*
Initialize a transformation given an input matrix 'm'
**/
  void set_matrix(MatrixXd m)
  {
    M = m;
  }

/*
Apply the transformation to all vertices stored in a matrix 'V' <br>

@param V  vector storing the input points
**/

  void transform(MatrixXd &V, RowVector3d &bar) {
      // Transform each points (one point=one row of V1) with the transformation T
      for(int i=0; i<V.rows(); i++){
          V.row(i) =  transform(bar - V.row(i)) + bar;
      }
  }

    /**
     * Apply the transformation to a 3d (row) vector 'v' <br>
   *
   * Remark: use homogeneous coordinates
   *
   * @return  the vector after transformation
     */
    RowVector3d transform(RowVector3d v) {

        std::cout << "Before the extra-dimension adding" << v << std::endl;
        // Modify the coordinates by adding an extra-dimension (x, y, z) ---> (x, y, z, 1)
        RowVector4d extent(v[0], v[1], v[2], 1.0);
        std::cout << "After the extra-dimension adding" << extent << std::endl;

        std::cout << "Transformation matrice T: "<< get_matrix() << std::endl;
        extent = extent * get_matrix();
        std::cout << "After transformation"<< extent << std::endl;

        // Return to homogeneous coordinates (x, y, z, w) --> (x/w, y/w, z/w)
        v = RowVector3d(extent[0]/extent[3], extent[1]/extent[3], extent[2]/extent[3]) ;
        std::cout << "Final coordinates:"<< v << std::endl;
        return v;
    }

    /**
     * Compose the current transformation with a transfomation 't': return a new transformation
     */
    MeshTransformation compose(MeshTransformation t) {
        MatrixXd m = t.get_matrix() * get_matrix();
        t.set_matrix(m);
        return t;
    }

    /**
     * Print the matrix transformation
     */
  friend std::ostream& operator<<(std::ostream &os, MeshTransformation& t) {
    return os << "matrix:\n" << t.get_matrix() << std::endl;
  }

private:
  MatrixXd M; // a 4x4 matrix representing a linear transformation
};
