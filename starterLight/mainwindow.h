#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    EdgeAttributes( OpenMesh::Attributes::Color );
    // vertex thickness
    VertexTraits{float thickness;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    float compute_area(MyMesh *_mesh);
    float get_min_area(MyMesh *_mesh);
    float get_max_area(MyMesh *_mesh);
    float compute_face_area(MyMesh *_mesh, int n_face);
    void displayMesh(MyMesh *_mesh);
    void resetAllColorsAndThickness(MyMesh *_mesh);
    bool containIsolated_points();
    MyMesh::Point getBarycenterFromFace(VertexHandle vh ,FaceHandle fh,MyMesh* _mesh);
    //MyMesh::Point getNormalFace (MyMesh* _mesh,VertexHandle vertexFromFace, float barycentre);
    std::vector<MyMesh::Point> getNormalFace (MyMesh* _mesh,VertexHandle vertexFromFace, float barycentre);
    void test();

    /**
     * @brief checkAllTriangularFace
     * @param _mesh
     * @return bool
     * @details vérifier que les fichiers contiennent des faces triangulaires.
     */
    bool checkAllTriangularFace(MyMesh* _mesh);

    /**
     * @brief checkOnlyPoint
     * @param _mesh
     * @return bool
     * @details vérifier que les fichiers contiennent seulement des points 3D.
     */
    bool checkOnlyPoint(MyMesh* _mesh);

    /**
     * @brief checkGlobalNeighbours
     * @param _mesh
     * @return bool
     * @details qu'il n'y a pas de faces sans voisines,
     * de point n'appartenant pas à une arête,
     * et qu'il n'y a pas d'arêtes n'appartenant pas à une face.
     */
    bool checkGlobalNeighbours(MyMesh* _mesh);

private slots:


    void on_pushButton_chargement_clicked();

    void on_pushButton_barycentre_clicked();

    void on_getInformation_clicked();

    void on_boundingBox_clicked();

    void on_pushButton_area_clicked();

    void on_triangleSurface_proportion_clicked();

    void on_meshIsValid_clicked();

private:

    bool modevoisinage;

    MyMesh mesh;

    int vertexSelection;
    int edgeSelection;
    int faceSelection;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
