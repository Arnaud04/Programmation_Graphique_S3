#include "mainwindow.h"
#include "ui_mainwindow.h"


/* **** début de la partie boutons et IHM **** */

void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}

/* **** fin de la partie boutons et IHM **** */


/* **** fonctions supplémentaires **** */

//return total mesh surface
float MainWindow::compute_area(MyMesh *_mesh)
{
    float area = 0.0;
    for(unsigned int i=0;i<_mesh->n_faces();i++)
    {
        area += compute_face_area(_mesh , i);
    }
    return area;
}

//return The smallest face area
float MainWindow::get_min_area(MyMesh *_mesh)
{
    float min = compute_face_area(_mesh,0);
    for(unsigned int i=0;i<_mesh->n_faces();i++)
    {
        if(min > compute_face_area(_mesh , i))
            min = compute_face_area(_mesh , i);

    }
    return min;
}

//return the biggest face area
float MainWindow::get_max_area(MyMesh *_mesh)
{
    float max = compute_face_area(_mesh,0);
    for(unsigned int i=0;i>_mesh->n_faces();i++)
    {
        if(max > compute_face_area(_mesh , i))
            max = compute_face_area(_mesh , i);

    }
    return max;
}

/*
 * retourne l'air d'une face triangulaire passé en paramètre
 * Ses paramètres sont :
 *      - le maillage chargé
 *      - La face dont on veut l'air
*/
float MainWindow::compute_face_area(MyMesh *_mesh, int n_face)
{
    QVector<float> vectors;
    float area = 0.0;

    FaceHandle fh = _mesh->face_handle(n_face);

    for(MyMesh::FaceVertexIter curVertex = _mesh->fv_iter(fh); curVertex.is_valid(); curVertex ++)
    {
        VertexHandle vh = *curVertex;
        vectors.append(_mesh->point(vh)[0]);
        vectors.append(_mesh->point(vh)[1]);
        vectors.append(_mesh->point(vh)[2]);
    }
    //on a nos 3 sommets de chaque face : vertex a = vector[0][1][2] vertex b = [3][4][5] vertex c = [6][7][8]

    QVector<float> vectAB;
    QVector<float> vectAC;
    QVector<float> produitVect;

    //On calcul les vecteur AB et AC
    for(int i=0; i<3; i++)
    {
        vectAB.append(vectors[i+3]-vectors[i]);
        vectAC.append(vectors[i+6]-vectors[i]);
    }
    // AB : (x1,y1,z1)    AC : (x2,y2,z2)
    produitVect.append(vectAB[1]*vectAC[2]-vectAB[2]*vectAC[1]); //produitVectoriel : vx
    produitVect.append(vectAB[0]*vectAC[2]-vectAB[2]*vectAC[0]); //produitVectoriel: vy
    produitVect.append(vectAB[0]*vectAC[1]-vectAB[1]*vectAC[0]); //produitVectoriel: vz

    //calcul area = 1/2||produiVect||
    area = 0.5*sqrt(produitVect[0]*produitVect[0]+produitVect[1]*produitVect[1]+produitVect[2]*produitVect[2]);

    return area;
}

std::vector<VertexHandle> getNormalFace (MyMesh* _mesh,VertexHandle v1, VertexHandle v2)
{
         std::vector<MyMesh::Point> normal(2);

         //On fait le produit vectoriel de nos 2 vetceurs
         float a,b,c;
         a = (_mesh->point(v1)[1] * _mesh->point(v2)[2])-(_mesh->point(v1)[2]*_mesh->point(v2)[1]);
         b = (_mesh->point(v1)[0] * _mesh->point(v2)[2])-(_mesh->point(v1)[2]*_mesh->point(v2)[0]);
         c = (_mesh->point(v1)[0] * _mesh->point(v2)[1])-(_mesh->point(v1)[1]*_mesh->point(v2)[0]);
        /*
        //===========================================

        VertexHandle vh = _mesh->vertex_handle(vertexID);
        FaceHandle fh = _mesh->face_handle(faceID);
        QVector<VertexHandle> listePoints;
        QVector<VertexHandle> listePointsOnFace;
        QVector<float> vectors;
        //On identifie les point voisins de vertexID qui appartiennent à faceID

        //tout les points voisins de vertexID
        for (MyMesh::VertexVertexIter curVertex = _mesh->vv_iter(vh); curVertex.is_valid(); curVertex ++)
        {
            VertexHandle v = *curVertex;
            listePoints.append(v);
        }
        //parmis ces points ceux qui appartiennent à la faceID

        for(MyMesh::FaceVertexIter curVertex = _mesh->fv_iter(fh); curVertex.is_valid(); curVertex ++)
        {
            VertexHandle v = *curVertex;
            if(listePoints.contains(v) && v.idx() != vertexID)
            {
                listePointsOnFace.append(v);
            }
        }

        //On créer des vecteurs a partir des point obtenu
        for(int i=0; i<listePointsOnFace.size();i++)
        {
            vectors.append((_mesh->point(listePointsOnFace[i])[0])-(_mesh->point(vh)[0]));
            vectors.append((_mesh->point(listePointsOnFace[i])[1])-(_mesh->point(vh)[1]));
            vectors.append((_mesh->point(listePointsOnFace[i])[2])-(_mesh->point(vh)[2]));

        }

        //on normalise les vecteurs obtenu
        QVector<float> norme;

        float tmp = vectors[0]*vectors[0]+vectors[1]*vectors[1]+vectors[2]*vectors[2];
        norme.append(sqrt(tmp)); //sqrt(x^2+y^2+z^2) : C'est la norme du premier vecteur

        tmp = vectors[3]*vectors[3]+vectors[4]*vectors[4]+vectors[5]*vectors[5];
        norme.append(sqrt(tmp)); //C'est la norme du second vecteurs

        //normalisation des vecteur :  chacun des vecteur / par la norme
        for(int i=0; i<listePointsOnFace.size(); i++) //2 vecteur : 2 iterations : 6points
        {
            vectors[i*3] = vectors[i*3]/norme[i];
            vectors[i*3+1] = vectors[i*3+1]/norme[i];
            vectors[i*3+2] = vectors[i*3+2]/norme[i];
        }

        float a = vectors[0]*vectors[0+3];
        float b = vectors[1]*vectors[1+3];
        float c = vectors[2]*vectors[2+3];

        float prodScal = a+b+c;
        abs(prodScal); //garantie que le produit scalaire soit positif

        float angle = acos(prodScal)/**180/M_PI*/;
       // return angle; //en radians

    //===========================================

    */
    return normal;

}


/*
 *  permet d'initialiser les couleurs et les épaisseurs des élements du maillage
*/
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
    MyMesh::ConstFaceVertexIter fvIt;
    int i = 0;
    for (; fIt!=fEnd; ++fIt)
    {
        fvIt = _mesh->cfv_iter(*fIt);
        triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
        triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
        triIndiceArray[i] = i;

        i++; ++fvIt;
        triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
        triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
        triIndiceArray[i] = i;

        i++; ++fvIt;
        triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
        triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
        triIndiceArray[i] = i;

        i++;
    }

    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }

    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    modevoisinage = false;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

/*
 * Affiche les coordonnées du barycentre dans la console
 * - peut aussi affiché les coordonnées de chacun des points
 *   du mesh pour s'en convaincre
*/
void MainWindow::on_pushButton_barycentre_clicked()
{

    //On parcourt toutes les faces du mesh

    std::vector<double> coordonnees;

    float x = 0;
    float y = 0;
    float z = 0;

    for(MyMesh::VertexIter vit = mesh.vertices_begin(); vit != mesh.vertices_end(); vit++)
    {

        VertexHandle vh = *vit;
        x = x + mesh.point(vh)[0];
        y = y + mesh.point(vh)[1];
        z = z + mesh.point(vh)[2];

        //visualisation des coordonnées de chacun des points du mesh
        /*
        qDebug()<< "============";
        qDebug() << "x " <<mesh.point(vh)[0] << " ";
        qDebug() << "y " <<mesh.point(vh)[1] << " ";
        qDebug() << "z " <<mesh.point(vh)[2] << " ";*/
    }

            float barycenter_x=0.0; float barycenter_y=0.0;float barycenter_z=0.0;
    if((x && y && z) != 0)
    {
        barycenter_x = x/(mesh.n_faces());
        barycenter_y = y/(mesh.n_faces());
        barycenter_z = z/(mesh.n_faces());
        qDebug()<<"coordonnées barycentriques : x= "<< barycenter_x << "y= " << barycenter_y << "z= " << barycenter_z;
    }

}

/*
 * THis function return the Barycentrique point of a face
 * It is used to compute normal face
*/
MyMesh::Point getBarycenterFromFace(VertexHandle vh ,FaceHandle fh, MyMesh* _mesh)
{
    //amelioration la rendre générique pour nimporte quelle face
    //pour ca remplacer 3 par noombre de sommet sur FaceCourante
    std::vector<double> coordonnees;

    float x = 0;
    float y = 0;
    float z = 0;

    x = x + _mesh->point(vh)[0];
    y = y + _mesh->point(vh)[1];
    z = z + _mesh->point(vh)[2];

    float barycenter_x=0.0; float barycenter_y=0.0;float barycenter_z=0.0;

    barycenter_x = x/3.0;
    barycenter_y = y/3.0;
    barycenter_z = z/3.0;

    MyMesh::Point barycentriquVertex(barycenter_x, barycenter_y, barycenter_z);
    return barycentriquVertex;

}

void MainWindow::on_getInformation_clicked()
{
     int nbVert = mesh.n_vertices();
     int nbFaces = mesh.n_faces();

     ui->verticesNumber->display(nbVert);
     ui->facesNumber->display(nbFaces);

     qDebug()<< "===vertex list ===";
     for(MyMesh::VertexIter vit = mesh.vertices_begin(); vit != mesh.vertices_end(); vit++)
     {
         VertexHandle vh = *vit;
         qDebug() << "id : " <<vh.idx();
     }

     qDebug()<< "===Face list ===";
     for(MyMesh::FaceIter f = mesh.faces_begin(); f != mesh.faces_end(); f++)
     {
         FaceHandle fh = *f;
         qDebug() << "id : " <<fh.idx();
     }


}

/*
 * Fonction qui calcul et affiche les valeur limite du maillage chargé
 */
void MainWindow::on_boundingBox_clicked()
{
    std::vector<double> minCoordonnees(3,DBL_MAX); //(x_min, y_min, z_min)
    std::vector<double> maxCoordonnees(3,DBL_MIN); //(x_max, y_max, z_max)

    for(MyMesh::VertexIter vit = mesh.vertices_begin(); vit != mesh.vertices_end(); vit++)
    {
        VertexHandle vh = *vit;
        if (minCoordonnees[0] > mesh.point(vh)[0])
            minCoordonnees[0] = mesh.point(vh)[0];

        if (minCoordonnees[1] > mesh.point(vh)[1])
            minCoordonnees[1] = mesh.point(vh)[1];

        if (minCoordonnees[2] > mesh.point(vh)[2])
            minCoordonnees[2] = mesh.point(vh)[2];

        if (maxCoordonnees[0] < mesh.point(vh)[0])
            maxCoordonnees[0] = mesh.point(vh)[0];

        if (maxCoordonnees[1] < mesh.point(vh)[1])
            maxCoordonnees[1] = mesh.point(vh)[1];

        if (maxCoordonnees[2] < mesh.point(vh)[2])
            maxCoordonnees[2] = mesh.point(vh)[2];

    }
    //On a nos coordonnées max et min respectivement pour x,y et z
    /*qDebug()<< maxCoordonnees[0] << maxCoordonnees[1] << maxCoordonnees[2];
    qDebug()<< minCoordonnees[0] << minCoordonnees[1] << minCoordonnees[2];*/

    /*Ainsi la boite englobante à 2^3 points et chaque coordonnées est le fruit d'une composition
     * des coordonnées précedement calculées.
    */


    qDebug()<<" Voici les différents points de la boites";
    qDebug()<<"p1" << minCoordonnees[0] <<","<<minCoordonnees[1] <<","<<minCoordonnees[2];
    qDebug()<<"p2" << maxCoordonnees[0] <<","<<minCoordonnees[1] <<","<<minCoordonnees[2];
    qDebug()<<"p3" << minCoordonnees[0] <<","<<maxCoordonnees[1] <<","<<minCoordonnees[2];
    qDebug()<<"p4" << maxCoordonnees[0] <<","<<maxCoordonnees[1] <<","<<minCoordonnees[2];
    qDebug()<<"p5" << minCoordonnees[0] <<","<<minCoordonnees[1] <<","<<maxCoordonnees[2];
    qDebug()<<"p6" << maxCoordonnees[0] <<","<<minCoordonnees[1] <<","<<maxCoordonnees[2];
    qDebug()<<"p7" << minCoordonnees[0] <<","<<maxCoordonnees[1] <<","<<maxCoordonnees[2];
    qDebug()<<"p8" << maxCoordonnees[0] <<","<<maxCoordonnees[1] <<","<<maxCoordonnees[2];

    //======k===

        MyMesh mesh;

        // on construit une liste de sommets
        MyMesh::VertexHandle sommets[8];
        sommets[0] = mesh.add_vertex(MyMesh::Point(minCoordonnees[0], minCoordonnees[1], minCoordonnees[2]));
        sommets[1] = mesh.add_vertex(MyMesh::Point(minCoordonnees[0], minCoordonnees[1], maxCoordonnees[2]));
        sommets[2] = mesh.add_vertex(MyMesh::Point(minCoordonnees[0], maxCoordonnees[1], minCoordonnees[2]));
        sommets[3] = mesh.add_vertex(MyMesh::Point(minCoordonnees[0], maxCoordonnees[1], maxCoordonnees[2]));
        sommets[4] = mesh.add_vertex(MyMesh::Point(maxCoordonnees[0], minCoordonnees[1], minCoordonnees[2]));
        sommets[5] = mesh.add_vertex(MyMesh::Point(maxCoordonnees[0], minCoordonnees[1], maxCoordonnees[2]));
        sommets[6] = mesh.add_vertex(MyMesh::Point(maxCoordonnees[0], maxCoordonnees[1], minCoordonnees[2]));
        sommets[7] = mesh.add_vertex(MyMesh::Point(maxCoordonnees[0], maxCoordonnees[1], maxCoordonnees[2]));


        // on construit des faces à partir des sommets

        std::vector<MyMesh::VertexHandle> uneNouvelleFace;


        uneNouvelleFace.clear();
        uneNouvelleFace.push_back(sommets[7]);
        uneNouvelleFace.push_back(sommets[6]);
        uneNouvelleFace.push_back(sommets[4]);
        uneNouvelleFace.push_back(sommets[5]);
        mesh.add_face(uneNouvelleFace);

        uneNouvelleFace.clear();
        uneNouvelleFace.push_back(sommets[5]);
        uneNouvelleFace.push_back(sommets[4]);
        uneNouvelleFace.push_back(sommets[0]);
        uneNouvelleFace.push_back(sommets[1]);
        mesh.add_face(uneNouvelleFace);

        uneNouvelleFace.clear();
        uneNouvelleFace.push_back(sommets[0]);
        uneNouvelleFace.push_back(sommets[2]);
        uneNouvelleFace.push_back(sommets[3]);
        uneNouvelleFace.push_back(sommets[1]);
        mesh.add_face(uneNouvelleFace);

        uneNouvelleFace.clear();
        uneNouvelleFace.push_back(sommets[2]);
        uneNouvelleFace.push_back(sommets[6]);
        uneNouvelleFace.push_back(sommets[7]);
        uneNouvelleFace.push_back(sommets[3]);
        mesh.add_face(uneNouvelleFace);

        uneNouvelleFace.clear();
        uneNouvelleFace.push_back(sommets[2]);
        uneNouvelleFace.push_back(sommets[0]);
        uneNouvelleFace.push_back(sommets[4]);
        uneNouvelleFace.push_back(sommets[6]);
        mesh.add_face(uneNouvelleFace);


        mesh.update_normals();
        this->mesh = mesh;
        // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
        resetAllColorsAndThickness(&mesh);

        // on affiche le maillage
        displayMesh(&mesh);


        /*   delete_vertices() delete_faces()  delete_edges()*/
        //delete_vertices();
        //delete_faces();


}

void MainWindow::on_pushButton_area_clicked()
{
    qDebug()<<"Surface total du mesh" <<compute_area(&mesh);
    qDebug()<<"Surface de la plus petite face du mesh" <<get_min_area(&mesh);
    qDebug()<<"Surface de la plus grande face du mesh" <<get_max_area(&mesh);
}



void MainWindow::on_triangleSurface_proportion_clicked()
{

    float fullMesh_area = compute_area(&mesh);
    std::vector<float> areaRepartition(10);

    for(MyMesh::FaceIter f = mesh.faces_begin(); f != mesh.faces_end(); f++)
    {
        FaceHandle fh = *f;

        float fh_area = compute_face_area(&mesh, fh.idx());
        float areaFrequency = fullMesh_area/10.0;


        for(int i=0;i<10;i++)
        {
            if((fh_area > areaFrequency*i)&&(fh_area <= areaFrequency*(i+1)))
                areaRepartition[i] = areaRepartition[i]+1;
        }


    }

    for(int i=0; i<10;i++)
    {
        qDebug()<< areaRepartition[i] << "face(s) on une taille compris entre"<< i*10 << "% et" <<(i+1)*10 <<"% de la surface total du mesh";
        qDebug()<< ((float)areaRepartition[i]*100.0)/(float)mesh.n_faces() <<"% de triangle \n";
    }

}

void MainWindow::test()
{
    qDebug()<<"TESTTTTTTTTTTTT";
}

bool MainWindow::containIsolated_points()
{
    bool hasIsolatePoint = false;

    //to be complated

    return hasIsolatePoint;
}




void MainWindow::on_meshIsValid_clicked()
{
    bool a,b;
    test();
    //a = containPoints(&mesh);
    //b = containTriangles(&mesh);
    //appel de au fonction précédent non fonctionnel

    // partie a mettre dans fonction containPoints()
    unsigned int numberNeighborsPoint;
    bool hasSinglePoints = false;

    for(MyMesh::VertexIter vit = mesh.vertices_begin(); vit != mesh.vertices_end(); vit ++)
    {
        VertexHandle vh = *vit;
        numberNeighborsPoint=0;
        for (MyMesh::VertexVertexIter curentVertex = mesh.vv_iter(vh); curentVertex.is_valid(); curentVertex ++)
        {
            numberNeighborsPoint ++;
        }
        if(numberNeighborsPoint==0)
        {
            hasSinglePoints = true;
        }

    }


    //========== partie a mettre dans fonction containTriangle()==============

    int counter = 0;
    bool hasTriangle = false;

    for (MyMesh::FaceIter fit = mesh.faces_begin(); fit != mesh.faces_end(); fit++)
    {
        counter = 0;
        for(MyMesh::VertexIter vit = mesh.vertices_begin(); vit != mesh.vertices_end(); vit ++)
        {
            counter ++;
        }
        if(counter == 3)
        {
            hasTriangle = true;

            break;
        }
    }

    //==========================

    if(hasTriangle && hasSinglePoints)
        qDebug()<<"Mesh ayant des points isolé et également des faces triangulaire";
    else
        qDebug()<<"Maillage conforme";
}
