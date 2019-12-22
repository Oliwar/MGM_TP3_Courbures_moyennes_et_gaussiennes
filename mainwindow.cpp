#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <cmath>

/* **** début de la partie à compléter **** */
/**
 * Calcul l'aire du triangle faceID
 * @brief MainWindow::faceArea
 * @param _mesh
 * @param faceID
 * @return l'aire du triangle faceID
 */
float MainWindow::faceArea(MyMesh* _mesh, int faceID)
{
    FaceHandle fh = _mesh->face_handle(faceID);

    HalfedgeHandle heh = _mesh->halfedge_handle(fh);

    VertexHandle vh_A = _mesh->to_vertex_handle(heh);

    heh = _mesh->next_halfedge_handle(heh);
    VertexHandle vh_B = _mesh->to_vertex_handle(heh);

    heh = _mesh->next_halfedge_handle(heh);
    VertexHandle vh_C = _mesh->to_vertex_handle(heh);

    return ( (_mesh->point(vh_B) - _mesh->point(vh_A)) % (_mesh->point(vh_C) - _mesh->point(vh_A)) ).norm() / 2;
}

/**
 * Calcul l'angle signé entre les normales des faces voisines à l'arête vertID1 - vertID0
 * @brief MainWindow::angleFF
 * @param _mesh
 * @param faceID0
 * @param faceID1
 * @param vertID0
 * @param vertID1
 * @return l'angle entre deux faces
 */
float MainWindow::angleFF(MyMesh* _mesh, int faceID0,  int faceID1, int vertID0, int vertID1)
{
    FaceHandle fh0 = _mesh->face_handle(faceID0);
    FaceHandle fh1 = _mesh->face_handle(faceID1);

    VertexHandle vh0 = _mesh->vertex_handle(vertID0);
    VertexHandle vh1 = _mesh->vertex_handle(vertID1);

    if(((_mesh->normal(fh0) % _mesh->normal(fh1)) | (_mesh->point(vh1) - _mesh->point(vh0))) < 0){
        return -std::acos(_mesh->normal(fh0) | _mesh->normal(fh1));
    }

    return std::acos(_mesh->normal(fh0) | _mesh->normal(fh1));
}

/**
 * Calcul l'angle absolu entre les deux arêtes de la face faceID depuis le sommet vertexID
 * @brief MainWindow::angleEE
 * @param _mesh
 * @param vertexID
 * @param faceID
 * @return l'angle entre deux arêtes
 */
float MainWindow::angleEE(MyMesh* _mesh, int vertexID,  int faceID)
{
    VertexHandle vh = _mesh->vertex_handle(vertexID);
    FaceHandle fh = _mesh->face_handle(faceID);

    HalfedgeHandle heh = _mesh->halfedge_handle(fh);

    VertexHandle vhA;
    VertexHandle vhB;
    VertexHandle vhC;

    for(int i = 0; i < 3; ++i){
        vhA = _mesh->to_vertex_handle(heh);
        if(vhA == vh){
            vhB = _mesh->from_vertex_handle(heh);
            heh = _mesh->next_halfedge_handle(heh);
            vhC = _mesh->to_vertex_handle(heh);
            break;
        }
        heh = _mesh->next_halfedge_handle(heh);
    }

    return std::acos(
                (_mesh->point(vhB) - _mesh->point(vhA)).normalize()
                |
                (_mesh->point(vhC) - _mesh->point(vhA)).normalize()
                );
}

/**
 * Calcul l'aire barycentrique autour du sommet vertexID
 * @brief MainWindow::areaBary
 * @param _mesh
 * @param vertexID
 * @return l'aire barycentrique du sommet vertexID
 */
float MainWindow::areaBary(MyMesh* _mesh, int vertexID){
    VertexHandle vh = _mesh->vertex_handle(vertexID);

    float aireBary = 0;
    FaceHandle fh;
    for (MyMesh::VertexFaceIter vf_it=_mesh->vf_iter(vh); vf_it.is_valid(); ++vf_it)
    {
        fh = *vf_it;
        aireBary += faceArea(_mesh, fh.idx());
    }

    return aireBary /= 3.0f;
}

/**
 * Calcul la somme des angles FF autour d'un sommet
 * @brief MainWindow::angleH_Curv
 * @param _mesh
 * @param vertexID
 * @return
 */
float MainWindow::angleH_Curv(MyMesh* _mesh, int vertexID){
    VertexHandle vh = _mesh->vertex_handle(vertexID);
    float angle = 0;
    HalfedgeHandle heh;
    FaceHandle fh0;
    FaceHandle fh1;
    VertexHandle vh0;
    VertexHandle vh1;
    for (MyMesh::VertexOHalfedgeIter voh_it=_mesh->voh_iter(vh); voh_it.is_valid(); ++voh_it)
    {
        heh = *voh_it;

        fh0 = _mesh->face_handle(heh);
        if(!fh0.is_valid()) continue;

        fh1 = _mesh->opposite_face_handle(heh);
        if(!fh1.is_valid()) continue;

        vh0 = _mesh->from_vertex_handle(heh);
        vh1 = _mesh->to_vertex_handle(heh);

        angle += angleFF(_mesh, fh0.idx(), fh1.idx(), vh0.idx(), vh1.idx()) * ( _mesh->point(vh1) - _mesh->point(vh0) ).norm();
    }
    return angle;
}

/**
 * Calcul la somme des angles EE autour d'un sommet
 * @brief MainWindow::angleK_Curv
 * @param _mesh
 * @param vertexID
 * @return
 */
float MainWindow::angleK_Curv(MyMesh* _mesh, int vertexID){
    VertexHandle vh = _mesh->vertex_handle(vertexID);
    float angle = 0;
    FaceHandle fh;
    for (MyMesh::VertexFaceIter vf_it=_mesh->vf_iter(vh); vf_it.is_valid(); ++vf_it)
    {
        fh = *vf_it;
        angle += angleEE(_mesh, vh.idx(), fh.idx());
    }
    return angle;
}

/**
 * Calcul la courbure moyenne
 * @brief MainWindow::H_Curv
 * @param _mesh
 */
void MainWindow::H_Curv(MyMesh* _mesh)
{
    float aireBary;
    float angle;
    VertexHandle vh_cur;
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert!=_mesh->vertices_end(); ++curVert)
    {
        vh_cur = *curVert;

        aireBary = areaBary(_mesh, vh_cur.idx());

        angle = angleH_Curv(_mesh, vh_cur.idx());

        _mesh->data(vh_cur).value = (1 / (4*aireBary) ) * angle;
    }
}

/**
 * Calcul la courbure Gaussienne
 * @brief MainWindow::K_Curv
 * @param _mesh
 */
void MainWindow::K_Curv(MyMesh* _mesh)
{
    float aireBary;
    float angle;
    VertexHandle vh_cur;
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert!=_mesh->vertices_end(); ++curVert)
    {
        vh_cur = *curVert;

        aireBary = areaBary(_mesh, vh_cur.idx());

        angle = angleK_Curv(_mesh, vh_cur.idx());

        _mesh->data(vh_cur).value = (1/aireBary)*(2*M_PI - angle);
    }
}
/* **** fin de la partie à compléter **** */



/* **** début de la partie boutons et IHM **** */
void MainWindow::on_pushButton_H_clicked()
{
    H_Curv(&mesh);
    displayMesh(&mesh, true); // true permet de passer en mode "carte de temperatures", avec une gestion automatique de la couleur (voir exemple)
}

void MainWindow::on_pushButton_K_clicked()
{
    K_Curv(&mesh);
    displayMesh(&mesh, true); // true permet de passer en mode "carte de temperatures", avec une gestion automatique de la couleur (voir exemple)
}

/*
    Cette fonction est à utiliser UNIQUEMENT avec le fichier testAngleArea.obj
    Elle est appelée par le bouton "Test angles/aires"

    Elle permet de vérifier les fonctions faceArea, angleFF et angleEE.
    Elle doit afficher :

    Aire de la face 0 : 2
    Aire de la face 1 : 2
    Angle entre les faces 0 et 1 : 1.5708
    Angle entre les faces 1 et 0 : -1.5708
    Angle au sommet 1 sur la face 0 : 0.785398
*/

void MainWindow::on_pushButton_angleArea_clicked()
{
    qDebug() << "Aire de la face 0 :" << faceArea(&mesh, 0);
    qDebug() << "Aire de la face 1 :" << faceArea(&mesh, 1);

    qDebug() << "Angle entre les faces 0 et 1 :" << angleFF(&mesh, 0, 1, 1, 2);
    qDebug() << "Angle entre les faces 1 et 0 :" << angleFF(&mesh, 1, 0, 1, 2);

    qDebug() << "Angle au sommet 1 sur la face 0 :" << angleEE(&mesh, 1, 0);
    qDebug() << "Angle au sommet 3 sur la face 1 :" << angleEE(&mesh, 3, 1);
}

void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}
/* **** fin de la partie boutons et IHM **** */

/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
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
void MainWindow::displayMesh(MyMesh* _mesh, bool isTemperatureMap, float mapRange)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(isTemperatureMap)
    {
        QVector<float> values;

        if(mapRange == -1)
        {
            for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                values.append(fabs(_mesh->data(*curVert).value));
            qSort(values);
            mapRange = values.at(values.size()*0.8);
            qDebug() << "mapRange" << mapRange;
        }

        float range = mapRange;
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
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
