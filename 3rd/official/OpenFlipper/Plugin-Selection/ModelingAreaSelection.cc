/*===========================================================================*\
 *                                                                           *
 *                              OpenFlipper                                  *
 *      Copyright (C) 2001-2010 by Computer Graphics Group, RWTH Aachen      *
 *                           www.openflipper.org                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is part of OpenFlipper.                                        *
 *                                                                           *
 *  OpenFlipper is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  OpenFlipper is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with OpenFlipper. If not,                                  *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/

/*===========================================================================*\
 *                                                                           *
 *   $Revision: 9595 $                                                         *
 *   $Author: moebius $                                                      *
 *   $Date: 2010-06-17 12:48:23 +0200 (Do, 17. Jun 2010) $                   *
 *                                                                           *
\*===========================================================================*/



#include "SelectionPlugin.hh"

#include <iostream>

#include <MeshTools/MeshSelectionT.hh>
#include <OpenFlipper/BasePlugin/PluginFunctions.hh>


//=========================================================
//==== Modeling Area selections
//=========================================================

void SelectionPlugin::selectModelingVertices( int objectId , IdList _vertexList ) {
  BaseObjectData* object;
  if ( ! PluginFunctions::getObject(objectId,object) ) {
    emit log(LOGERR,tr("selectModelingVertices : unable to get object") );
    return;
  }

  if ( _vertexList.size() == 0 )
    return;

  if ( object->dataType() == DATA_TRIANGLE_MESH ){
      MeshSelection::setArea(PluginFunctions::triMesh(object) , _vertexList , AREA, true);
      update_regions( PluginFunctions::triMesh(object) );
  } else if ( object->dataType() == DATA_POLY_MESH ){
      MeshSelection::setArea(PluginFunctions::polyMesh(object) , _vertexList , AREA, true);
      update_regions( PluginFunctions::polyMesh(object) );
#ifdef ENABLE_TSPLINEMESH_SUPPORT
  } else if ( object->dataType() == DATA_TSPLINE_MESH ){
      MeshSelection::setArea(PluginFunctions::tsplineMesh(object) , _vertexList , AREA, true);
      update_regions( PluginFunctions::tsplineMesh(object) );
#endif
  } else {
      emit log(LOGERR,tr("selectModelingVertices : Unsupported object Type") );
      return;
  }

  QString selection = "selectModelingVertices( ObjectId , [ " + QString::number(_vertexList[0]);

  for ( uint i = 1 ; i < _vertexList.size(); ++i) {
    selection +=  " , " + QString::number(_vertexList[i]);
  }

  selection += " ] )";

  emit updatedObject(object->id(), UPDATE_ALL);
  emit scriptInfo( selection );
}

//=========================================================

void SelectionPlugin::unselectModelingVertices( int objectId , IdList _vertexList ) {
  BaseObjectData* object;
  if ( ! PluginFunctions::getObject(objectId,object) ) {
    emit log(LOGERR,tr("unselectModelingVertices : unable to get object") );
    return;
  }

  if ( _vertexList.size() == 0 )
    return;

  if ( object->dataType() == DATA_TRIANGLE_MESH )
      MeshSelection::setArea(PluginFunctions::triMesh(object) , _vertexList , AREA, false);
  else if ( object->dataType() == DATA_POLY_MESH )
      MeshSelection::setArea(PluginFunctions::polyMesh(object) , _vertexList , AREA, false);
#ifdef ENABLE_TSPLINEMESH_SUPPORT
  else if ( object->dataType() == DATA_TSPLINE_MESH )
      MeshSelection::setArea(PluginFunctions::tsplineMesh(object) , _vertexList , AREA, false);
#endif
  else{
      emit log(LOGERR,tr("unselectModelingVertices : Unsupported object Type") );
      return;
  }

  QString selection = "unselectModelingVertices( ObjectId , [ " + QString::number(_vertexList[0]);

  for ( uint i = 1 ; i < _vertexList.size(); ++i) {
    selection +=  " , " + QString::number(_vertexList[i]);
  }

  selection += " ] )";

  emit updatedObject(object->id(), UPDATE_ALL);
  emit scriptInfo( selection );
}

//=========================================================

void SelectionPlugin::clearModelingVertices( int objectId ) {
  BaseObjectData* object;
  if ( ! PluginFunctions::getObject(objectId,object) ) {
    emit log(LOGERR,tr("clearModelingVertices : unable to get object") );
    return;
  }

  if ( object->dataType() == DATA_TRIANGLE_MESH )
      MeshSelection::setArea(PluginFunctions::triMesh(object) , AREA, false);
  else if ( object->dataType() == DATA_POLY_MESH )
      MeshSelection::setArea(PluginFunctions::polyMesh(object) , AREA, false);
#ifdef ENABLE_TSPLINEMESH_SUPPORT
  else if ( object->dataType() == DATA_TSPLINE_MESH )
      MeshSelection::setArea(PluginFunctions::tsplineMesh(object) , AREA, false);
#endif
  else{
      emit log(LOGERR,tr("clearModelingVertices : Unsupported object Type") );
      return;
  }

  emit updatedObject(object->id(), UPDATE_ALL);
  emit scriptInfo( "clearModelingVertices( ObjectId )" );
}

//=========================================================

void SelectionPlugin::setAllModelingVertices( int objectId  ) {
  BaseObjectData* object;
  if ( ! PluginFunctions::getObject(objectId,object) ) {
    emit log(LOGERR,tr("setAllModelingVertices : unable to get object") );
    return;
  }

  if ( object->dataType() == DATA_TRIANGLE_MESH )
      MeshSelection::setArea(PluginFunctions::triMesh(object) , AREA, true);
  else if ( object->dataType() == DATA_POLY_MESH )
      MeshSelection::setArea(PluginFunctions::polyMesh(object) , AREA, true);
#ifdef ENABLE_TSPLINEMESH_SUPPORT
  else if ( object->dataType() == DATA_TSPLINE_MESH )
      MeshSelection::setArea(PluginFunctions::tsplineMesh(object) , AREA, true);
#endif
  else{
      emit log(LOGERR,tr("setAllModelingVertices : Unsupported object Type") );
      return;
  }

  emit updatedObject(object->id(), UPDATE_ALL);
  emit scriptInfo( "setAllModelingVertices( ObjectId )" );
}

//=========================================================
IdList SelectionPlugin::getModelingVertices( int objectId  ) {
  BaseObjectData* object;
  if ( ! PluginFunctions::getObject(objectId,object) ) {
    emit log(LOGERR,tr("getModelingVertices : unable to get object") );
    return IdList(0);
  }

  emit scriptInfo( "getModelingVertices( ObjectId )" );

  if ( object->dataType() == DATA_TRIANGLE_MESH )
      return MeshSelection::getArea(PluginFunctions::triMesh(object) , AREA);
  else if ( object->dataType() == DATA_POLY_MESH )
      return MeshSelection::getArea(PluginFunctions::polyMesh(object) , AREA);
#ifdef ENABLE_TSPLINEMESH_SUPPORT
  else if ( object->dataType() == DATA_TSPLINE_MESH )
      return MeshSelection::getArea(PluginFunctions::tsplineMesh(object) , AREA);
#endif
  else{
      emit log(LOGERR,tr("getModelingVertices : Unsupported object Type") );
      return IdList(0);
  }

  return IdList(0);
}

//=========================================================

void SelectionPlugin::loadFlipperModelingSelection( int _objectId , QString _filename ) {
  QFile file(_filename);

  if ( ! file.exists() ) {
    emit log(LOGERR,tr("Unable to find file : ") + _filename );
    return;
  }

   if (file.open(QFile::ReadOnly)) {
    QTextStream input(&file);

    QString header = input.readLine();

    if ( !header.contains("Selection") ) {
       emit log(LOGERR,tr("Wrong file header! should be Selection but is ") + header );
       return;
    }

    header = input.readLine();

    bool ok = false;

    uint vertexCount = header.toUInt(&ok);

    if ( !ok ) {
       emit log(LOGERR,tr("Unable to parse header. Cant get vertex count from string : ") + header );
       return;
    }

    //compare VertexCount
    BaseObjectData* object;
    if ( ! PluginFunctions::getObject(_objectId, object) ) {
      emit log(LOGERR,tr("loadSelection : unable to get object") );
      return;
    }

    if ( object->dataType() == DATA_TRIANGLE_MESH ){
        if ( PluginFunctions::triMesh(object)->n_vertices() != vertexCount )
          return;
    } else if ( object->dataType() == DATA_POLY_MESH ){
        if ( PluginFunctions::polyMesh(object)->n_vertices() != vertexCount )
          return;
#ifdef ENABLE_TSPLINEMESH_SUPPORT
    } else if ( object->dataType() == DATA_TSPLINE_MESH ){
        if ( PluginFunctions::tsplineMesh(object)->n_vertices() != vertexCount )
          return;
#endif
    } else {
      return;
    }

    //get selected handles
    IdList handleVertices;
    IdList modelingVertices;

    uint vertexId = 0;

    do {
      // Split into two substrings
      QStringList inputList = input.readLine().split(" ");

      if ( inputList.size() != 2 ) {
        emit log(LOGERR,tr("Unable to parse entry at vertex index ") + QString::number( vertexId ) );
        return;
      }

      if ( inputList[0] == "1" )
        modelingVertices.push_back(vertexId);

      if ( inputList[1] == "1" )
        handleVertices.push_back(vertexId);

      ++vertexId;

    } while (!input.atEnd());

    clearModelingVertices(_objectId);
    selectModelingVertices(_objectId,modelingVertices);

    clearHandleVertices(_objectId);
    selectHandleVertices(_objectId,handleVertices);

  } else
    emit log(LOGERR,tr("Unable to open selection file!"));
}

//=========================================================

void SelectionPlugin::saveFlipperModelingSelection( int _objectId , QString _filename ) {
  QFile file(_filename);

  if (file.open(QFile::WriteOnly)) {
    QTextStream input(&file);

    //get the object
    BaseObjectData* object;
    if ( ! PluginFunctions::getObject(_objectId,object) ) {
      emit log(LOGERR,tr("saveFlipperModelingSelection : unable to get object") );
      return;
    }

    if ( object->dataType() == DATA_TRIANGLE_MESH ) {
      TriMesh* mesh = PluginFunctions::triMesh(object);

      //header
      input << "Selection" << endl;
      input << mesh->n_vertices() << endl;

      std::vector< int > modelingVertices = MeshSelection::getArea(mesh, AREA);
      std::vector< int > handleVertices   = MeshSelection::getArea(mesh, HANDLEAREA);

      std::vector< bool > modelingAll(mesh->n_vertices(), false);
      std::vector< bool > handleAll(mesh->n_vertices(), false);

      for (uint i=0; i < modelingVertices.size(); i++)
        modelingAll[ modelingVertices[i] ] = true;

      for (uint i=0; i < handleVertices.size(); i++)
        handleAll[ handleVertices[i] ] = true;

      for (uint i=0; i < mesh->n_vertices(); i++)
        input << (int) modelingAll[i] << " " << (int) handleAll[i] << endl;

    } else if ( object->dataType() == DATA_POLY_MESH){

      PolyMesh* mesh = PluginFunctions::polyMesh(object);

      //header
      input << "Selection" << endl;
      input << mesh->n_vertices() << endl;

      std::vector< int > modelingVertices = MeshSelection::getArea(mesh, AREA);
      std::vector< int > handleVertices   = MeshSelection::getArea(mesh, HANDLEAREA);

      std::vector< bool > modelingAll(mesh->n_vertices(), false);
      std::vector< bool > handleAll(mesh->n_vertices(), false);

      for (uint i=0; i < modelingVertices.size(); i++)
        modelingAll[ modelingVertices[i] ] = true;

      for (uint i=0; i < handleVertices.size(); i++)
        handleAll[ handleVertices[i] ] = true;

      for (uint i=0; i < mesh->n_vertices(); i++)
        input << (int) modelingAll[i] << " " << (int) handleAll[i] << endl;

#ifdef ENABLE_TSPLINEMESH_SUPPORT
    } else if ( object->dataType() == DATA_TSPLINE_MESH){

      TSplineMesh* mesh = PluginFunctions::tsplineMesh(object);

      //header
      input << "Selection" << endl;
      input << mesh->n_vertices() << endl;

      std::vector< int > modelingVertices = MeshSelection::getArea(mesh, AREA);
      std::vector< int > handleVertices   = MeshSelection::getArea(mesh, HANDLEAREA);

      std::vector< bool > modelingAll(mesh->n_vertices(), false);
      std::vector< bool > handleAll(mesh->n_vertices(), false);

      for (uint i=0; i < modelingVertices.size(); i++)
        modelingAll[ modelingVertices[i] ] = true;

      for (uint i=0; i < handleVertices.size(); i++)
        handleAll[ handleVertices[i] ] = true;

      for (uint i=0; i < mesh->n_vertices(); i++)
        input << (int) modelingAll[i] << " " << (int) handleAll[i] << endl;

#endif
    } else {
      emit log(LOGERR, tr("saveFlipperModelingSelection : Unsupported Type."));
    }
  } else
    emit log(LOGERR,tr("Unable to open selection file!"));
}


