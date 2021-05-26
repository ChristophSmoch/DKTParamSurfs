#ifndef __INCLUDESVTK_H
#define __INCLUDESVTK_H


// #ifdef PESOPT_WITH_VTK

#include "vtkAutoInit.h" 
VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
VTK_MODULE_INIT(vtkInteractionStyle);

#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkAppendPolyData.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkCellTypes.h>
#include <vtkClipPolyData.h>
#include <vtkColorTransferFunction.h>
#include <vtkContourFilter.h>
#include <vtkContourTriangulator.h>
#include <vtkDataReader.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetReader.h>
#include <vtkDecimatePro.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkFloatArray.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkGenericDataObjectWriter.h>
#include <vtkImageData.h>
#include <vtkImageMapper.h>
#include <vtkImageSlice.h>
#include <vtkImageToStructuredPoints.h>
#include <vtkImageViewer2.h>
#include <vtkInteractorStyleImage.h>
#include <vtkKdTreePointLocator.h>
#include <vtkLight.h>
#include <vtkLightKit.h>
#include <vtkLookupTable.h>
#include <vtkLoopSubdivisionFilter.h>
#include <vtkMapper.h>
#include <vtkMarchingCubes.h>
#include <vtkMath.h>
#include <vtkNamedColors.h>
#include <vtkPiecewiseFunction.h>
#include <vtkPLYReader.h>
#include <vtkPNGReader.h>
#include <vtkPNGWriter.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkProperty.h>
#include <vtkPostScriptWriter.h>
#include <vtkQuadricDecimation.h>
#include <vtkRectilinearGrid.h>
#include <vtkRectilinearGridWriter.h> 
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkScalarBarActor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridVolumeMapper.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkVersion.h>
#include <vtkVoxelContoursToSurfaceFilter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkXMLCompositeDataReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLReader.h>
#include <vtkXMLRectilinearGridReader.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>




//! data can either belong to vertices or to faces
enum VTKDataSupp { VERTEX_DATA, FACE_DATA };
//! vector-valued data can be saved as vectors, normals or texture coordinates (the file format also supports color scalars and lookup tables, which are not use here).
enum VTKVectorSpec { VECTORS, NORMALS };

enum VTKDATATYPE {
     VTKUNSTRUCTUREDGRID = 0,
     VTKSTRUCTUREDGRID = 1,
     VTKPOLYDATA = 2,
     VTKRECTILINEARGRID = 3,
     VTKSTRUCTUREDPOINTS = 4,
     VTKIMAGEDATA = 5
};


#endif
