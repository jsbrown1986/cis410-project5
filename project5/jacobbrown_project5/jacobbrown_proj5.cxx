/*=========================================================================

				Project 5 for CIS 410 (W18)
				Isoline code and Interpolation function
				implemented by Jacob Brown 2/9/2018

===========================================================================*/
/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#define ISOVALUE 3.2f

// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}


class SegmentList
{
   public:
                   SegmentList() { maxSegments = 10000; segmentIdx = 0; pts = new float[4*maxSegments]; };
     virtual      ~SegmentList() { delete [] pts; };

     void          AddSegment(float X1, float Y1, float X2, float Y2);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxSegments;
     int           segmentIdx;
};

void
SegmentList::AddSegment(float X1, float Y1, float X2, float Y2)
{
    pts[4*segmentIdx+0] = X1;
    pts[4*segmentIdx+1] = Y1;
    pts[4*segmentIdx+2] = X2;
    pts[4*segmentIdx+3] = Y2;
    segmentIdx++;
}

vtkPolyData *
SegmentList::MakePolyData(void)
{
    int nsegments = segmentIdx;
    int numPoints = 2*(nsegments);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nsegments ; i++)
    {
        double pt[3];
        pt[0] = pts[4*i];
        pt[1] = pts[4*i+1];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[4*i+2];
        pt[1] = pts[4*i+3];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx+1, pt);
        vtkIdType ids[2] = { ptIdx, ptIdx+1 };
        lines->InsertNextCell(2, ids);
        ptIdx += 2;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetLines(lines);
    lines->Delete();
    vtk_pts->Delete();

    return pd;
}

// Interpolation Function
float Interpolation(float value, float a, float b, float fa, float fb) {
	return (b - a) != 0 ? fa + ((value - a) / (b - a)) * (fb - fa) : 0;
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj5.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    // Add 4 segments that put a frame around your isolines.  This also
    // documents how to use "AddSegment".
    SegmentList sl;
    sl.AddSegment(-10, -10, +10, -10); // Add segment (-10,-10) -> (+10, -10)
    sl.AddSegment(-10, +10, +10, +10);
    sl.AddSegment(-10, -10, -10, +10);
    sl.AddSegment(+10, -10, +10, +10);

	// **Beginning of my code**

	int no_cells = GetNumberOfCells(dims);
	int no_segments, case_id, v0, v1, v2, v3;

	int bottom_left_log[2], bottom_right_log[2], top_left_log[2], top_right_log[2];

	float bottom_left_f, bottom_right_f, top_left_f, top_right_f;

	int numSegments[16] = { 0, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 0 };

	// Note: same output as screenshot.png from the website.
	// For same output as screenshot2.png, swap cases #6
	// and #9
	int lookup[16][4] = {
		{-1, -1, -1, -1}, //0
		{+0, +3, -1, -1}, //1
		{+0, +1, -1, -1}, //2
		{+1, +3, -1, -1}, //3
		{+2, +3, -1, -1}, //4
		{+0, +2, -1, -1}, //5
		{+0, +3, +1, +2}, //6
		{+1, +2, -1, -1}, //7
		{+1, +2, -1, -1}, //8
		{+0, +1, +2, +3}, //9
		{+0, +2, -1, -1}, //10
		{+2, +3, -1, -1}, //11
		{+1, +3, -1, -1}, //12
		{+0, +1, -1, -1}, //13
		{+0, +3, -1, -1}, //14
		{-1, -1, -1, -1} //15
	};

	// Runs through every cell
	for (int i = 0; i < no_cells; i++) {

		// Initialize case ID
		case_id = 0;

		// Determines the logical cell indices based on
		// current cell number
		GetLogicalCellIndex(bottom_left_log, i, dims);

		bottom_right_log[0] = bottom_left_log[0] + 1;
		bottom_right_log[1] = bottom_left_log[1];

		top_left_log[0] = bottom_left_log[0];
		top_left_log[1] = bottom_left_log[1] + 1;

		top_right_log[0] = bottom_left_log[0] + 1;
		top_right_log[1] = bottom_left_log[1] + 1;

		// Grabs the field values at every index
		bottom_left_f = F[GetPointIndex(bottom_left_log, dims)];
		bottom_right_f = F[GetPointIndex(bottom_right_log, dims)];
		top_left_f = F[GetPointIndex(top_left_log, dims)];
		top_right_f = F[GetPointIndex(top_right_log, dims)];

		// Compares field values to isovalue to determine
		// which case to use
		bottom_left_f < ISOVALUE ? v0 = 0 : v0 = 1;
		bottom_right_f < ISOVALUE ? v1 = 0 : v1 = 1;
		top_left_f < ISOVALUE ? v2 = 0 : v2 = 1;
		top_right_f < ISOVALUE ? v3 = 0 : v3 = 1;
		
		// "Binary" to find case number
		v0 == 1 ? case_id += 1 : case_id += 0;
		v1 == 1 ? case_id += 2 : case_id += 0;
		v2 == 1 ? case_id += 4 : case_id += 0;
		v3 == 1 ? case_id += 8 : case_id += 0;

		// Grabs number of segments based on the case
		no_segments = numSegments[case_id];

		// For each segment in a specified case...
		for (int j = 0; j < no_segments; j++) {

			// Grabs the edge number based on the case number
			int edge1 = lookup[case_id][2 * j];
			int edge2 = lookup[case_id][(2 * j) + 1];

			float pt1[2], pt2[2];

			// Interpolates first point based on edge location
			if (edge1 == 0) {
				pt1[0] = Interpolation(ISOVALUE, bottom_left_f, bottom_right_f, X[bottom_left_log[0]], X[bottom_right_log[0]]);
				pt1[1] = Y[bottom_left_log[1]];
			} else if (edge1 == 1) {
				pt1[0] = X[bottom_right_log[0]];
				pt1[1] = Interpolation(ISOVALUE, bottom_right_f, top_right_f, Y[bottom_right_log[1]], Y[top_right_log[1]]);
			} else if (edge1 == 2) {
				pt1[0] = Interpolation(ISOVALUE, top_left_f, top_right_f, X[top_left_log[0]], X[top_right_log[0]]);
				pt1[1] = Y[top_left_log[1]];
			} else if (edge1 == 3) {
				pt1[0] = X[bottom_left_log[0]];
				pt1[1] = Interpolation(ISOVALUE, bottom_left_f, top_left_f, Y[bottom_left_log[1]], Y[top_left_log[1]]);
			}

			// Interpolates second point based on edge location
			if (edge2 == 0) {
				pt2[0] = Interpolation(ISOVALUE, bottom_left_f, bottom_right_f, X[bottom_left_log[0]], X[bottom_right_log[0]]);
				pt2[1] = Y[bottom_left_log[1]];
			} else if (edge2 == 1) {
				pt2[0] = X[bottom_right_log[0]];
				pt2[1] = Interpolation(ISOVALUE, bottom_right_f, top_right_f, Y[bottom_right_log[1]], Y[top_right_log[1]]);
			} else if (edge2 == 2) {
				pt2[0] = Interpolation(ISOVALUE, top_left_f, top_right_f, X[top_left_log[0]], X[top_right_log[0]]);
				pt2[1] = Y[top_left_log[1]];
			} else if (edge2 == 3) {
				pt2[0] = X[bottom_left_log[0]];
				pt2[1] = Interpolation(ISOVALUE, bottom_left_f, top_left_f, Y[bottom_left_log[1]], Y[top_left_log[1]]);
			}

			// Adds line segment to be drawn
			sl.AddSegment(pt1[0], pt1[1], pt2[0], pt2[1]);
		}
	}

	// **End of my code**

    vtkPolyData *pd = sl.MakePolyData();

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
