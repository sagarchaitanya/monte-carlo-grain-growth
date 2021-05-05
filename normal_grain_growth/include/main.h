#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <vector>
#include <random>
#include <unistd.h>

#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkCellData.h>
#include <vtkHexahedron.h>
#include <vtkSmartPointer.h>
#include <vtkIntArray.h>
#include <vtkCellArray.h>

#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetMapper.h>
#include <vtkCamera.h>
#include <vtkPolyDataMapper.h>
#include <vtkOutlineFilter.h>
#include <vtkProperty.h>

#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>

using namespace std;

// [[ MEMBER FUNCTIONS ]]
void preprocessor();
void run();
void postprocessor();

// VTK grid:
vtkSmartPointer<vtkUnstructuredGrid> allocate_vtk_grid(int NX, int NY, int NZ);

// Vector arrays:
vector<vector<vector<int>>> cell_state;		// state of the cell (grain orientation)
vector<vector<vector<bool>>> boundary;		// if true == boundary, if false == grain interior


// [[ MAPPERS ]]
vtkSmartPointer<vtkDataSetMapper> mapper_grid;
vtkSmartPointer<vtkDataSetMapper> mapper_selected_state;
vtkSmartPointer<vtkDataSetMapper> mapper_boundary;
vtkSmartPointer<vtkPolyDataMapper> mapper_outline;

// [[ ACTORS ]]
vtkSmartPointer<vtkActor> actor_grid;
vtkSmartPointer<vtkActor> actor_selected_state;
vtkSmartPointer<vtkActor> actor_boundary;
vtkSmartPointer<vtkActor> actor_outline;

vtkSmartPointer<vtkRenderer> renderer;
vtkSmartPointer<vtkRenderWindow> renderWindow;
vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;

// [[ PNG WRITER ]]
//vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter;
//vtkSmartPointer<vtkPNGWriter> writer;

// [[ SCALARS ]]
int nCellsX;
int nCellsY;
int nCellsZ;

int numberOfStates;

std::default_random_engine generator;

int itr = 0;
int frame = 0;

double rotation_rate = 0.00001;

#endif
