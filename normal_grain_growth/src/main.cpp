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

#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetMapper.h>
#include <vtkCamera.h>

#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>

vtkSmartPointer<vtkUnstructuredGrid> allocate_vtk_grid(int NX, int NY, int NZ);

using namespace std;

int main(int argc, char *argv[])
{
	cout << "Main program has started ..." << endl;

	int numberOfStates = 1000000;

	// Define the VTK grid
	// Number of cells = number of cellular automata cells, I won't be plotting point data
	int nCellsX = 200;
	int nCellsY = 200;
	int nCellsZ = 200;
	
	// Allocate vector
	vector<vector<vector<int>>> cell_state;
	cell_state.resize(nCellsX);
	for (int i=0; i<nCellsX; i++)			
	{
		cell_state[i].resize(nCellsY);
		for (int j=0; j<nCellsY; j++)
			cell_state[i][j].resize(nCellsZ);
	}

	// Initialize vector
	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0,numberOfStates-1);
	std::uniform_int_distribution<int> distribution_i(1,nCellsX-2);	//exclude boundaries
	std::uniform_int_distribution<int> distribution_j(1,nCellsY-2);
	std::uniform_int_distribution<int> distribution_k(1,nCellsZ-2);

	// Initialize with random noise
	int l = 0;
	for (int i=0; i<nCellsX; i++)
	for (int j=0; j<nCellsY; j++)
	for (int k=0; k<nCellsZ; k++)
	{
		cell_state[i][j][k] = distribution(generator);
	}

	// March in time
	int frame = 0;
	for (int t=0; t<100000001; t++)
	{
		cout << "t = " << t << endl;

		// Pick a random lattice
		int i = distribution_i(generator);
		int j = distribution_j(generator);
		int k = distribution_j(generator);
		cout << "[i,j,k] = " << i << ", " << j << ", " << k << endl;

		// Identify the nearest neighbours
		std::vector<int> q;
		std::vector<int> q_unique;

		q.push_back(cell_state[i-1][j-1][k-1]);
		q.push_back(cell_state[i-1][j-1][k]);
		q.push_back(cell_state[i-1][j-1][k+1]);
		q.push_back(cell_state[i-1][j][k-1]);
		q.push_back(cell_state[i-1][j][k]);
		q.push_back(cell_state[i-1][j][k+1]);
		q.push_back(cell_state[i-1][j+1][k-1]);
		q.push_back(cell_state[i-1][j+1][k]);
		q.push_back(cell_state[i-1][j+1][k+1]);
		q.push_back(cell_state[i][j-1][k-1]);
		q.push_back(cell_state[i][j-1][k]);
		q.push_back(cell_state[i][j-1][k+1]);
		q.push_back(cell_state[i][j][k-1]);
//		q.push_back(cell_state[i][j][k]);	// exclude the state itself
		q.push_back(cell_state[i][j][k+1]);
		q.push_back(cell_state[i][j+1][k-1]);
		q.push_back(cell_state[i][j+1][k]);
		q.push_back(cell_state[i][j+1][k+1]);
		q.push_back(cell_state[i+1][j-1][k-1]);
		q.push_back(cell_state[i+1][j-1][k]);
		q.push_back(cell_state[i+1][j-1][k+1]);
		q.push_back(cell_state[i+1][j][k-1]);
		q.push_back(cell_state[i+1][j][k]);
		q.push_back(cell_state[i+1][j][k+1]);
		q.push_back(cell_state[i+1][j+1][k-1]);
		q.push_back(cell_state[i+1][j+1][k]);
		q.push_back(cell_state[i+1][j+1][k+1]);


		// Make q unique, name it q_unique
		for (int qq=0; qq<q.size(); qq++)
		{
			int this_state = q[qq];
			bool flag_contains = false;
			for (int qqu=0; qqu<q_unique.size(); qqu++)
			{
				if (q_unique[qqu] == this_state)
				flag_contains = true;
				continue;
			}
	
			if (flag_contains == false)
				q_unique.push_back(this_state);			
		}
		
/*
		// Print q and q unique (for testing)
		cout << "Print results (testing) ..." << endl;
		for (int qq=0; qq<q.size(); qq++)
			cout << "q[ " << qq << " ] = " << q[qq] << endl; 
		for (int qq=0; qq<q_unique.size(); qq++)
			cout << "q_unique[ " << qq << " ] = " << q_unique[qq] << endl; 
*/
		
		// Calculate the current energy:
		int E = 0;
		for (int qq=0; qq<q.size(); qq++)
		if (cell_state[i][j][k] == q[qq])
		{
			// do nothing
		}
		else
		{
			// add energy
			E++;
		}


		// Randomly assign one of the unlike neighbours
		std::uniform_int_distribution<int> distribution_neighbour(0,q_unique.size()-1);
		int new_state = q_unique[distribution_neighbour(generator)];
//		cout << "new_state = " << new_state << endl;

		// Calculate the energy change:
		int E_new = 0;
		for (int qq=0; qq<q.size(); qq++)
		if (new_state == q[qq])
		{
			// do nothing
		}
		else
		{
			// add energy
			E_new++;
		}
		int dE = E_new - E;

		cout << "E = " << E << endl;
		cout << "E_new = " << E_new << endl;
		cout << "dE = " << dE << endl;

		// Calculate the probability of acquiring this state
		double P = 0.5*(1.0-tanh(double(dE)));
		cout << "P = " << P << endl;
		std::vector<bool> pp;
		for (int ppp=0; ppp<100; ppp++)
			pp.push_back(true);

		if (P == 0.0)
		{
			for (int ppp=0; ppp<100; ppp++)
			pp[ppp] =  false;
		}
		else
		{
			for (int ppp=int(P*100); ppp<100; ppp++)
			pp[ppp] = false;
		}

		std::uniform_int_distribution<int> distribution_change(0,99);
		bool changeTakesPlace = pp[distribution_change(generator)];
		cout << "changeTakesPlace = " << changeTakesPlace << endl;

		if (changeTakesPlace)
			cell_state[i][j][k] = new_state;

		// Periodic boundary conditions
		for (int i=0; i<nCellsX; i++)
		for (int j=0; j<nCellsY; j++)
		{
			cell_state[i][j][0] = cell_state[i][j][nCellsZ-2];
			cell_state[i][j][nCellsZ-1] = cell_state[i][j][1];
		}
		for (int j=0; j<nCellsY; j++)
		for (int k=0; k<nCellsZ; k++)
		{
			cell_state[0][j][k] = cell_state[nCellsX-2][j][k];
			cell_state[nCellsX-1][j][k] = cell_state[1][j][k];
		}
		for (int i=0; i<nCellsX; i++)
		for (int k=0; k<nCellsZ; k++)
		{
			cell_state[i][0][k] = cell_state[i][nCellsY-2][k];
			cell_state[i][nCellsY-1][k] = cell_state[i][1][k];
		}

		// [[ POST-PROCESSING ]]
		if (t%100000 == 0)
		{
			cout << "post-processing ... " << endl;
			vtkSmartPointer<vtkUnstructuredGrid> grid = allocate_vtk_grid(nCellsX, nCellsY, nCellsZ);
			vtkSmartPointer<vtkIntArray> vtk_cell_state = vtkSmartPointer<vtkIntArray>::New();
			vtk_cell_state->SetName("cell state");
			vtk_cell_state->SetNumberOfComponents(1);
			vtk_cell_state->SetNumberOfTuples(nCellsX*nCellsY*nCellsZ);

			int l = 0;
			for (int k=0; k<nCellsZ; k++)
			for (int j=0; j<nCellsY; j++)
			for (int i=0; i<nCellsX; i++)
			{
				vtk_cell_state->SetTuple1(l,cell_state[i][j][k]);
				l++;
			}

			grid->GetCellData()->AddArray(vtk_cell_state);
			grid->GetCellData()->SetActiveScalars("cell state");

			vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
			mapper->SetInputData(grid);
			mapper->ScalarVisibilityOn();
			mapper->SetScalarModeToUseCellData();
			mapper->SetScalarRange(0, numberOfStates-1);
		//	mapper->SetLookupTable(lookupTable);

			vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);
			actor->RotateX(25.0);
			actor->RotateY(-25.0);

			vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
			renderer->AddActor(actor);
			renderer->ResetCamera();
			renderer->GetActiveCamera()->Zoom(1.5);
			
			vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
			renderWindow->AddRenderer(renderer);

			vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
			renderWindowInteractor->SetRenderWindow(renderWindow);

			renderWindow->Render();

			
			cout << "Write into a file ..." << endl;
			vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
			windowToImageFilter->SetInput(renderWindow);
			windowToImageFilter->SetInputBufferTypeToRGBA();
			windowToImageFilter->ReadFrontBufferOff();
			windowToImageFilter->SetScale(4,4);
			windowToImageFilter->Update();

			std::string std_frame; 
			std_frame = to_string(frame);
			int std_frame_size = std_frame.size();
			for (int s=0; s<4-std_frame_size; s++)
				std_frame.insert(0, "0");


			vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
			writer->SetFileName(("file-"+std_frame+".png").c_str());
			writer->SetInputConnection(windowToImageFilter->GetOutputPort());
			writer->Write();
			frame++;
		}
	//	renderWindowInteractor->Start();

	}

//	usleep(5000000);


	cout << "Main program has successfully finished!" << endl;
	return EXIT_SUCCESS;
}

vtkSmartPointer<vtkUnstructuredGrid> allocate_vtk_grid(int NX, int NY, int NZ)
{
	vtkSmartPointer<vtkUnstructuredGrid> vtkGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkGrid->Allocate(NX*NY*NZ);

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (int k=0;k<NZ+1;k++){
	for (int j=0;j<NY+1;j++){
	for (int i=0;i<NX+1;i++){
			double point_x = double(i)/double(NX);
			double point_y = double(j)/double(NY);
			double point_z = double(k)/double(NZ);

			points->InsertNextPoint(point_x,point_y,point_z);
	}
	}
	}
    vtkGrid->SetPoints(points);

	vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
	for (int k=0;k<NZ;k++){
	for (int j=0;j<NY;j++){
	for (int i=0;i<NX;i++){
	
		hex->GetPointIds()->SetId(0,i+(NX+1)*j+(NX+1)*(NY+1)*k);
		hex->GetPointIds()->SetId(1,i+(NX+1)*j+(NX+1)*(NY+1)*k + 1);
		hex->GetPointIds()->SetId(2,i+(NX+1)*j+(NX+1)*(NY+1)*k + (NX+1)+1);
		hex->GetPointIds()->SetId(3,i+(NX+1)*j+(NX+1)*(NY+1)*k + (NX+1));
		hex->GetPointIds()->SetId(4,i+(NX+1)*j+(NX+1)*(NY+1)*k + (NX+1)*(NY+1));
		hex->GetPointIds()->SetId(5,i+(NX+1)*j+(NX+1)*(NY+1)*k + (NX+1)*(NY+1)+1);
		hex->GetPointIds()->SetId(6,i+(NX+1)*j+(NX+1)*(NY+1)*k + (NX+1)*(NY+1) + (NX+1)+1);
		hex->GetPointIds()->SetId(7,i+(NX+1)*j+(NX+1)*(NY+1)*k + (NX+1)*(NY+1) + (NX+1));

		vtkGrid->InsertNextCell(hex->GetCellType(),hex->GetPointIds());
	}
	}
	}

return vtkGrid;
}

