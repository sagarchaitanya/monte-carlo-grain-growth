#include "main.h"

int main(int argc, char *argv[])
{
	cout << "Main program has started ..." << endl;

	preprocessor();
	run();

	cout << "Main program has successfully finished!" << endl;
	return EXIT_SUCCESS;
}

void run()
{
	std::uniform_int_distribution<int> distribution_i(1,nCellsX-2);	//exclude boundaries
	std::uniform_int_distribution<int> distribution_j(1,nCellsY-2);
	std::uniform_int_distribution<int> distribution_k(1,nCellsZ-2);

	// [[ MAIN LOOP IN TIME ]]
	// March in time
	itr = 0;
	frame = 0;
	while(1)
	{
		cout << "itr = " << itr << endl;

		// Pick a random lattice
		int i = distribution_i(generator);
		int j = distribution_j(generator);
		int k = distribution_k(generator);
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
		if (itr%100000 == 0)
		{
			postprocessor();	
		}
		itr++;
	}

	return;
}

void preprocessor()
{
	numberOfStates = 1000;

	// Define the VTK grid
	// Number of cells = number of cellular automata cells, I won't be plotting point data
	nCellsX = 50;
	nCellsY = 50;
	nCellsZ = 50;
	
	// Allocate vector
	cell_state.resize(nCellsX);
	boundary.resize(nCellsX);
	for (int i=0; i<nCellsX; i++)			
	{
		cell_state[i].resize(nCellsY);
		boundary[i].resize(nCellsY);
		for (int j=0; j<nCellsY; j++)
		{
			cell_state[i][j].resize(nCellsZ);
			boundary[i][j].resize(nCellsZ);
		}
	}

	// Initialize vector
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
		boundary[i][j][k] = false;
	}

	return;
}

void postprocessor()
{
	cout << "post-processing ... " << endl;

	for (int i=1; i<nCellsX-1; i++)
	for (int j=1; j<nCellsY-1; j++)
	for (int k=1; k<nCellsZ-1; k++)
	{
		int thisState = cell_state[i][j][k];
		std::vector<int> q;
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

		boundary[i][j][k] = false;
		for (int qq=0; qq<q.size(); qq++)
			if (thisState != q[qq])
			{
				boundary[i][j][k] = true;
				break;
			}
	}
	// Periodic boundary conditions
	for (int i=0; i<nCellsX; i++)
	for (int j=0; j<nCellsY; j++)
	{
		boundary[i][j][0] = boundary[i][j][nCellsZ-2];
		boundary[i][j][nCellsZ-1] = boundary[i][j][1];
	}
	for (int j=0; j<nCellsY; j++)
	for (int k=0; k<nCellsZ; k++)
	{
		boundary[0][j][k] = boundary[nCellsX-2][j][k];
		boundary[nCellsX-1][j][k] = boundary[1][j][k];
	}
	for (int i=0; i<nCellsX; i++)
	for (int k=0; k<nCellsZ; k++)
	{
		boundary[i][0][k] = boundary[i][nCellsY-2][k];
		boundary[i][nCellsY-1][k] = boundary[i][1][k];
	}

	vtkSmartPointer<vtkUnstructuredGrid> grid = allocate_vtk_grid(nCellsX, nCellsY, nCellsZ);

	vtkSmartPointer<vtkIntArray> vtk_cell_state = vtkSmartPointer<vtkIntArray>::New();
	vtk_cell_state->SetName("cell state");
	vtk_cell_state->SetNumberOfComponents(1);
	vtk_cell_state->SetNumberOfTuples(nCellsX*nCellsY*nCellsZ);

	// Unstructure grid of one cell state
	vtkSmartPointer<vtkUnstructuredGrid> grid_cell = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkCellArray> cell_array = vtkSmartPointer<vtkCellArray>::New();			

	vtkSmartPointer<vtkUnstructuredGrid> grid_boundary = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkCellArray> cell_array_boundary = vtkSmartPointer<vtkCellArray>::New();

	grid_cell->SetPoints(grid->GetPoints());
	grid_boundary->SetPoints(grid->GetPoints());

	// Find the dominant state (the most number of cells)
	std::vector<int> numberOfCells;
	numberOfCells.resize(numberOfStates);
	for (int i=0; i<numberOfStates; i++)
		numberOfCells[i] = 0;

	for (int i=0; i<nCellsX; i++)
	for (int j=0; j<nCellsY; j++)
	for (int k=0; k<nCellsZ; k++)
	{
		numberOfCells[ cell_state[i][j][k] ]++;
	}

	int positionOfMax = 0;
	int maxNumberOfCells = 0;
	for (int i=0; i<numberOfStates; i++)
		if (numberOfCells[i] > maxNumberOfCells)
		{
			maxNumberOfCells = numberOfCells[i];
			positionOfMax = i;
		}

	int selected_state = positionOfMax;
	for (int i=0; i<nCellsX; i++)
	for (int j=0; j<nCellsY; j++)
	for (int k=0; k<nCellsZ; k++)
		if (cell_state[i][j][k] == selected_state)
		{
			int thisCellId = nCellsX*nCellsY*k + nCellsY*j +i;
			vtkCell *cell = grid->GetCell(thisCellId);
			cell_array->InsertNextCell(cell);
		}	
	grid_cell->SetCells(VTK_HEXAHEDRON, cell_array);


	for (int i=0; i<nCellsX; i++)
	for (int j=0; j<nCellsY; j++)
	for (int k=0; k<nCellsZ; k++)
	{
		if (boundary[i][j][k])
		{
			int thisCellId = nCellsX*nCellsY*k + nCellsY*j +i;
			vtkCell *cell = grid->GetCell(thisCellId);
			cell_array_boundary->InsertNextCell(cell);
		}
	}
	grid_boundary->SetCells(VTK_HEXAHEDRON, cell_array_boundary);

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


	mapper_grid = vtkSmartPointer<vtkDataSetMapper>::New();
	mapper_grid->SetInputData(grid);
	mapper_grid->ScalarVisibilityOn();
	mapper_grid->SetScalarModeToUseCellData();
	mapper_grid->SetScalarRange(0, numberOfStates-1);
//	mapper_grid->SetLookupTable(lookupTable);

	mapper_selected_state = vtkSmartPointer<vtkDataSetMapper>::New();
	mapper_selected_state->SetInputData(grid_cell);
	mapper_selected_state->ScalarVisibilityOn();
	mapper_selected_state->SetScalarModeToUseCellData();
	mapper_selected_state->SetScalarRange(0, numberOfStates-1);

	mapper_boundary = vtkSmartPointer<vtkDataSetMapper>::New();
	mapper_boundary->SetInputData(grid_boundary);
	mapper_boundary->ScalarVisibilityOn();
	mapper_boundary->SetScalarModeToUseCellData();
	mapper_boundary->SetScalarRange(0, numberOfStates-1);

	vtkSmartPointer<vtkActor> actor_grid = vtkSmartPointer<vtkActor>::New();
	actor_grid->SetMapper(mapper_grid);
	actor_grid->RotateX(25.0);
	actor_grid->RotateY(-25.0+0.000001*double(itr));
	actor_grid->GetProperty()->SetOpacity(0.5);

	actor_selected_state = vtkSmartPointer<vtkActor>::New();
	actor_selected_state->SetMapper(mapper_selected_state);
	actor_selected_state->RotateX(25.0);
	actor_selected_state->RotateY(-25.0+0.000001*double(itr));
	actor_selected_state->GetProperty()->SetColor(0.34,0.17,0.94);

	actor_boundary = vtkSmartPointer<vtkActor>::New();
	actor_boundary->SetMapper(mapper_boundary);
	actor_boundary->RotateX(25.0);
	actor_boundary->RotateY(-25.0+0.000001*double(itr));
	actor_boundary->GetProperty()->SetOpacity(0.02);
	actor_boundary->GetProperty()->SetColor(0.34,0.17,0.94);

	renderer = vtkSmartPointer<vtkRenderer>::New();

	mapper_outline = vtkSmartPointer<vtkPolyDataMapper>::New();	

	vtkSmartPointer<vtkOutlineFilter> outline = vtkSmartPointer<vtkOutlineFilter>::New();
	outline->SetInputData(grid);
	mapper_outline->SetInputConnection(outline->GetOutputPort());

	actor_outline = vtkSmartPointer<vtkActor>::New();
	actor_outline->SetMapper(mapper_outline);
	actor_outline->RotateX(25.0);
	actor_outline->RotateY(-25.0+0.000001*double(itr));

	// Decide what is shown in plot : 
//			renderer->AddActor(actor_grid);					// plot the full grid (all states)
	renderer->AddActor(actor_selected_state);			// plot only cells with one selected state
	renderer->AddActor(actor_boundary);			// plot grain boundary
	renderer->AddActor(actor_outline);			// plot outline
	renderer->ResetCamera();
	renderer->GetActiveCamera()->Zoom(1.5);
	
	renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);

	renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	renderWindow->Render();
	renderWindow->WaitForCompletion();

	cout << "Write into a file ..." << endl;
	windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
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

	writer = vtkSmartPointer<vtkPNGWriter>::New();
	writer->SetFileName(("../screenshots/file-"+std_frame+".png").c_str());
	writer->SetInputConnection(windowToImageFilter->GetOutputPort());
	writer->Write();
	frame++;

	renderWindow->Finalize();

	return;
}


vtkSmartPointer<vtkUnstructuredGrid> allocate_vtk_grid(int NX, int NY, int NZ)
{
	vtkSmartPointer<vtkUnstructuredGrid> vtkGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkGrid->Allocate(NX*NY*NZ);

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (int k=0;k<NZ+1;k++){
	for (int j=0;j<NY+1;j++){
	for (int i=0;i<NX+1;i++){
			double point_x = double(i)/double(1);
			double point_y = double(j)/double(1);
			double point_z = double(k)/double(1);

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

