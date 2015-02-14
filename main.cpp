/*=========================================================================

  Visualization Term Project
  Instructor:		Dr. Xiao
  Name		:		Adithya Addanki
  UA Net ID	:		aa207
  
==========================================================================*/
//
// Oct-Tree Preprocessing- Scattered Data Visualization
//
/*

Step 1. Preprocessing Data using Oct-Tree analysis.

The whole grid ranging co-ordinate values for x,y,z from [0,1] is considered as a cube.
This grid is split into octants based on the granularity required by the user, and then sent through shepard method for interpolation into scientific data

Step 2. Modeling.

Use the Shepard’s method (with alpha index = 2) to interpolate the scattered data into a uniform grid and output it into a
text file in the same format as the head data. Allow the user to specify the dimensions (nx,ny,nz) of the grid at the start of your program.
The xyzv ranges are [0.0,1.0][0.0,1.0][0.0,1.0][0.0,1.0].

Step 3. Rendering.

Follow the head example (below) to read the above grid data. 
Create and display isosurfaces using a contour filter.  
You should allow the user to use the “i” key to increase the contour value and the “d” key to decrease the contour value. 
In a separate viewport, display the operation instructions.

*/


// Include the required header files for the VTK classes that are used in the program
#include "stdlib.h"
#include "vtkActor.h"
#include "vtkActorCollection.h"
#include "vtkBoxWidget.h"
#include "vtkCamera.h"
#include "vtkColor.h"
#include "vtkProperty.h"
#include "vtkCommand.h"
#include "vtkContourFilter.h"
#include "vtkCubeSource.h"
#include "vtkDelimitedTextReader.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkLight.h"
#include "vtkTextActor.h"
#include "vtkOutlineFilter.h"
#include "vtkLookupTable.h"
#include "vtkParticleReader.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkProperty2D.h"
#include "vtkRendererCollection.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include <vtkObjectFactory.h>
#include "vtkStripper.h"
#include "vtkStructuredPointsReader.h"
#include "vtkTable.h"
#include "vtkTransform.h"
#include "vtkVolume16Reader.h"
#include <deque>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <vector>
#include <vtkConeSource.h>
#include <vtkContourWidget.h>
#include <vtkOrientedGlyphContourRepresentation.h>
#include <vtkWidgetEvent.h>
#include <vtkWidgetEventTranslator.h>
#include <windows.h>

using namespace std;

struct sc_pt 
{
	double x,y,z,value;
};

float contVal=0.1;
vtkLookupTable *lut=vtkLookupTable::New();
	

class vtkMyCallback : public vtkCommand
{
public:
  static vtkMyCallback *New() 
    { return new vtkMyCallback; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkTransform *t = vtkTransform::New();
      vtkBoxWidget *widget = reinterpret_cast<vtkBoxWidget*>(caller);
      widget->GetTransform(t);
      widget->GetProp3D()->SetUserTransform(t);
      t->Delete();
    }
};

class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
  public:
    static KeyPressInteractorStyle* New();
    vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);
 
    virtual void OnKeyPress() 
    {
      // Get the keypress
      vtkRenderWindowInteractor *rwi = this->Interactor;
      std::string key = rwi->GetKeySym();
	  vtkRenderWindow *rw=rwi->GetRenderWindow();
	  vtkRendererCollection *rc=rw->GetRenderers();
	  rc->InitTraversal();
	  vtkRenderer *r=rc->GetNextItem();
	  vtkActorCollection *vac=r->GetActors();
	  vac->InitTraversal();
	  vtkActor *va=vac->GetNextActor();

      // Output the key that was pressed
      //std::cout << "Pressed " << key << std::endl;
 
      // Handle i
      if(key == "i")
        {
			std::cout << "i was pressed." << std::endl;
			contVal=contVal+0.001;
		}
 
      // Handle d
      if(key == "d")
        {
			std::cout << "d was pressed." << std::endl;
			contVal=contVal-0.001;
	    }
		vtkContourFilter *cF = vtkContourFilter::New();
		
		vtkStructuredPointsReader *v16 = vtkStructuredPointsReader::New();
		v16->SetFileName("processedData.vtk");

		cF->SetInputConnection(v16->GetOutputPort());
		cF->SetValue(1,contVal);

		vtkPolyDataNormals *cN = vtkPolyDataNormals::New();
		cN->SetInputConnection(cF->GetOutputPort());
		cN->SetFeatureAngle(60.0);

		vtkStripper *cS = vtkStripper::New();
		cS->SetInputConnection(cN->GetOutputPort());

		vtkPolyDataMapper *cM = vtkPolyDataMapper::New();
		cM->SetLookupTable(lut);
		cM->SetColorModeToDefault();
		cM->SetInputConnection(cS->GetOutputPort());
		va->SetMapper(cM);
		r->ResetCamera ();
		rw->Render();
      // Forward events
      vtkInteractorStyleTrackballCamera::OnKeyPress();
    }
};
vtkStandardNewMacro(KeyPressInteractorStyle);

double calcShepardVal(vector<sc_pt>, double,double,double,int,int);

int main()
{
	const char *filename="data.txt";
	cout<<"Reading Scattered Data from data.txt\n";
	vtkCollection *vc=vtkCollection::New();
	string line;
	ifstream myfile (filename);
	std::vector<sc_pt> vec;
	std::vector<double> sgValues;
	int count=0,nx,ny,nz;
	double dx=0,dy=0,dz=0;
	if (myfile.is_open())
	{    
		while ( getline (myfile,line) )
		{
			double x=0,y=0,z=0,v=0;
			istringstream sin(line);
			sc_pt sp;
			sin>>sp.x>>sp.y>>sp.z>>sp.value;
			vec.emplace_back(sp);
			++count;
		}
		myfile.close();
	}
	cout<<"Enter Nx,Ny,Nz:\n";
	cin>>nx>>ny>>nz;

	dx=(double)1/nx;
	dy=(double)1/ny;
	dz=(double)1/nz;
	int quadLevel=0;
	cout<<"\nPlease enter the oct tree level for constructing the grid values ranging 0-n\n";
	cout<<"   0-consider range 0-1 as a single cube; divide into octants in other cases\n";
	cin>>quadLevel;
	cout<<"Preparing data for constructing a grid; using Shepard Method\n";
	for(float k=0;k<1;k=k+dz)
	{
		for(float j=0;j<1;j=j+dy)
		{
			for(float i=0;i<1;i=i+dx)
			{
				sgValues.emplace_back(calcShepardVal(vec,i,j,k,2,quadLevel));
			}
		}
	}
	cout<<"Writing to file........\n";
	std::ofstream ofile;
	ofile.open("processedData.vtk");
	ofile<<"# vtk DataFile Version 3.0"<<endl;
	ofile<<"vtk output"<<endl;
	ofile<<"ASCII"<<endl;
	ofile<<"DATASET STRUCTURED_POINTS"<<endl;
	ofile<<"DIMENSIONS "<<nx<<" "<<ny<<" "<<nz<<endl;
	ofile<<"SPACING "<<dx<<" "<<dy<<" "<<dz<<endl;
	ofile<<"ORIGIN 0 0 0"<<endl;
	ofile<<"POINT_DATA "<<nx*ny*nz<<endl;
	ofile<<"SCALARS scalars float"<<endl;
	ofile<<"LOOKUP_TABLE default"<<endl;
	int ind=0;
	for(int k=0;k<nz;k++)
	{
		for(int j=0;j<ny;j++)
		{
			for(int i=0;i<nx;i++)
			{
				ofile<<(sgValues[ind])<<"  ";
				//cout<< (sgValues[ind])<<"  ";
				ind++;
			}
			ofile<<endl;
			//cout<<endl;
		}
	}
	ofile.close();
	cout<<"Reading the generated file; processedData.vtk\n";
	
	vtkContourFilter *cF = vtkContourFilter::New();
	vtkRenderWindow *renWin = vtkRenderWindow::New();
	vtkRenderer *aRenderer = vtkRenderer::New();
	aRenderer->SetViewport(0,0,1,0.6);
	vtkRenderer *iRenderer = vtkRenderer::New();
	iRenderer->SetViewport(0,0.6,1,1);
	vtkTextActor *vta3=vtkTextActor::New();
	vta3->SetInput("Visualization\nTerm Project: Oct Tree Analysis & Preprocessing for Scattered Data Visualization\nInstructor : Dr.Xiao\nName : Adithya Addanki\nUA Net ID: aa207 \n\nDocumentation:\n----------------------\n1. Read the file ; data.txt[x,y,z,v] \n2. Set dimensions of the grid to user specified dimensions (Uniform Grid)\n3. Oct tree level [Depth of Oct Tree] for constructing the grid values ranging 0-n\n            0-consider range 0-1 as a single cube; divide into octants in other cases\n4. Create a grid using the dimensions specified; Shepard Method\n5. Created a contour surface \n6. User could interact; press \"i\" to increase contour value and \"d\" to decrease it\n7. Grid renders itself with the new surfaces\n");
	iRenderer->AddActor(vta3);

	
	vtkStructuredPointsReader *v16 = vtkStructuredPointsReader::New();
	v16->SetFileName("processedData.vtk");

    cF->SetInputConnection(v16->GetOutputPort());
    cF->SetValue(0, 0.1);

	vtkPolyDataNormals *cN = vtkPolyDataNormals::New();
    cN->SetInputConnection(cF->GetOutputPort());

	vtkStripper *cS = vtkStripper::New();
    cS->SetInputConnection(cN->GetOutputPort());

	vtkPolyDataMapper *cM = vtkPolyDataMapper::New();
    cM->SetInputConnection(cS->GetOutputPort());
	cM->SetLookupTable(lut);
	cM->SetColorModeToDefault();

	vtkOutlineFilter *outlineData = vtkOutlineFilter::New();
    outlineData->SetInputConnection(v16->GetOutputPort());

	vtkPolyDataMapper *mapOutline = vtkPolyDataMapper::New();
    mapOutline->SetInputConnection(outlineData->GetOutputPort());

	vtkActor *outline = vtkActor::New();
    outline->SetMapper(mapOutline);
    outline->GetProperty()->SetColor(0,0,0);


	vtkActor *con = vtkActor::New();
    con->SetMapper(cM);
	aRenderer->AddActor(con);
	aRenderer->AddActor(outline);
    aRenderer->ResetCamera();
	aRenderer->SetBackground(255,255,255);


	vtkTextActor *vtaX=vtkTextActor::New();
	vtaX->SetInput("X-Axis->");
	vtaX->GetProperty()->SetColor(0,0,0);
	vtaX->SetDisplayPosition(0,0);
	
	vtkTextActor *vtaY=vtkTextActor::New();
	vtaY->SetInput("Y-Axis\n^\n|\n");
	vtaY->GetProperty()->SetColor(0,0,0);
	vtaY->SetDisplayPosition(0,5);
		
	//aRenderer->AddActor(vtaX);
	//aRenderer->AddActor(vtaY);

	//aRenderer->GetActiveCamera()->Elevation(90);
	renWin->AddRenderer(aRenderer);
	renWin->AddRenderer(iRenderer);
	renWin->SetSize(480, 570);

	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);
	KeyPressInteractorStyle *style = KeyPressInteractorStyle::New();
	iren->SetInteractorStyle(style);

	cout<<"Rendering the grid\n";
	renWin->Render();
	
	iren->Initialize();
    iren->Start(); 
	v16->Delete();
    cF->Delete();
    cN->Delete();
    cS->Delete();
    cM->Delete();
    con->Delete();
    aRenderer->Delete();

    renWin->Delete();
    iren->Delete();
 return 0;
}

double calcShepardVal(vector<sc_pt> vect, double x,double y,double z, int alpha, int quadLvl)
{
	double weightedDis=0;
	double distPow=0;
	double nm=0,dm=0;
	double dist=0;
	double rangeInt=(double)1/pow((double)2,quadLvl); 
	double tempx=0,tempy=0,tempz=0;
	int xt=0,yt=0,zt=0;
	double xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1;

	// finding octant range for x
	tempx=x;
	while (tempx >= rangeInt)
	{
		tempx-=rangeInt;
		++xt;
	}
	xmin=rangeInt*(xt);
	xmax=rangeInt*(xt+1);
	
	// finding octant range for y
	tempy=y;
	while (tempy >= rangeInt)
	{
		tempy-=rangeInt;
		++yt;
	}
	ymin=rangeInt*(yt);
	ymax=rangeInt*(yt+1);
	
	// finding octant range for z
	tempz=z;
	while (tempz >= rangeInt)
	{
		tempz-=rangeInt;
		++zt;
	}
	zmin=rangeInt*(zt);
	zmax=rangeInt*(zt+1);

	//cout<< "xmin: " << xmin << "xmax: "<< xmax << "\n"; 
	//cout<< "ymin: " << ymin << "ymax: "<< ymax << "\n"; 
	//cout<< "zmin: " << zmin << "zmax: "<< zmax << "\n"; 

	//cout<<"Distance range to be checked is : "<<rangeInt<<"\n";
	double shVal=0;
	if(quadLvl==0){
		for(int i=0;i<vect.size(); i++){
			dist=sqrt(pow((x-vect[i].x),2)+pow((y-vect[i].y),2)+pow((z-vect[i].z),2));
			nm+=(vect[i].value/pow(dist,alpha));
			dm+=(1/pow(dist,alpha));
		}
		if (dm!=0)
			shVal= nm/dm;
	}
	else if (quadLvl > 0)
	{	
		for(int i=0;i<vect.size(); i++){
		// checking the range if the point lies near by
		if(vect[i].x >= xmin && vect[i].x <= xmax && vect[i].y >= ymin && vect[i].y <= ymax && vect[i].z >= zmin && vect[i].z <= zmax){
			dist=sqrt(pow((x-vect[i].x),2)+pow((y-vect[i].y),2)+pow((z-vect[i].z),2));
			nm+=(vect[i].value/pow(dist,alpha));
			dm+=(1/pow(dist,alpha));
			}
		}
		if(dm!=0)
			shVal= nm/dm;
	}
	return shVal;
}