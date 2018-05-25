/*=========================================================================

  Program:   SpaceVoyage
  Module:    SpaceVoyage.cxx

  Author: Samuel G. Deitle
=========================================================================*/

#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyleTrackballActor.h>
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkPNGReader.h"
#include "vtkImageData.h"
#include "vtkRendererCollection.h"
#include "vtkTexture.h"
#include "vtkMatrix4x4.h"
#include "vtkIndent.h"

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>

#include <vector>
#include <math.h>

class Triangle
{
  public:
	  double         X[3];
	  double         Y[3];
	  double         Z[3];
	  double 		 Tu[3];
	  double 		 Tv[3];
};

class Matrix
{
public:
	double A[16];

	Matrix(){;}

	Matrix(double left[4], double up[4], double fore[4], double pos[4])
	{
		for (int i = 0; i < 4; ++i)
		{
			A[0+i]   = left[i];
			A[4+i] = up[i];
			A[8+i] = fore[i];
			A[12+i] = pos[i];
		}
	}

	void SetA(double left[4], double up[4], double fore[4], double pos[4])
	{
		for (int i = 0; i < 4; ++i)
		{
			A[0+i]  = left[i];
			A[4+i]  = up[i];
			A[8+i]  = fore[i];
			A[12+i] = pos[i];
		}
	}

	void SetColumn(int i, double a[4])
	{
		for (int j = 0; j < 4; ++j){A[4*i + j] = a[j];}
	}

	void PrintA()
	{
		for (int i = 0; i < 4; i ++)
			cerr << "(" << A[i] << ",  " << A[4+i] << ",  " << A[8+i] << ",  " << A[12+i] <<")" << endl;
	}
	
	void Pitch(int theta)
	{
		double c  = cos(M_PI*theta/180.0);
		double s  = sin(M_PI*theta/180.0);
		double yx = A[4];
		double yy = A[5];
		double yz = A[6];
		double zx = A[8];
		double zy = A[9];
		double zz = A[10];

		A[4] = c*yx + s*zx;
		A[5] = c*yy + s*zy;
		A[6] = c*yz + s*zz;
		A[8] = c*zx - s*yx;
		A[9] = c*zy - s*yy;
		A[10]= c*zz - s*yz;
	}

	void Yaw(int theta)
	{
		double c  = cos(M_PI*theta/180.0);
		double s  = sin(M_PI*theta/180.0);
		double xx = A[0];
		double xy = A[1];
		double xz = A[2];
		double zx = A[8];
		double zy = A[9];
		double zz = A[10];

		A[0] = c*xx - s*zx;
		A[1] = c*xy - s*zy;
		A[2] = c*xz - s*zz;
		A[8] = c*zx + s*xx;
		A[9] = c*zy + s*xy;
		A[10]= c*zz + s*xz;
	}

	void Roll(int theta)
	{
		double c  = cos(M_PI*theta/180.0);
		double s  = sin(M_PI*theta/180.0);
		double xx = A[0];
		double xy = A[1];
		double xz = A[2];
		double yx = A[4];
		double yy = A[5];
		double yz = A[6];

		A[0] = c*xx + s*yx;
		A[5] = c*xy + s*yy;
		A[6] = c*xz + s*yz;
		A[8] = c*yx - s*xx;
		A[9] = c*yy - s*xy;
		A[10]= c*yz - s*xz;
	}
	
};


void Normalize(double *U)
{
	double length = sqrt( (U[0]*U[0]) + (U[1]*U[1]) + (U[2]*U[2]) );

	for (int i = 0; i < 3; i++)
	{
		U[i] = (U[i]/length);
	}
}

void Normalize(float *U)
{
	double length = sqrt( (U[0]*U[0]) + (U[1]*U[1]) + (U[2]*U[2]) );

	for (int i = 0; i < 3; i++)
	{
		U[i] = (U[i]/length);
	}
}

double DotProduct(double *U, double *V)
{
	double product = 0;
	for (int i = 0; i < 3; i++)
	{product += (U[i]*V[i]);}
		
	return product;
}

double DotProduct(float *U, float *V)
{
	double product = 0;
	for (int i = 0; i < 3; i++)
	{product += (U[i]*V[i]);}
		
	return product;
}


void CrossProduct(double *U, double *V, double *N)
{
// cerr << "In CrossProduct..." << endl;
// cerr << "U: (" << U[0] << ", " << U[1] << ", " << U[2] << ")" << endl;
// cerr << "V: (" << V[0] << ", " << V[1] << ", " << V[2] << ")" << endl;
	N[0] = (U[1]*V[2]) - (U[2]*V[1]);
	N[1] = (U[2]*V[0]) - (U[0]*V[2]);
	N[2] = (U[0]*V[1]) - (U[1]*V[0]);
// cerr << "N: (" << N[0] << ", " << N[1] << ", " << N[2] << ")" << endl;
}

void CrossProduct(float *U, float *V, float *N)
{
// cerr << "In CrossProduct..." << endl;
// cerr << "U: (" << U[0] << ", " << U[1] << ", " << U[2] << ")" << endl;
// cerr << "V: (" << V[0] << ", " << V[1] << ", " << V[2] << ")" << endl;
	N[0] = (U[1]*V[2]) - (U[2]*V[1]);
	N[1] = (U[2]*V[0]) - (U[0]*V[2]);
	N[2] = (U[0]*V[1]) - (U[1]*V[0]);
// cerr << "N: (" << N[0] << ", " << N[1] << ", " << N[2] << ")" << endl;
}

class Rocket
{
public:
	double direction[4]; //direction of forward-axis
	double position[4];  //xyzw-coordinates
	double rotation[3];  //pitch, yaw, roll
	double up[4];		//direction of up-axis
	double uxd[4];		//direction of left-axis (Per right handed coordinates)

	Matrix        axes;

//	bool   RThrustEngaged;
//	bool   LThrustEngaged;
//	bool   TThrustEngaged;
//	bool   BThrustEngaged;

	bool   ThrustersEngaged;

	double length;
	double width;
	double lwratio;

	Rocket(void)
	{
		direction[0] = 0;
		direction[1] = 0;
		direction[2] = -1;
		direction[3] = 0;

		position[0]  = 0;
		position[1]  = 0;
		position[2]  = 0;
		position[3]  = 1;

		rotation[0]  = 0;  //pitch
		rotation[1]  = 0;   //yaw
		rotation[2]  = 0;   //roll

		up[0]		 = 0;
		up[1]		 = 1;
		up[2]		 = 0;
		up[3]		 = 0;

		CrossProduct(up, direction, uxd);
		axes.SetA(uxd, up, direction, position);
		// axes.PrintA();

     // RThrustEngaged = false;
	 //    LThrustEngaged = false;
	 //    TThrustEngaged = false;
	 //    BThrustEngaged = false;

		width		 = 3.1;
		lwratio		 = 3;
		length		 = width*lwratio;


		ThrustersEngaged = false;
	}


	void UpdateDirection()
	{
		direction[0] = axes.A[8];
		direction[1] = axes.A[9];
		direction[2] = axes.A[10];
	}
};

Rocket 						 rocket;
vtkSmartPointer<vtkRenderer> ren1;

// class FollowCam : public vtkCamera
// {
// public:
// 	static FollowCam *New();
// 	void AddSubject(Rocket *rocket)
// 	{
// 		subject = rocket;
// 	}

// 	void RealignView(void)
// 	{
// 		double *pos = subject->position;

// 		double  CamPos[3] = {pos[0],
//   					         pos[1] + 5*rocket.width,
//   					   		 pos[2] + 10*rocket.length};
//   		SetFocalPoint(pos[0],pos[1],pos[2]);
//   		SetPosition(CamPos[0], CamPos[1], CamPos[2]);
//   		OrthogonalizeViewUp();
// 	}

// protected:
// 	Rocket *subject;
// };

// vtkStandardNewMacro(FollowCam);

class vtk441InteractorStyle : public vtkInteractorStyleTrackballCamera
{
	public:
			static vtk441InteractorStyle *New();

			vtk441InteractorStyle()
			{
					shouldPick = false;
					FollowCam  = false;
					cam        = ren1->GetActiveCamera();
			}

			void RealignCamera(double pos[3])
			{
				double  CamPos[3] = {pos[0],
				  				     pos[1] + 5*rocket.width,
				  					 pos[2] + 10*rocket.length};
				cam->SetFocalPoint(pos[0],pos[1],pos[2]-2);
				cam->SetPosition(CamPos[0], CamPos[1], CamPos[2]);
				cam->OrthogonalizeViewUp();
			}

/*			virtual void OnChar()
			{
				vtkRenderWindowInteractor *rwi = this->Interactor;
				char pressedChar = rwi->GetKeyCode();
				cerr << "Pressed " << pressedChar << endl;
				if (pressedChar == 'p')
				{
						cerr << "Should pick!" << endl;
						shouldPick = true;
				}

				vtkInteractorStyleTrackballCamera::OnChar();
			};
*/


			virtual void OnKeyDown()
			{
				rocket.ThrustersEngaged = true;
				vtkRenderWindowInteractor *rwi = this->Interactor;
				vtkRenderWindow *rw            = rwi->GetRenderWindow();

				double *pos = rocket.position;
				double *dir = rocket.direction;
				double *rot = rocket.rotation;

				std::string key = rwi->GetKeySym();
				// cout << "Pressed " << key << endl;

				if (key.compare("w") == 0)
				{
					// double speed = 3.0;
					double p[3];
					cam->GetPosition(p);

					for(int i = 1; i < 4; ++i)
					{

						pos[0] -= dir[0]/i;
						pos[1] -= dir[1]/i;
						pos[2] -= dir[2]/i;
						p[0]   -= dir[0]/i;
						p[1]   -= dir[1]/i;
						p[2]   -= dir[2]/i;
						rocket.axes.SetColumn(3, pos);
						cam->SetPosition(p);
						cam->SetFocalPoint(pos);
						rw->Render();
					}
				}

				else if (key.compare("s") == 0)
				{
					for(int i = 0; i < 3; ++i){pos[i] -= dir[i];}
					//RealignCamera(pos);
				}

				else if (key.compare("f") == 0)
				{
					if(FollowCam){FollowCam = false;}
					else{FollowCam = true;}
				}


				else
				{
				int   theta = 5;
//Left, Right control yaw
// 		=> Left = Rotate CCW, Right = Rotate CW

					if (key.compare("Left") == 0)
					{
						int StartTheta = rot[1];
						int EndTheta   = StartTheta + theta;
						for(int i = StartTheta; i < EndTheta; ++i)
						{rocket.axes.Yaw(2); if(FollowCam){cam->Azimuth(2);}rw->Render();}


					}
					else if (key.compare("Right") == 0)
					{
						int StartTheta = rot[1];
						int EndTheta   = StartTheta + theta;
						for(int i = StartTheta; i < EndTheta; ++i)
						{rocket.axes.Yaw(-2); if(FollowCam){cam->Azimuth(-2);}rw->Render();}
						

					}

//Up, Down control pitch
//		=> Up = Rotate CCW, Down = Rotate CW
					else if (key.compare("Up") == 0)
					{
						int StartTheta = rot[1];
						int EndTheta   = StartTheta + theta;
						for(int i = StartTheta; i < EndTheta; ++i)
						{rocket.axes.Pitch(2); if(FollowCam){cam->Elevation(-2);}rw->Render();}
						

					}
					else if (key.compare("Down") == 0)
					{
						int StartTheta = rot[1];
						int EndTheta   = StartTheta + theta;
						for(int i = StartTheta; i < EndTheta; ++i)
						{rocket.axes.Pitch(-2); if(FollowCam){cam->Elevation(2);}rw->Render();}
						
					}
					
					else if (key.compare("z") == 0)
					{
						int StartTheta = rot[2];
						int EndTheta   = StartTheta + theta;
						for(int i = StartTheta; i < EndTheta; ++i)
						{rocket.axes.Roll(2); rw->Render();}

					}
					else if (key.compare("slash") == 0)
					{
						int StartTheta = rot[2];
						int EndTheta   = StartTheta + theta;
						for(int i = StartTheta; i < EndTheta; ++i)
						{rocket.axes.Roll(-2); rw->Render();}
					}

					if ((int)rot[0] % 360 == 0)
					{rot[0] = 0; dir[1] = 0;}
					else if ((int)rot[1] % 360 == 0)
					{rot[1] = 0; dir[0] = 0;}
					else if ((int)rot[2] % 360 == 0)
					{rot[2] = 0;}
					if ((rot[0] == 0) && (rot[1] == 0))
					{dir[2] = -1;}
					rocket.UpdateDirection();
				}
				rocket.ThrustersEngaged = false;
				rw->Render();

				// rocket.axes.PrintA();
				// Normalize(dir);
				// cerr << "Position:  (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << endl;
				// cerr << "Direction: (" << dir[0] << ", " << dir[1] << ", " << dir[2] << ")" << endl;
				// cerr << "Rotation:  (" << rot[0] << ", " << rot[1] << ", " << rot[2] << ")" << endl;
	  // forward events
				vtkInteractorStyleTrackballCamera::OnKeyPress();
				// Normalize(dir);
			}


			virtual void OnLeftButtonDown()
			{

				vtkRenderWindowInteractor *rwi = this->Interactor;                
				vtkCamera *cam =
				  rwi->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera();
				double *pos = cam->GetDirectionOfProjection();
				cerr << "Camera Direction: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << endl;
				pos = cam->GetPosition();
				cerr << "Camera Position: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << endl;

				if (shouldPick)
				{

					vtkRenderWindow *rw = rwi->GetRenderWindow();
					int *size = rw->GetSize();
					int x = this->Interactor->GetEventPosition()[0];
					int y = this->Interactor->GetEventPosition()[1];
					vtkRenderer *ren = rwi->FindPokedRenderer(x, y);
					double pos[3];
					pos[0] = 2.0*((double)x/(double)size[0])-1;
					pos[1] = 2.0*((double)y/(double)size[1])-1;
					pos[2] = ren->GetZ(x, y);
					cerr << "Picked on " << x << ", " << y << endl;
					cerr << " = " << pos[0] << ", " << pos[1] << ", " << pos[2] << endl;
					ren->ViewToWorld(pos[0], pos[1], pos[2]);
					cerr << " converted to " << pos[0] << ", " << pos[1] << ", " << pos[2] << endl;

					shouldPick = false;
				}

				vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
			};

	private:
			bool shouldPick;
			bool FollowCam;
			vtkCamera *cam;
};

vtkStandardNewMacro(vtk441InteractorStyle);

class vtk441Mapper : public vtkOpenGLPolyDataMapper
{
  protected:
   GLuint displayList;
   bool   initialized;
   double rlength, rwidth;

  public:
   vtk441Mapper()
   {
	 initialized = false;
	 rlength = rocket.length;
	 rwidth  = rocket.width;
   }
	
   void
   RemoveVTKOpenGLStateSideEffects()
   {
	 float Info[4] = { 0, 0, 0, 1 };
	 glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Info);
	 float ambient[4] = { 1,1, 1, 1.0 };
	 glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
	 float diffuse[4] = { 1, 1, 1, 1.0 };
	 glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
	 float specular[4] = { 1, 1, 1, 1.0 };
	 glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
   }


void SetupLight(void)
   {
       glEnable(GL_LIGHTING);
       glEnable(GL_LIGHT0);
       GLfloat diffuse0[4] = { 0.6, 0.6, 0.6, 1 };
       GLfloat ambient0[4] = { 0.2, 0.2, 0.2, 1 };
       GLfloat specular0[4] = { 0.0, 0.0, 0.0, 1 };
       GLfloat pos0[4] = { 0, .707, 0.707, 0 };
       glLightfv(GL_LIGHT0, GL_POSITION, pos0);
       glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse0);
       glLightfv(GL_LIGHT0, GL_AMBIENT, ambient0);
       glLightfv(GL_LIGHT0, GL_SPECULAR, specular0);
       glEnable(GL_LIGHT1);
       GLfloat pos1[4] = { .707, -.707, 0, 0 };
       glLightfv(GL_LIGHT1, GL_POSITION, pos1);
       glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse0);
       glLightfv(GL_LIGHT1, GL_AMBIENT, ambient0);
       glLightfv(GL_LIGHT1, GL_SPECULAR, specular0);
       glDisable(GL_LIGHT2);
       glDisable(GL_LIGHT3);
       glDisable(GL_LIGHT5);
       glDisable(GL_LIGHT6);
       glDisable(GL_LIGHT7);
       // glEnable(GL_LIGHT2);
       // double dpos2[3];
       // ren1->GetActiveCamera()->GetDirectionOfProjection(dpos2);
       // float *pos2 = (float *)dpos2;
       // GLfloat ambient2[4] = { 1, 1, 1, 1 };
       // glLightfv(GL_LIGHT2, GL_POSITION, pos2);
       // glLightfv(GL_LIGHT2, GL_DIFFUSE, diffuse0);
       // glLightfv(GL_LIGHT2, GL_AMBIENT, ambient2);
       // glLightfv(GL_LIGHT2, GL_SPECULAR, specular0);

   }

  
};

class vtk441MapperPart3 : public vtk441Mapper
{
 public:
   static vtk441MapperPart3 *New();
   
   GLuint texture1;
   GLuint flames;
   GLuint texture3;
   GLuint vortex;
   GLuint texture5;
   GLuint texture6;
   GLuint SpaceBox;
   double *pos;
   double *dir;
   double *rot;
   double *up;
   double *uxd;

   vtk441MapperPart3()
   {
	 initialized = false;
	 pos = rocket.position;
     dir = rocket.direction;
     up  = rocket.up;
     uxd = rocket.uxd;
     rot = rocket.rotation;
   }

   void SetUpTexture()
   {

		glGenTextures(1, &texture1);
		glBindTexture(GL_TEXTURE_2D, texture1);

		vtkJPEGReader *rdr = vtkJPEGReader::New();
		rdr->SetFileName("stars1.jpg");
		rdr->Update();
		vtkImageData *img = rdr->GetOutput();
		int dims[3];
		img->GetDimensions(dims);
		unsigned char *buffer = (unsigned char *) img->GetScalarPointer(0,0,0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, dims[0], dims[1], 0, GL_RGB, GL_UNSIGNED_BYTE, buffer);

		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		// glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);


		glGenTextures(1, &flames);
		glBindTexture(GL_TEXTURE_2D, flames);
		rdr->SetFileName("flames.jpg");
		rdr->Update();
		img = rdr->GetOutput();
		img->GetDimensions(dims);
		buffer = (unsigned char *) img->GetScalarPointer(0,0,0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, dims[0], dims[1], 0, GL_RGB, GL_UNSIGNED_BYTE, buffer);

		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

		glGenTextures(1, &texture3);
		glBindTexture(GL_TEXTURE_2D, texture3);
		vtkPNGReader *prdr = vtkPNGReader::New();
		prdr->SetFileName("manyPixel1.png");
		prdr->Update();
		img = prdr->GetOutput();
		img->GetDimensions(dims);
		buffer = (unsigned char *) img->GetScalarPointer(0,0,0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, dims[0], dims[1], 0, GL_RGB, GL_UNSIGNED_BYTE, buffer);

		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		glGenTextures(1, &vortex);
		glBindTexture(GL_TEXTURE_2D, vortex);
		rdr->SetFileName("vortex.jpg");
		rdr->Update();
		img = rdr->GetOutput();
		img->GetDimensions(dims);
		buffer = (unsigned char *) img->GetScalarPointer(0,0,0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, dims[0], dims[1], 0, GL_RGB, GL_UNSIGNED_BYTE, buffer);

		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		glGenTextures(1, &texture5);
		glBindTexture(GL_TEXTURE_2D, texture5);
		rdr->SetFileName("img-thing.jpg");
		rdr->Update();
		img = rdr->GetOutput();
		img->GetDimensions(dims);
		buffer = (unsigned char *) img->GetScalarPointer(0,0,0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, dims[0], dims[1], 0, GL_RGB, GL_UNSIGNED_BYTE, buffer);

		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		glGenTextures(1, &texture6);
		glBindTexture(GL_TEXTURE_2D, texture6);
		rdr->SetFileName("y_o.jpg");
		rdr->Update();
		img = rdr->GetOutput();
		img->GetDimensions(dims);
		buffer = (unsigned char *) img->GetScalarPointer(0,0,0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, dims[0], dims[1], 0, GL_RGB, GL_UNSIGNED_BYTE, buffer);

		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

   }

   void SetupSpaceBox()
   {

		int d_f = 150;
		float i = 3;
		float ambient[3] = {1,1,1};

		SpaceBox = glGenLists(1);
		glNewList(SpaceBox, GL_COMPILE);

		glBindTexture(GL_TEXTURE_2D, texture1);
		glEnable(GL_TEXTURE_2D);

		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
		glBegin(GL_QUADS);

//Back Face
		glNormal3f(0, 0, 1);
		glTexCoord2f(0, 0);
		glVertex3i( d_f,-d_f,-d_f);
		glTexCoord2f(0, i);
		glVertex3i( d_f, d_f,-d_f);
		glTexCoord2f(i, i);
		glVertex3i(-d_f, d_f,-d_f);
		glTexCoord2f(i, 0);
		glVertex3i(-d_f,-d_f,-d_f);

//Right Face
		glNormal3f(-1, 0, 0);
		glTexCoord2f(0, 0);
		glVertex3i( d_f,-d_f, d_f);
		glTexCoord2f(0, i);
		glVertex3i( d_f, d_f, d_f);
		glTexCoord2f(i, i);
		glVertex3i( d_f, d_f,-d_f);
		glTexCoord2f(i, 0);
		glVertex3i( d_f,-d_f,-d_f);

//Front Face (Inward)
		glNormal3f(0, 0, -1);
		glTexCoord2f(0, 0);
		glVertex3i(-d_f,-d_f, d_f);
		glTexCoord2f(0, i);
		glVertex3i(-d_f, d_f, d_f);
		glTexCoord2f(i, i);
		glVertex3i( d_f, d_f, d_f);
		glTexCoord2f(i, 0);
		glVertex3i( d_f,-d_f, d_f);

//Left Face
		glNormal3f(1, 0, 0);
		glTexCoord2f(0, 0);
		glVertex3i(-d_f,-d_f,-d_f);
		glTexCoord2f(0, i);
		glVertex3i(-d_f, d_f,-d_f);
		glTexCoord2f(i, i);
		glVertex3i(-d_f, d_f, d_f);
		glTexCoord2f(i, 0);
		glVertex3i(-d_f,-d_f, d_f);

//Bottom Face
		glNormal3f(0, 1, 0);
		glTexCoord2f(0, 0);
		glVertex3i(-d_f,-d_f,-d_f);
		glTexCoord2f(0, i);
		glVertex3i(-d_f,-d_f, d_f);
		glTexCoord2f(i, i);
		glVertex3i( d_f,-d_f, d_f);
		glTexCoord2f(i, 0);
		glVertex3i( d_f,-d_f,-d_f);        

//Top Face
		glNormal3f(0, -1, 0);
		glTexCoord2f(0, 0);
		glVertex3i( d_f, d_f,-d_f);
		glTexCoord2f(0, i);
		glVertex3i( d_f, d_f, d_f);
		glTexCoord2f(i, i);
		glVertex3i(-d_f, d_f, d_f);
		glTexCoord2f(i, 0);
		glVertex3i(-d_f, d_f,-d_f);

		glEnd();
		glDisable(GL_TEXTURE_2D);

		glEndList();
 
   }

   void DrawCylinder()
   {
	   int nfacets = 30;
	   glBegin(GL_TRIANGLES);
	   for (int i = 0 ; i < nfacets ; i++)
	   {
		   double angle = 3.14159*2.0*i/nfacets;
		   double nextAngle = (i == nfacets-1 ? 0 : 3.14159*2.0*(i+1)/nfacets);
		   glNormal3f(0,0,1);
		   glVertex3f(0, 0, 1);
		   glVertex3f(cos(angle), sin(angle), 1);
		   glVertex3f(cos(nextAngle), sin(nextAngle), 1);
		   glNormal3f(0,0,-1);
		   glVertex3f(0, 0, 0);
		   glVertex3f(cos(angle), sin(angle), 0);
		   glVertex3f(cos(nextAngle), sin(nextAngle), 0);
	   }
	   glEnd();
	   glBegin(GL_QUADS);
	   for (int i = 0 ; i < nfacets ; i++)
	   {
		   double angle = 3.14159*2.0*i/nfacets;
		   double nextAngle = (i == nfacets-1 ? 0 : 3.14159*2.0*(i+1)/nfacets);
		   glNormal3f(cos(angle), sin(angle), 0);
		   glVertex3f(cos(angle), sin(angle), 1);
		   glVertex3f(cos(angle), sin(angle), 0);
		   glNormal3f(cos(nextAngle), sin(nextAngle), 0);
		   glVertex3f(cos(nextAngle), sin(nextAngle), 0);
		   glVertex3f(cos(nextAngle), sin(nextAngle), 1);
	   }
	   glEnd();
   }

   void DrawTexCylinder(GLuint *texture)
   {
   	   glBindTexture(GL_TEXTURE_2D, *texture);
   	   glEnable(GL_TEXTURE_2D);
	   int nfacets = 30;
	   glBegin(GL_TRIANGLES);
	   for (int i = 0 ; i < nfacets ; i++)
	   {
		   double angle = 3.14159*2.0*i/nfacets;
		   double nextAngle = (i == nfacets-1 ? 0 : 3.14159*2.0*(i+1)/nfacets);

		   const double ca  = cos(angle);
		   const double sa  = sin(angle);
		   const double cnA = cos(nextAngle);
		   const double snA = sin(nextAngle);

		   const double tX  = (1 + ca)/2.0;
		   const double tY  = (1 + sa)/2.0;
		   const double tnX = (1 + cnA)/2.0;
		   const double tnY = (1 + snA)/2.0;

		   glNormal3f(0,0,1);
		   glTexCoord2f(.5, .5);
		   glVertex3f(0, 0, 1);
		   glTexCoord2f(tX, tY);
		   glVertex3f(ca, sa, 1);
		   glTexCoord2f(tnX, tnY);
		   glVertex3f(cnA, snA, 1);

		   glNormal3f(0,0,-1);
		   glTexCoord2f(.5, .5);
		   glVertex3f(0, 0, 0);
		   glTexCoord2f(tX, tY);
		   glVertex3f(ca, sa, 0);
		   glTexCoord2f(tnX, tnY);
		   glVertex3f(cnA, snA, 0);
	   }
	   glEnd();
	   glBegin(GL_QUADS);
	   for (int i = 0 ; i < nfacets ; i++)
	   {
		   double angle = 3.14159*2.0*i/nfacets;
		   double nextAngle = (i == nfacets-1 ? 0 : 3.14159*2.0*(i+1)/nfacets);

		   const double ca  = cos(angle);
		   const double sa  = sin(angle);
		   const double cnA = cos(nextAngle);
		   const double snA = sin(nextAngle);
		   const double tex = (double)i/(double)nfacets;
		   const double ntex= (double)(i+1)/(double)nfacets;

		   glNormal3f(ca, sa, 0);
		   glTexCoord2f(tex, 1);
		   glVertex3f(ca, sa, 1);
		   glTexCoord2f(tex, 0);
		   glVertex3f(ca, sa, 0);
		   glNormal3f(cnA, snA, 0);
		   glTexCoord2f(ntex, 0);
		   glVertex3f(cnA, snA, 0);
		   glTexCoord2f(ntex, 1);
		   glVertex3f(cnA, snA, 1);
	   }
	   glEnd();
	   glDisable(GL_TEXTURE_2D);
   }

   std::vector<Triangle> SplitTriangle(std::vector<Triangle> &list)
   {
	   std::vector<Triangle> output(4*list.size());
	   for (unsigned int i = 0 ; i < list.size() ; i++)
	   {
		   double mid1[5], mid2[5], mid3[5];
		   mid1[0] = (list[i].X[0]+list[i].X[1])/2;
		   mid1[1] = (list[i].Y[0]+list[i].Y[1])/2;
		   mid1[2] = (list[i].Z[0]+list[i].Z[1])/2;
		   mid1[3] = (list[i].Tu[0]+list[i].Tu[1])/2;
		   mid1[4] = (list[i].Tv[0]+list[i].Tv[1])/2;

		   mid2[0] = (list[i].X[1]+list[i].X[2])/2;
		   mid2[1] = (list[i].Y[1]+list[i].Y[2])/2;
		   mid2[2] = (list[i].Z[1]+list[i].Z[2])/2;
		   mid2[3] = (list[i].Tu[1]+list[i].Tu[2])/2;
		   mid2[4] = (list[i].Tv[1]+list[i].Tv[2])/2;
		   
		   mid3[0] = (list[i].X[0]+list[i].X[2])/2;
		   mid3[1] = (list[i].Y[0]+list[i].Y[2])/2;
		   mid3[2] = (list[i].Z[0]+list[i].Z[2])/2;
		   mid3[3] = (list[i].Tu[0]+list[i].Tu[2])/2;
		   mid3[4] = (list[i].Tv[0]+list[i].Tv[2])/2;
		   
		   output[4*i+0].X[0]  = list[i].X[0];
		   output[4*i+0].Y[0]  = list[i].Y[0];
		   output[4*i+0].Z[0]  = list[i].Z[0];
		   output[4*i+0].Tu[0] = list[i].Tu[0];
		   output[4*i+0].Tv[0] = list[i].Tv[0];
		   
		   output[4*i+0].X[1]  = mid1[0];
		   output[4*i+0].Y[1]  = mid1[1];
		   output[4*i+0].Z[1]  = mid1[2];
		   output[4*i+0].Tu[1] = mid1[3];
		   output[4*i+0].Tv[1] = mid1[4];
		   
		   output[4*i+0].X[2]  = mid3[0];
		   output[4*i+0].Y[2]  = mid3[1];
		   output[4*i+0].Z[2]  = mid3[2];
		   output[4*i+0].Tu[2] = mid3[3];
		   output[4*i+0].Tv[2] = mid3[4];
		   
		   output[4*i+1].X[0]  = list[i].X[1];
		   output[4*i+1].Y[0]  = list[i].Y[1];
		   output[4*i+1].Z[0]  = list[i].Z[1];
		   output[4*i+1].Tu[0] = list[i].Tu[1];
		   output[4*i+1].Tv[0] = list[i].Tv[1];
		   
		   output[4*i+1].X[1]  = mid2[0];
		   output[4*i+1].Y[1]  = mid2[1];
		   output[4*i+1].Z[1]  = mid2[2];
		   output[4*i+1].Tu[1] = mid2[3];
		   output[4*i+1].Tv[1] = mid2[4];
		   
		   output[4*i+1].X[2]  = mid1[0];
		   output[4*i+1].Y[2]  = mid1[1];
		   output[4*i+1].Z[2]  = mid1[2];
		   output[4*i+1].Tu[2] = mid1[3];
		   output[4*i+1].Tv[2] = mid1[4];
		   
		   output[4*i+2].X[0]  = list[i].X[2];
		   output[4*i+2].Y[0]  = list[i].Y[2];
		   output[4*i+2].Z[0]  = list[i].Z[2];
		   output[4*i+2].Tu[0] = list[i].Tu[2];
		   output[4*i+2].Tv[0] = list[i].Tv[2];
		   
		   output[4*i+2].X[1]  = mid3[0];
		   output[4*i+2].Y[1]  = mid3[1];
		   output[4*i+2].Z[1]  = mid3[2];
		   output[4*i+2].Tu[1] = mid3[3];
		   output[4*i+2].Tv[1] = mid3[4];		   
		   
		   output[4*i+2].X[2]  = mid2[0];
		   output[4*i+2].Y[2]  = mid2[1];
		   output[4*i+2].Z[2]  = mid2[2];
		   output[4*i+2].Tu[2] = mid2[3];
		   output[4*i+2].Tv[2] = mid2[4];
		   
		   output[4*i+3].X[0]  = mid1[0];
		   output[4*i+3].Y[0]  = mid1[1];
		   output[4*i+3].Z[0]  = mid1[2];
		   output[4*i+3].Tu[0] = mid1[3];
		   output[4*i+3].Tv[0] = mid1[4];
		   
		   output[4*i+3].X[1]  = mid2[0];
		   output[4*i+3].Y[1]  = mid2[1];
		   output[4*i+3].Z[1]  = mid2[2];
		   output[4*i+3].Tu[1] = mid2[3];
		   output[4*i+3].Tv[1] = mid2[4];
		   
		   output[4*i+3].X[2]  = mid3[0];
		   output[4*i+3].Y[2]  = mid3[1];
		   output[4*i+3].Z[2]  = mid3[2];
		   output[4*i+3].Tu[2] = mid3[3];
		   output[4*i+3].Tv[2] = mid3[4];
	   }
	   return output;
   }

   void DrawSphere()
   {
	   int recursionLevel = 3;
	   Triangle t;
	   t.X[0]  = 1;
	   t.Y[0]  = 0;
	   t.Z[0]  = 0;
	   t.Tu[0] = 1;
	   t.Tv[0] = 0;
	   
	   t.X[1]  = 0;
	   t.Y[1]  = 1;
	   t.Z[1]  = 0;
	   t.Tu[1] = 1;
	   t.Tv[1] = 0;

	   t.X[2]  = 0;
	   t.Y[2]  = 0;
	   t.Z[2]  = 1;
	   t.Tu[2] = 1;
	   t.Tv[2] = 1;
	   std::vector<Triangle> list;
	   list.push_back(t);
	   for (int r = 0 ; r < recursionLevel ; r++)
	   {
		   list = SplitTriangle(list);
	   }

	   // really draw `
	   for (int octent = 0 ; octent < 8 ; octent++)
	   {
		   glPushMatrix();
		   glRotatef(90*(octent%4), 1, 0, 0);
		   if (octent >= 4)
			   glRotatef(180, 0, 0, 1);
		   glBegin(GL_TRIANGLES);
		   for (unsigned int i = 0 ; i < list.size() ; i++)
		   {
			   for (int j = 0 ; j < 3 ; j++)
			   {
				   double ptMag = sqrt(list[i].X[j]*list[i].X[j]+
									   list[i].Y[j]*list[i].Y[j]+
									   list[i].Z[j]*list[i].Z[j]);

				   glNormal3f(list[i].X[j]/ptMag, list[i].Y[j]/ptMag, list[i].Z[j]/ptMag);
				   glVertex3f(list[i].X[j]/ptMag, list[i].Y[j]/ptMag, list[i].Z[j]/ptMag);
			   }
		   }
		   glEnd();
		   glPopMatrix();
	   }
	}
	void DrawTexSphere(GLuint *texture)
   	{
   	   glBindTexture(GL_TEXTURE_2D, *texture);
       glEnable(GL_TEXTURE_2D);
	   float ambient[3] = {.5,.5,.5};
	   glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, ambient);
	   
	   int recursionLevel = 4;
	   Triangle t;
	   t.X[0]  = 1;
	   t.Y[0]  = 0;
	   t.Z[0]  = 0;
	   t.Tu[0] = 0;
	   t.Tv[0] = 0;

	   t.X[1]  = 0;
	   t.Y[1]  = 1;
	   t.Z[1]  = 0;
	   t.Tu[1] = .5;
	   t.Tv[1] = 1;

	   t.X[2]  = 0;
	   t.Y[2]  = 0;
	   t.Z[2]  = 1;
	   t.Tu[2] = 1;
	   t.Tv[2] = 0;
	   std::vector<Triangle> list;
	   list.push_back(t);
	   for (int r = 0 ; r < recursionLevel ; r++)
	   {
		   list = SplitTriangle(list);
	   }

	   // really draw `
	   for (int octent = 0 ; octent < 8 ; octent++)
	   {
		   glPushMatrix();
		   glRotatef(90*(octent%4), 0, 1, 0);
		   if (octent >= 4)
			   glRotatef(180, 0, 0, 1);
		   glBegin(GL_TRIANGLES);
		   for (unsigned int i = 0 ; i < list.size() ; i++)
		   {
			   for (int j = 0 ; j < 3 ; j++)
			   {
				   double ptMag = sqrt(list[i].X[j]*list[i].X[j]+
									   list[i].Y[j]*list[i].Y[j]+
									   list[i].Z[j]*list[i].Z[j]);

				   glTexCoord2f(list[i].Tu[j]/ptMag, list[i].Tv[j]/ptMag);
				   glNormal3f(list[i].X[j]/ptMag, list[i].Y[j]/ptMag, list[i].Z[j]/ptMag);
				   glVertex3f(list[i].X[j]/ptMag, list[i].Y[j]/ptMag, list[i].Z[j]/ptMag);
			   }
		   }
		   glEnd();
		   glPopMatrix();
	   }
	   glDisable(GL_TEXTURE_2D);
   	}



   void Brown(void)      { glColor3ub(205, 133, 63); };
   void LightBrown(void) { glColor3ub(245, 222, 179); };
   void DarkBrown(void)  { glColor3ub(162, 82, 45); };
   void Pink(void)       { glColor3ub(250, 128, 114); };
   void White(void)      { glColor3ub(255, 255, 255); };
   void Black(void)      { glColor3ub(0, 0, 0); };
   void BlueGrey(void)   { glColor3ub(128, 145, 183); };
   void Red(void)        { glColor3ub(255, 0, 0); };

   void DrawRocket()
   {
   		glEnable(GL_COLOR_MATERIAL);
		glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

		//Place rocket at current position and orientation;
		glMatrixMode(GL_MODELVIEW);
		// glLoadIdentity();
		glMultMatrixd(rocket.axes.A);
		// cerr << endl << endl;
		// glTranslatef(pos[0], pos[1], pos[2]);
		// glRotatef(rot[0], uxd[0], uxd[1], uxd[2]);
		// glRotatef(rot[1], up[0], up[1], up[2]);
		// glRotatef(rot[2], dir[0], dir[1], dir[2]);

		// glGetDoublev(GL_MODELVIEW_MATRIX, rocket.axes.A);
		// rocket.axes.PrintA();




		//Draw Body
		float ambient[3]  = {1, 1, 1};
		float specular[3] = {.2, .2, .2};
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
		glPushMatrix();
		glScalef(rwidth, rwidth, rlength);
		White();
		DrawSphere();
		glPopMatrix();

		//Draw Thrusters
		Red();
		glPushMatrix();

		glTranslatef(0, 0, .5*rlength);
		glScalef(1, 1.3, 1.2*rlength);
		float diffuse[4] = {1, 1, 1, 1};
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
		for (int i = 1; i < 7; ++i)
		{
			glPushMatrix();
			glRotatef(i*60 - 45, 0, 0, 1);
			glTranslatef(.75*rwidth, .75*rwidth, 0);
			DrawCylinder();
			if(rocket.ThrustersEngaged)
			{
				glPushMatrix();
				glTranslatef(0, 0, 1);
				glScalef(1.1, 1, .5);
				DrawTexCylinder(&flames);
				glPopMatrix();
			}
			glScalef(1, 1, 1./(1.2*rlength));
			DrawSphere();
			glPopMatrix();
		}
		Black();
		glPushMatrix();
		glTranslatef(0, 0, -.01);
		glScalef(1.1*rwidth, rwidth, 3./rlength);
		DrawCylinder();
		glPopMatrix();

		glPopMatrix();
   }
  

   virtual void RenderPiece(vtkRenderer *ren, vtkActor *act)
   {
   		if(!initialized)
   		{
   			SetUpTexture();
   			// SetupSpaceBox();
   			initialized = true;
   		}

		RemoveVTKOpenGLStateSideEffects();
		SetupLight();
		glEnable(GL_COLOR_MATERIAL);
	
		//glBindTexture(GL_TEXTURE_2D, texture1);
		glMatrixMode(GL_MODELVIEW);

		glCallList(SpaceBox);
	
		glPushMatrix();
		DrawRocket();
		glPopMatrix();

		glPushMatrix();
		glTranslatef(40, 40, 40);
		glScalef(15, 15, 15);
		DrawTexSphere(&texture5);
		glPopMatrix();


		glPushMatrix();
		glTranslatef(40, 0, 100);
		glScalef(15, 15, 15);
		DrawTexSphere(&texture6);
		glPopMatrix();

		glPushMatrix();
		glTranslatef(0, 30, 149);
		glScalef(20, 20, 1);
		DrawTexCylinder(&vortex);

		glPopMatrix();


		// glEndList();
	}
};

vtkStandardNewMacro(vtk441MapperPart3);


int main()
{
  // Dummy input so VTK pipeline mojo is happy.
  //
  vtkSmartPointer<vtkSphereSource> sphere =
	vtkSmartPointer<vtkSphereSource>::New();
  sphere->SetThetaResolution(100);
  sphere->SetPhiResolution(50);

  vtkSmartPointer<vtk441MapperPart3> win3Mapper =
	vtkSmartPointer<vtk441MapperPart3>::New();
  win3Mapper->SetInputConnection(sphere->GetOutputPort());

  vtkSmartPointer<vtkActor> win3Actor =
	vtkSmartPointer<vtkActor>::New();
  win3Actor->SetMapper(win3Mapper);

  ren1 = vtkSmartPointer<vtkRenderer>::New();

  vtkSmartPointer<vtkRenderWindow> renWin =
	vtkSmartPointer<vtkRenderWindow>::New();
  renWin->AddRenderer(ren1);
  ren1->SetViewport(0, 0, 1, 1);

  vtkSmartPointer<vtkRenderWindowInteractor> iren =
	vtkSmartPointer<vtkRenderWindowInteractor>::New();
  iren->SetRenderWindow(renWin);

  vtk441InteractorStyle *style = vtk441InteractorStyle::New();
  iren->SetInteractorStyle(style);
  style->SetInteractor(iren);

  // Add the actors to the renderer, set the background and size.
  //
  bool doWindow3 = true;
  if (doWindow3)
	  ren1->AddActor(win3Actor);

  vtkJPEGReader *rdr =
	vtkJPEGReader::New();
  rdr->SetFileName("stars1.jpg");
  rdr->Update();

  vtkTexture *texture = vtkTexture::New();
  texture->SetInputConnection(rdr->GetOutputPort());

  ren1->TexturedBackgroundOn();
  ren1->SetBackgroundTexture(texture);
  // ren1->SetBackground(.3, .3, .3);
  renWin->SetSize(1000, 1000);
  renWin->SetWindowName("Space Voyage");

  // Set up the lighting.
  //
  vtkSmartPointer<vtkLight> light =
	vtkSmartPointer<vtkLight>::New();
  light->SetFocalPoint(1.875,0.6125,0);
  light->SetPosition(0.875,1.6125,1);
  light->SetPosition(0,20,20);
  ren1->AddLight(light);

//***Ideally this would set up the FollowCam***\\

  vtkSmartPointer<FollowCam> cam =
	vtkSmartPointer<FollowCam>::New();
  cam->AddSubject(&rocket);
  ren1->SetActiveCamera(cam);

  double *pos = rocket.position;

  
  double  CamPos[3] = {pos[0],
  					   pos[1] + 6*rocket.width,
  					   pos[2] - 10*rocket.length};
  ren1->GetActiveCamera()->SetFocalPoint(pos[0],pos[1],pos[2]);
  ren1->GetActiveCamera()->SetPosition(CamPos[0], CamPos[1], CamPos[2]);
  ren1->GetActiveCamera()->OrthogonalizeViewUp();
  ren1->GetActiveCamera()->SetClippingRange(10,400);
  // ren1->GetActiveCamera()->SetDistance(20);
  
  // This starts the event loop and invokes an initial render.
  //
  ((vtkInteractorStyle *)iren->GetInteractorStyle())->SetAutoAdjustCameraClippingRange(0);
  renWin->Render();
  iren->Initialize();
  iren->Start();

  return EXIT_SUCCESS;
}