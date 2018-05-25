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
#include "vtkImageData.h"
#include "vtkRendererCollection.h"
#include "vtkTexture.h"

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
};

void Normalize(double *U)
{
	double length = sqrt( (U[0]*U[0]) + (U[1]*U[1]) + (U[2]*U[2]) );

	for (int i = 0; i < 3; i++)
	{
		U[i] = (U[i]/length);
	}
}

class Rocket
{
public:
	double direction[3]; //direction of roll-axis vector
	double position[3];  //xyz-coordinates
	double rotation[3];  //pitch, yaw, roll

	double length;
	double width;
	double lwratio;

	Rocket(void)
	{
		direction[0] = 0;
		direction[1] = 0;
		direction[2] = -1;

		position[0]  = 0;
		position[1]  = 0;
		position[2]  = 0;

		rotation[0]  = 0;  //pitch
		rotation[1]  = 0;   //yaw
		rotation[2]  = 0;   //roll

		width		 = 3.1;
		lwratio		 = 3.1;
		length		 = width*lwratio;
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
			}

			void RealignCamera(double pos[3])
			{
				vtkCamera *cam 	  = ren1->GetActiveCamera();

				double  CamPos[3] = {pos[0],
				  				     pos[1] + 5*rocket.width,
				  					 pos[2] + 10*rocket.length};
				cam->SetFocalPoint(pos[0],pos[1],pos[2] - 2);
				cam->SetPosition(CamPos[0], CamPos[1], CamPos[2]);
				cam->OrthogonalizeViewUp();
				cam->SetClippingRange(10,250);
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
				vtkRenderWindowInteractor *rwi = this->Interactor;
				vtkRenderWindow *rw            = rwi->GetRenderWindow();
				vtkCamera *cam 				   = ren1->GetActiveCamera();

				double *pos = rocket.position;
				double *dir = rocket.direction;
				double *rot = rocket.rotation;

				std::string key = rwi->GetKeySym();
				cout << "Pressed " << key << endl;

				if (key.compare("w") == 0)
				{
					double p[3];
					cam->GetPosition(p);
					for(int i = 0; i < 3; ++i){pos[i] += dir[i]; p[i] += dir[i];}
					cam->SetPosition(p);
					cam->SetFocalPoint(pos);
				}

				else if (key.compare("s") == 0)
				{
					double p[3];
					cam->GetPosition(p);
					for(int i = 0; i < 3; ++i){pos[i] -= dir[i]; p[i] -= dir[i];}
					cam->SetPosition(p);
					cam->SetFocalPoint(pos);
				}

				else
				{
				int   theta = 5;
				double t  = (theta/180.0)*3.14159;
				double ct = cos(t);
				double st = sin(t);

//Left, Right control yaw
// 		=> Left = Rotate CCW, Right = Rotate CW

				if (key.compare("Left") == 0)
				{
					rot[1] += theta;

					if ((int)rot[1] % 360 == 0)
					{
						rot[1]  = 0;
						dir[0]  = 0;
						dir[2] *= ct;
					}
					else
					{
						double dir_x = ct*dir[0] + st*dir[2];
						double dir_z = ct*dir[2] - st*dir[0];

						dir[0]  = dir_x;
						dir[2]  = dir_z;
					}
					cam->Azimuth(theta);
				}
				else if (key.compare("Right") == 0)
				{
					rot[1] -= theta;
					if ((int)rot[1] % 360 == 0)
					{
						rot[1]  = 0;
						dir[0]  = 0;
						dir[2] *= ct;
					}
					else
					{
						double dir_x = ct*dir[0] - st*dir[2];
						double dir_z = ct*dir[2] + st*dir[0];

						dir[0]  = dir_x;
						dir[2]  = dir_z;
					}
					cam->Azimuth(-theta);
				}

//Up, Down control pitch
//		=> Up = Rotate CCW, Down = Rotate CW
				else if (key.compare("Up") == 0)
				{
					rot[0] -= theta;
					if ((int)rot[0] % 360 == 0)
					{
						rot[0]  = 0;
						dir[1]  = 0;
						dir[2] *= ct;
					}
					else
					{
						double dir_y = ct*dir[1] + st*dir[2];
						double dir_z = ct*dir[2] - st*dir[1];

						dir[1]  = dir_y;
						dir[2]  = dir_z;
					}
					cam->Elevation(theta);					
				}
				else if (key.compare("Down") == 0)
				{
					rot[0] += theta;

					if ((int)rot[0] % 360 == 0)
					{
						rot[0]  = 0;
						dir[1]  = 0;
						dir[2] *= ct;
					}
					else
					{
						double dir_y = ct*dir[1] - st*dir[2];
						double dir_z = ct*dir[2] + st*dir[1];

						dir[1]  = dir_y;
						dir[2]  = dir_z;
					}
					cam->Elevation(-theta);			
				}
				if ((rot[0] == 0) && (rot[1] == 0))
				{
					dir[2] = -1;
					//RealignCamera(pos);
				}
				}				

				cam->OrthogonalizeViewUp();
				rw->Render();
				// Normalize(dir);
				cerr << "Position:  (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << endl;
				cerr << "Direction: (" << dir[0] << ", " << dir[1] << ", " << dir[2] << ")" << endl;
				cerr << "Rotation:  (" << rot[0] << ", " << rot[1] << ", " << rot[2] << ")" << endl;

	  // forward events
				vtkInteractorStyleTrackballCamera::OnKeyPress();
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
   GLuint displayList;
   double *pos;
   double *dir;
   double *rot;

//   bool   initialized;
   

   vtk441MapperPart3()
   {
	 initialized = false;
	 pos = rocket.position;
     dir = rocket.direction;
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
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

		initialized = true;
   }

   void DrawSpacebox()
   {

		int d_f = 90;

		

		glEnable(GL_TEXTURE_2D);
		float ambient[3] = {1,1,1};
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
		glBegin(GL_QUADS);

//Back Face
		glNormal3f(0, 0, 1);
		glTexCoord2f(0, 0);
		glVertex3i( d_f,-d_f,-d_f);
		glTexCoord2f(0, 1);
		glVertex3i( d_f, d_f,-d_f);
		glTexCoord2f(1, 1);
		glVertex3i(-d_f, d_f,-d_f);
		glTexCoord2f(1, 0);
		glVertex3i(-d_f,-d_f,-d_f);

//Right Face
		glNormal3f(-1, 0, 0);
		glTexCoord2f(0, 0);
		glVertex3i( d_f,-d_f, d_f);
		glTexCoord2f(0, 1);
		glVertex3i( d_f, d_f, d_f);
		glTexCoord2f(1, 1);
		glVertex3i( d_f, d_f,-d_f);
		glTexCoord2f(1, 0);
		glVertex3i( d_f,-d_f,-d_f);

//Front Face (Inward)
		glNormal3f(0, 0, -1);
		glTexCoord2f(0, 0);
		glVertex3i(-d_f,-d_f, d_f);
		glTexCoord2f(0, 1);
		glVertex3i(-d_f, d_f, d_f);
		glTexCoord2f(1, 1);
		glVertex3i( d_f, d_f, d_f);
		glTexCoord2f(1, 0);
		glVertex3i( d_f,-d_f, d_f);

//Left Face
		glNormal3f(1, 0, 0);
		glTexCoord2f(0, 0);
		glVertex3i(-d_f,-d_f,-d_f);
		glTexCoord2f(0, 1);
		glVertex3i(-d_f, d_f,-d_f);
		glTexCoord2f(1, 1);
		glVertex3i(-d_f, d_f, d_f);
		glTexCoord2f(1, 0);
		glVertex3i(-d_f,-d_f, d_f);

//Bottom Face
		glNormal3f(0, 1, 0);
		glTexCoord2f(0, 0);
		glVertex3i(-d_f,-d_f,-d_f);
		glTexCoord2f(0, 1);
		glVertex3i(-d_f,-d_f, d_f);
		glTexCoord2f(1, 1);
		glVertex3i( d_f,-d_f, d_f);
		glTexCoord2f(1, 0);
		glVertex3i( d_f,-d_f,-d_f);        

//Top Face
		glNormal3f(0, -1, 0);
		glTexCoord2f(0, 0);
		glVertex3i( d_f, d_f,-d_f);
		glTexCoord2f(0, 1);
		glVertex3i( d_f, d_f, d_f);
		glTexCoord2f(1, 1);
		glVertex3i(-d_f, d_f, d_f);
		glTexCoord2f(1, 0);
		glVertex3i(-d_f, d_f,-d_f);

		glEnd();
		glDisable(GL_TEXTURE_2D);
 
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
   std::vector<Triangle> SplitTriangle(std::vector<Triangle> &list)
   {
	   std::vector<Triangle> output(4*list.size());
	   for (unsigned int i = 0 ; i < list.size() ; i++)
	   {
		   double mid1[3], mid2[3], mid3[3];
		   mid1[0] = (list[i].X[0]+list[i].X[1])/2;
		   mid1[1] = (list[i].Y[0]+list[i].Y[1])/2;
		   mid1[2] = (list[i].Z[0]+list[i].Z[1])/2;
		   mid2[0] = (list[i].X[1]+list[i].X[2])/2;
		   mid2[1] = (list[i].Y[1]+list[i].Y[2])/2;
		   mid2[2] = (list[i].Z[1]+list[i].Z[2])/2;
		   mid3[0] = (list[i].X[0]+list[i].X[2])/2;
		   mid3[1] = (list[i].Y[0]+list[i].Y[2])/2;
		   mid3[2] = (list[i].Z[0]+list[i].Z[2])/2;
		   output[4*i+0].X[0] = list[i].X[0];
		   output[4*i+0].Y[0] = list[i].Y[0];
		   output[4*i+0].Z[0] = list[i].Z[0];
		   output[4*i+0].X[1] = mid1[0];
		   output[4*i+0].Y[1] = mid1[1];
		   output[4*i+0].Z[1] = mid1[2];
		   output[4*i+0].X[2] = mid3[0];
		   output[4*i+0].Y[2] = mid3[1];
		   output[4*i+0].Z[2] = mid3[2];
		   output[4*i+1].X[0] = list[i].X[1];
		   output[4*i+1].Y[0] = list[i].Y[1];
		   output[4*i+1].Z[0] = list[i].Z[1];
		   output[4*i+1].X[1] = mid2[0];
		   output[4*i+1].Y[1] = mid2[1];
		   output[4*i+1].Z[1] = mid2[2];
		   output[4*i+1].X[2] = mid1[0];
		   output[4*i+1].Y[2] = mid1[1];
		   output[4*i+1].Z[2] = mid1[2];
		   output[4*i+2].X[0] = list[i].X[2];
		   output[4*i+2].Y[0] = list[i].Y[2];
		   output[4*i+2].Z[0] = list[i].Z[2];
		   output[4*i+2].X[1] = mid3[0];
		   output[4*i+2].Y[1] = mid3[1];
		   output[4*i+2].Z[1] = mid3[2];
		   output[4*i+2].X[2] = mid2[0];
		   output[4*i+2].Y[2] = mid2[1];
		   output[4*i+2].Z[2] = mid2[2];
		   output[4*i+3].X[0] = mid1[0];
		   output[4*i+3].Y[0] = mid1[1];
		   output[4*i+3].Z[0] = mid1[2];
		   output[4*i+3].X[1] = mid2[0];
		   output[4*i+3].Y[1] = mid2[1];
		   output[4*i+3].Z[1] = mid2[2];
		   output[4*i+3].X[2] = mid3[0];
		   output[4*i+3].Y[2] = mid3[1];
		   output[4*i+3].Z[2] = mid3[2];
	   }
	   return output;
   }

   void DrawSphere()
   {
	   int recursionLevel = 3;
	   Triangle t;
	   t.X[0] = 1;
	   t.Y[0] = 0;
	   t.Z[0] = 0;
	   t.X[1] = 0;
	   t.Y[1] = 1;
	   t.Z[1] = 0;
	   t.X[2] = 0;
	   t.Y[2] = 0;
	   t.Z[2] = 1;
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



   void Brown(void)      { glColor3ub(205, 133, 63); };
   void LightBrown(void) { glColor3ub(245, 222, 179); };
   void DarkBrown(void)  { glColor3ub(162, 82, 45); };
   void Pink(void)       { glColor3ub(250, 128, 114); };
   void White(void)      { glColor3ub(255, 255, 255); };
   void Black(void)      { glColor3ub(0, 0, 0); };
   void BlueGrey(void)   { glColor3ub(128, 145, 183); };
   void Red(void)        { glColor3ub(255, 0, 0); };
  

   virtual void RenderPiece(vtkRenderer *ren, vtkActor *act)
   {
   		if(!initialized)
   		{
   			SetUpTexture();
   		}

		RemoveVTKOpenGLStateSideEffects();
		SetupLight();
		// displayList = glGenLists(1);
		// glNewList(displayList, GL_COMPILE_AND_EXECUTE);
		glEnable(GL_COLOR_MATERIAL);
	
		glBindTexture(GL_TEXTURE_2D, texture1);
		glMatrixMode(GL_MODELVIEW);

//		DrawSpacebox();
		glPushMatrix();
		glEnable(GL_COLOR_MATERIAL);
		glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

		//Place rocket at current position and direction
		glTranslatef(pos[0], pos[1], pos[2]);
		glRotatef(rot[0], 1, 0, 0);
		glRotatef(rot[1], 0, 1, 0);
		glRotatef(rot[2], 0, 0, 1);


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
		glScalef(1, 1.3, 1.3*rlength);
		float diffuse[4] = {1, 1, 1, 1};
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
		for (int i = 1; i < 7; ++i)
		{
			glPushMatrix();
			glRotatef(i*60 - 5, 0, 0, 1);
			glTranslatef(.75*rwidth, .75*rwidth, 0);
			DrawCylinder();
			glScalef(1, 1, 1./(1.3*rlength));
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
//
 //  vtkSmartPointer<FollowCam> cam =
	// vtkSmartPointer<FollowCam>::New();
 //  cam->AddSubject(&rocket);
 //  ren1->SetActiveCamera(cam);

  double *pos = rocket.position;

  
  double  CamPos[3] = {pos[0],
  					   pos[1] + 6*rocket.width,
  					   pos[2] + 10*rocket.length};
  ren1->GetActiveCamera()->SetFocalPoint(pos[0],pos[1],pos[2] - 2);
  ren1->GetActiveCamera()->SetPosition(CamPos[0], CamPos[1], CamPos[2]);
  ren1->GetActiveCamera()->OrthogonalizeViewUp();
  ren1->GetActiveCamera()->SetClippingRange(10,250);
//  ren1->GetActiveCamera()->SetDistance(20);
  
  // This starts the event loop and invokes an initial render.
  //
  ((vtkInteractorStyle *)iren->GetInteractorStyle())->SetAutoAdjustCameraClippingRange(0);
//  renWin->Render();
  iren->Initialize();
  iren->Start();


  return EXIT_SUCCESS;
}




