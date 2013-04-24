/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "showCube.h"
#include <vector>

int pointMap(int side, int i, int j)
{
  int r;

  switch (side)
  {
  case 1: //[i][j][0] bottom face
    r = 64 * i + 8 * j;
    break;
  case 6: //[i][j][7] top face
    r = 64 * i + 8 * j + 7;
    break;
  case 2: //[i][0][j] front face
    r = 64 * i + j;
    break;
  case 5: //[i][7][j] back face
    r = 64 * i + 56 + j;
    break;
  case 3: //[0][i][j] left face
    r = 8 * i + j;
    break;
  case 4: //[7][i][j] right face
    r = 448 + 8 * i + j;
    break;
  }

  return r;
}

void showCube(struct world * jello)
{
  int i,j,k,ip,jp,kp;
  point r1,r2,r3; // aux variables
  
  /* normals buffer and counter for Gourad shading*/
  struct point normal[8][8];
  int counter[8][8];

  int face;
  double faceFactor, length;

  if (fabs(jello->p[0][0][0].x) > 10)
  {
    printf ("Your cube somehow escaped way out of the box.\n");
    exit(0);
  }

  
  #define NODE(face,i,j) (*((struct point * )(jello->p) + pointMap((face),(i),(j))))

  
  #define PROCESS_NEIGHBOUR(di,dj,dk) \
    ip=i+(di);\
    jp=j+(dj);\
    kp=k+(dk);\
    if\
    (!( (ip>7) || (ip<0) ||\
      (jp>7) || (jp<0) ||\
    (kp>7) || (kp<0) ) && ((i==0) || (i==7) || (j==0) || (j==7) || (k==0) || (k==7))\
       && ((ip==0) || (ip==7) || (jp==0) || (jp==7) || (kp==0) || (kp==7))) \
    {\
      glVertex3f(jello->p[i][j][k].x,jello->p[i][j][k].y,jello->p[i][j][k].z);\
      glVertex3f(jello->p[ip][jp][kp].x,jello->p[ip][jp][kp].y,jello->p[ip][jp][kp].z);\
    }\

 
  if (viewingMode==0) // render wireframe
  {
    glLineWidth(1);
    glPointSize(5);
    glDisable(GL_LIGHTING);
    for (i=0; i<=7; i++)
      for (j=0; j<=7; j++)
        for (k=0; k<=7; k++)
        {
          if (i*j*k*(7-i)*(7-j)*(7-k) != 0) // not surface point
            continue;

          glBegin(GL_POINTS); // draw point
            glColor4f(0,0,0,0);  
            glVertex3f(jello->p[i][j][k].x,jello->p[i][j][k].y,jello->p[i][j][k].z);        
          glEnd();

          //
          //if ((i!=7) || (j!=7) || (k!=7))
          //  continue;

          glBegin(GL_LINES);      
          // structural
          if (structural == 1)
          {
            glColor4f(0,0,1,1);
            PROCESS_NEIGHBOUR(1,0,0);
            PROCESS_NEIGHBOUR(0,1,0);
            PROCESS_NEIGHBOUR(0,0,1);
            PROCESS_NEIGHBOUR(-1,0,0);
            PROCESS_NEIGHBOUR(0,-1,0);
            PROCESS_NEIGHBOUR(0,0,-1);
          }
          
          // shear
          if (shear == 1)
          {
            glColor4f(0,1,0,1);
            PROCESS_NEIGHBOUR(1,1,0);
            PROCESS_NEIGHBOUR(-1,1,0);
            PROCESS_NEIGHBOUR(-1,-1,0);
            PROCESS_NEIGHBOUR(1,-1,0);
            PROCESS_NEIGHBOUR(0,1,1);
            PROCESS_NEIGHBOUR(0,-1,1);
            PROCESS_NEIGHBOUR(0,-1,-1);
            PROCESS_NEIGHBOUR(0,1,-1);
            PROCESS_NEIGHBOUR(1,0,1);
            PROCESS_NEIGHBOUR(-1,0,1);
            PROCESS_NEIGHBOUR(-1,0,-1);
            PROCESS_NEIGHBOUR(1,0,-1);

            PROCESS_NEIGHBOUR(1,1,1)
            PROCESS_NEIGHBOUR(-1,1,1)
            PROCESS_NEIGHBOUR(-1,-1,1)
            PROCESS_NEIGHBOUR(1,-1,1)
            PROCESS_NEIGHBOUR(1,1,-1)
            PROCESS_NEIGHBOUR(-1,1,-1)
            PROCESS_NEIGHBOUR(-1,-1,-1)
            PROCESS_NEIGHBOUR(1,-1,-1)
          }
          
          // bend
          if (bend == 1)
          {
            glColor4f(1,0,0,1);
            PROCESS_NEIGHBOUR(2,0,0);
            PROCESS_NEIGHBOUR(0,2,0);
            PROCESS_NEIGHBOUR(0,0,2);
            PROCESS_NEIGHBOUR(-2,0,0);
            PROCESS_NEIGHBOUR(0,-2,0);
            PROCESS_NEIGHBOUR(0,0,-2);
          }           
          glEnd();
        }
    glEnable(GL_LIGHTING);
  }
  
  else
  {
    glPolygonMode(GL_FRONT, GL_FILL); 
    
    for (face=1; face <= 6; face++) 
      // face == face of a cube
      // 1 = bottom, 2 = front, 3 = left, 4 = right, 5 = far, 6 = top
    {
      
      if ((face==1) || (face==3) || (face==5))
        faceFactor=-1; // flip orientation
      else
        faceFactor=1;
      

      for (i=0; i <= 7; i++) // reset buffers
        for (j=0; j <= 7; j++)
        {
          normal[i][j].x=0;normal[i][j].y=0;normal[i][j].z=0;
          counter[i][j]=0;
        }

      /* process triangles, accumulate normals for Gourad shading */
  
      for (i=0; i <= 6; i++)
        for (j=0; j <= 6; j++) // process block (i,j)
        {
          pDIFFERENCE(NODE(face,i+1,j),NODE(face,i,j),r1); // first triangle
          pDIFFERENCE(NODE(face,i,j+1),NODE(face,i,j),r2);
          CROSSPRODUCTp(r1,r2,r3); pMULTIPLY(r3,faceFactor,r3);
          pNORMALIZE(r3);
          pSUM(normal[i+1][j],r3,normal[i+1][j]);
          counter[i+1][j]++;
          pSUM(normal[i][j+1],r3,normal[i][j+1]);
          counter[i][j+1]++;
          pSUM(normal[i][j],r3,normal[i][j]);
          counter[i][j]++;

          pDIFFERENCE(NODE(face,i,j+1),NODE(face,i+1,j+1),r1); // second triangle
          pDIFFERENCE(NODE(face,i+1,j),NODE(face,i+1,j+1),r2);
          CROSSPRODUCTp(r1,r2,r3); pMULTIPLY(r3,faceFactor,r3);
          pNORMALIZE(r3);
          pSUM(normal[i+1][j],r3,normal[i+1][j]);
          counter[i+1][j]++;
          pSUM(normal[i][j+1],r3,normal[i][j+1]);
          counter[i][j+1]++;
          pSUM(normal[i+1][j+1],r3,normal[i+1][j+1]);
          counter[i+1][j+1]++;
        }

      
        /* the actual rendering */
        for (j=1; j<=7; j++) 
        {

          if (faceFactor  > 0)
            glFrontFace(GL_CCW); // the usual definition of front face
          else
            glFrontFace(GL_CW); // flip definition of orientation
         
          glBegin(GL_TRIANGLE_STRIP);
          for (i=0; i<=7; i++)
          {
            glNormal3f(normal[i][j].x / counter[i][j],normal[i][j].y / counter[i][j],
              normal[i][j].z / counter[i][j]);
            glVertex3f(NODE(face,i,j).x, NODE(face,i,j).y, NODE(face,i,j).z);
            glNormal3f(normal[i][j-1].x / counter[i][j-1],normal[i][j-1].y/ counter[i][j-1],
              normal[i][j-1].z / counter[i][j-1]);
            glVertex3f(NODE(face,i,j-1).x, NODE(face,i,j-1).y, NODE(face,i,j-1).z);
          }
          glEnd();
        }
        
        
    }  
  } // end for loop over faces
  glFrontFace(GL_CCW);
}

void showBoundingBox()
{
  int i,j;

  glColor4f(0.6,0.6,0.6,0);

  glBegin(GL_LINES);

  // front face
  for(i=-2; i<=2; i++)
  {
    glVertex3f(i,-2,-2);
    glVertex3f(i,-2,2);
  }
  for(j=-2; j<=2; j++)
  {
    glVertex3f(-2,-2,j);
    glVertex3f(2,-2,j);
  }

  // back face
  for(i=-2; i<=2; i++)
  {
    glVertex3f(i,2,-2);
    glVertex3f(i,2,2);
  }
  for(j=-2; j<=2; j++)
  {
    glVertex3f(-2,2,j);
    glVertex3f(2,2,j);
  }

  // left face
  for(i=-2; i<=2; i++)
  {
    glVertex3f(-2,i,-2);
    glVertex3f(-2,i,2);
  }
  for(j=-2; j<=2; j++)
  {
    glVertex3f(-2,-2,j);
    glVertex3f(-2,2,j);
  }

  // right face
  for(i=-2; i<=2; i++)
  {
    glVertex3f(2,i,-2);
    glVertex3f(2,i,2);
  }
  for(j=-2; j<=2; j++)
  {
    glVertex3f(2,-2,j);
    glVertex3f(2,2,j);
  }
  
  glEnd();

  return;
}

void showPlane(struct world * jello)
{
	// define 8 points
	// front surface
	double frontLeft;
	double frontRight;
	// back surface
	double backLeft;
	double backRight;
	// up surface
	double upFront;
	double upBack;
	// bottom surface
	double bottomFront;
	double bottomBack;
	point p[8];
	bool intersectPoint[8];
	
	for (int i = 0; i < 8; i++)
	{
		intersectPoint[i] = false;
	}

	// check if the front surface has intersections
	// front left
	p[0].x = -2;
	frontLeft = -((jello->a * (-2.0) + jello->c * (-2.0) + jello->d) / jello->b);
	p[0].z = -2;
	if ((frontLeft >= -2 && frontLeft <= 2))
	{
		intersectPoint[0] = true;
		p[0].y = frontLeft;
	}
	// front right
	p[1].x = 2;
	frontRight = -((jello->a * (2.0) + jello->c * (-2.0) + jello->d) / jello->b);
	p[1].z = -2;
	if ((frontRight >= -2 && frontRight <= 2))
	{
		intersectPoint[1] = true;
		p[1].y = frontRight;
	}

	// check the back surface has intersections
	// back left
	p[2].x = -2;
	backLeft = -((jello->a * (-2.0) + jello->c * (2.0) + jello->d) / jello->b);
	p[2].z = 2;
	if ((backLeft >= -2 && backLeft <= 2))
	{
		intersectPoint[2] = true;
		p[2].y = backLeft;
	}
	// back right
	p[3].x = 2;
	backRight = -((jello->a * (2.0) + jello->c * (2.0) + jello->d) / jello->b);
	p[3].z = 2;
	if ((backRight >= -2 && backRight <= 2))
	{
		intersectPoint[3] = true;
		p[3].y = backRight;
	}

	// check the up surface has intersections
	// up front
	upFront = -((jello->b * (2.0) + jello->c * (-2.0) + jello->d) / jello->a);
	p[4].y = 2;
	p[4].z = -2;
	if ((upFront > -2 && upFront < 2))
	{
		intersectPoint[4] = true;
		p[4].x = upFront;
	}
	// up back
	upBack = -((jello->b * (2.0) + jello->c * (2.0) + jello->d) / jello->a);
	p[5].y = 2;
	p[5].z = 2;
	if ((upBack > -2 && upBack < 2))
	{
		intersectPoint[5] = true;
		p[5].x = upBack;
	}

	// check the bottom surface has intersections
	// bottom front
	bottomFront = -((jello->b * (-2.0) + jello->c * (-2.0) + jello->d) / jello->a);
	p[6].y = -2;
	p[6].z = -2;
	if ((bottomFront > -2 && bottomFront < 2))
	{
		intersectPoint[6] = true;
		p[6].x = bottomFront;
	}
	// bottom back
	bottomBack = -((jello->b * (-2.0) + jello->c * (2.0) + jello->d) / jello->a);
	p[7].y = -2;
	p[7].z = 2;
	if ((bottomBack > -2 && bottomBack < 2))
	{
		intersectPoint[7] = true;
		p[7].x = bottomBack;
	}

	// check how many points intersect with the bounding box
	int pointCount = 0;
	std::vector<point> drawPoint;
	for (int i = 0; i < 8; i++)
	{
		if (intersectPoint[i])
		{
			pointCount++;
			drawPoint.push_back(p[i]);
		}
	}

	// draw triangle
	if (pointCount == 3)
	{
		glBegin(GL_POLYGON);
		glColor4f(0, 0, 1, 1);
		glVertex3f(drawPoint[2].x, drawPoint[2].y, drawPoint[2].z); // 1st vertex
		glColor4f(0, 1, 0, 1);
		glVertex3f(drawPoint[1].x, drawPoint[1].y, drawPoint[1].z); // 2nd vertex
		glColor4f(1, 0, 0, 1);
		glVertex3f(drawPoint[0].x, drawPoint[0].y, drawPoint[0].z); // last vertex
		glEnd();
	}
	
	// otherwise, draw quad
	else
	{
		glBegin(GL_QUADS);
		glColor4f(0, 0, 1, 1);
		glVertex3f(p[3].x, backRight, p[3].z);
		glVertex3f(p[0].x, frontLeft, p[0].z);
		glColor4f(0, 1, 0, 1);
		glVertex3f(p[1].x, frontRight, p[1].z);
		glVertex3f(p[2].x, backLeft, p[2].z);
		glEnd();
	}
}