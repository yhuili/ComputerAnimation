#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <algorithm>
#include "motion.h"
#include "interpolator.h"
#include "types.h"

Interpolator::Interpolator()
{
  //Set default interpolation type
  m_InterpolationType = LINEAR;

  //set default angle representation to use for interpolation
  m_AngleRepresentation = EULER;
}

Interpolator::~Interpolator()
{
}

//Create interpolated motion
void Interpolator::Interpolate(Motion * pInputMotion, Motion ** pOutputMotion, int N) 
{
  //Allocate new motion
  *pOutputMotion = new Motion(pInputMotion->GetNumFrames(), pInputMotion->GetSkeleton()); 

  //Perform the interpolation
  if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == EULER))
    LinearInterpolationEuler(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == QUATERNION))
    LinearInterpolationQuaternion(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == EULER))
    BezierInterpolationEuler(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == QUATERNION))
    BezierInterpolationQuaternion(pInputMotion, *pOutputMotion, N);
  else
  {
    printf("Error: unknown interpolation / angle representation type.\n");
    exit(1);
  }
}

void Interpolator::LinearInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

  int startKeyframe = 0;
  while (startKeyframe + N + 1 < inputLength)
  {
    int endKeyframe = startKeyframe + N + 1;

    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    // interpolate in between
    for(int frame=1; frame<=N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N+1);

      // interpolate root position
      interpolatedPosture.root_pos = startPosture->root_pos * (1-t) + endPosture->root_pos * t;

      // interpolate bone rotations
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
        interpolatedPosture.bone_rotation[bone] = startPosture->bone_rotation[bone] * (1-t) + endPosture->bone_rotation[bone] * t;

      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;
  }

  for(int frame=startKeyframe+1; frame<inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::Rotation2Euler(double R[9], double angles[3])
{
  double cy = sqrt(R[0]*R[0] + R[3]*R[3]);

  if (cy > 16*DBL_EPSILON) 
  {
    angles[0] = atan2(R[7], R[8]);
    angles[1] = atan2(-R[6], cy);
    angles[2] = atan2(R[3], R[0]);
  } 
  else 
  {
    angles[0] = atan2(-R[5], R[4]);
    angles[1] = atan2(-R[6], cy);
    angles[2] = 0;
  }

  for(int i=0; i<3; i++)
    angles[i] *= 180 / M_PI;
}

double Interpolator::degree2Radian(double degree)
{
	return (double)(degree / 180) * M_PI;
}

void Interpolator::matrixMul(const RotMatrix matrix_top, const RotMatrix matrix_new, RotMatrix matrix_out)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j =  0; j < 3; j++)
		{
			matrix_out[i][j] = 0;
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j =  0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				matrix_out[i][j] = matrix_out[i][j] + matrix_top[i][k] * matrix_new[k][j];
			}
		}
	}
}

void Interpolator::Euler2Rotation(double angles[3], double R[9])
{
	//trans degree to radian
	double radian0 = degree2Radian(angles[0]);
	double radian1 = degree2Radian(angles[1]);
	double radian2 = degree2Radian(angles[2]);

	// rotation matrix
	RotMatrix rot0 = { 
		1.0,	0.0,	0.0,
		0.0,	cos(radian0),	-sin(radian0),
		0.0,	sin(radian0),	cos(radian0)
	};

	RotMatrix rot1 = { 
		cos(radian1),	0.0,	sin(radian1),
		0.0,	1.0,	0.0, 
		-sin(radian1),	0.0,	cos(radian1)
	};

    RotMatrix rot2 = { 
		cos(radian2),	-sin(radian2),	0.0,
		sin(radian2),	cos(radian2),	0.0,
		0.0,	0.0,	1.0
    };

	// matrix multiplication
	RotMatrix temp;
	RotMatrix out;
	
	matrixMul(rot2, rot1, temp); // r2 * r1
	matrixMul(temp, rot0, out); // r2 * r1 * r0
	for (int i = 0; i < 3; i++)
	{
		for (int j =  0; j < 3; j++)
		{
			// write array out to R
			R[j + 3 * i] = out[i][j];
		}
	}
}

void Interpolator::BezierInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
	// Bezier Interpolation Euler implementation
	int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

	int startKeyframe = 0;
	// key points
	vector p0;
	vector p1;
	vector p2;
	vector p3;
	// control points
	vector a1;
	vector b2;
	// bone orientation key points
	vector bp0;
	vector bp1;
	vector bp2;
	vector bp3;
	// bone orientation control points
	vector ba1;
	vector bb2;

	while (startKeyframe + N + 1 < inputLength)
	{
		int endKeyframe = startKeyframe + N + 1;
		int previousKeyframe = startKeyframe - N - 1;
		int nextEndKeyframe = endKeyframe + N + 1;

		Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
		Posture * endPosture = pInputMotion->GetPosture(endKeyframe);
		Posture * previousPosture = NULL;
		Posture * nextPosture = NULL;
		// copy start and end keyframe
		pOutputMotion->SetPosture(startKeyframe, *startPosture);
		pOutputMotion->SetPosture(endKeyframe, *endPosture);

		if (previousKeyframe >= 0)
		{
			previousPosture = pInputMotion->GetPosture(previousKeyframe);
		}

		if (nextEndKeyframe < inputLength)
		{
			nextPosture = pInputMotion->GetPosture(nextEndKeyframe);
			pOutputMotion->SetPosture(nextEndKeyframe, *nextPosture);
		}

		// interpolate in between
		for (int frame = 1; frame <= N; frame++)
		{
			Posture interpolatedPosture;
			double t = 1.0 * frame / (N + 1);

			// interpolate root position
			// the root key points position p1 and p2
			p1 = startPosture->root_pos;
			p2 = endPosture->root_pos;
			// position temp, temp1, temp2 are used to calculate control points
			vector temp;
			vector a2;
			vector temp1;
			vector temp2;
			// the first frame
			if (startKeyframe == 0)
			{
				// the first frame does not have the previous posture
				p3 = nextPosture->root_pos;
				// calculate the first control point a1
				a1 = p2 - p3 + p2;
				// set the control point to 1/3 of the original length
				a1 = (a1 - p1) * (1.0 / 3) + p1;
				
				// calculate the second control point b2
				temp = p2 - p1 + p2;
				a2 = (temp + p3) * 0.5;
				// set the control point to 1/3 of the original length
				a2 = (a2 - p2) * (1.0 / 3) + p2;
				b2 = p2 - a2 + p2;
			}
			// the last frame
			else if (nextEndKeyframe > inputLength)
			{
				// the last frame does not have the next posture
				p0 = previousPosture->root_pos;
				// calculate the first control point a1
				temp = p1 - p0 + p1;
				a1 = (temp + p2) * 0.5;
				// set the control point to 1/3 of the original length
				a1 = (a1 - p1) * (1.0 / 3) + p1;

				// calculate the second control point b2
				a2 = p1 - p0 + p1;
				// set the control point to 1/3 of the original length
				b2 = (a2 - p2) * (1.0 / 3) + p2;
			}
			// the intermediate frame
			else
			{
				p0 = previousPosture->root_pos;
				p3 = nextPosture->root_pos;
				// calculate the first control point a1
				temp1 = p1 - p0 + p1;
				a1 = (temp1 + p2) * 0.5;
				// set the control point to 1/3 of the original length
				a1 = (a1 - p1) * (1.0 / 3) + p1;

				// calculate the second control point b2
				temp2 = p2 - p1 + p2;
				a2 = (temp2 + p3) * 0.5;
				// set the control point to 1/3 of the original length
				a2 = (a2 - p2) * (1.0 / 3) + p2;
				b2 = p2 - a2 + p2;
			}
			interpolatedPosture.root_pos = DeCasteljauEuler(t, p1, a1, b2, p2);

			// interpolate bone rotations
			for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
			{
				// the key points bp1 and bp2 based on bone rotation
				bp1 = startPosture->bone_rotation[bone];
				bp2 = endPosture->bone_rotation[bone];
				// btemp, btemp1, btemp2 are used to calculate control points
				vector btemp;
				vector ba2;
				vector btemp1;
				vector btemp2;
				// the first frame
				if (startKeyframe == 0)
				{
					// the first frame does not have the previous posture
					bp3 = nextPosture->bone_rotation[bone];
					// calculate the first control point ba1 based on the orientation
					ba1 = bp2 - bp3 + bp2;
					// set the control point to 1/3 of the original length
					ba1 = (ba1 - bp1) * (1.0 / 3) + bp1;
				
					// calculate the second control point bb2 based on the orientation
					btemp = bp2 - bp1 + bp2;
					ba2 = (btemp + bp3) * 0.5;
					// set the control point to 1/3 of the original length
					ba2 = (ba2 - bp2) * (1.0 / 3) + bp2;
					bb2 = bp2 - ba2 + bp2;
				}
				// the last frame
				else if (nextEndKeyframe > inputLength)
				{
					// the last frame does not have the next posture
					bp0 = previousPosture->bone_rotation[bone];
					// calculate the first control point ba1 based on the orientation
					btemp = bp1 - bp0 + bp1;
					ba1 = (btemp + bp2) * 0.5;
					// set the control point to 1/3 of the original length
					ba1 = (ba1 - bp1) * (1.0 / 3) + bp1;

					// calculate the second control point bb2 based on the orientation
					ba2 = bp1 - bp0 + bp1;
					// set the control point to 1/3 of the original length
					bb2 = (ba2 - bp2) * (1.0 / 3) + bp2;
				}
				// the intermediate frame
				else
				{
					bp0 = previousPosture->bone_rotation[bone];
					bp3 = nextPosture->bone_rotation[bone];
					// calculate the first control point ba1 based on the orientation
					btemp1 = bp1 - bp0 + bp1;
					ba1 = (btemp1 + bp2) * 0.5;
					// set the control point to 1/3 of the original length
					ba1 = (ba1 - bp1) * (1.0 / 3) + bp1;

					// calculate the second control point bb2 based on the orientation
					btemp2 = bp2 - bp1 + bp2;
					ba2 = (btemp2 + bp3) * 0.5;
					// set the control point to 1/3 of the original length
					ba2 = (ba2 - bp2) * (1.0 / 3) + bp2;
					bb2 = bp2 - ba2 + bp2;
				}
				interpolatedPosture.bone_rotation[bone] = DeCasteljauEuler(t, bp1, ba1, bb2, bp2);
			}
			pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
		}
		startKeyframe = endKeyframe;
	}

	for( int frame = startKeyframe + 1; frame < inputLength; frame++)
	{
		pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
	}
}

void Interpolator::LinearInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
	// Linear Interpolation Quaternion implementation
	int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

	int startKeyframe = 0;
	while (startKeyframe + N + 1 < inputLength)
	{
		int endKeyframe = startKeyframe + N + 1;

		Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
		Posture * endPosture = pInputMotion->GetPosture(endKeyframe);

		// copy start and end keyframe
		pOutputMotion->SetPosture(startKeyframe, *startPosture);
		pOutputMotion->SetPosture(endKeyframe, *endPosture);

		// interpolate in between
		for (int frame = 1; frame <= N; frame++)
		{
			Posture interpolatedPosture;
			double t = 1.0 * frame / (N + 1);

			// interpolate root position
			interpolatedPosture.root_pos = startPosture->root_pos * (1 - t) + endPosture->root_pos * t;

			// initialize the start, end and interpolate orientation on every bone in the form of quaternion
			Quaternion<double> bStart, bEnd, bInterpolate;
			// interpolate bone rotations
			for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
			{
				// convert from euler angles to quaternion
				Euler2Quaternion(startPosture->bone_rotation[bone].p, bStart);
				Euler2Quaternion(endPosture->bone_rotation[bone].p, bEnd);
				// interpolate bone orientation in the form of quaternion
				bInterpolate = Slerp(t, bStart, bEnd);
				// convert from quaternion back to euler angles
				Quaternion2Euler(bInterpolate, interpolatedPosture.bone_rotation[bone].p);
			}
			pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
		}
		startKeyframe = endKeyframe;
	}

	for (int frame = startKeyframe + 1; frame < inputLength; frame++)
	{
		pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
	}
}

void Interpolator::BezierInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
	// Bezier Interpolation Quaternion implementation
	int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

	int startKeyframe = 0;
	// key points
	vector p0;
	vector p1;
	vector p2;
	vector p3;
	// control points
	vector a1;
	vector b2;
	// bone orientation key points on unit sphere in forms of quaternion
	Quaternion<double> bp0;
	Quaternion<double> bp1;
	Quaternion<double> bp2;
	Quaternion<double> bp3;
	// bone orientation control points on unit sphere in forms of quaternion
	Quaternion<double> ba1;
	Quaternion<double> bb2;
	// interpolated orientation on unit sphere in forms of quaternion
	Quaternion<double> interpolate;

	while (startKeyframe + N + 1 < inputLength)
	{
		int endKeyframe = startKeyframe + N + 1;
		int previousKeyframe = startKeyframe - N - 1;
		int nextEndKeyframe = endKeyframe + N + 1;

		Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
		Posture * endPosture = pInputMotion->GetPosture(endKeyframe);
		Posture * previousPosture = NULL;
		Posture * nextPosture = NULL;
		// copy start and end keyframe
		pOutputMotion->SetPosture(startKeyframe, *startPosture);
		pOutputMotion->SetPosture(endKeyframe, *endPosture);

		if (previousKeyframe >= 0)
		{
			previousPosture = pInputMotion->GetPosture(previousKeyframe);
		}

		if (nextEndKeyframe <= inputLength)
		{
			nextPosture = pInputMotion->GetPosture(nextEndKeyframe);
			pOutputMotion->SetPosture(nextEndKeyframe, *nextPosture);
		}

		// interpolate in between
		for (int frame = 1; frame <= N; frame++)
		{
			Posture interpolatedPosture;
			double t = 1.0 * frame / (N+1);

			// interpolate root position
			// the key points position p1 and p2
			p1 = startPosture->root_pos;
			p2 = endPosture->root_pos;
			// position temp, temp1, temp2 are used to calculate control points
			vector temp;
			vector a2;
			vector temp1;
			vector temp2;

			// the first frame
			if (startKeyframe == 0)
			{
				// the first frame does not have the previous posture
				p3 = nextPosture->root_pos;
				// calculate the first control point a1
				a1 = p2 - p3 + p2;
				// set the control point to 1/3 of the original length
				a1 = (a1 - p1) * (1.0 / 3.0) + p1;

				// calculate the second control point b2
				temp = p2 - p1 + p2;
				a2 = (temp + p3) * 0.5;
				// set the control point to 1/3 of the original length
				a2 = (a2 - p2) * (1.0 / 3.0) + p2;
				b2 = p2 - a2 + p2;
			}
			// the last frame
			else if (nextEndKeyframe > inputLength)
			{
				// the last frame does not have the next posture
				p0 = previousPosture->root_pos;
				// calculate the first control point a1
				temp = p1 - p0 + p1;
				a1 = (temp + p2) * 0.5;
				// set the control point to 1/3 of the original length
				a1 = (a1 - p1) * (1.0 / 3.0) + p1;

				// calculate the second control point b2
				a2 = p1 - p0 + p1;
				// set the control point to 1/3 of the original length
				b2 = (a2 - p2) * (1.0 / 3.0) + p2;
			}
			// the intermediate frame
			else
			{
				p0 = previousPosture->root_pos;
				p3 = nextPosture->root_pos;
				// calculate the first control point a1
				temp1 = p1 - p0 + p1;
				a1 = (temp1 + p2) * 0.5;
				// set the control point to 1/3 of the original length
				a1 = (a1 - p1) * (1.0 / 3.0) + p1;

				// calculate the second control point b2
				temp2 = p2 - p1 + p2;
				a2 = (temp2 + p3) * 0.5;
				// set the control point to 1/3 of the original length
				a2 = (a2 - p2) * (1.0 / 3.0) + p2;
				b2 = p2 - a2 + p2;
			}
			interpolatedPosture.root_pos = DeCasteljauEuler(t, p1, a1, b2, p2);

			// interpolate bone rotations
			for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
			{
				// the key points bp1 and bp2 based on bone rotation in forms of quaternion
				Euler2Quaternion(startPosture->bone_rotation[bone].p, bp1);
				Euler2Quaternion(endPosture->bone_rotation[bone].p, bp2);
				// btemp, btemp1, btemp2 are used to calculate control points
				Quaternion<double> btemp;
				Quaternion<double> ba2;
				Quaternion<double> btemp1;
				Quaternion<double> btemp2;

				// the first frame
				if (startKeyframe == 0)
				{
					// the first frame does not have the previous posture
					Euler2Quaternion(nextPosture->bone_rotation[bone].p, bp3);
					// calculate the first control point ba1 based on the orientation in forms of quaternion
					ba1 = Double(bp3, bp2);
					// set the control point to 1/3 of the original length
					ba1 = Slerp((1.0 / 3.0), bp1, ba1);
				
					// calculate the second control point bb2 based on the orientation in forms of quaternion
					btemp = Double(bp1, bp2);
					ba2 = Slerp(0.5, btemp, bp3);
					// set the control point to 1/3 of the original length
					bb2 = Slerp((-1.0 / 3.0), bp2, ba2);
				}
				// the last frame
				else if (nextEndKeyframe > inputLength)
				{
					// the last frame does not have the next posture
					Euler2Quaternion(previousPosture->bone_rotation[bone].p, bp0);
					// calculate the first control point ba1 based on the orientation in forms of quaternion
					btemp = Double(bp0, bp1);
					ba1 = Slerp(0.5, btemp, bp2);
					// set the control point to 1/3 of the original length
					ba1 = Slerp((1.0 / 3.0), bp1, ba1);

					// calculate the second control point bb2 based on the orientation in forms of quaternion
					ba2 = Double(bp0, bp1);
					// set the control point to 1/3 of the original length
					bb2 = Slerp((1.0 / 3.0), bp2, ba2);
				}
				// the intermediate frame
				else
				{
					Euler2Quaternion(previousPosture->bone_rotation[bone].p, bp0);
					Euler2Quaternion(nextPosture->bone_rotation[bone].p, bp3);
					// calculate the first control point ba1 based on the orientation in forms of quaternion
					btemp1 = Double(bp0, bp1);
					ba1 = Slerp(0.5, btemp1, bp2);
					// set the control point to 1/3 of the original length
					ba1 = Slerp((1.0 / 3.0), bp1, ba1);

					// calculate the second control point bb2 based on the orientation in forms of quaternion
					btemp2 = Double(bp1, bp2);
					ba2 = Slerp(0.5, btemp2, bp3);
					// set the control point to 1/3 of the original length
					bb2 = Slerp((-1.0 / 3.0), bp2, ba2);
				}
				interpolate = DeCasteljauQuaternion(t, bp1, ba1, bb2, bp2);
				Quaternion2Euler(interpolate, interpolatedPosture.bone_rotation[bone].p);
			}
			pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
		}
		startKeyframe = endKeyframe;
	}

	for( int frame = startKeyframe + 1; frame < inputLength; frame++)
	{
		pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
	}
}

void Interpolator::Euler2Quaternion(double angles[3], Quaternion<double> & q) 
{
	// first, convert from euler angles to rotation matrix
	double R[9];
	Euler2Rotation(angles, R);
	// second, convert from rotation matrix to quaternion
	q = Quaternion<double>::Matrix2Quaternion(R);
	q.Normalize();
}

void Interpolator::Quaternion2Euler(Quaternion<double> & q, double angles[3]) 
{
	// first, convert from quaternion to rotation matrix
	double R[9];
	q.Quaternion2Matrix(R);
	// second, convert from rotation matrix to euler angles
	Rotation2Euler(R, angles);
}

Quaternion<double> Interpolator::Slerp(double t, Quaternion<double> & qStart, Quaternion<double> & qEnd)
{
	Quaternion<double> result;
	// cos = q1 ¡¤ q2 = s1 * s2 + x1 * x2 + y1 * y2 + z1 * z2
	double cosTheta = qStart.Gets() * qEnd.Gets() + qStart.Getx() * qEnd.Getx() + qStart.Gety() * qEnd.Gety() + qStart.Getz() * qEnd.Getz();
	double angle;
	double sinTheta;

	// we always pick the shortest path
	if (cosTheta < 0)
	{
		cosTheta = cosTheta * (-1.0);
		qEnd = qEnd * (-1.0);
	}

	// calculate the angle between qStart and qEnd
	angle = acos(cosTheta);
	// sin(theta)
	sinTheta = sin(angle); 
	
	// if the start point and the end point are the same point, we do not need to interpolate, just return this point
	if(fabs(sinTheta) == 0.0)
	{
		result = qStart;
		return result;
	}

	// if the angle is too small, we implement linear interpolation instead of slerp
	if(fabs(sinTheta) < 1E-10)
	{
		result = (1 - t) * qStart + t * qEnd;
		result.Normalize();
		return result;
	}

	result = (sin((1 - t) * angle) / sin(angle)) * qStart + (sin(t * angle) / sin(angle)) * qEnd;
	result.Normalize();
	return result;
}

Quaternion<double> Interpolator::Double(Quaternion<double> p, Quaternion<double> q)
{
  Quaternion<double> result;
  result = Slerp(2.0, p, q);
  return result;
}

vector Interpolator::DeCasteljauEuler(double t, vector p0, vector p1, vector p2, vector p3)
{
	// initialize the control points
	vector result;
	vector q0, q1, q2;
	vector r0, r1;
	// apply DeCasteljau construction using euler angles
	q0 = p0 * (1 - t) + p1 * t;
	q1 = p1 * (1 - t) + p2 * t;
	q2 = p2 * (1 - t) + p3 * t;
	r0 = q0 * (1 - t) + q1 * t;
	r1 = q1 * (1 - t) + q2 * t;
	result = r0 * (1 - t) + r1 * t;
	return result;
}

Quaternion<double> Interpolator::DeCasteljauQuaternion(double t, Quaternion<double> p0, Quaternion<double> p1, Quaternion<double> p2, Quaternion<double> p3)
{
	// initialize the control points
	Quaternion<double> result;
	Quaternion<double> q0, q1, q2;
	Quaternion<double> r0, r1;
	// apply DeCasteljau construction using quaternion
	q0 = Slerp(t, p0, p1);
	q1 = Slerp(t, p1, p2);
	q2 = Slerp(t, p2, p3);
	r0 = Slerp(t, q0, q1);
	r1 = Slerp(t, q1, q2);
	result = Slerp(t, r0, r1);
	return result;
}

